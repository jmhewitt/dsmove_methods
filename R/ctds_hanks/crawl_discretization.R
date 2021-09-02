crawl_discretization = function(
  xlocs, ylocs, semi_majors = NULL, semi_minors = NULL, orientations = NULL,
  times, domain, nimputed, crawl_proj, obs_proj, 
  pred_times = NULL, delta = 30, discrete_coords, discrete_lons, discrete_lats
) {
  # Approximate posterior distribution for the 2 parameter CTDS model using the 
  # crawl model as an approximate imputation distribution (AID) for multiple 
  # imputation strategy.
  #
  # Parameters:
  #
  #  xlocs - longitudes or (projected) x-axis coordinates
  #  ylocs - latitudes or (projected) y-axis coordinates
  #  semi_majors - semi-major axis lengths for error ellipses, or NULL if 
  #    xlocs and ylocs are observed without error information
  #  semi_minors - semi-minor axis lengths for error ellipses, or NULL if 
  #    xlocs and ylocs are observed without error information
  #  orientations - orientations for error ellipses, or NULL if 
  #    xlocs and ylocs are observed without error information
  #  times - times at which xlocs and ylocs are observed
  #  domain - Raster object specifying CTDS domain
  #  crawl_proj - coordinate projection for CTDS model fitting
  #  obs_proj - coordinate projection for xlocs and ylocs
  #  nimputed - number of complete trajectories to impute via AID
  #  pred_times - times at which locations should be imputed for
  #  delta - timestep to be used for building a complete trajectory
  #  discrete_coords - lon/lat coordinates for discretized spatial domain
  #  discrete_lons - unique longitude points in grid
  #  discrete_lats - unique latitude points in grid
  # 
  #  params - list containing initial model parameters
  #  niter - number of MCMC samples to draw
  #  priors - list containing specification of model priors
  #  states - matrix of observation states/coordinates
  #  times - times at which observations are made
  #  dims - dimensions of spatial domain
  
  #
  # 1. fit AID
  #
  
  # transform crs strings into crs objects, if needed
  if(inherits(obs_proj, 'character')) { obs_proj = sp::CRS(obs_proj) }
  if(inherits(crawl_proj, 'character')) { crawl_proj = sp::CRS(crawl_proj) }
  
  # project observation coordinates
  obs_projected = spTransform(
    x = SpatialPoints(
      coords = cbind(lon = xlocs, lat = ylocs), 
      proj4string = obs_proj
    ), 
    CRSobj = crawl_proj
  )
  
  # crawl-format for error ellipse information
  error_df = argosDiag2Cov(
    Major = semi_majors, Minor = semi_minors, Orientation = orientations
  )
  
  # munge observation coordinates and times with error ellipse information
  crawl_spdf = SpatialPointsDataFrame(
    coords = obs_projected, 
    data = cbind(error_df, time = times)
  )
  
  # weak priors for crawl parameters
  prior <- function(p) {
    sum(dnorm(x = p, mean = 0, sd = 1e2, log = TRUE))
  }
  
  # MAP for crawl model parameters
  fit = crwMLE(
    mov.model = ~ 1,
    err.model = list(
      x = ~ ln.sd.x - 1,
      y = ~ ln.sd.y - 1,
      rho = ~ error.corr
    ),
    activity = NULL,
    drift = FALSE,
    data = crawl_spdf,
    Time.name = 'time',
    time.scale = 'sec',
    prior = prior,
    # fix error parameters, only estimate movement parameters
    fixPar = c(1, 1, NA, NA),
    attempts = 8,
    control = list(trace = 0),
    initialSANN = list(maxit = 1500, trace = 0)
  )
  
  #
  # 2. sample and discretize multiple imputations from AID at specific times
  #
  
  # build an object to sample from crawl imputation posterior
  simulator = crwSimulator(
    object.crwFit = fit, 
    predTime = sort(unique(c(
      pred_times, seq(from = times[1], to = tail(times,1), by = delta)
    )))
  )
  
  message('Sampling trajectories ...')
  
  # draw samples
  imputations = replicate(
    n = nimputed, 
    expr = list(crwPostIS(object.sim = simulator))
  )
  
  message('Discretizing trajectories ...')
  
  # project grid coordinates into coordinate space for crawl
  grid_projected = spTransform(
    x = SpatialPoints(coords = discrete_coords, proj4string = obs_proj), 
    CRSobj = crawl_proj
  )
  
  # push imputed trajectories to CTDS domain coordinates
  imputations_domain = lapply(imputations, function(x) {

    # map imputation points to nearest grid cell in projected space
    imputation_gridded_projected = nn2(
      query = x$alpha.sim[, c('mu.x', 'mu.y')],
      data = coordinates(grid_projected),
      k = 1
    )
    
    # back-transform coordinates
    res = data.frame(
      discrete_coords[imputation_gridded_projected$nn.idx,]
    )
    colnames(res) = c('mapped_lon', 'mapped_lat')
    
    # enrich results
    res$lon_ind = sapply(
      res$mapped_lon, function(x) which(x == discrete_lons)
    )
    res$lat_ind = sapply(
      res$mapped_lat, function(x) which(x == discrete_lats)
    )
    res$time = as.POSIXct(x = x$time, tz = 'UTC', 
                          origin = '1970-01-01 00:00.00 UTC')
    res
  })
  
  message('Summarizing trajectories ...')
  
  # extract imputation at specific timepoints
  imputations_marginal = lapply(imputations_domain, function(x) {
    # find time indices of imputed locations
    pred_inds = sapply(pred_times, function(s) which(s == x$time))
    # extract imputed data points
    x[pred_inds,]
  })
  
  # build complete trajectories from imputations
  imputations_complete = lapply(imputations_domain, function(x) {
    # L1 distance between imputation points
    coord_diffs = abs(diff(x$lon_ind)) + abs(diff(x$lat_ind))
    # notification if we cannot strictly build a complete trajectory
    if(any(coord_diffs > 1)) {
      warning(paste('Imputation has ', sum(coord_diffs > 1), 
                    ' unobserved transitions', sep = ''))
    }
    # indices of timepoints where a new location is transitioned to
    tx_inds = c(1, which(coord_diffs > 0) + 1)
    # subset down to the complete trajectory
    res = x[tx_inds, ]
    # enrich with residence time (sec) and transition direction information
    res$duration = c(diff(as.numeric(res$time)), NA)
    diff(res$lon_ind)
    
    # compute direction of movement taken to arrive at each location
    dlon = diff(res$lon_ind)
    dlat = diff(res$lat_ind)
    res$direction_of_movement = factor(c(NA, ifelse(
      dlon < 0, 'W', 
      ifelse(dlon > 0, 'E', 
             ifelse(dlat < 0, 'S', 'N')
             )
    )), levels = c('N','E','S','W'))
    # return complete trajectory
    res
  })
  
  # package results
  list(
    imputations_at_times = imputations_marginal,
    complete_trajectories = imputations_complete
  )
}