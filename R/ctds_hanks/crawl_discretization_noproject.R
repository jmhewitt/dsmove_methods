crawl_discretization_noproject = function(
  xlocs, ylocs, times, nimputed, delta = 30
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
  
  # weak priors for crawl parameters
  prior <- function(p) {
    sum(dnorm(x = p, mean = 0, sd = 1e2, log = TRUE))
  }
  
  # MAP for crawl model parameters
  fit = crwMLE(
    mov.model = ~ 1,
    err.model = NULL,
    activity = NULL,
    drift = FALSE,
    data = data.frame(x = xlocs, y = ylocs, time = times),
    Time.name = 'time',
    time.scale = 'sec',
    prior = prior,
    # fix error parameters, only estimate movement parameters
    fixPar = c(NA, NA),
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
    predTime = seq(from = times[1], to = tail(times,1), by = delta)
  )
  
  message('Sampling trajectories ...')
  
  # draw samples
  imputations = replicate(
    n = nimputed, 
    expr = list(crwPostIS(object.sim = simulator))
  )
  
  message('Discretizing trajectories ...')
  
  # push imputed trajectories to CTDS domain coordinates
  imputations_domain = lapply(imputations, function(x) {
    data.frame(round(x$alpha.sim[, c('mu.x', 'mu.y')]),
               time = x$time)
  })
  
  message('Summarizing trajectories ...')
  
  
  # build complete trajectories from imputations
  imputations_complete = lapply(imputations_domain, function(x) {
    # L1 distance between imputation points
    coord_diffs = abs(diff(x$mu.x)) + abs(diff(x$mu.y))
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
    diff(res$mu.x)
    
    # compute direction of movement taken to arrive at each location
    dlon = diff(res$mu.x)
    dlat = diff(res$mu.y)
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
    complete_trajectories = imputations_complete
  )
}