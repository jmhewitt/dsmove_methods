marginal_distn = function(
  pkg, times, speed, cell_size, betaAR
) {
  # Marginal distribution for location given model parameters and data
  # 
  # Parameters:
  #   pkg - flattened dataset containing data to analyze (see flatten_data.R)
  #   times - vector of times to compute  marginal distribution for.  times 
  #     must align with the time discretization to be valid.
  #   speed - speed parameter for CTDS model
  #   cell_size - average distance traveled per grid cell transition
  #   betaAR - directional persistence parameter for CTDS model
  
  # verify all marginal distn'. times are well defined wrt time discretization
  if(!all(times %in% pkg$mapped_times)) {
    stop('All time values must appear in the time discretization from pkg.')
  }
  
  # time-discretization indices at which marginal locations are computed
  pred_inds = sapply(times, function(t) which(t == pkg$mapped_times))
  
  # model parameters
  theta = c(1 - speed / cell_size * pkg$delta, betaAR)
  
  # forward-based predictive distributions for marginal locations
  pred_dists = dsmovetools2d:::SattagPredDist(
    init_dsts = matrix(
      pkg$mapped_data[1, c('lon_ind', 'lat_ind')] - 1,
      nrow = 4, ncol = 2,
      byrow = TRUE
    ),
    init_srcs = matrix(
      c(pkg$mapped_data[1, 'lon_ind'], pkg$mapped_data[1, 'lat_ind'] - 1,
        pkg$mapped_data[1, 'lon_ind'], pkg$mapped_data[1, 'lat_ind'] + 1,
        pkg$mapped_data[1, 'lon_ind'] - 1, pkg$mapped_data[1, 'lat_ind'],
        pkg$mapped_data[1, 'lon_ind'] + 1, pkg$mapped_data[1, 'lat_ind']) - 1,
      nrow = 4, ncol = 2, 
      byrow = TRUE
    ), 
    init_log_probs = rep(log(1/4),4), 
    gps_trunc_alpha = .05, 
    obs_lons = pkg$mapped_data[,'lon'],
    obs_lats = pkg$mapped_data[,'lat'], 
    obs_semi_majors = pkg$mapped_data[,'semi_major'], 
    obs_semi_minors = pkg$mapped_data[,'semi_minor'], 
    obs_orientations = pkg$mapped_data[,'orientation'], 
    obs_depths = pkg$mapped_data[,'depth'], 
    lon_gridvals = pkg$grid$lons, 
    lat_gridvals = pkg$grid$lats, 
    surface_heights = pkg$bathymetry_matrix, 
    min_elevation = -1e4,
    max_elevation = 0, 
    log_self_tx = log(theta[1]), 
    betaAR = theta[2], 
    pred_steps = pred_inds - 1
  )
  
  # backward information filters
  ind = nrow(pkg$mapped_data)
  last_coords = pkg$mapped_data[nrow(pkg$mapped_data), c('lon_ind', 'lat_ind')]
  bif_dists = dsmovetools2d:::BackInfoFilteringDist(
    init_dsts = matrix(
      last_coords - 1,
      nrow = 4, ncol = 2,
      byrow = TRUE
    ),
    init_srcs = matrix(
      c(last_coords[1], last_coords[2] - 1,
        last_coords[1], last_coords[2] + 1,
        last_coords[1] - 1, last_coords[2],
        last_coords[1] + 1, last_coords[2]) - 1,
      nrow = 4, ncol = 2, 
      byrow = TRUE
    ),  
    init_log_probs = rep(log(1/4),4), 
    gps_trunc_alpha = .05, 
    obs_lons = pkg$mapped_data[,'lon'],
    obs_lats = pkg$mapped_data[,'lat'], 
    obs_semi_majors = pkg$mapped_data[,'semi_major'], 
    obs_semi_minors = pkg$mapped_data[,'semi_minor'], 
    obs_orientations = pkg$mapped_data[,'orientation'], 
    obs_depths = pkg$mapped_data[,'depth'], 
    lon_gridvals = pkg$grid$lons, 
    lat_gridvals = pkg$grid$lats, 
    surface_heights = pkg$bathymetry_matrix, 
    min_elevation = -1e4,
    max_elevation = 0, 
    log_self_tx = log(theta[1]), 
    betaAR = theta[2], 
    pred_steps = pred_inds - 1
  )
  
  # convert c++ output to 1-based indexing for locations
  pred_dists = lapply(pred_dists, function(pd) {
    pd[,c('lon_from_ind', 'lat_from_ind', 'lon_to_ind', 'lat_to_ind')] = 
      pd[,c('lon_from_ind', 'lat_from_ind', 'lon_to_ind', 'lat_to_ind')] + 1
    pd
  })
  bif_dists = lapply(bif_dists, function(pd) {
    pd[,c('lon_from_ind', 'lat_from_ind', 'lon_to_ind', 'lat_to_ind')] = 
      pd[,c('lon_from_ind', 'lat_from_ind', 'lon_to_ind', 'lat_to_ind')] + 1
    pd
  })
  
  # merge distributions, yielding marginal distribution wrt CTMC states
  res = mapply(function(pred, bif) {
    inner_join(
      x = data.frame(pred),
      y = data.frame(bif),   
      by = c('lon_from_ind', 'lat_from_ind', 'lon_to_ind', 'lat_to_ind')
    ) %>% 
      mutate(
        lp = log_prob.x + log_prob.y,
        lp = lp - dsmovetools2d:::log_sum_c(lp)
      ) %>% 
      select(lon_from_ind, lat_from_ind, lon_to_ind, lat_to_ind, lp)
  }, pred_dists, bif_dists, SIMPLIFY = FALSE)
  
  # collapse distributions according to the observable location
  marginal_location = lapply(res, function(d) {
    d %>% 
      group_by(lon_to_ind, lat_to_ind) %>% 
      summarise(lp = dsmovetools2d:::log_sum_c(lp), .groups = 'keep') %>% 
      ungroup()
  })
  
  # package results
  list(
    marginal_location = marginal_location,
    params = list(speed = speed, cell_size = cell_size, betaAR = betaAR),
    times = times
  ) 
}