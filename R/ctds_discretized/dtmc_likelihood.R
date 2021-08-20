dtmc_likelihood = function(
  pkg, speed, cell_size, betaAR
) {
  # Marginal distribution for location given model parameters and data
  # 
  # Parameters:
  #   pkg - flattened dataset containing data to analyze (see flatten_data.R)
  #   speed - speed parameter for CTDS model
  #   cell_size - average distance traveled per grid cell transition
  #   betaAR - directional persistence parameter for CTDS model
  
  # model parameters
  theta = c(1 - speed / cell_size * pkg$delta, betaAR)
  
  ll = dsmovetools2d:::SattagFilteredLL(
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
    betaAR = theta[2]
  )
  
  # package results
  data.frame(
    speed = speed, 
    betaAR = betaAR, 
    cell_size = cell_size, 
    delta = pkg$delta, 
    ll = ll
  )
}