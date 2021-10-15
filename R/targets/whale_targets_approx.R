whale_targets_approx = list(
  
  
  # pre-process and package data for analysis
  tar_target(
    name = whale_pkg_approx,
    command = {
      # load data sets
      x = readRDS('data/BeakedWhale/zctag095_gps.rds')
      depths_raw = read.csv(
        'data/BeakedWhale/ZcTag095_DUML_series_20200703.csv'
      )
      depths_raw$Date = as.POSIXct(x = depths_raw$Date, tz = 'UTC', 
                                   origin = '1970-01-01 00:00.00 UTC')
      # subset data rules
      if(whale_pkg_subset == 'with_depth') {
        depth_data = -depths_raw$Depth
      } else {
        depth_data = NULL
      }
      # extract data
      res = sampleFromEllipses(
        obs_lon = x$Longitude, obs_lat = x$Latitude, obs_time = x$date_time, 
        bathymetry = whale_domain, delta = whale_ll_delta, max_speed_filter = 6, 
        max_speed_iter = 2,
        obs_depth = depth_data, 
        obs_depth_time = depths_raw$Date, 
        obs_semi_major = x$Error.Semi.major.axis, 
        obs_semi_minor = x$Error.Semi.minor.axis, 
        obs_orientation = x$Error.Ellipse.orientation, 
        analysis_window = cee_analysis_window
      )
      # enrich and return data
      res$subset = whale_pkg_subset
      res$tagID = 'Zc095'
      list(res)
    }, 
    pattern = map(whale_pkg_subset)
  ),
  
  tar_target(
    name = whale_imputed_approx_post,
    command = {
      browser()
      
      pkg = whale_pkg_approx[[1]]
      
      betaAR = 1
      speed = 1.6
      cell_size = 623
      
      theta = c(1 - speed / cell_size * pkg$delta, betaAR)
      
      
      dsmovetools2d:::ExactSattagFilteredLL(
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
        obs_lons = pkg$mapped_data[,'lon'],
        obs_lats = pkg$mapped_data[,'lat'], 
        obs_depths = pkg$mapped_data[,'depth'],
        lon_gridvals = pkg$grid$lons, 
        lat_gridvals = pkg$grid$lats, 
        surface_heights = pkg$bathymetry_matrix, 
        min_elevation = -1e4,
        max_elevation = 0, 
        log_self_tx = log(theta[1]), 
        betaAR = theta[2], lptrunc = Inf
      )
      head(pkg$mapped_data)
      
      
    },
    pattern = map(whale_pkg_approx)
  )
  
)
