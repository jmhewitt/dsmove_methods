whale_targets = list(
  
  #
  # gridded likelihood evaluation
  #
  
  # time-discretization step size (seconds)
  tar_target(whale_ll_delta, 30),
  
  tar_target(
    name = cee_start, 
    command = strptime(
      x = '8/19/2019  7:11:00 PM UTC', 
      format = '%m/%d/%Y %I:%M:%S %p',
      tz = 'UTC'
    )
  ),
  
  # subset of depth/location data to analyze
  tar_target(
    name = cee_analysis_window, 
    command = cee_start + duration(num = c(-1,1), units = 'days')
  ),
  
  # times at which to compute posterior distribution for location
  tar_target(
    name = cee_posterior_location_times,
    command = {
      pred_window = cee_start + duration(num = c(0,1), units = 'hours')
      pred_times = whale_pkg[[1]]$mapped_times[
        (pred_window[1] <= whale_pkg[[1]]$mapped_times) &
        (whale_pkg[[1]]$mapped_times <= pred_window[2]) &
        (is.finite(whale_pkg[[1]]$mapped_data[,'depth']))
      ]
    }
  ),
  
  # discrete spatial domain for analysis
  tar_target(
    name = whale_domain,
    command = {
      # load and crop bathymetry file
      bathy = crop(
        x = raster(file.path('data', 'STRM_30_PLUS', 'w100n40.nc')), 
        y = extent(c(min_lon = -77, max_lon = -72, min_lat = 33, max_lat = 37))
      )
      # set CRS
      crs(bathy) = 
        '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
      # return raster object
      bathy
    }
  ),
  
  # projection used to develop AID approximation to CTDS posterior
  tar_target(
    name = brs_crs,
    command = paste("+proj=aea +lat_1=27.33333333333333",
                    "+lat_2=40.66666666666666 +lat_0=34 +lon_0=-78",
                    "+x_0=0 +y_0=0",
                    "+ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  ),
  
  # labels for data subsets to analyze
  tar_target(whale_pkg_subset, c('with_depth', 'without_depth')),
  
  # pre-process and package data for analysis
  tar_target(
    name = whale_pkg,
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
      res = flatten_data(
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
    name = whale_ll_grid, 
    command = expand.grid(
      speed = c(.01, seq(from = .1, to = 6, by = .1)), # m/s
      betaAR = seq(from = 0, to = 5, length.out = 15)
    )
  ),
  
  tar_target(
    name = whale_ll_approx,
    command = {
      tick = proc.time()[3]
      res = dtmc_likelihood(
        pkg = whale_pkg[[1]], 
        speed = unlist(whale_ll_grid['speed']),
        betaAR = unlist(whale_ll_grid['betaAR']),
        cell_size = 623.7801
      )
      tock = proc.time()[3]
      res$delta = whale_pkg[[1]]$delta
      res$subset = whale_pkg[[1]]$subset
      res$computation_time = tock - tick
      res$computation_date = Sys.time()
      res$computation_node = Sys.info()['nodename']
      res
    } ,
    pattern = cross(whale_ll_grid, map(whale_pkg)),
    deployment = 'worker',
    storage = 'worker',
    retrieval = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(
    name = whale_marginal_approx,
    command = {
      # compute posterior
      tick = proc.time()[3]
      res = marginal_distn(
        pkg = whale_pkg[[1]], 
        times = cee_posterior_location_times,
        speed = unlist(whale_ll_grid['speed']),
        betaAR = unlist(whale_ll_grid['betaAR']),
        cell_size = 623.7801
      )
      tock = proc.time()[3]
      res$subset = whale_pkg[[1]]$subset
      res$delta = whale_pkg[[1]]$delta
      res$computation_time = tock - tick
      res$computation_date = Sys.time()
      res$computation_node = Sys.info()['nodename']
      # save output, return file name
      d = file.path('output', 'posterior_marginal_location', res$subset)
      dir.create(path = d, showWarnings = FALSE, recursive = TRUE)
      f = file.path(d, paste(tar_name(), '.rds', sep = ''))
      saveRDS(res, file = f)
      f
    } ,
    pattern = cross(whale_ll_grid, map(whale_pkg)),
    deployment = 'worker',
    storage = 'worker',
    retrieval = 'worker',
    memory = 'transient',
    error = 'continue'
  )
  
)
