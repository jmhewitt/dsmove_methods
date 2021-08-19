flatten_data = function(
  obs_lon, obs_lat, obs_time, bathymetry, delta,
  max_speed_filter = NULL, max_speed_iter = 2, obs_depth = NULL, 
  obs_depth_time = NULL, obs_semi_major, obs_semi_minor, obs_orientation, 
  analysis_window = NULL
) {
  # Discretize and filter data in preparation for analysis.  Observations will 
  # be mapped to a discrete grid and filtered such that 1) observations will 
  # fall within the specified time window, and 2) pointwise maximum-speed 
  # filtering rules are satisfied.  Depth data will be filtered so that it is 
  # contained within the range of lon/lat observation times.
  # 
  # Parameters:
  #  obs_lon - vector of observed longitudes
  #  obs_lat - vector of observed latitudes
  #  obs_time - vector of lon/lat observation times (POSIXct format)
  #  max_speed_filter - if not NULL, then the maximum speed above which 
  #    observations should be removed from analysis
  #  max_speed_iter - max_speed_filter will be applied iteratively, this many
  #    times 
  #  obs_depth - vector of observed depths, or NULL if depth data will not 
  #    be incorporated into analysis
  #  obs_depth_time - vector of depth observation times (POSIXct format)
  #  obs_semi_major - vector of semi-major values for lon/lat error ellipses
  #  obs_semi_minor - vector of semi-minor values for lon/lat error ellipses
  #  obs_orientation - vector of lon/lat error ellipse orientations
  #  analysis_window - if not NULL, then lon/lat observations outside this 
  #    date range (start, end) will be filtered out of the analysis data
  #  bathymetry - raster object containing spatial domain to which movement
  #    is discretized
  #  delta - time discretization to use (sec)
  #
  # Return:
  #  a list containing the reduced dataset.  Indices of mapped lon/lat 
  #  coordinates are not 0-indexed (i.e., their range is 1...nlons).  This must 
  #  be done elsewhere if this is important (i.e., for using as C++ indices).
  
  # consolidate location data
  x = data.frame(
    lon = obs_lon,
    lat = obs_lat,
    semi_major = obs_semi_major,
    semi_minor = obs_semi_minor,
    orientation = obs_orientation,
    time = obs_time
  )
  
  # filter out observations that are missing (Error) information
  x = x[complete.cases(x),]
  
  # filter out observations beyond analysis window
  if(!is.null(analysis_window)) {
    x = x %>% filter(analysis_window[1] <= time, time <= analysis_window[2])
  }
  
  # extract discretized coordinate grid
  coords = coordinates(bathymetry)
  lons = unique(coords[,1])
  lats = unique(coords[,2])
  
  # map the observed coordinates to the bathymetry grid
  mapped_coords = map_coords(
    lon = x$lon, lat = x$lat, lon_grid = lons, lat_grid = lats, 
    coord_grid = coords
  )
  
  # extract the discretized coordinates
  x$mapped_lon_ind = mapped_coords$lon_ind
  x$mapped_lat_ind = mapped_coords$lat_ind
  x$mapped_lon = lons[mapped_coords$lon_ind]
  x$mapped_lat = lats[mapped_coords$lat_ind]
  
  # crude max-speed filtering loop
  if(all(!is.null(max_speed_filter), max_speed_iter > 0)) {
    for(i in 1:max_speed_iter) {
      # compute empirical speed between mapped lon/lat observations (m/s)
      x$mapped_speed = c(0, sapply(2:nrow(x), function(ind) {
        rdist.earth(x1 = as.matrix(x[ind-1, c('mapped_lon', 'mapped_lat')]),
                    x2 = as.matrix(x[ind, c('mapped_lon', 'mapped_lat')]),
                    miles = FALSE) * 1e3 /
          diff(as.numeric(x$time[ind + (-1:0)]))
      }))
      # filter out observations associated with large empirical speeds
      x = x %>% filter(mapped_speed <= max_speed_filter)
    }
  }
  
  #
  # flatten data for likelihood approximation
  #
  
  # "total" observation error, for assisting in choosing between multiple 
  # observation at the same discrete timestep
  x$sqerror = x$semi_major^2 + x$semi_minor^2
  
  if(!inherits(bathymetry, 'RasterLayer')) {
    stop('bathymetry must be a RasterLayer to ensure accurate data extraction.')
  }
  
  # munge bathymetry to matrix format
  zsurf = matrix(
    getValues(bathymetry), 
    nrow = length(lats), 
    ncol = length(lons),
    byrow = TRUE
  )
  
  # discretized timesteps
  tseq = seq(
    from = min(x$time),
    to = max(x$time),
    by = duration(delta, 'seconds')
  )
  
  # container for discretized data
  flattened_coords = matrix(NA, nrow = length(tseq), ncol = 9)
  colnames(flattened_coords) = c('lon_ind', 'lat_ind', 'sqerror', 'lon', 'lat',
                                 'semi_major', 'semi_minor', 'orientation', 
                                 'depth')
  
  # map lon/lat observations
  flattened_coords[, 'sqerror'] = Inf
  for(ind in 1:nrow(x)) {
    tgt_ind = which.min(abs(x$time[ind] - tseq))
    # transfer observation if it has less observation error
    if(flattened_coords[tgt_ind, 'sqerror'] > x$sqerror[ind]) {
      flattened_coords[tgt_ind, 'lon_ind'] = x$mapped_lon_ind[ind]
      flattened_coords[tgt_ind, 'lat_ind'] = x$mapped_lat_ind[ind]
      flattened_coords[tgt_ind, 'sqerror'] = x$sqerror[ind]
      flattened_coords[tgt_ind, 'lon'] = x$mapped_lon[ind]
      flattened_coords[tgt_ind, 'lat'] = x$mapped_lat[ind]
      flattened_coords[tgt_ind, 'semi_major'] = x$semi_major[ind]
      flattened_coords[tgt_ind, 'semi_minor'] = x$semi_minor[ind]
      flattened_coords[tgt_ind, 'orientation'] = x$orientation[ind]
    }
  }
  
  #
  # flatten depth data
  #
  
  # only process depth data if it is available
  if(all(!is.null(obs_depth), !is.null(obs_depth_time))) {
    
    # consolidate depth data
    x_depth = data.frame(
      depth = obs_depth,
      time = obs_depth_time
    )
    
    # filter depth data to lie within the location data date range
    obs_window = range(x$time)
    x_depth = x_depth %>% filter(obs_window[1] <= time, time <= obs_window[2])
    
    # map depth observations
    for(ind in 1:nrow(x_depth)) {
      tgt_ind = which.min(abs(x_depth$time[ind] - tseq))
      flattened_coords[tgt_ind, 'depth'] = x_depth$depth[ind]
    }
    
  } else {
    x_depth = NULL
  }
  
  #
  # package data
  #
  
  list(
    data = list(loc = x, depth = x_depth),
    grid = list(lons = lons, lats = lats, coords = coords),
    analysis_window = analysis_window,
    bathymetry_matrix = zsurf,
    mapped_data = flattened_coords[
      1:nrow(flattened_coords), setdiff(colnames(flattened_coords), 'sqerror')
    ]
  )
}