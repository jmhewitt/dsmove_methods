whale_targets = list(
  
  #
  # gridded likelihood evaluation
  #
  
  tar_target(whale_ll_delta, 300),
  tar_target(whale_ll_delta_srtm30, 60),
  
  tar_target(
    name = whale_ll_grid, 
    command = expand.grid(
      # m/s
      speed = seq(from = .01, to = 3, length.out = 20),
      betaAR = seq(from = -3, to = 3, length.out = 20)
    )
  ),
  
  tar_target(
    name = whale_ll,
    command = {
      
      # raw bathymetry
      bathy = read.csv('data/BeakedWhale/bathyEC.nvo', header = FALSE, 
                       sep =  '')
      colnames(bathy) = c('lat', 'lon', 'depth')
      
      # conversion to gridded sp object
      g = SpatialGridDataFrame(
        grid = GridTopology(
          cellcentre.offset = c(min(bathy[,'lon']), min(bathy[,'lat'])), 
          cellsize = rep(1/3, 2), 
          cells.dim = c(length(unique(bathy[,'lon'])), 
                        length(unique(bathy[,'lat'])))
        ), 
        data = data.frame(depth = bathy[,'depth']), 
        proj4string = '+proj=longlat +datum=WGS84 +no_defs'
      )
      
      # zc095 gps observations
      x = readRDS('data/BeakedWhale/zctag095_gps.rds')
      
      # filter out observations that are missing (Error) information
      x = x[complete.cases(x),]
      
      #
      # map data to grid
      #
      
      # map the observed lat-lons to the bathymetry lon/lat grid
      # (needed to set initial filtering locations)
      coords = coordinates(g)
      lons = unique(coords[,1])
      lats = unique(coords[,2])
      mapped_coords = do.call(rbind, lapply(1:nrow(x), function(ind) {
        coord_ind = which.min(
          spDists(x = x[ind, c('Longitude','Latitude')], y = coords)
        )
        lon_lat = coordinates(g)[coord_ind,]
        c(which(lon_lat[2] == lats), which(lon_lat[1] == lons)) - 1
      }))
      colnames(mapped_coords) = c('lat_ind', 'lon_ind')
      
      #
      # likelihood approximation
      #
      
      # append dummy depth variable and dimension information
      states = cbind(mapped_coords, 0)
      dims = c(length(lats), length(lons), 1)
      # height of domain surface
      zsurf = matrix(
        -g$depth, 
        nrow = g@grid@cells.dim[1], 
        ncol = g@grid@cells.dim[2]
      )
      # dummy height of vertical layer
      zval = 0
      
      # discretized timesteps
      delta = whale_ll_delta
      tseq = seq(
        from = min(x$date_time),
        to = max(x$date_time),
        by = duration(delta, 'seconds')
      )
      
      # map observations to discretized timesteps
      obs_coords = matrix(NA, nrow = length(tseq), ncol = 3)
      for(ind in 1:nrow(x)) {
        obs_coords[which.min(abs(x$date_time[ind] - tseq)),] = states[ind,]
      }
      
      # map gps observations to discretized timesteps
      obs_gps = data.frame(
        lon = rep(NA, length(tseq)),
        lat = rep(NA, length(tseq)),
        semi_major = rep(NA, length(tseq)),
        semi_minor = rep(NA, length(tseq)),
        orientation = rep(NA, length(tseq))
      )
      for(ind in 1:nrow(x)) {
        dat_ind = which.min(abs(x$date_time[ind] - tseq))
        obs_gps$lon[dat_ind] = lons[1 + mapped_coords[ind, 'lon_ind']]
        obs_gps$lat[dat_ind] = lats[1 + mapped_coords[ind, 'lat_ind']]
        obs_gps$semi_major[dat_ind] = x$Error.Semi.major.axis[ind]
        obs_gps$semi_minor[dat_ind] = x$Error.Semi.minor.axis[ind]
        obs_gps$orientation[dat_ind] = x$Error.Ellipse.orientation[ind]
      }
      
      # map depth observations to discretized timesteps
      obs_depths = rep(NA, length(tseq))
      template_bins = read.csv('data/BeakedWhale/template_bins.csv')
      depths_raw = read.csv('data/BeakedWhale/ZcTag095_DUML_series_20200703.csv')
      depths_raw$Date = as.POSIXct(x = depths_raw$Date, tz = 'UTC', 
                                   origin = '1970-01-01 00:00.00 UTC')
      for(ind in 1:nrow(depths_raw)) {
        obs_depths[which.min(abs(depths_raw$Date[ind] - tseq))] = which.min(abs(
          depths_raw$Depth[ind] - template_bins$center
        )) - 1
      }
      
      # uniform mass over neighborhood structure for initial state
      a0_prev_coords = dsmovetools:::TestZConstrainedRookNeighborhood(
        dims = dims, 
        x = states[1,], 
        zfield = zsurf, 
        zvals = zval
      )
      log_a0val = rep(1, nrow(a0_prev_coords))
      log_a0val = log(log_a0val / sum(log_a0val))
      
      # log-likelihood for observations
      ll = function(speed, betaAR) {
        theta = c(1 - speed / 37107.96 * delta, betaAR)
        dsmovetools:::SattagFilteredLL(
          a0 = matrix(states[1,], 
                      nrow = nrow(a0_prev_coords), 
                      ncol = ncol(states), byrow = TRUE),
          a0_prev_coords = a0_prev_coords, 
          gps_trunc_alpha = .05, 
          obs_lons = obs_gps$lon, 
          obs_lats = obs_gps$lat,
          obs_semi_majors = obs_gps$semi_major,
          obs_semi_minors = obs_gps$semi_minor, 
          obs_orientations = obs_gps$orientation, 
          obs_depth_bins = obs_depths, 
          lon_gridvals = lons, 
          lat_gridvals = lats, 
          log_a0val = log_a0val, 
          dims = dims, 
          surface_heights = zsurf, domain_heights = -template_bins$center, 
          log_self_tx = log(theta[1]), betaAR = theta[2]
        ) 
      }
      
      res = data.frame(
        speed = whale_ll_grid['speed'],
        betaAR = whale_ll_grid['betaAR'],
        ll = ll(speed = unlist(whale_ll_grid['speed']),
                betaAR = unlist(whale_ll_grid['betaAR'])),
        delta = delta
      )
      
      res
    },
    pattern = map(whale_ll_grid),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(
    name = whale_ll_plot,
    command = {
      
      whale_ll = readRDS('whale_ll.rds')
      
      ggplot(whale_ll, aes(x = speed, y = betaAR, fill = ll, z = ll)) + 
        geom_raster() + 
        scale_fill_viridis(direction = -1) + 
        theme_few() + 
        stat_contour2(col = 'grey90') + 
        geom_text_contour(col = 'grey90')
      
      whale_ll$lp = whale_ll$ll - log_sum(whale_ll$ll)
      
      sum(whale_ll$speed * exp(whale_ll$lp))
      sum(whale_ll$betaAR * exp(whale_ll$lp))
      
      marginal_post_speed = whale_ll %>% 
        group_by(speed) %>%
        summarise(lp = log_sum(lp))
      
      plot(exp(lp) ~ speed, marginal_post_speed, type = 'l')
      
      # TODO: use smoothing splines, like INLA, to interpolate the posterior
      #  between grid locations
      
      # TODO: add a dx component to the smooth so that we get a proper density
    }
  ),
  
  tar_target(
    name = whale_ll_srtm30,
    command = {
      
      # set spatial extents to extract
      extent_vec = c(min_lon = -77, max_lon = -72, min_lat = 33, max_lat = 37)
      
      domain.extent = extent(extent_vec)
      
      # load, crop, and merge bathymetry files
      bathy = crop(raster('data/STRM_30_PLUS/w100n40.nc'), domain.extent)
      
      # set CRS
      crs(bathy) = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
      
      # zc095 gps observations
      x = readRDS('data/BeakedWhale/zctag095_gps.rds')
      
      # filter out observations that are missing (Error) information
      x = x[complete.cases(x),]
      
      #
      # map data to grid
      #
      
      # map the observed lat-lons to the bathymetry lon/lat grid
      # (needed to set initial filtering locations)
      coords = coordinates(bathy)
      lons = unique(coords[,1])
      lats = unique(coords[,2])
      mapped_coords = do.call(rbind, lapply(1:nrow(x), function(ind) {
        coord_ind = which.min(
          spDists(x = x[ind, c('Longitude','Latitude')], y = coords)
        )
        lon_lat = coordinates(bathy)[coord_ind,]
        c(which(lon_lat[2] == lats), which(lon_lat[1] == lons)) - 1
      }))
      colnames(mapped_coords) = c('lat_ind', 'lon_ind')
      
      # load depth data
      template_bins = read.csv('data/BeakedWhale/template_bins.csv')
      depths_raw = read.csv('data/BeakedWhale/ZcTag095_DUML_series_20200703.csv')
      depths_raw$Date = as.POSIXct(x = depths_raw$Date, tz = 'UTC', 
                                   origin = '1970-01-01 00:00.00 UTC')
      
      #
      # likelihood approximation
      #
      
      # append dummy depth variable and dimension information
      states = cbind(mapped_coords, 0)
      dims = c(length(lats), length(lons), 1)
      # height of domain surface
      zsurf = matrix(
        getValues(t(bathy)), 
        nrow = nrow(bathy), 
        ncol = ncol(bathy)
      )
      # dummy height of vertical layer
      zval = 0
      
      # discretized timesteps
      delta = whale_ll_delta_srtm30
      tseq = seq(
        from = min(x$date_time),
        to = max(x$date_time),
        by = duration(delta, 'seconds')
      )
      
      # map observations to discretized timesteps
      obs_coords = matrix(NA, nrow = length(tseq), ncol = 3)
      for(ind in 1:nrow(x)) {
        obs_coords[which.min(abs(x$date_time[ind] - tseq)),] = states[ind,]
      }
      
      # map gps observations to discretized timesteps
      obs_gps = data.frame(
        lon = rep(NA, length(tseq)),
        lat = rep(NA, length(tseq)),
        semi_major = rep(NA, length(tseq)),
        semi_minor = rep(NA, length(tseq)),
        orientation = rep(NA, length(tseq))
      )
      for(ind in 1:nrow(x)) {
        dat_ind = which.min(abs(x$date_time[ind] - tseq))
        obs_gps$lon[dat_ind] = lons[1 + mapped_coords[ind, 'lon_ind']]
        obs_gps$lat[dat_ind] = lats[1 + mapped_coords[ind, 'lat_ind']]
        obs_gps$semi_major[dat_ind] = x$Error.Semi.major.axis[ind]
        obs_gps$semi_minor[dat_ind] = x$Error.Semi.minor.axis[ind]
        obs_gps$orientation[dat_ind] = x$Error.Ellipse.orientation[ind]
      }
      
      # map depth observations to discretized timesteps
      obs_depths = rep(NA, length(tseq))
      for(ind in 1:nrow(depths_raw)) {
        obs_depths[which.min(abs(depths_raw$Date[ind] - tseq))] = which.min(abs(
          depths_raw$Depth[ind] - template_bins$center
        )) - 1
      }
      
      # uniform mass over neighborhood structure for initial state
      a0_prev_coords = dsmovetools:::TestZConstrainedRookNeighborhood(
        dims = dims, 
        x = states[1,], 
        zfield = zsurf, 
        zvals = zval
      )
      log_a0val = rep(1, nrow(a0_prev_coords))
      log_a0val = log(log_a0val / sum(log_a0val))
      
      # log-likelihood for observations
      ll = function(speed, betaAR) {
        theta = c(1 - speed / 1187.295 * delta, betaAR)
        dsmovetools:::SattagFilteredLL(
          a0 = matrix(states[1,], 
                      nrow = nrow(a0_prev_coords), 
                      ncol = ncol(states), byrow = TRUE),
          a0_prev_coords = a0_prev_coords, 
          gps_trunc_alpha = .05, 
          obs_lons = obs_gps$lon, 
          obs_lats = obs_gps$lat,
          obs_semi_majors = obs_gps$semi_major,
          obs_semi_minors = obs_gps$semi_minor, 
          obs_orientations = obs_gps$orientation, 
          obs_depth_bins = obs_depths, 
          lon_gridvals = lons, 
          lat_gridvals = lats, 
          log_a0val = log_a0val, 
          dims = dims, 
          surface_heights = zsurf, domain_heights = -template_bins$center, 
          log_self_tx = log(theta[1]), betaAR = theta[2]
        ) 
      }
      
      res = data.frame(
        speed = whale_ll_grid['speed'],
        betaAR = whale_ll_grid['betaAR'],
        ll = ll(speed = unlist(whale_ll_grid['speed']),
                betaAR = unlist(whale_ll_grid['betaAR'])),
        delta = delta
      )
      
      res
    },
    pattern = map(whale_ll_grid),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  )
  
)
