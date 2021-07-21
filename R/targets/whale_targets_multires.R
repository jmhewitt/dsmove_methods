whale_targets_multires = list(
  
  
  tar_target(whale_res_factor, 2:4),
  
  
  tar_target(
    name = whale_ll_srtm30_filtered_scaled,
    command = {
      
      # set spatial extents to extract
      extent_vec = c(min_lon = -77, max_lon = -72, min_lat = 33, max_lat = 37)
      
      domain.extent = extent(extent_vec)
      
      # load, crop, and merge bathymetry files
      bathy = crop(raster('data/STRM_30_PLUS/w100n40.nc'), domain.extent)
      
      # set CRS
      crs(bathy) = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
      
      bathy = aggregate(x = bathy, fact = whale_res_factor, fun = min)
      
      # scale bathymetry
      
      # zc095 gps observations
      x = readRDS('data/BeakedWhale/zctag095_gps.rds')
      
      # map the observed lat-lons to the bathymetry lon/lat grid
      # (needed to set initial filtering locations)
      coords = coordinates(bathy)
      lons = unique(coords[,1])
      lats = unique(coords[,2])
      x_mapped_coords = map_coords(
        lon = x$Longitude, 
        lat = x$Latitude, 
        lon_grid = lons, 
        lat_grid = lats, 
        coord_grid = coords
      )
      
      x$mapped_lon = lons[x_mapped_coords$lon_ind]
      x$mapped_lat = lats[x_mapped_coords$lat_ind]
      
      # filter out observations that are missing (Error) information
      x = x[complete.cases(x),]
      
      # x$mapped_filter = sdafilter(
      #   lat = x$mapped_lat, lon = x$mapped_lon, dtime = x$date_time, 
      #   lc = x$Quality, vmax = 5
      # )
      # 
      # x = x[x$mapped_filter != 'removed',]
      
      # append with empirical speed estimates (m/s)
      x$empirical_speed = c(0, sapply(2:nrow(x), function(ind) {
        rdist.earth(x1 = as.matrix(x[ind-1, c('Longitude', 'Latitude')]),
                    x2 = as.matrix(x[ind, c('Longitude', 'Latitude')]), 
                    miles = FALSE) * 1e3 / 
          diff(as.numeric(x$date_time[ind + (-1:0)]))
      }))
      
      # empirical speed estimates wrt mapped coordinates
      x$mapped_speed = c(0, sapply(2:nrow(x), function(ind) {
        rdist.earth(x1 = as.matrix(x[ind-1, c('mapped_lon', 'mapped_lat')]),
                    x2 = as.matrix(x[ind, c('mapped_lon', 'mapped_lat')]), 
                    miles = FALSE) * 1e3 / 
          diff(as.numeric(x$date_time[ind + (-1:0)]))
      }))
      
      # (rounded) minutes between observations
      x$time_diff = round(c(0, diff(as.numeric(x$date_time)))/60)
      
      # filter out observations associated with large empirical speeds
      x = x %>% filter(mapped_speed <= 6)

      # recompute empirical speed estimates (m/s)
      x$empirical_speed = c(0, sapply(2:nrow(x), function(ind) {
        rdist.earth(x1 = as.matrix(x[ind-1, c('Longitude', 'Latitude')]),
                    x2 = as.matrix(x[ind, c('Longitude', 'Latitude')]),
                    miles = FALSE) * 1e3 /
          diff(as.numeric(x$date_time[ind + (-1:0)]))
      }))

      x$mapped_speed = c(0, sapply(2:nrow(x), function(ind) {
        rdist.earth(x1 = as.matrix(x[ind-1, c('mapped_lon', 'mapped_lat')]),
                    x2 = as.matrix(x[ind, c('mapped_lon', 'mapped_lat')]),
                    miles = FALSE) * 1e3 /
          diff(as.numeric(x$date_time[ind + (-1:0)]))
      }))
      
      # filter out observations associated with large empirical speeds
      x = x %>% filter(mapped_speed <= 6)
      
      # recompute empirical speed estimates (m/s)
      x$empirical_speed = c(0, sapply(2:nrow(x), function(ind) {
        rdist.earth(x1 = as.matrix(x[ind-1, c('Longitude', 'Latitude')]),
                    x2 = as.matrix(x[ind, c('Longitude', 'Latitude')]),
                    miles = FALSE) * 1e3 /
          diff(as.numeric(x$date_time[ind + (-1:0)]))
      }))
      
      x$mapped_speed = c(0, sapply(2:nrow(x), function(ind) {
        rdist.earth(x1 = as.matrix(x[ind-1, c('mapped_lon', 'mapped_lat')]),
                    x2 = as.matrix(x[ind, c('mapped_lon', 'mapped_lat')]),
                    miles = FALSE) * 1e3 /
          diff(as.numeric(x$date_time[ind + (-1:0)]))
      }))

      # (rounded) minutes between observations
      x$time_diff = round(c(0, diff(as.numeric(x$date_time)))/60)
      
      x = x %>% filter(
        date(date_time) %in% (date(cee_start) + duration(0:1, units= 'days'))
      )

      
      #
      # map data to grid
      #
      
      # map the observed lat-lons to the bathymetry lon/lat grid
      # (needed to set initial filtering locations)
      coords = coordinates(bathy)
      lons = unique(coords[,1])
      lats = unique(coords[,2])
      mapped_coords = map_coords(
        lon = x$Longitude, 
        lat = x$Latitude, 
        lon_grid = lons, 
        lat_grid = lats, 
        coord_grid = coords
      )
      
      x$coord_diffs = c(0, sapply(2:nrow(mapped_coords), function(ind) {
        sum(abs(mapped_coords[ind - 1,] - mapped_coords[ind,]))
      }))
      
      # load depth data
      template_bins = read.csv('data/BeakedWhale/template_bins.csv')
      depths_raw = read.csv('data/BeakedWhale/ZcTag095_DUML_series_20200703.csv')
      depths_raw$Date = as.POSIXct(x = depths_raw$Date, tz = 'UTC', 
                                   origin = '1970-01-01 00:00.00 UTC')
      
      # filter depth data to lie within the GPS date range
      depths_raw = depths_raw %>% filter(
        min(x$date_time) <= Date, Date <= max(x$date_time)
      )
      
      #
      # likelihood approximation
      #
      
      # append dummy depth variable and dimension information
      states = cbind(
        as.matrix(mapped_coords[,c('lon_ind','lat_ind')] - 1),
        sqerror = x$Error.Semi.major.axis^2 + x$Error.Semi.minor.axis^2,
        x_ind = 1:nrow(mapped_coords)
      )
      # height of domain surface
      zsurf = matrix(
        getValues(bathy), 
        nrow = length(lats), 
        ncol = length(lons),
        byrow = TRUE
      )
      
      # discretized timesteps
      delta = 30 #whale_ll_delta_srtm30
      tseq = seq(
        from = min(x$date_time),
        to = max(x$date_time),
        by = duration(delta, 'seconds')
      )
      
      # map observations to discretized timesteps
      obs_coords = matrix(NA, nrow = length(tseq), ncol = ncol(states))
      colnames(obs_coords) = colnames(states)
      for(ind in 1:nrow(x)) {
        tgt_ind = which.min(abs(x$date_time[ind] - tseq))
        if(any(is.finite(obs_coords[tgt_ind,]))) {
          if(obs_coords[tgt_ind,'sqerror'] > states[ind,'sqerror']) {
            # transfer observation with less observation error
            obs_coords[tgt_ind,] = states[ind,]
          }
          
        } else {
          # transfer observation
          obs_coords[tgt_ind,] = states[ind,]
        }
      }
      
      # obs_steps = complete.cases(obs_coords)
      # coord_diffs = sapply(2:nrow(obs_coords[obs_steps,]), function(ind) {
      #   sum(abs(
      #     obs_coords[obs_steps,][ind,c('lon_ind','lat_ind')] - 
      #     obs_coords[obs_steps,][ind - 1,c('lon_ind','lat_ind')]
      #   ))
      # })
      # 
      # tsteps = diff(which(obs_steps))
      # which(!(coord_diffs <= tsteps))
      
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
        obs_gps$lon[dat_ind] = lons[mapped_coords[ind, 'lon_ind']]
        obs_gps$lat[dat_ind] = lats[mapped_coords[ind, 'lat_ind']]
        obs_gps$semi_major[dat_ind] = x$Error.Semi.major.axis[ind]
        obs_gps$semi_minor[dat_ind] = x$Error.Semi.minor.axis[ind]
        obs_gps$orientation[dat_ind] = x$Error.Ellipse.orientation[ind]
      }
      
      # map depth observations to discretized timesteps
      obs_depths = rep(NA, length(tseq))
      for(ind in 1:nrow(depths_raw)) {
        obs_depths[which.min(abs(depths_raw$Date[ind] - tseq))] = 
          obs_depths = -depths_raw$Depth[ind]
      }
      
      grid_res = rdist.earth(x1 = coords[1, , drop = FALSE], 
                             x2 = coords[2, , drop = FALSE], 
                             miles = FALSE) * 1e3
      
      # log-likelihood for observations
      ll = function(speed, betaAR, start_state = states[1,]) {
        
        theta = c(1 - speed / grid_res * delta, betaAR)
        
        browser()
        dsmovetools2d:::SattagFilteredLL(
          init_dsts = matrix(
            c(states[1,'lon_ind'], states[1,'lat_ind']), 
            nrow = 4, ncol = 2,
            byrow = TRUE
          ),
          init_srcs = matrix(
            c(states[1,'lon_ind'], states[1,'lat_ind'] - 1,
              states[1,'lon_ind'], states[1,'lat_ind'] + 1,
              states[1,'lon_ind'] - 1, states[1,'lat_ind'],
              states[1,'lon_ind'] + 1, states[1,'lat_ind']
              ), 
            nrow = 4, ncol = 2,
            byrow = TRUE
          ), 
          init_log_probs = rep(log(1/4),4), 
          gps_trunc_alpha = .05, 
          obs_lons = obs_gps$lon, 
          obs_lats = obs_gps$lat, 
          obs_semi_majors = obs_gps$semi_major, 
          obs_semi_minors = obs_gps$semi_minor, 
          obs_orientations = obs_gps$orientation, 
          obs_depths = obs_depths, 
          lon_gridvals = lons, 
          lat_gridvals = lats, 
          surface_heights = zsurf, 
          min_elevation = -10000, 
          max_elevation = 0, 
          log_self_tx = log(theta[1]), 
          betaAR = theta[2]
        ) 
      }
      
      tick = proc.time()[3]
      res = data.frame(
        speed = whale_ll_grid['speed'],
        betaAR = whale_ll_grid['betaAR'],
        ll = ll(speed = unlist(whale_ll_grid['speed']),
                betaAR = unlist(whale_ll_grid['betaAR']),
                start_state = states[1,]),
        delta = delta,
        grid_res = grid_res
      )
      tock = proc.time()[3]
      
      res$computation_time = tock-tick
      res$computation_date = Sys.time()
      res$computation_node = Sys.info()['nodename']
      
      print(res)
      
      res
    },
    pattern = cross(map(whale_ll_grid), whale_res_factor),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  )
  
)
