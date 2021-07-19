whale_targets = list(
  
  #
  # gridded likelihood evaluation
  #
  
  tar_target(whale_ll_delta, 300),
  tar_target(whale_ll_delta_srtm30, 60),
  
  tar_target(
    name = cee_start, 
    command = strptime(
      x = '8/19/2019  7:11:00 PM UTC', 
      format = '%m/%d/%Y %I:%M:%S %p',
      tz = 'UTC'
    )
  ),
  
  tar_target(
    name = whale_ll_grid, 
    command = expand.grid(
      # m/s
      speed = seq(from = .01, to = 4, length.out = 41),
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
          cellsize = c(
            diff(unique(bathy[,'lon']))[1],
            diff(unique(bathy[,'lat']))[1]
          ),
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
      
      pl = ggplot(whale_ll, aes(x = speed, y = betaAR, fill = ll, z = ll)) + 
        geom_raster() + 
        scale_fill_viridis(direction = -1) + 
        theme_few() + 
        stat_contour2(col = 'grey90') + 
        geom_text_contour(col = 'grey90')
      
      ggsave(pl, filename = 'whale_ll_plot.png')
      
      whale_ll$lp = whale_ll$ll - log_sum(whale_ll$ll)
      
      sum(whale_ll$speed * exp(whale_ll$lp))
      sum(whale_ll$betaAR * exp(whale_ll$lp))
      
      marginal_post_speed = whale_ll %>% 
        group_by(speed) %>%
        summarise(lp = log_sum(lp))
      
      plot(exp(lp) ~ speed, marginal_post_speed, type = 'l')
      
      pl = ggplot(marginal_post_speed, aes(x = speed, y = exp(lp))) + 
        geom_line() +
        theme_few()
      
      ggsave(pl, filename = 'marginal_post_speed.png')
      
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
      
      # append with empirical speed estimates (m/s)
      x$empirical_speed = c(0, sapply(2:nrow(x), function(ind) {
        rdist.earth(x1 = as.matrix(x[ind-1, c('Longitude', 'Latitude')]),
                    x2 = as.matrix(x[ind, c('Longitude', 'Latitude')]), 
                    miles = FALSE) * 1e3 / 
          diff(as.numeric(x$date_time[ind + (-1:0)]))
      }))
      
      x$time_diff = round(c(0, diff(as.numeric(x$date_time)))/60)
      
      #
      # data filtering
      #
      
      # deconflict entries with similar timestamps
      x_filtered = do.call(rbind, lapply(1:nrow(x), function(ind) {
        # identify all observations within 5min of current observation
        obs_inds = which(
          abs(as.numeric(x$date_time[ind]) - as.numeric(x$date_time)) <= 300
        )
        # export the observation with the smallest observation error
        export_ind = obs_inds[which.min(
          x$Error.Semi.major.axis[obs_inds]^2 + 
            x$Error.Semi.minor.axis[obs_inds]^2
        )]
        x[export_ind,]
      }))
      x_filtered$time_diff = NULL
      x_filtered$empirical_speed = NULL
      x_filtered = unique(x_filtered)
      
      # remove data points with large observation error
      x_filtered = x_filtered %>% 
        filter(
          Error.Semi.major.axis < 50000/2,
          Error.Semi.minor.axis < 50000/2
        )
      
      # recompute empirical speeds
      x_filtered$time_diff = c(0, diff(as.numeric(x_filtered$date_time)))
      x_filtered$time_steps = round(x_filtered$time_diff/60)
      x_filtered$empirical_speed = c(0, sapply(2:nrow(x_filtered), function(ind) {
        rdist.earth(x1 = as.matrix(x_filtered[ind-1, c('Longitude', 'Latitude')]),
                    x2 = as.matrix(x_filtered[ind, c('Longitude', 'Latitude')]), 
                    miles = FALSE) * 1e3 / 
          diff(as.numeric(x_filtered$date_time[ind + (-1:0)]))
      }))
      
      # remove observations contributing to large speeds
      x_filtered = x_filtered %>% filter(empirical_speed < 10)
      
      # recompute empirical speeds after filtering
      x_filtered$time_diff = c(0, diff(as.numeric(x_filtered$date_time)))
      x_filtered$time_steps = round(x_filtered$time_diff/60)
      x_filtered$empirical_speed = c(0, sapply(2:nrow(x_filtered), function(ind) {
        rdist.earth(x1 = as.matrix(x_filtered[ind-1, c('Longitude', 'Latitude')]),
                    x2 = as.matrix(x_filtered[ind, c('Longitude', 'Latitude')]), 
                    miles = FALSE) * 1e3 / 
          diff(as.numeric(x_filtered$date_time[ind + (-1:0)]))
      }))
      
      # replace raw dataset with filtered dataset
      x = x_filtered
      
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
      
      coord_diffs = c(0, sapply(2:nrow(mapped_coords), function(ind) {
        sum(abs(mapped_coords[ind - 1,] - mapped_coords[ind,]))
      }))
      
      delta_diffs = round(c(0, diff(as.numeric(x$date_time))/delta))
      
      # View(data.frame(
      #   mapped_coords,
      #   coord_diffs,
      #   delta_diffs,
      #   infeasible = delta_diffs < coord_diffs,
      #   slack = abs(delta_diffs - coord_diffs)
      # ))
      
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
  ),
  
  tar_target(
    name = whale_ll_srtm30_plot,
    command = {
      
      # set spatial extents to extract
      extent_vec = c(min_lon = -75.6, max_lon = -74, min_lat = 34.5, max_lat = 36)
      
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
      mapped_coords = map_coords(
        lon = x$Longitude, 
        lat = x$Latitude, 
        lon_grid = lons, 
        lat_grid = lats, 
        coord_grid = coords
      )
      
      # load depth data
      template_bins = read.csv('data/BeakedWhale/template_bins.csv')
      depths_raw = read.csv('data/BeakedWhale/ZcTag095_DUML_series_20200703.csv')
      depths_raw$Date = as.POSIXct(x = depths_raw$Date, tz = 'UTC', 
                                   origin = '1970-01-01 00:00.00 UTC')
      
      # plotting projection
      plot_crs = '+proj=aea +lat_1=34.5 +lat_2=36 +lon_0=-75'
      
      obs_ind = 1
      obs_ind = 225
      
      # append dummy depth variable and dimension information
      states = cbind(mapped_coords - 1, 0)
      dims = c(length(lats), length(lons), 1)
      # height of domain surface
      zsurf = matrix(
        getValues(t(bathy)), 
        nrow = nrow(bathy), 
        ncol = ncol(bathy)
      )
      # dummy height of vertical layer
      zval = 0
      
      # uniform mass over neighborhood structure for initial state
      a0_prev_coords = dsmovetools:::TestZConstrainedRookNeighborhood(
        dims = dims, 
        x = as.matrix(states[obs_ind,])[,c(2,1,3)], 
        zfield = zsurf, 
        zvals = zval
      )
      log_a0val = rep(1, nrow(a0_prev_coords))
      log_a0val = log(log_a0val / sum(log_a0val))
      
      expanded_nbhd = dsmovetools:::SattagExpandNeighborhood(
        a0 = matrix(unlist(states[obs_ind,c(2,1,3)]), 
                    nrow = nrow(a0_prev_coords), 
                    ncol = ncol(states),
                    byrow = TRUE),
        a0_prev_coords = a0_prev_coords, 
        log_a0val = log_a0val, 
        gps_trunc_alpha = .05, 
        obs_lons = lons[mapped_coords[obs_ind,'lon_ind']],
        obs_lats = lats[mapped_coords[obs_ind,'lat_ind']],
        obs_semi_majors = x$Error.Semi.major.axis[obs_ind],
        obs_semi_minors = x$Error.Semi.minor.axis[obs_ind], 
        obs_orientations = x$Error.Ellipse.orientation[obs_ind], 
        obs_depth_bins = NA, 
        lon_gridvals = lons, 
        lat_gridvals = lats, 
        dims = dims, 
        surface_heights = zsurf, domain_heights = -template_bins$center, 
        obs_ind = 0, 
        nsteps = 50
      )
      
      expanded_nbhd = unique(expanded_nbhd[,1:2])
      expanded_nbhd_coords = cbind(
        lon = lons[expanded_nbhd[,2] + 1],
        lat = lats[expanded_nbhd[,1] + 1]
      )
      
      expanded_nbhd_coords.unprojected = SpatialPoints(
        coords = expanded_nbhd_coords,
        proj4string = crs(bathy)
      )
      
      expanded_nbhd_coords.projected = spTransform(
        expanded_nbhd_coords.unprojected, CRS(plot_crs)
      )
      
      animalobs.unprojected = SpatialPoints(
        coords = x[obs_ind, c('Longitude', 'Latitude')], 
        proj4string = crs(bathy)
      )
      
      obserrors.projected = do.call(rbind,apply(
        x[obs_ind, c('Longitude', 'Latitude', 'Error.Semi.major.axis', 
              'Error.Semi.minor.axis', 'Error.Ellipse.orientation')], 1,
        function(r) {
          # project lon/lat GPS reading to a local coordinate system
          # so that the coordinates are compatible with the error scale (m)
          ctr = SpatialPoints(coords = matrix(r[c('Longitude', 'Latitude')], 
                                              nrow = 1), 
                              proj4string = crs(bathy))
          ctr.projected = spTransform(ctr, CRS(plot_crs))
          # points for error ellipse
          res = data.frame(
            ellipsePoints(center = as.numeric(coordinates(ctr.projected)), 
                          major = as.numeric(r['Error.Semi.major.axis']), 
                          minor = as.numeric(r['Error.Semi.minor.axis']), segments = 50,
                          angle = (90 - as.numeric(r['Error.Ellipse.orientation'])) / 180 * pi)
          )
          res$id = 1
          res
        }))
      
      # projected GPS coordinates, for plotting
      animalobs.projected = spTransform(animalobs.unprojected, CRS(plot_crs))
      
      # subsample data, for faster plotting
      bathy.downsampled = sampleRegular(x = bathy, size = 1e6, asRaster = TRUE)
      
      # projected bathymetry coordinates, for plotting
      bathy.projected = projectRaster(from = bathy.downsampled, crs = plot_crs)
      
      # separate land from bathymetry
      land.mask = bathy.projected
      land.mask[land.mask <= 0] = NA
      land.mask[land.mask > 0] = 1
      bathy.projected[bathy.projected > 0] = NA
      bathy.projected = bathy.projected * -1
      
      # munge raster to tidy format
      df = data.frame(coordinates(bathy.projected), 
                      depth = raster::values(bathy.projected))
      df.land = data.frame(coordinates(land.mask), 
                           mask = raster::values(land.mask))
      
      # plot gps track against bathymetry and coastline
      ggplot(mapping = aes(x = x, y = y)) + 
        # bathymetry
        geom_tile(mapping = aes(fill = depth), data = df, col = 1) +
        # geom_point(data = df, pch = 20, size = .1) +
        scale_fill_distiller('Depth (m)', palette = 'Blues', direction = -1,
                             na.value = 'transparent', trans = 'reverse') +
        # land mask
        new_scale_fill() + 
        guides(fill = 'none') + 
        geom_raster(mapping = aes(fill = mask), data = df.land) + 
        scale_fill_viridis_c(na.value = 'transparent') +
        # observation error ellipses
        geom_polygon(mapping = aes(x = X1, y = X2, group = id),
                     data = obserrors.projected, inherit.aes = FALSE,
                     fill = NA, col = 'black') +
        # expanded neighborhood points
        geom_point(
          data = data.frame(coordinates(expanded_nbhd_coords.projected)), 
          mapping = aes(x = lon, y = lat), inherit.aes = FALSE, col = 2,
          size = .5
        ) + 
        # observations
        # geom_path(mapping = aes(x = Longitude, y = Latitude), 
        #           inherit.aes = FALSE, alpha = .3, lwd = .5,
        #           data = data.frame(animalobs.projected)) +
        geom_point(mapping = aes(x = Longitude, y = Latitude), 
                   inherit.aes = FALSE,
                   data = data.frame(animalobs.projected), size = .5) +
        # map extent labels
        scale_x_continuous(
          breaks = range(coordinates(obserrors.projected)[,1]) + 50e3 * c(-1,1),
          limits = range(coordinates(obserrors.projected)[,1]) + 50e3 * c(-1,1),
          oob = scales::oob_keep,
          labels = telefit:::lon_trans()$format(
            extent_vec[c('min_lon', 'max_lon')]
          )) +
        scale_y_continuous(
          breaks = range(coordinates(obserrors.projected)[,2]) + 1e3 * c(-1,1),
          limits = range(coordinates(obserrors.projected)[,2]) + 1e3 * c(-1,1),
          oob = scales::oob_keep,
          labels = telefit:::lat_trans()$format(
            extent_vec[c('min_lat', 'max_lat')]
          )) +
        # formatting
        coord_equal() + 
        theme_few() + 
        theme(panel.border = element_blank(), 
              axis.title = element_blank(), 
              axis.ticks = element_blank())
      
      ggplot(
        x %>% pivot_longer(cols = c('Latitude', 'Longitude'), 
                           names_to = 'Coordinate', values_to = 'Value'),
        aes(x = date_time, y = Value)
      ) + 
        geom_line() + 
        geom_point() + 
        geom_vline(xintercept = cee_start, lty = 3) + 
        facet_wrap(~Coordinate, nrow = 2, scales = 'free_y') + 
        theme_few()
      
      
      
    },
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(
    name = whale_ll_srtm30_filtered,
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
      
      # append with empirical speed estimates (m/s)
      x$empirical_speed = c(0, sapply(2:nrow(x), function(ind) {
        rdist.earth(x1 = as.matrix(x[ind-1, c('Longitude', 'Latitude')]),
                    x2 = as.matrix(x[ind, c('Longitude', 'Latitude')]), 
                    miles = FALSE) * 1e3 / 
          diff(as.numeric(x$date_time[ind + (-1:0)]))
      }))
      
      x$time_diff = round(c(0, diff(as.numeric(x$date_time)))/60)
      
      
      #
      # data filtering
      #
      
      # deconflict entries with similar timestamps
      x_filtered = do.call(rbind, lapply(1:nrow(x), function(ind) {
        # identify all observations within 5min of current observation
        obs_inds = which(
          abs(as.numeric(x$date_time[ind]) - as.numeric(x$date_time)) <= 300
        )
        # export the observation with the smallest observation error
        export_ind = obs_inds[which.min(
          x$Error.Semi.major.axis[obs_inds]^2 + 
            x$Error.Semi.minor.axis[obs_inds]^2
        )]
        x[export_ind,]
      }))
      x_filtered$time_diff = NULL
      x_filtered$empirical_speed = NULL
      x_filtered = unique(x_filtered)
      
      # remove data points with large observation error
      x_filtered = x_filtered %>% 
        filter(
          Error.Semi.major.axis < 50000/2,
          Error.Semi.minor.axis < 50000/2
        )
      
      # recompute empirical speeds
      x_filtered$time_diff = c(0, diff(as.numeric(x_filtered$date_time)))
      x_filtered$time_steps = round(x_filtered$time_diff/60)
      x_filtered$empirical_speed = c(0, sapply(2:nrow(x_filtered), function(ind) {
        rdist.earth(x1 = as.matrix(x_filtered[ind-1, c('Longitude', 'Latitude')]),
                    x2 = as.matrix(x_filtered[ind, c('Longitude', 'Latitude')]), 
                    miles = FALSE) * 1e3 / 
          diff(as.numeric(x_filtered$date_time[ind + (-1:0)]))
      }))
      
      # remove observations contributing to large speeds
      x_filtered = x_filtered %>% filter(empirical_speed < 10)
      
      # recompute empirical speeds after filtering
      x_filtered$time_diff = c(0, diff(as.numeric(x_filtered$date_time)))
      x_filtered$time_steps = round(x_filtered$time_diff/60)
      x_filtered$empirical_speed = c(0, sapply(2:nrow(x_filtered), function(ind) {
        rdist.earth(x1 = as.matrix(x_filtered[ind-1, c('Longitude', 'Latitude')]),
                    x2 = as.matrix(x_filtered[ind, c('Longitude', 'Latitude')]), 
                    miles = FALSE) * 1e3 / 
          diff(as.numeric(x_filtered$date_time[ind + (-1:0)]))
      }))
      
      # replace raw dataset with filtered dataset
      x = x_filtered
      
      # study segments of data with bounded time between observations
      gaps = rle(c(0, diff(as.numeric(x_filtered$date_time))) / 3600 <= 5)
      # total number of observations across data sections
      sum(gaps$lengths[gaps$values == TRUE])
      sum(gaps$lengths[gaps$values == TRUE][gaps$lengths[gaps$values == TRUE] > 1])
      # lengths of contiguous data sections
      gaps$lengths[gaps$values == TRUE]
      
      x_filtered$gap = x_filtered$time_diff / 3600 > 5
      
      table(date(x_filtered$date_time), x_filtered$Source, x_filtered$gap)
      
      # restrict dataset to a more targeted time range
      x_filtered = x_filtered %>% filter(
        date(date_time) %in% (date(cee_start) + duration(0:1, units= 'days'))
      )
      
      x = x_filtered
      
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
      states = as.matrix(cbind(mapped_coords[,c('lat_ind','lon_ind')] - 1, 0))
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
      
      coord_diffs = c(0, sapply(2:nrow(mapped_coords), function(ind) {
        sum(abs(mapped_coords[ind - 1,] - mapped_coords[ind,]))
      }))
      
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
        obs_gps$lon[dat_ind] = lons[mapped_coords[ind, 'lon_ind'] - 1]
        obs_gps$lat[dat_ind] = lats[mapped_coords[ind, 'lat_ind'] - 1]
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
      
      # log-likelihood for observations
      ll = function(speed, betaAR, start_ind = 1, end_ind = nrow(obs_gps),
                    start_state = states[1,]) {
        
        # uniform mass over neighborhood structure for initial state
        a0_prev_coords = dsmovetools:::TestZConstrainedRookNeighborhood(
          dims = dims, 
          x = start_state, 
          zfield = zsurf, 
          zvals = zval
        )
        log_a0val = rep(1, nrow(a0_prev_coords))
        log_a0val = log(log_a0val / sum(log_a0val))
        
        # expand neighborhood!
        
        expanded_nbhd = dsmovetools:::SattagExpandNeighborhood(
          a0 = matrix(start_state, 
                      nrow = nrow(a0_prev_coords), 
                      ncol = ncol(states), byrow = TRUE),
          a0_prev_coords = a0_prev_coords, 
          gps_trunc_alpha = .05, 
          obs_lons = obs_gps$lon[start_ind], 
          obs_lats = obs_gps$lat[start_ind],
          obs_semi_majors = obs_gps$semi_major[start_ind],
          obs_semi_minors = obs_gps$semi_minor[start_ind], 
          obs_orientations = obs_gps$orientation[start_ind], 
          obs_depth_bins = obs_depths[start_ind], 
          lon_gridvals = lons, 
          lat_gridvals = lats, 
          log_a0val = log_a0val, 
          dims = dims, 
          surface_heights = zsurf, domain_heights = -template_bins$center, 
          obs_ind = 0, 
          nsteps = 50
        ) 
        
        theta = c(1 - speed / 1187.295 * delta, betaAR)
        
        dsmovetools:::SattagFilteredLL(
          a0 = expanded_nbhd[,1:3],
          a0_prev_coords = expanded_nbhd[,4:6], 
          gps_trunc_alpha = .05, 
          obs_lons = obs_gps$lon[start_ind:end_ind], 
          obs_lats = obs_gps$lat[start_ind:end_ind],
          obs_semi_majors = obs_gps$semi_major[start_ind:end_ind],
          obs_semi_minors = obs_gps$semi_minor[start_ind:end_ind], 
          obs_orientations = obs_gps$orientation[start_ind:end_ind], 
          obs_depth_bins = obs_depths[start_ind:end_ind], 
          lon_gridvals = lons, 
          lat_gridvals = lats, 
          log_a0val = rep(-log(nrow(expanded_nbhd)), nrow(expanded_nbhd)), 
          dims = dims, 
          surface_heights = zsurf, domain_heights = -template_bins$center, 
          log_self_tx = log(theta[1]), betaAR = theta[2]
        ) 
      }
      
      # flag indices with observations; annotate with gap between observations
      obs_runs = rle(is.na(obs_gps$lon))
      obs_summary = cbind(
        obs_ind = cumsum(obs_runs$lengths)[obs_runs$values == FALSE],
        dist_to_next_obs = c(obs_runs$lengths[obs_runs$values == TRUE], 0) + 1
      )
      
      # identify segments of observations that are relatively close together
      run_chunks = rle(obs_summary[,'dist_to_next_obs'] <= 150)
      run_chunks_summary = cbind(
        obs_summary_ind = cumsum(run_chunks$lengths)[run_chunks$values == FALSE]
      )
      
      # start/stop indices of likelihood data for each contiguous segment
      start_inds = obs_summary[
        cumsum(c(1, run_chunks$lengths))[c(run_chunks$values, FALSE) == TRUE],
        'obs_ind'
      ]
      end_inds = rbind(obs_summary, c(nrow(obs_gps), 0))[
        cumsum(c(1, run_chunks$lengths))[c(run_chunks$values, FALSE) == FALSE],
        'obs_ind'
      ]
      
      # grid-indices from which diffusions should be initialized
      state_inds = (1:nrow(obs_summary))[
        cumsum(c(1, run_chunks$lengths))[c(run_chunks$values, FALSE) == TRUE]
      ]
      
      res = data.frame(
        speed = whale_ll_grid['speed'],
        betaAR = whale_ll_grid['betaAR'],
        ll = sum(sapply(1:length(start_inds), function(seg_ind) {
          ll(speed = unlist(whale_ll_grid['speed']),
             betaAR = unlist(whale_ll_grid['betaAR']), 
             start_ind = start_inds[seg_ind], 
             end_ind = end_inds[seg_ind],
             start_state = states[state_inds[seg_ind],])
        })),
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
    name = whale_ll_srtm30_filtered_alt,
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
        ncol = length(lons)
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
      
      # log-likelihood for observations
      ll = function(speed, betaAR, start_state = states[1,]) {
        
        theta = c(1 - speed / 1187.295 * delta, betaAR)
        
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
        delta = delta
      )
      tock = proc.time()[3]
      
      res$computation_time = tock-tick
      
      res
    },
    pattern = map(whale_ll_grid),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  )
  
)
