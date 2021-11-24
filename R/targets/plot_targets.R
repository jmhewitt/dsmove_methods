plot_targets = list(
  
  # 3-panel plot comparing posterior locations between different approximations
  tar_target(
    name = zc095_post_loc_comparison, 
    command = {
      
      #
      # load data
      #
      
      # load data
      pkg = whale_pkg[[1]]
      # load posterior locations from CTDS model
      post_marginal = readRDS('marginal_locations.rds')
      # load posterior locations from crawl
      post_crawl = readRDS('crawl_discretized_post_locations100.rds')
      
      #
      # subset data
      #
      
      # extract prediction times
      pred_times = as.POSIXct(
        x = sapply(post_marginal, function(r) r$meta$time), 
        tz = 'UTC', origin = '1970-01-01 00:00.00 UTC'
      )
      # target prediction time
      offset = 2
      ptime_without_depth = 18 + offset
      ptime_with_depth = 6 + offset
      # CTDS-based posterior location
      x_with_depth = post_marginal[[ptime_with_depth]]$posterior
      x_without_depth = post_marginal[[ptime_without_depth]]$posterior
      
      #
      # prepare data
      #
      
      # normalize posterior
      x_with_depth$lp = x_with_depth$lp - dsmovetools2d:::log_sum_c(x_with_depth$lp)
      x_without_depth$lp = x_without_depth$lp - 
        dsmovetools2d:::log_sum_c(x_without_depth$lp)
      # crawl-based posterior location
      x_crawl = do.call(rbind, lapply(post_crawl, function(imputation) {
        ind = which.min(abs(imputation$time - pred_times[ptime_without_depth]))
        imputation[ind, c('mapped_lon', 'mapped_lat')]
      }))
      # condense to empirical probability distribution
      x_crawl = x_crawl %>% 
        mutate(lp = log(1 /n())) %>%
        group_by(mapped_lon, mapped_lat) %>%
        summarise(lp = dsmovetools2d:::log_sum_c(lp)) %>% 
        ungroup() %>% 
        select(lon = mapped_lon, lat = mapped_lat, lp)
      # map the observed lat-lons to the bathymetry lon/lat grid
      coords = pkg$grid$coords
      lons = pkg$grid$lons
      lats = pkg$grid$lats
      # extract coordinates from CTDS via time discretization
      x_with_depth$lon = lons[x_with_depth$lon_to_ind]
      x_with_depth$lat = lats[x_with_depth$lat_to_ind]
      x_without_depth$lon = lons[x_without_depth$lon_to_ind]
      x_without_depth$lat = lats[x_without_depth$lat_to_ind]
      
      #
      # project data
      #
      
      # plotting projection
      plot_crs = '+proj=aea +lat_1=34.5 +lat_2=36 +lon_0=-75'
      
      # projected bathymetry coordinates, for plotting
      bathy.projected = projectRaster(from = whale_domain, crs = plot_crs)
      
      # project posteriors for marginal location
      loc_with_depth = project_to_grid(
        data = x_with_depth, 
        data_crs = crs(whale_domain), 
        rasterObj = bathy.projected, 
        plot_crs = plot_crs
      )
      loc_without_depth = project_to_grid(
        data = x_without_depth,  
        data_crs = crs(whale_domain), 
        rasterObj = bathy.projected, 
        plot_crs = plot_crs
      )
      loc_crawl = project_to_grid(
        data = x_crawl, 
        data_crs = crs(whale_domain), 
        rasterObj = bathy.projected, 
        plot_crs = plot_crs
      )
      
      #
      # build plot
      #
      
      # separate land from bathymetry
      land.mask = bathy.projected
      land.mask[land.mask <= 0] = NA
      land.mask[land.mask > 0] = 1
      
      # re-label depth data
      bathy.projected[bathy.projected > 0] = NA
      bathy.projected = bathy.projected * -1
      
      # munge land/depth rasters to tidy format
      bathy_df = data.frame(coordinates(bathy.projected), 
                            depth = raster::values(bathy.projected))
      land_df = data.frame(coordinates(land.mask), 
                           mask = raster::values(land.mask))
      
      # trim plotting data
      bathy_df$depth[bathy_df$depth > 2000] = NA
      bathy_contour = bathy_df
      
      # combine data for plotting
      df = rbind(
        loc_with_depth %>% 
          select(lon.projected, lat.projected, lp) %>% 
          filter(lp >= log(.005)) %>% 
          mutate(Model = 'State space (w/Depth)'), 
        loc_without_depth %>% 
          select(lon.projected, lat.projected, lp) %>% 
          filter(lp >= log(.005)) %>%
          mutate(Model = 'State space'),
        loc_crawl %>% 
          select(lon.projected, lat.projected, lp) %>% 
          mutate(Model = 'CTCRW-AID')
      ) %>% 
        mutate(Model = factor(Model, levels = c('CTCRW-AID', 
                                                'State space',
                                                'State space (w/Depth)')))
      
      pl = ggplot(mapping = aes(x = x, y = y)) + 
        # bathymetry
        geom_raster(mapping = aes(fill = depth), data = bathy_df) + 
        scale_fill_distiller('Depth (m)', palette = 'Blues', direction = -1,
                             na.value = 'transparent', trans = 'reverse') + 
        # land mask
        new_scale_fill() + 
        guides(fill = 'none') + 
        geom_raster(mapping = aes(fill = mask), data = land_df) + 
        scale_fill_viridis_c(na.value = 'transparent') + 
        # marginal posterior densities
        geom_point(mapping = aes(x = lon.projected, 
                                 y = lat.projected, 
                                 col = exp(lp)), 
                   data = df, 
                   inherit.aes = FALSE,
                   size = 2) + 
        scale_color_viridis_c('Prob.', direction = -1, trans = 'log', 
                              breaks = c(.5, .25, .125, .05, .01)) + 
        # bathymetry contours
        geom_contour2(mapping = aes(z = depth), data = bathy_contour, 
                      col = 'grey50', breaks = seq(500,2000,by=100)) + 
        geom_text_contour(mapping = aes(z = depth), data = bathy_contour, 
                          breaks = seq(500,2000,by=100)) + 
        # formatting
        facet_wrap(~Model) + 
        xlim(15000 + 1e3, 28000 - 1e3) + 
        ylim(3765000 - 3e3, 3775000 + 3e3) + 
        coord_equal() + 
        theme_few() + 
        theme(panel.border = element_blank()) +
        xlab('Easting (m)') + 
        ylab('Northing (m)') 
      
      f = file.path('output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '.png', sep = ''))
      ggsave(pl, filename = f, width = 18/2, height = 9/2, dpi = 'print')
      
      #
      # plot meta data and summary
      #
      
      res = list(
        animal = pkg$tagID,
        observation_time = pred_times[ptime_with_depth],
        depth_at_observation = pkg$mapped_data[
          which(pkg$mapped_times == pred_times[ptime_with_depth]), 'depth'
        ],
        path = f
      )
      
      res
    }
  )
)
