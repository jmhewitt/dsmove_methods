plot_targets = list(
  
  # 3-panel plot comparing posterior locations between different approximations
  # at times when we have depth data but no location data
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
          mutate(Model = 'HMM (w/Depth)'), 
        loc_without_depth %>% 
          select(lon.projected, lat.projected, lp) %>% 
          filter(lp >= log(.005)) %>%
          mutate(Model = 'HMM'),
        loc_crawl %>% 
          select(lon.projected, lat.projected, lp) %>% 
          mutate(Model = 'CTCRW-AID')
      ) %>% 
        mutate(Model = factor(Model, levels = c('CTCRW-AID', 
                                                'HMM',
                                                'HMM (w/Depth)')))
      
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
  ),
  
  # comparing posterior locations between different approximations at times when
  # we have low-quality location data (i.e., from ARGOS)
  tar_target(
    name = zc095_post_uncertainty_comparison, 
    command = {
      
      #
      # load data
      #
      
      # load data
      pkg = whale_pkg[[1]]
      # load posterior locations from CTDS model
      post_marginal = readRDS('whale_marginal_additional_location_post2.rds')
      # load posterior locations from crawl
      post_crawl = readRDS(
        'crawl_discretized_post_locations_additional2_100.rds'
      )
      
      # extract prediction info
      manifest = do.call(rbind, lapply(post_marginal, function(x) x$meta)) %>% 
        mutate(ind = 1:n())
      
      # extract prediction times
      pred_times = unique(manifest$time)
      
      # plotting projection
      plot_crs = '+proj=aea +lat_1=34.5 +lat_2=36 +lon_0=-75'
      
      # projected bathymetry coordinates, for plotting
      bathy.projected = projectRaster(from = whale_domain, crs = plot_crs)
      
      
      #
      # compute uncertainty per timepoint
      #
      
      df_uncertainty = do.call(rbind, lapply(
        1:length(pred_times), function(tind) {
         # target prediction time
         ptime_without_depth = manifest %>% 
           filter(time == pred_times[tind],
                  data_subset == 'without_depth',
                  model_specification == 'speed_parameterization') %>% 
           select(ind) %>% 
           as.numeric()
         ptime_with_depth = manifest %>% 
           filter(time == pred_times[tind],
                  data_subset == 'with_depth',
                  model_specification == 'speed_parameterization') %>% 
           select(ind) %>% 
           as.numeric()
         # CTDS-based posterior location
         x_with_depth = post_marginal[[ptime_with_depth]]$posterior
         x_without_depth = post_marginal[[ptime_without_depth]]$posterior
         # normalize posteriors
         x_with_depth$lp = x_with_depth$lp - 
           dsmovetools2d:::log_sum_c(x_with_depth$lp)
         x_without_depth$lp = x_without_depth$lp - 
           dsmovetools2d:::log_sum_c(x_without_depth$lp)
         # crawl-based posterior location
         x_crawl = do.call(rbind, lapply(post_crawl, function(imputation) {
           ind = which.min(abs(imputation$time - pred_times[tind]))
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
         # plotting projection
         plot_crs = '+proj=aea +lat_1=34.5 +lat_2=36 +lon_0=-75'
         # project crawl posterior for marginal location
         loc_crawl = project_to_grid(
           data = x_crawl, 
           data_crs = crs(whale_domain), 
           rasterObj = bathy.projected, 
           plot_crs = plot_crs
         )
         # convex hull of crawl points
         crawl_pattern = ppp(
           x = loc_crawl$lon.projected, 
           y = loc_crawl$lat.projected, 
           window = owin(xrange = range(loc_crawl$lon.projected), 
                         yrange = range(loc_crawl$lat.projected))
         )
         crawl_pattern_hull = convexhull(crawl_pattern)
         # area of other HPD regions
         hpdregion_area = function(lp, scale, prob = 0.95) {
           # Parameters:
           #  lp - log-probability value for point pattern
           #  prob - probability threshold
           #  scale - area of each cell in the point pattern
           lp = sort(x = lp, decreasing = TRUE)
           lprob = log(prob)
           lcdf = lp[1]
           ncells = 1
           len = length(lp)
           while(lcdf <= lprob) {
             if(ncells >= len) { break }
             ncells = ncells + 1
             lcdf = dsmovetools2d:::log_sum_c(c(lcdf, lp[ncells]))
           }
           ncells * scale
         }
         
         hpdregion = function(locs, lp, prob = .95) {
           # Parameters:
           #  locs - lon/lat locations associated with lp values
           #  lp - log-probability value for point pattern
           #  prob - probability threshold
           o = order(x = lp, decreasing = TRUE)
           lprob = log(prob)
           lcdf = lp[o[1]]
           ncells = 1
           len = length(lp)
           while(lcdf <= lprob) {
             if(ncells >= len) { break }
             ncells = ncells + 1
             lcdf = dsmovetools2d:::log_sum_c(c(lcdf, lp[o[ncells]]))
           }
           # locations that comprise the hpd region
           locs[o[1:ncells],c('lon', 'lat')]
         }
         
         hpdregion_with_depth = hpdregion(
           locs = x_with_depth[,c('lon','lat')],
           lp = x_with_depth$lp
         )
         
         hpdregion_without_depth = hpdregion(
           locs = x_without_depth[,c('lon','lat')],
           lp = x_without_depth$lp, 
         )
         
         # project crawl posterior for marginal location
         hpdregion_with_depth_projected = project_to_grid(
           data = hpdregion_with_depth, 
           data_crs = crs(whale_domain), 
           rasterObj = bathy.projected, 
           plot_crs = plot_crs
         )
         
         # project crawl posterior for marginal location
         hpdregion_without_depth_projected = project_to_grid(
           data = hpdregion_without_depth, 
           data_crs = crs(whale_domain), 
           rasterObj = bathy.projected, 
           plot_crs = plot_crs
         )
         
         # package uncertainty results
         data.frame(
           time = pred_times[tind],
           area = c(
             # area of convex hull in km^2
             area(crawl_pattern_hull) / 1e6,
             # area in km^2
             hpdregion_area(lp = x_without_depth$lp, scale = (1187.295/1e3)^2),
             hpdregion_area(lp = x_with_depth$lp, scale = (1187.295/1e3)^2)
           ),
           prop_ss_w_depth_contained_within_region = c(
             mean(sp::point.in.polygon(
               point.x = hpdregion_with_depth_projected$lon.projected,
               point.y = hpdregion_with_depth_projected$lat.projected,
               pol.x = crawl_pattern_hull$bdry[[1]]$x,
               pol.y = crawl_pattern_hull$bdry[[1]]$y
             ) != 0),
             nrow(plyr::match_df(x = hpdregion_without_depth,
                                 y = hpdregion_with_depth)) / 
               nrow(hpdregion_with_depth),
             1
           ),
           prop_ss_wout_depth_contained_within_region = c(
             mean(sp::point.in.polygon(
               point.x = hpdregion_without_depth_projected$lon.projected,
               point.y = hpdregion_without_depth_projected$lat.projected,
               pol.x = crawl_pattern_hull$bdry[[1]]$x,
               pol.y = crawl_pattern_hull$bdry[[1]]$y
             ) != 0),
             1,
             nrow(plyr::match_df(x = hpdregion_with_depth,
                                 y = hpdregion_without_depth)) / 
               nrow(hpdregion_with_depth)
           ),
           method = c('CTCRW-AID', 'HMM', 'HMM (w/Depth)')
         )
      }))
      
      # times at which location observations are available
      obs_times = pkg$mapped_times[is.finite(pkg$mapped_data[,'lon_ind'])]
      obs_times = obs_times[
        (min(df_uncertainty$time) <= obs_times) & 
        (obs_times <= max(df_uncertainty$time))
      ]
      
      
      
      pl = ggplot(df_uncertainty, aes(x = time, y = area, col = method, 
                                      lty = method, shape = method)) + 
        # times of location observations
        geom_vline(xintercept = obs_times, lty = 3) + 
        # posterior uncertainty in location
        geom_line() + 
        geom_point() + 
        # formatting
        scale_color_brewer(type = 'qual', palette = 'Dark2') + 
        ylab(expression('95% HPD region'~(km^2))) + 
        xlab('Time (UTC)') + 
        theme_few() + 
        theme(legend.title = element_blank())
      
      f = file.path('output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '_uncertainty.png', sep = ''))
      ggsave(pl, filename = f, width = 18/2, height = 9/2, dpi = 'print')
      
      
      pl = ggplot(df_uncertainty %>% 
                    filter(method %in% c('CTCRW-AID', 'HMM')), 
                  aes(x = time, y = prop_ss_w_depth_contained_within_region, 
                      col = method)) + 
        # times of location observations
        geom_vline(xintercept = obs_times, lty = 3) + 
        # posterior uncertainty in location
        geom_line() + 
        geom_point() + 
        # formatting
        scale_color_brewer(
          '$\\mathcal M_{a}$', type = 'qual', palette = 'Dark2'
        ) + 
        ylab('$\\mathcal O(\\mathcal M_a, t;\\mathcal M_0)$') +
        xlab('Time (UTC)') + 
        theme_few() + 
        theme(axis.title.y = element_text(angle = 0, vjust = .5),
              text = element_text(family = 'sans'))
      
      f = file.path('output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '_overlap.tex', sep = ''))

      tikz(
        file = f, width = 8, height = 9/2, standAlone = TRUE,
        packages = c(
          "\\usepackage{tikz}",
          "\\usepackage[active,tightpage,psfixbb]{preview}",
          "\\PreviewEnvironment{pgfpicture}",
          "\\setlength\\PreviewBorder{0pt}",
          "\\usepackage{amssymb}",
          # set font to helvetica
          "\\tikzset{every picture/.style={/utils/exec={\\fontfamily{phv}}}}"
        )
      )
      print(pl)
      dev.off()
      
      cwd = getwd()
      setwd(dirname(f))
      tools::texi2dvi(file = basename(f), pdf = TRUE, clean = TRUE)
      setwd(cwd)
      
      
      pl = ggplot(df_uncertainty %>% filter(method %in% c('CTCRW-AID')), 
                  aes(x = time, y = prop_ss_wout_depth_contained_within_region, 
                      col = method)) + 
        # times of location observations
        geom_vline(xintercept = obs_times, lty = 3) + 
        # posterior uncertainty in location
        geom_line() + 
        geom_point() + 
        # formatting
        scale_color_brewer(
          '$\\mathcal M_{a}$', type = 'qual', palette = 'Dark2'
        ) + 
        ylab('$\\mathcal O(\\mathcal M_a, t;\\mathcal M_{a\'})$') +
        xlab('Time (UTC)') + 
        theme_few() + 
        theme(axis.title.y = element_text(angle = 0, vjust = .5),
              text = element_text(family = 'sans'))
      
      f = file.path('output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '_overlap_ss_without_depth.tex', 
                             sep = ''))
      
      tikz(
        file = f, width = 8, height = 9/2, standAlone = TRUE,
        packages = c(
          "\\usepackage{tikz}",
          "\\usepackage[active,tightpage,psfixbb]{preview}",
          "\\PreviewEnvironment{pgfpicture}",
          "\\setlength\\PreviewBorder{0pt}",
          "\\usepackage{amssymb}",
          # set font to helvetica
          "\\tikzset{every picture/.style={/utils/exec={\\fontfamily{phv}}}}"
        )
      )
      print(pl)
      dev.off()
      
      cwd = getwd()
      setwd(dirname(f))
      tools::texi2dvi(file = basename(f), pdf = TRUE, clean = TRUE)
      setwd(cwd)
      
      
      # averages!
      df_uncertainty %>% 
        filter(method %in% c('CTCRW-AID')) %>% 
        select(prop_ss_wout_depth_contained_within_region) %>% 
        unlist() %>% 
        mean()
      
      df_uncertainty %>% 
        filter(method %in% c('CTCRW-AID', 'HMM'))
      
      #
      # panel plot at last predicted timepoint
      #
      
      #
      # subset data
      #
      
      
      # target prediction time
      tind = 13
      ptime_without_depth = manifest %>% 
        filter(time == pred_times[tind],
               data_subset == 'without_depth',
               model_specification == 'speed_parameterization') %>% 
        select(ind) %>% 
        as.numeric()
      ptime_with_depth = manifest %>% 
        filter(time == pred_times[tind],
               data_subset == 'with_depth',
               model_specification == 'speed_parameterization') %>% 
        select(ind) %>% 
        as.numeric()
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
        ind = which.min(abs(imputation$time - pred_times[tind]))
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
      
      
      # #
      # # build comparison with raw data
      # #
      # 
      # pkg_ind = which(pkg$mapped_times == pred_times[tind])
      # x_obs = data.frame(coords)
      # colnames(x_obs) = c('lon', 'lat')
      # # evaluate raw likelihood for observation across grid
      # x_obs$lp = dsmovetools2d:::GpsLikEvalGroup(
      #   obs_lons = pkg$mapped_data[,'lon'], 
      #   obs_lats = pkg$mapped_data[,'lat'], 
      #   semi_majors = pkg$mapped_data[,'semi_major'], 
      #   semi_minors = pkg$mapped_data[,'semi_minor'], 
      #   orientations = pkg$mapped_data[,'orientation'], 
      #   alpha = .05, 
      #   test_lon = x_obs[,'lon'], 
      #   test_lat = x_obs[,'lat'], 
      #   ind = pkg_ind - 1
      # )
      # # standardize distribution
      # x_obs$lp = x_obs$lp - dsmovetools2d:::log_sum_c(x_obs$lp)
      
      
      #
      # project data
      #
      
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
      # loc_raw = project_to_grid(
      #   data = x_obs, 
      #   data_crs = crs(whale_domain), 
      #   rasterObj = bathy.projected, 
      #   plot_crs = plot_crs
      # )
      
      
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
          mutate(Model = 'HMM (w/Depth)'), 
        loc_without_depth %>% 
          select(lon.projected, lat.projected, lp) %>% 
          filter(lp >= log(.005)) %>%
          mutate(Model = 'HMM'),
        loc_crawl %>% 
          select(lon.projected, lat.projected, lp) %>% 
          mutate(Model = 'CTCRW-AID')
      ) %>% 
        mutate(Model = factor(Model, levels = c('CTCRW-AID', 
                                                'HMM',
                                                'HMM (w/Depth)')))
      
      xrange = c(15000 + 1e3, 28000 - 1e3)
      yrange = c(3765000 - 4e3, 3775000 + 3e3) - 6e3
      
      pl = ggplot(mapping = aes(x = x, y = y)) + 
        # bathymetry
        geom_raster(mapping = aes(fill = depth), 
                    data = bathy_df %>% 
                      filter(xrange[1] <= x, x <= xrange[2],
                             yrange[1] <= y, y <= yrange[2])) + 
        scale_fill_distiller('Depth (m)', palette = 'Blues', direction = -1,
                             na.value = 'transparent', trans = 'reverse') + 
        # land mask
        new_scale_fill() + 
        guides(fill = 'none') + 
        geom_raster(mapping = aes(fill = mask), 
                    data = land_df %>% 
                      filter(xrange[1] <= x, x <= xrange[2],
                             yrange[1] <= y, y <= yrange[2])) + 
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
        geom_contour2(mapping = aes(z = depth), 
                      data = bathy_contour %>% 
                        filter(xrange[1] <= x, x <= xrange[2],
                               yrange[1] <= y, y <= yrange[2]), 
                      col = 'grey50', breaks = seq(500,2000,by=100)) + 
        geom_text_contour(mapping = aes(z = depth), 
                          data = bathy_contour %>% 
                            filter(xrange[1] <= x, x <= xrange[2],
                                   yrange[1] <= y, y <= yrange[2]), 
                          breaks = seq(500,2000,by=100)) + 
        # formatting
        facet_wrap(~Model, nrow = 1) + 
        xlim(xrange) + 
        ylim(yrange) + 
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
        observation_time = pred_times[tind],
        depth_at_observation = pkg$mapped_data[
          which(pkg$mapped_times == pred_times[tind]), 'depth'
        ],
        path = f
      )
      
      res
    }
  )
)
