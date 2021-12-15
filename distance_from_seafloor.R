library(raster)
library(sp)
library(metR)
library(ggplot2)
library(ggthemes)
library(ggnewscale)
library(lubridate)
library(dplyr)
library(tidyr)

source('R/utils/log_posteriors.R')

# load data
targets::tar_load(whale_pkg)
pkg = whale_pkg$whale_pkg_5bc05e7c

# load posteriors
post_marginal = readRDS('marginal_locations.rds')
post_crawl = readRDS('crawl_discretized_post_locations100.rds')

# extract prediction times
pred_times = as.POSIXct(
  x = sapply(post_marginal, function(r) r$meta$time), 
  tz = 'UTC', origin = '1970-01-01 00:00.00 UTC'
)

targets::tar_load(cee_start)
targets::tar_load(whale_domain)

# map the observed lat-lons to the bathymetry lon/lat grid
# (needed to set initial filtering locations)
coords = pkg$grid$coords
lons = pkg$grid$lons
lats = pkg$grid$lats


# manifest of marginal distributions
depth_times = do.call(rbind, lapply(post_marginal, function(r) r$meta)) %>% 
  mutate(pred_index = 1:n()) %>% 
  filter(model_specification == 'speed_parameterization')

seafloor_summaries = do.call(rbind, apply(depth_times, 1, function(r) {
  
  pred_index = as.numeric(r['pred_index'])
  
  # CTDS-based posterior location
  ctds_post = post_marginal[[pred_index]]$posterior
  
  # normalize posterior
  ctds_post$lp = ctds_post$lp - dsmovetools2d:::log_sum_c(ctds_post$lp)
  
  # crawl-based posterior location
  crawl_post = do.call(rbind, lapply(post_crawl, function(imputation) {
    ind = which.min(abs(imputation$time - depth_times$time[pred_index]))
    imputation[ind, c('mapped_lon', 'mapped_lat', 'lon_ind', 'lat_ind')]
  }))
  colnames(crawl_post)[3:4] = c('lon_to_ind', 'lat_to_ind')
  
  # condense to empirical probability distribution
  crawl_post = crawl_post %>% 
    mutate(lp = log(1 /n())) %>%
    group_by(lon_to_ind, lat_to_ind) %>%
    summarise(lp = dsmovetools2d:::log_sum_c(lp)) %>% 
    ungroup()
  
  # extract seafloor depths for locations
  ctds_post$seafloor_depth = apply(ctds_post, 1, function(r) {
    pkg$bathymetry_matrix[r['lat_to_ind'], r['lon_to_ind']]
  })
  crawl_post$seafloor_depth = apply(crawl_post, 1, function(r) {
    pkg$bathymetry_matrix[r['lat_to_ind'], r['lon_to_ind']]
  })
  
  observed_depth = pkg$mapped_data[
    which(pkg$mapped_times == depth_times$time[pred_index]), 'depth'
  ]
  
  # marginalize distributions for seafloor depth
  ctds_post = ctds_post %>% 
    group_by(seafloor_depth) %>%
    summarize(lp = dsmovetools2d:::log_sum_c(lp)) %>% 
    ungroup() %>% 
    mutate(distance_from_seafloor = observed_depth - seafloor_depth)
  crawl_post = crawl_post %>% 
    group_by(seafloor_depth) %>%
    summarize(lp = dsmovetools2d:::log_sum_c(lp)) %>% 
    ungroup() %>% 
    mutate(distance_from_seafloor = observed_depth - seafloor_depth)
  
  # combine datasets
  df = rbind(
    cbind(ctds_post, Model = r['data_subset']),
    cbind(crawl_post, Model = 'CRAWL as AID')
  ) %>% 
    pivot_longer(cols = c('seafloor_depth', 'distance_from_seafloor'), 
                 names_to = 'Metric', values_to = 'value') %>%
    group_by(Model, Metric) %>%
    summarise(lwr = hpd.marginal(x = value, p = exp(lp))[1],
              upr = hpd.marginal(x = value, p = exp(lp))[2],
              mean = mean.marginal(x = value, p = exp(lp))) %>% 
    ungroup() %>% 
    mutate(time = depth_times$time[pred_index],
           subset = r['data_subset'],
           observed_depth = observed_depth)
  
  df
}))


pl = ggplot(seafloor_summaries %>% 
              filter(Metric == 'seafloor_depth') %>% 
              mutate(Model = gsub('CRAWL as AID', 'CTCRW-AID', Model),
                     Model = gsub('with_depth', 'State space (w/Depth)', Model),
                     Model = gsub('without_depth', 'State space', Model)), 
       aes(x = time, y = abs(mean), ymin = abs(lwr), ymax = abs(upr))) + 
  geom_pointrange(col = '#a4b9b5') + 
  xlab('Time') + 
  scale_y_reverse('Seafloor depth (m)') + 
  geom_line(mapping = aes(x = time, y = abs(observed_depth)), 
            inherit.aes = FALSE) +
  geom_point(mapping = aes(x = time, y = abs(observed_depth)), 
             inherit.aes = FALSE) + 
  facet_grid(~Model, scales = 'free') + 
  scale_x_datetime() + 
  theme_few()

ggsave(pl, filename = 'seafloor_depth.png', width = 8, height = 4, 
       dpi = 'print')
