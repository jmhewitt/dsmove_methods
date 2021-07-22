library(metR)
library(ggplot2)
library(ggthemes)
library(viridis)
library(dplyr)

# tar_load(whale_ll_srtm30_filtered)

# whale_ll_srtm30_filtered = readRDS('whale_ll_srtm30_filtered.rds')

whale_ll_srtm30_filtered = readRDS('whale_ll_srtm30_filtered_alt_byrow_bathy.rds')

source('R/utils/log_add.R')

priors = list(
  speed = c('shape' = 2.25, 'rate' = 1.5)
)

lC = log_sum(
  whale_ll_srtm30_filtered$ll + 
  dgamma(x = whale_ll_srtm30_filtered$speed,
         shape = priors$speed['shape'],
         rate = priors$speed['rate'],
         log = TRUE) +
  log(diff(unique(whale_ll_srtm30_filtered$speed))[1]) + 
  log(diff(unique(whale_ll_srtm30_filtered$betaAR))[1])
) 


ggplot(whale_ll_srtm30_filtered %>% 
         filter( is.finite(ll)  ),
       aes(x = speed, y = betaAR, fill = ll, z = ll)) + 
  geom_raster() + 
  geom_contour2(col = 'grey50') +
  geom_text_contour(col = 'grey90') +
  theme_few() + 
  scale_fill_viridis(direction = -1)


whale_ll_srtm30_filtered %>%
  group_by(speed) %>% 
  summarise(
    lp = log_sum(
      ll +
      dgamma(x = speed, shape = priors$speed['shape'],
             rate = priors$speed['rate'], log = TRUE) +
      log(diff(unique(whale_ll_srtm30_filtered$betaAR))[1]) -
      lC
      ) 
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = speed, y = exp(lp))) + 
  stat_function(
    fun = function(x) { dgamma(x = x, shape = priors$speed['shape'],
                               rate = priors$speed['rate'], log = FALSE)},
                lty = 3, geom = 'line') +
  geom_line() + 
  theme_few() 


whale_ll_srtm30_filtered %>% 
  group_by(betaAR) %>% 
  summarise(
    lp = log_sum(ll) + dgamma(x = speed, shape = 9, rate = 3, log = TRUE)
  ) %>% 
  # summarise(lp = log_sum(lp)  + dgamma(x = speed, shape = 2, rate = 1, log = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = betaAR, y = exp(lp))) + 
  geom_line() + 
  theme_few() 



inf_nodes = unique(whale_ll_srtm30_filtered %>% filter(is.infinite(ll)) %>% select(computation_node) %>% unlist())
non_inf_nodes = setdiff(unique(whale_ll_srtm30_filtered$computation_node), inf_nodes)
