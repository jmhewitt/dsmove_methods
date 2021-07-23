library(metR)
library(ggplot2)
library(ggthemes)
library(viridis)
library(dplyr)

# tar_load(whale_ll_srtm30_filtered)

# whale_ll_srtm30_filtered = readRDS('whale_ll_srtm30_filtered.rds')

whale_ll_srtm30_filtered = readRDS('whale_ll_srtm30_filtered_alt_no5422.rds')

source('R/utils/log_add.R')
source('R/utils/distribution_parameters.R')


priors = list(
  speed = gamma.param(mu = .5, sd = .5)
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
         filter( is.finite(ll), ll > -700  ),
       aes(x = speed, y = betaAR, fill = ll, z = ll)) + 
  geom_raster() + 
  geom_contour2(col = 'grey50') +
  geom_text_contour(col = 'grey90') +
  theme_few() + 
  scale_fill_viridis(direction = -1)


speed_marginal = whale_ll_srtm30_filtered %>%
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
  ungroup()

mean.marginal = function(x, p) {
  # Mean when given values and probabilities
  # 
  # Parameters:
  #  x - grid of values
  #  p - probability of value
  
  sum(x * p)
}

hpd.marginal = function(x, p, level = .95) {
  # HPD interval when given values and probabilities
  # 
  # Parameters:
  #  x - grid of values
  #  p - probability of value
  
  # ensure input is ordered
  o = order(x)
  x = x[o]
  p = p[o]
  
  cdf = cumsum(p)
  
  lwr.ind = 1
  upr.ind = min(which(cdf - cdf[lwr.ind] >= level))
  len = x[upr.ind] - x[lwr.ind]
  for(i in 2:length(x)) {
    valid = cdf - cdf[i] >= level
    if(any(valid)) {
      upr.ind.prop = min(which(valid))
      len.prop = x[upr.ind.prop] - x[i]
      if(len.prop < len) {
        lwr.ind = i
        upr.ind = upr.ind.prop
      }
    }
  }
  
  x[c(lwr.ind, upr.ind)]
}

hpd_speed = hpd.marginal(
  x = speed_marginal$speed,
  p = exp(speed_marginal$lp) * diff(speed_marginal$speed)[1]
)

mean_speed = mean.marginal(
  x = speed_marginal$speed,
  p = exp(speed_marginal$lp) * diff(speed_marginal$speed)[1]
)

ggplot(speed_marginal, aes(x = speed, y = exp(lp))) + 
  stat_function(
    fun = function(x) { dgamma(x = x, shape = priors$speed['shape'],
                               rate = priors$speed['rate'], log = FALSE)},
                lty = 3, geom = 'line') +
  geom_line() + 
  geom_vline(xintercept = hpd_speed, lty = 3) + 
  theme_few() 


whale_ll_srtm30_filtered %>% 
  group_by(betaAR) %>% 
  summarise(
    lp = ll +
      dgamma(x = speed, shape = priors$speed['shape'],
             rate = priors$speed['rate'], log = TRUE) +
      log(diff(unique(whale_ll_srtm30_filtered$speed))[1]) -
      lC
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = betaAR, y = exp(lp))) + 
  geom_line() + 
  theme_few() 
