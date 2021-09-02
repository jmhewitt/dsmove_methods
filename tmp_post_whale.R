library(metR)
library(ggplot2)
library(ggthemes)
library(viridis)
library(dplyr)

# tar_load(whale_ll_srtm30_filtered)

# whale_ll_srtm30_filtered = readRDS('whale_ll_srtm30_filtered.rds')

# whale_ll_srtm30_filtered = readRDS('whale_ll_srtm30_filtered_alt_no5422.rds')
whale_ll_srtm30_filtered = readRDS('dcc/whale_ll_srtm30_filtered_alt_orientation.rds')

# rescale speed grid to align with average meters traveled per cell transition
whale_ll_srtm30_filtered$speed = whale_ll_srtm30_filtered$speed * 623.7801 / 1187.295

# shrink beta0 toward the speed there should be if betaAR = 0. ideally it might 
# make more sense to be able to shrink toward a known speed range, independently
# of making assumptions about betaAR
log(1.2) - log(4 * 623.7801)


source('R/utils/log_add.R')
source('R/utils/distribution_parameters.R')


priors = list(
  # speed = gamma.param(mu = .5, sd = .5),
  speed = gamma.param(mu = 1, sd = 1),
  # speed = c(shape = 2, rate = 1)
  betaAR = c(mean = 0, sd = 1e2),
  beta0_hanks = c(mean = 0, sd = 1e2)
)

# one way to potentially help set the speed prior
qgamma(.99, shape = priors$speed['shape'], rate = priors$speed['rate'])

#
# back-transform hanks model parameterizations
#

# append induced hanks location parameter
whale_ll_srtm30_filtered$beta0_hanks = 
  log(whale_ll_srtm30_filtered$speed) -
  log(623.7801) -
  log(2 + exp(-whale_ll_srtm30_filtered$betaAR) + 
        exp(whale_ll_srtm30_filtered$betaAR) )

whale_ll_srtm30_filtered$lp_hanks = 
  # likelihood
  whale_ll_srtm30_filtered$ll + 
  # prior
  dnorm(x = whale_ll_srtm30_filtered$beta0_hanks, 
        mean = priors$beta0_hanks['mean'], 
        sd = priors$beta0_hanks['sd'], log = TRUE) +
  dnorm(x = whale_ll_srtm30_filtered$betaAR,
        mean = priors$betaAR['mean'], 
        sd = priors$betaAR['sd'], log = TRUE) -
  # jacobian for transformation
  log(whale_ll_srtm30_filtered$speed)

whale_ll_srtm30_filtered$lp = 
  # likelihood
  whale_ll_srtm30_filtered$ll + 
  # prior
  dgamma(x = whale_ll_srtm30_filtered$speed,
         shape = priors$speed['shape'],
         rate = priors$speed['rate'],
         log = TRUE) + 
  dnorm(x = whale_ll_srtm30_filtered$betaAR,
        mean = priors$betaAR['mean'], 
        sd = priors$betaAR['sd'], log = TRUE)
  

lC = log_sum(
  whale_ll_srtm30_filtered$lp + 
  log(diff(unique(whale_ll_srtm30_filtered$speed))[1]) + 
  log(diff(unique(whale_ll_srtm30_filtered$betaAR))[1])
) 

lC_hanks = log_sum(
  whale_ll_srtm30_filtered$lp_hanks + 
    log(diff(unique(whale_ll_srtm30_filtered$speed))[1]) + 
    log(diff(unique(whale_ll_srtm30_filtered$betaAR))[1])
) 


pl = ggplot(whale_ll_srtm30_filtered %>% 
              mutate(ll = lp_hanks - lC_hanks
              ) %>% 
          filter( ll > -20  ),
       aes(x = speed, y = betaAR, fill = ll, z = ll)) + 
  # geom_tile(height = 0.3157895, width = 0.3139539) +
  geom_raster() + 
  geom_contour2(col = 'grey60') + 
  geom_text_contour(col = 'grey90') + 
  theme_few() + 
  scale_fill_viridis(direction = -1) + 
  ggtitle('Zc093 prior mean = 1, sd = 1')

pl

ggsave(pl, filename = 'Zc093_mu_1_sd_1.pdf')

speed_marginal = whale_ll_srtm30_filtered %>%
  group_by(speed) %>% 
  summarise(
    lp = log_sum(
      lp + 
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

var.marginal = function(x, p) {
  # Variance when given values and probabilities
  # 
  # Parameters:
  #  x - grid of values
  #  p - probability of value
  
  mu = mean.marginal(x = x, p = p)
  
  sum((x-mu)^2 * p)
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

sd_speed = sqrt(var.marginal(
  x = speed_marginal$speed,
  p = exp(speed_marginal$lp) * diff(speed_marginal$speed)[1]
))


print(hpd_speed)
print(mean_speed)
print(sd_speed)

pl = ggplot(speed_marginal, aes(x = speed, y = exp(lp))) + 
  stat_function(
    fun = function(x) { dgamma(x = x, shape = priors$speed['shape'],
                               rate = priors$speed['rate'], log = FALSE)},
    # fun = function(x) {
    #   # prior assumes betaAR = 0, which is a bit of a kludge for the transformation
    #   beta0 = log(x) - log(623.7801) - log(4)
    #   dnorm(x = beta0, mean = priors$beta0_hanks['mean'],
    #         sd = priors$beta0_hanks['sd'], log = FALSE) / x
    # },
                lty = 3, geom = 'line') +
  geom_line() + 
  theme_few() + 
  ggtitle('Zc093 prior mean = 1, sd = 1; post. mean = 1.8, sd = .5')

pl

ggsave(pl, filename = 'Zc093_mu_1_sd_1_speed.pdf')

pl = whale_ll_srtm30_filtered %>% 
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

pl

ggsave(pl, filename = 'Zc073_fastloc_mu_3_sd_1_betaAR.pdf')


# CONCLUSION for comparison of 2 and 4, in speed space
#  - for weak priors, the posteriors between the two parameterizations are 
#    relatively similar (wrt the speed parameter), but it is much easier to 
#    interpret the (mu,\beta1) prior.