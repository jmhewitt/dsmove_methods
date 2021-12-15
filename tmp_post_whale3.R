library(metR)
library(ggplot2)
library(ggthemes)
library(viridis)
library(dplyr)
library(tidyr)

source('R/utils/log_add.R')
source('R/utils/distribution_parameters.R')
source('R/utils/log_posteriors.R')

# load likelihood surfaces
res = readRDS('whale_ll_subsets.rds')

# define prior distributions
priors = list(
  # speed = gamma.param(mu = .5, sd = .5),
  speed = gamma.param(mu = 1, sd = 1),
  # speed = c(shape = 2, rate = 1)
  betaAR = c(mean = 0, sd = 1e2),
  beta0 = c(mean = 0, sd = 1e2)
  # beta0 = c(mean = -10, sd = 1)
)

# one way to potentially help set the speed prior
qgamma(.99, shape = priors$speed['shape'], rate = priors$speed['rate'])

# compute joint posteriors for different parameterizations by data subset
joint_posteriors_by_data = lapply(split(res, res$subset), function(r) {

  # extract grid dimensions
  dx = diff(sort(unique(r$speed))[2:3])
  dy = diff(sort(unique(r$betaAR))[2:3])
  
  # reparameterized model parameter
  beta0 = log(r$speed) - log(r$cell_size) - 
    log(2 + exp(-r$betaAR) + exp(r$betaAR))
  
  # log-prior for reparameterized model, normalized to grid
  lprior_reparam = log_density_gridded(
    x = r$speed, y = r$betaAR, dx = dx, dy = dy,
    ld = # prior for reparameterized parameters
      dnorm(x = beta0, mean = priors$beta0['mean'], 
            sd = priors$beta0['sd'], log = TRUE) + 
      dnorm(x = r$betaAR, mean = priors$betaAR['mean'], 
            sd = priors$betaAR['sd'], log = TRUE) - 
      # jacobian for transformation
      log(r$speed) 
  )
  
  # log-prior for model, normalized to grid
  lprior = log_density_gridded(
    x = r$speed, y = r$betaAR, dx = dx, dy = dy,
    ld = # prior for parameters
      dgamma(x = r$speed, shape = priors$speed['shape'], 
             rate = priors$speed['rate'], log = TRUE) + 
      dnorm(x = r$betaAR, mean = priors$betaAR['mean'], 
            sd = priors$betaAR['sd'], log = TRUE)
  )
  
  # posterior surfaces 
  list(
    # direct transition rate model
    speed_parameterization = data.frame(
      data_subset = r$subset,
      speed = r$speed,
      betaAR = r$betaAR,
      beta0 = beta0,
      lprior = lprior,
      parameterization = 'speed',
      lpost = log_density_gridded(
        x = r$speed, y = r$betaAR, ld = r$ll + lprior, dx = dx, dy = dy
      )
    ),
    # original transition rate model
    beta0_parameterization = data.frame(
      data_subset = r$subset,
      speed = r$speed,
      betaAR = r$betaAR,
      beta0 = beta0,
      lprior = lprior_reparam,
      parameterization = 'beta0',
      lpost = log_density_gridded(
        x = r$speed, y = r$betaAR, ld = r$ll + lprior_reparam, dx = dx, dy = dy
      )
    )
  )
})

# compute marginal distributions for different parameterizations by data subset
marginal_posteriors_by_data = lapply(joint_posteriors_by_data, function(param) {
  # process each parameterization separately
  lapply(param, function(r) {

    # extract grid dimensions
    dx = diff(sort(unique(r$speed))[2:3])
    dy = diff(sort(unique(r$betaAR))[2:3])

    # marginal prior and posterior for speed
    speed_marginal = rbind(
      log_marginal_x(x = r$speed, y = r$betaAR, ld = r$lpost, dy = dy) %>% 
        mutate(density = 'posterior', 
               parameter = 'speed',
               parameterization = r$parameterization[1],
               data_subset = r$data_subset[1]),
      log_marginal_x(x = r$speed, y = r$betaAR, ld = r$lprior, dy = dy) %>% 
        mutate(density = 'prior',
               parameter = 'speed',
               parameterization = r$parameterization[1],
               data_subset = r$data_subset[1])
    )
    
    # marginal prior and posterior for directional persistence
    persistence_marginal = rbind(
      log_marginal_x(x = r$betaAR, y = r$speed, ld = r$lpost, dy = dx) %>% 
        mutate(density = 'posterior', 
               parameter = 'betaAR',
               parameterization = r$parameterization[1],
               data_subset = r$data_subset[1]),
      log_marginal_x(x = r$betaAR, y = r$speed, ld = r$lprior, dy = dx) %>% 
        mutate(density = 'prior',
               parameter = 'betaAR',
               parameterization = r$parameterization[1],
               data_subset = r$data_subset[1])
    )
    
    # package results
    list(
      speed = speed_marginal,
      persistence = persistence_marginal
    )
  })
})
  
marginal_summaries = lapply(marginal_posteriors_by_data, function(subset) {
  lapply(subset, function(parameterization) {
    lapply(parameterization, function(param) {
      # standardize marginal distribution
      
      param %>% 
        group_by(density) %>% 
        summarise(
          mean = mean.marginal(x = x, p = exp(ld) * diff(x)[2]),
          var = var.marginal(x = x, p = exp(ld) * diff(x)[2]),
          lwr = hpd.marginal(x = x, p = exp(ld) * diff(x)[2])[1],
          upr = hpd.marginal(x = x, p = exp(ld) * diff(x)[2])[2],
        ) %>% 
        ungroup() %>% 
        mutate(
          parameterization = param$parameterization[1],
          parameter = param$parameter[1],
          data_subset = param$data_subset[1]
        )
    })
  })
})

# reducing uncertainty in movement speed and decreasing movement speed est.
rbind(
  marginal_summaries$with_depth$speed_parameterization$speed[1,],
  marginal_summaries$without_depth$speed_parameterization$speed[1,]
)

# finding evidence for directional persistence
rbind(
  marginal_summaries$with_depth$speed_parameterization$persistence[1,],
  marginal_summaries$without_depth$speed_parameterization$persistence[1,]
)

# plot joint posteriors
pl = ggplot(do.call(rbind, unlist(joint_posteriors_by_data, recursive = FALSE)),
            aes(x = speed, y = betaAR, fill = exp(lpost), z = exp(lpost))) + 
  # geom_tile(height = 0.3157895, width = 0.3139539) +
  # geom_raster() + 
  geom_contour2(col = 'grey60') + 
  geom_text_contour(col = 'grey20') + 
  facet_grid(data_subset~parameterization) + 
  theme_few() + 
  scale_fill_viridis(direction = -1) + 
  ggtitle('Posterior density surface comparisons')

pl

ggsave(pl, filename = 'posterior_surface_comparisons.pdf')


df = do.call(rbind, unlist(
    unlist(marginal_posteriors_by_data, recursive = FALSE), 
recursive = FALSE))

ggplot(df, aes(x = x, y = exp(ld), col = interaction(density, data_subset))) + 
  geom_line() + 
  ylim(0,1) + 
  facet_grid(parameterization~parameter, scales = 'free')

mean.marginal(
  x = marginal_posteriors_by_data$with_depth$speed_parameterization$speed$x,
  p = exp(marginal_posteriors_by_data$with_depth$speed_parameterization$speed$ld) * 
      diff(marginal_posteriors_by_data$with_depth$speed_parameterization$speed$x)[2]
)

mean.marginal(
  x = marginal_posteriors_by_data$with_depth$speed_parameterization$persistence$x,
  p = exp(marginal_posteriors_by_data$with_depth$speed_parameterization$persistence$ld) * 
    diff(marginal_posteriors_by_data$with_depth$speed_parameterization$persistence$x)[2]
)

marginal_summaries = speed_marginal %>% 
  group_by(subset, lp_type) %>% 
  summarise(
    mean = mean.marginal(x = speed, p = exp(lp) * dx),
    sd = sqrt(var.marginal(x = speed, p = exp(lp) * dx)),
    hpd_lwr = hpd.marginal(x = speed, p = exp(lp) * dx)[1],
    hpd_upr = hpd.marginal(x = speed, p = exp(lp) * dx)[2]
  )

marginal_summaries

pl = ggplot(speed_marginal %>% filter(lp_type %in% c('lp', 'lp_hanks')), 
            aes(x = speed, y = exp(lp), 
                col = interaction(subset, lp_type))) + 
  # priors
  geom_line(mapping = aes(x = speed, y = exp(lp), lty = lp_type), 
            data = speed_marginal %>% 
              filter(lp_type %in% c('log_prior_hanks', 'log_prior_speed')), 
            inherit.aes = FALSE, col = 'grey70') + 
  # posteriors
  geom_line() + 
  scale_color_brewer(type = 'qual', palette = 'Dark2') + 
  # formatting
  ylim(0, 1.25) + 
  theme_few() 

pl

ggsave(pl, filename = 'speed_comparisons.pdf')


pl = ggplot(persistence_marginal %>% filter(lp_type %in% c('lp', 'lp_hanks')), 
            aes(x = betaAR, y = exp(lp), 
                col = interaction(subset, lp_type))) + 
  # priors
  geom_line(mapping = aes(x = betaAR, y = exp(lp), lty = lp_type), 
            data = persistence_marginal %>% 
              filter(lp_type %in% c('log_prior_hanks', 'log_prior_speed')), 
            inherit.aes = FALSE, col = 'grey70') + 
  # posteriors
  geom_line() + 
  scale_color_brewer(type = 'qual', palette = 'Dark2') + 
  # formatting
  ylim(0, 1.25) + 
  theme_few() + 
  ggtitle('Zc093 prior mean = 1, sd = 1; post. mean = 1.8, sd = .5')


pl

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


#
# compute correlations
# 

# compute marginal distributions for different parameterizations by data subset
joint_summaries = lapply(joint_posteriors_by_data, function(param) {
  # process each parameterization separately
  lapply(param, function(r) {
    
    # extract grid dimensions
    dx = diff(sort(unique(r$speed))[2:3])
    dy = diff(sort(unique(r$betaAR))[2:3])
    
    # package results
    data.frame(
      exy = exp(dsmovetools2d:::log_sum_c(
        log(r$speed) + log(r$betaAR) + r$lpost + log(dx) + log(dy)
      )),
      data_subset = r$data_subset[1],
      parameterization = r$parameterization[1]
    )
  })
})


(joint_summaries$with_depth$speed_parameterization$exy - 
  marginal_summaries$with_depth$speed_parameterization$speed$mean[1] * 
  marginal_summaries$with_depth$speed_parameterization$persistence$mean[1]
) /
(
  sqrt(marginal_summaries$with_depth$speed_parameterization$speed$var[1]) *
  sqrt(marginal_summaries$with_depth$speed_parameterization$persistence$var[1])
)


(joint_summaries$without_depth$speed_parameterization$exy - 
    marginal_summaries$without_depth$speed_parameterization$speed$mean[1] * 
    marginal_summaries$without_depth$speed_parameterization$persistence$mean[1]
) /
  (
    sqrt(marginal_summaries$without_depth$speed_parameterization$speed$var[1]) *
      sqrt(marginal_summaries$without_depth$speed_parameterization$persistence$var[1])
  )
