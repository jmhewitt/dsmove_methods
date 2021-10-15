library(metR)
library(ggplot2)
library(ggthemes)
library(viridis)
library(dplyr)
library(tidyr)

source('R/utils/log_add.R')
source('R/utils/distribution_parameters.R')

res = readRDS('whale_ll_subsets.rds')

priors = list(
  # speed = gamma.param(mu = .5, sd = .5),
  speed = gamma.param(mu = 1, sd = 1),
  # speed = c(shape = 2, rate = 1)
  betaAR = c(mean = 0, sd = 1e2),
  beta0_hanks = c(mean = 0, sd = 1e2)
  # beta0_hanks = c(mean = -10, sd = 1)
)

# one way to potentially help set the speed prior
qgamma(.99, shape = priors$speed['shape'], rate = priors$speed['rate'])

#
# back-transform hanks model parameterizations
#

# append induced hanks location parameter
res = res %>% mutate(
  beta0_hanks = log(speed) - log(cell_size) - 
    log(2 + exp(-betaAR) + exp(betaAR))
)

res = res %>% mutate(
  log_prior_hanks = 
    dnorm(x = beta0_hanks, mean = priors$beta0_hanks['mean'], 
          sd = priors$beta0_hanks['sd'], log = TRUE) + 
    dnorm(x = betaAR, mean = priors$betaAR['mean'], 
          sd = priors$betaAR['sd'], log = TRUE) - 
    # jacobian 
    log(speed),
  lp_hanks = ll + log_prior_hanks,
  lp = ll + 
    # prior
    dgamma(x = speed, shape = priors$speed['shape'], 
           rate = priors$speed['rate'], log = TRUE) + 
    dnorm(x = betaAR, mean = priors$betaAR['mean'], 
          sd = priors$betaAR['sd'], log = TRUE),
  log_prior_speed =
    dgamma(x = speed, shape = priors$speed['shape'], 
           rate = priors$speed['rate'], log = TRUE) + 
    dnorm(x = betaAR, mean = priors$betaAR['mean'], 
          sd = priors$betaAR['sd'], log = TRUE)
)

dxy = log(diff(unique(res$speed))[1]) + log(diff(unique(res$betaAR))[1])

log_consts = res %>% 
  group_by(subset) %>% 
  summarise(
    lC = log_sum(lp + dxy),
    lC_hanks = log_sum(lp_hanks + dxy),
    lC_hanks_prior = log_sum(log_prior_hanks + dxy),
    lC_speed_prior = log_sum(log_prior_speed + dxy)
  )

df = res %>% 
  left_join(log_consts, by = 'subset') %>% 
  mutate(lp = lp - lC, lp_hanks = lp_hanks - lC_hanks,
         log_prior_hanks = log_prior_hanks - lC_hanks_prior,
         log_prior_speed = log_prior_speed - lC_speed_prior) %>% 
  pivot_longer(cols = c(lp, lp_hanks, log_prior_hanks, log_prior_speed), 
               names_to = 'lp_type', values_to = 'lp_value') %>% 
  mutate(lp_type = factor(lp_type),
         subset = factor(subset))

pl = ggplot(df %>% filter(lp_type %in% c('lp', 'lp_hanks')),
            aes(x = speed, y = betaAR, fill = exp(lp_value), 
                z = exp(lp_value))) + 
  # geom_tile(height = 0.3157895, width = 0.3139539) +
  # geom_raster() + 
  geom_contour2(col = 'grey60') + 
  geom_text_contour(col = 'grey20') + 
  facet_grid(subset~lp_type) + 
  theme_few() + 
  scale_fill_viridis(direction = -1) + 
  ggtitle('Posterior density surface comparisons')

pl

ggsave(pl, filename = 'posterior_surface_comparisons.pdf')



dx = log(diff(unique(res$betaAR))[1])
dx2 = log(diff(unique(res$speed))[2])

speed_marginal = df %>%
  group_by(speed, subset, lp_type) %>% 
  summarise(
    lp = log_sum(
      lp_value + dx
    ) 
  ) %>% 
  ungroup()

persistence_marginal = df %>%
  group_by(betaAR, subset, lp_type) %>% 
  summarise(
    lp = log_sum(
      lp_value + dx2
    ) 
  ) %>% 
  ungroup()



dx = diff(unique(speed_marginal$speed))[2]

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