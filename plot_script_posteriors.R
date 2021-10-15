# additional posterior comparisons

library(raster)
library(sp)
library(metR)
library(ggplot2)
library(ggthemes)
library(ggnewscale)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggpubr)

crawl_samples = readRDS('ctds_crawl_posterior_samples_zc95.rds')

targets::tar_load(whale_post_marginal_parameters)

#
# compare AID vs time discretization
#

# posterior for persistence
pl1 = ggplot(
  data.frame(x = crawl_samples[,'betaAR'], Density = 'CTCRW as AID posterior'),
  aes(x = x, col = Density)
) + 
  # MCMC samples
  stat_density(geom = 'line') + 
  # discretization approximation
  geom_line(
    data = rbind(
      cbind(
        whale_post_marginal_parameters[[
          'with_depth']][['speed_parameterization']][['persistence']] %>% 
          filter(density == 'posterior'),
        Density = 'Time discretization posterior (Surface + Depth)'
      ),
      cbind(
        whale_post_marginal_parameters[[
          'with_depth']][['speed_parameterization']][['persistence']] %>% 
          filter(density == 'prior'),
        Density = 'Prior'
      )
    ),
    mapping = aes(x = x, y = exp(ld), col = Density),
    inherit.aes = FALSE
  ) + 
  scale_color_brewer(type = 'qual', palette = 'Dark2') + 
  xlab(expression(rho)) + 
  ylab(expression(f(rho))) + 
  theme_few()

# posterior for speed
pl2 = ggplot(
  data.frame(x = exp(crawl_samples[,'log_speed']), 
             Density = 'CTCRW as AID posterior'),
  aes(x = x, col = Density)
) + 
  # MCMC samples
  stat_density(geom = 'line') + 
  # discretization approximation
  geom_line(
    data = rbind(
      cbind(
        whale_post_marginal_parameters[[
          'with_depth']][['speed_parameterization']][['speed']] %>% 
          filter(density == 'posterior'),
        Density = 'Time discretization posterior (Surface + Depth)'
      ),
      cbind(
        whale_post_marginal_parameters[[
          'with_depth']][['speed_parameterization']][['speed']] %>% 
          filter(density == 'prior'),
        Density = 'Prior'
      )
    ),
    mapping = aes(x = x, y = exp(ld), col = Density),
    inherit.aes = FALSE
  ) + 
  scale_color_brewer(type = 'qual', palette = 'Dark2') + 
  xlab('Speed (m/s)') + 
  ylab(expression(f(Speed))) + 
  theme_few()

# combined plot, as figure
ggsave(
  ggarrange(pl1, pl2, common.legend = TRUE, legend = 'bottom'),
  file = 'posterior_comparison_with_depth.png', dpi = 'print',
  width = 8, height = 5
)


# joint posterior
targets::tar_load(whale_post_joint_parameters)

pl = cbind(
  whale_post_joint_parameters$with_depth$speed_parameterization,
  Method = 'Discretization'
  ) %>% 
  filter(lpost > -50) %>% 
  mutate(lpost = exp(lpost)) %>% 
  ggplot(aes(x = speed, y = betaAR, fill = lpost,  z = lpost, 
             col = Method)) + 
  # geom_raster() + 
  geom_contour() + 
  geom_density_2d(
    mapping = aes(x = speed, y = betaAR, col = Method),
    data = data.frame(crawl_samples, Method = 'CRAWL AID') %>% 
      mutate(speed = exp(log_speed)),
    inherit.aes = FALSE
  ) + 
  theme_few() + 
  scale_color_brewer(type = 'qual', palette = 'Dark2') + 
  # ggtitle('Joint posterior associated with Fig. 2')
  # ggtitle('Posterior comparison') + 
  xlab('Speed (m/s)') + 
  ylab('Persistence')

pl


ggsave(pl, filename = 'joint_posterior_comparison_with_depth.png',
       dpi = 'print')

#
# compare time discretization with and without depth data
#

pl1 = ggplot(
  rbind(
    cbind(
      whale_post_marginal_parameters$with_depth$speed_parameterization$persistence,
      Dataset = 'Surface + Depth'
    ),
    cbind(
      whale_post_marginal_parameters$without_depth$speed_parameterization$persistence,
      Dataset = 'Surface only'
    )
  ) %>% filter(density == 'posterior'),
  aes(x = x, y = exp(ld), col = Dataset)
) +
  geom_line() + 
  xlab(expression(rho)) + 
  ylab(expression(f(rho))) + 
  theme_few()
  
pl2 = ggplot(
  rbind(
    cbind(
      whale_post_marginal_parameters$with_depth$speed_parameterization$speed,
      Dataset = 'Surface + Depth'
    ),
    cbind(
      whale_post_marginal_parameters$without_depth$speed_parameterization$speed,
      Dataset = 'Surface only'
    )
  ) %>% filter(density == 'posterior'),
  aes(x = x, y = exp(ld), col = Dataset)
) +
  geom_line() + 
  xlab('Speed (m/s)') + 
  ylab(expression(f(Speed))) + 
  theme_few()

# combined plot, as figure
ggsave(
  ggarrange(pl1, pl2, common.legend = TRUE, legend = 'bottom'),
  file = 'posterior_comparison_by_dataset.png', dpi = 'print',
  width = 8, height = 5
)


#
# parameterization comparison
#

pl1 = ggplot(
  rbind(
    cbind(
      whale_post_marginal_parameters$with_depth$speed_parameterization$persistence,
      Parameterization = 'Direct'
    ),
    cbind(
      whale_post_marginal_parameters$with_depth$beta0_parameterization$persistence,
      Parameterization = 'Original'
    )
  ) %>% filter(density == 'posterior'),
  aes(x = x, y = exp(ld), col = Parameterization)
) +
  geom_line() + 
  xlab(expression(rho)) + 
  ylab(expression(f(rho))) + 
  theme_few()

pl2 = ggplot(
  rbind(
    cbind(
      whale_post_marginal_parameters$with_depth$speed_parameterization$speed,
      Parameterization = 'Direct'
    ),
    cbind(
      whale_post_marginal_parameters$with_depth$beta0_parameterization$speed,
      Parameterization = 'Original'
    )
  ) %>% filter(density == 'posterior'),
  aes(x = x, y = exp(ld), col = Parameterization)
) +
  geom_line() + 
  xlab('Speed (m/s)') + 
  ylab(expression(f(Speed))) + 
  theme_few()

# combined plot, as figure
ggsave(
  ggarrange(pl1, pl2, common.legend = TRUE, legend = 'bottom'),
  file = 'posterior_comparison_by_parameteriation.png', dpi = 'print',
  width = 8, height = 5
)
s