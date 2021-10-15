targets::tar_load(sim_priors)
targets::tar_load(sim_obs)

library(sp)
library(crawl)
library(RANN)
library(bisque)
library(coda)
library(ggplot2)
library(ggthemes)

source('R/ctds_hanks/fit_hanks_imputed.R')
source('R/ctds_hanks/crawl_discretization_noproject.R')

obs = sim_obs[[1]]$obs

complete_imputations = crawl_discretization_noproject(
  xlocs = obs$states[,1],
  ylocs = obs$states[,2],
  times = obs$times,
  nimputed = 100,
  delta = .005
)

post_samples = fit_hanks_imputed(
  imputed_trajectories = complete_imputations$complete_trajectories, 
  niter = 1e4, cell_size = 1, 
  prior = function(betaAR, speed) {
    dnorm(x = betaAR, sd = 1e2, log = TRUE) + 
      dgamma(x = speed, shape = 1, rate = 1, log = TRUE)
  }
)

param_samples = do.call(rbind, lapply(post_samples, function(x) x$param_vec))

effectiveSize(mcmc(param_samples))

plot(mcmc(param_samples))
plot(mcmc(exp(param_samples[,'log_speed'])))

HPDinterval(mcmc(exp(param_samples[,'log_speed'])))

param_samples_thinned = param_samples[
  seq(from = 1, to = nrow(param_samples), length.out = 1e4),
]


plot(mcmc(param_samples_thinned))
plot(mcmc(exp(param_samples_thinned[,'log_speed'])))

HPDinterval(mcmc(exp(param_samples_thinned[,'log_speed'])))
mean(mcmc(exp(param_samples_thinned[,'log_speed'])))



ggplot(data.frame(param_samples_thinned), 
       aes(x = exp(log_speed), y = betaAR)) + 
  geom_density_2d(col = 'black') +
  theme_few() 
