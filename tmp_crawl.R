targets::tar_load(whale_pkg)
targets::tar_load(whale_domain)
targets::tar_load(brs_crs)

library(sp)
library(crawl)
library(RANN)
library(bisque)
library(coda)
library(ggplot2)
library(ggthemes)

whale_pkg = whale_pkg$whale_pkg_5bc05e7c

source('R/ctds_hanks/fit_hanks_imputed.R')
source('R/ctds_hanks/crawl_discretization.R')

pmmap = readRDS('post_marginal_map.rds')

pred_times = strptime(
  x = gsub(pattern = 'Prediction time: ', replacement = '', 
           x = names(pmmap$marginal_location)),
  format = '%Y-%m-%d %H:%M:%S', 
  tz = 'GMT'
)

complete_imputations = crawl_discretization(
  xlocs = whale_pkg$data$loc$mapped_lon,
  ylocs = whale_pkg$data$loc$mapped_lat, 
  semi_majors = whale_pkg$data$loc$semi_major, 
  semi_minors = whale_pkg$data$loc$semi_minor, 
  orientations = whale_pkg$data$loc$orientation, 
  times = whale_pkg$data$loc$time,
  domain = whale_domain, 
  nimputed = 100,
  # CRS used to project data for analysis
  crawl_proj = brs_crs,
  obs_proj = '+proj=longlat +datum=WGS84 +no_defs',
  pred_times = pred_times, 
  delta = 5, 
  discrete_coords = whale_pkg$grid$coords, 
  discrete_lons = whale_pkg$grid$lons,
  discrete_lats = whale_pkg$grid$lats
)

saveRDS(complete_imputations$imputations_at_times, 
        file = 'crawl_discretized_post_locations100.rds')

post_samples = fit_hanks_imputed(
  imputed_trajectories = complete_imputations$complete_trajectories, 
  niter = 1e4, cell_size = 623.7801, 
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

saveRDS(param_samples_thinned, 'ctds_crawl_posterior_samples_zc95.rds')

plot(mcmc(param_samples_thinned))
plot(mcmc(exp(param_samples_thinned[,'log_speed'])))

HPDinterval(mcmc(exp(param_samples_thinned[,'log_speed'])))
mean(mcmc(exp(param_samples_thinned[,'log_speed'])))



ggplot(data.frame(param_samples_thinned), 
       aes(x = exp(log_speed), y = betaAR)) + 
  geom_density_2d(col = 'black') +
  theme_few() 
