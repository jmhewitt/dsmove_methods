exact_targets = list(
  
  # parameters and initial condition for a CTRW trajectory 
  tar_target(
    name = sim_rw_params, 
    command = list(
      t0 = 0,         # start time for simulation
      tf = 500,       # end time for simulation
      theta = 1,      # transition rate for simulation
      x0 = rep(0,2)   # location from which simulated trajectory begins
    )
  ),
  
  # realization of a CTRW trajectory
  tar_target(
    name = sim_rw,
    command = {
      # number of transitions in realization
      n = rpois(
        n = 1, 
        lambda = sim_rw_params$theta * 
          (sim_rw_params$tf - sim_rw_params$t0)
      )
      # transition times
      times = c(0, sort(runif(
        n = n, 
        min = sim_rw_params$t0, 
        max = sim_rw_params$tf
      )))
      # valid state transitions
      valid_tx = list(c(1,0), c(-1,0), c(0,1), c(0,-1))
      # states visited
      states = matrix(sim_rw_params$x0, nrow = 1)
      for(i in 1:n) {
        loc = states[nrow(states),]
        states = rbind(states, loc  + sample(x = valid_tx, size = 1)[[1]])
      }
      # package results
      list(
        states = states,
        times = times,
        durations = diff(times)
      )
    }
  ),
  
  # observe CTRW trajectory under increasing temporal resolution
  tar_target(
    name = sim_rw_obs,
    command = list(list(
      obs = ctds.observe(
        states = sim_rw$states,
        times = sim_rw$times,
        t.obs = seq(
          from = sim_rw$times[1],
          to = tail(sim_rw$times, 1) + sim_rw_obs_interval,
          by = sim_rw_obs_interval
        )
      ),
      obs_interval = sim_rw_obs_interval
    )),
    pattern = map(sim_rw_obs_interval)
  ),
  
  tar_target(sim_rw_obs_interval, c(5, 1, .5, .25)),
  
  # priors used for sensitivity study
  tar_target(rw_priors, data.frame(
    name = c('uninformative', 'semi_weak', 'weak', 'strong'),
    shape = c(.001, 1, 2, 100),
    rate = c(.001, 1, 2, 100)
  )),
  
  # posterior approximation parameters
  tar_target(mc_rw_params, list(
    niter = 1e3,                   # number of MC samples to draw
    nseq = seq(from = 0, to = 100) # finite support for imputed path lengths
  )),
  
  # MC samples for exact posterior approximations
  tar_target(
    name = rw_post_samples, 
    command = {
      # extract observations
      obs = sim_rw_obs[[1]]$obs
      m = nrow(obs$states) - 1
      # generate and package results
      list(list(
        post_samples = replicate(n = mc_rw_params$niter, expr = {
          # draw path lengths
          N_imputed = sapply(1:m, function(obs_ind) {
            dens = dN(N = mc_rw_params$nseq, s0 = obs$states[obs_ind,], 
                      sf = obs$states[obs_ind + 1,], t0 = obs$times[obs_ind], 
                      tf = obs$times[obs_ind + 1], shape = rw_priors$shape, 
                      rate = rw_priors$rate, log = F) 
            sample(x = mc_rw_params$nseq, size = 1, prob = dens)
          })
          # draw and return speed parameter
          rgamma(n = 1, shape = rw_priors$shape + sum(N_imputed), 
                 rate = rw_priors$rate + diff(range(obs$times)))
        }),
        priors = rw_priors,
        obs_interval = sim_rw_obs[[1]]$obs_interval
      ))
    },
    pattern = cross(map(rw_priors), sim_rw_obs),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient'
  ),
      
  # parameters for most-likely path posterior approximations
  tar_target(
    name = rw_post_mostlikely,
    command = {
      # extract observations
      obs = sim_rw_obs[[1]]$obs
      m = nrow(obs$states) - 1
      # determine least number of steps to connect observations
      N_likely = sapply(1:m, function(obs_ind) {
        dens = dN(N = mc_rw_params$nseq, s0 = obs$states[obs_ind,], 
                  sf = obs$states[obs_ind + 1,], t0 = obs$times[obs_ind], 
                  tf = obs$times[obs_ind + 1], shape = rw_priors$shape, 
                  rate = rw_priors$rate, log = T) 
        mc_rw_params$nseq[which.max(dens)]
        
      })
      # compute posterior (Gamma) parameters and package results
      list(list(
        parameters = c(shape = rw_priors$shape + sum(N_likely),
                       rate = rw_priors$rate + diff(range(obs$times))),
        priors = rw_priors,
        obs_interval = sim_rw_obs[[1]]$obs_interval
      ))
    },
    pattern = cross(map(rw_priors), sim_rw_obs),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient'
  ),
  
  # parameters for shortest-path posterior approximations
  tar_target(
    name = rw_post_shortest,
    command = {
      # extract observations
      obs = sim_rw_obs[[1]]$obs
      m = nrow(obs$states) - 1
      # determine least number of steps to connect observations
      N_shortest = sapply(1:m, function(obs_ind) {
        sum(abs(obs$states[obs_ind + 1,] - obs$states[obs_ind,]))
      })
      # compute posterior (Gamma) parameters and package results
      list(list(
        parameters = c(shape = rw_priors$shape + sum(N_shortest),
                       rate = rw_priors$rate + diff(range(obs$times))),
        priors = rw_priors,
        obs_interval = sim_rw_obs[[1]]$obs_interval
      ))
    },
    pattern = cross(map(rw_priors), sim_rw_obs),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient'
  ),
  
  # MC samples for exact posterior approximations
  tar_target(
    name = rw_post_hanks, 
    command = {
      # extract observations
      obs = sim_rw_obs[[1]]$obs
      m = nrow(obs$states) - 1
      # generate and package results
      list(list(
        post_samples = hanks_imputation(
          states = obs$states, times = obs$times, reps = 50, 
          samples_per_rep = 1e3, priors = rw_priors
        ),
        priors = rw_priors,
        obs_interval = sim_rw_obs[[1]]$obs_interval
      ))
    },
    pattern = cross(map(rw_priors), sim_rw_obs),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient'
  ),
  
  # gaussian approximation to posterior for exact posterior approximations
  tar_target(
    name = rw_post_dtmc, 
    command = {
      # extract observations
      obs = sim_rw_obs[[1]]$obs
      m = nrow(obs$states) - 1
      # generate and package results
      list(list(
        posterior = dtmc_approximation(
          states = obs$states, times = obs$times, delta = .125, 
          priors = rw_priors, niter = mc_rw_params$niter
        ),
        priors = rw_priors,
        obs_interval = sim_rw_obs[[1]]$obs_interval
      ))
    },
    pattern = cross(map(rw_priors), sim_rw_obs),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient'
  ),
  
  # batch execution plan for importance sampling
  tar_target(n_is_batches, 5),
  tar_target(importance_sample_batch, 1:n_is_batches),
  
  # MC samples for exact posterior approximations
  tar_target(
    name = rw_post_dtmc_samples, 
    command = {
      # extract observations
      obs = sim_rw_obs[[1]]$obs
      m = nrow(obs$states) - 1
      # generate and package results
      list(list(
        samples = dtmc_approximation_is(
          states = obs$states, times = obs$times, delta = .125, 
          priors = rw_priors, niter = mc_rw_params$niter / n_is_batches, 
          gapprox = rw_post_dtmc[[1]]$posterior
        ),
        priors = rw_priors,
        obs_interval = sim_rw_obs[[1]]$obs_interval,
        batch_id = importance_sample_batch
      ))
    },
    pattern = cross(
      map(cross(map(rw_priors), sim_rw_obs), rw_post_dtmc),
      importance_sample_batch
    ),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient'
  ),
  
  # MC samples for exact posterior approximations
  tar_target(
    name = rw_post_hanks_uninformative, 
    command = {
      # extract observations
      obs = sim_rw_obs[[1]]$obs
      m = nrow(obs$states) - 1
      # generate and package results
      list(list(
        post_samples = hanks_imputation(
          states = obs$states, times = obs$times, reps = 50, 
          samples_per_rep = 1e3, priors = rw_priors,
          hanks.priors = list(
            a = .001, b = .001, r = .001, q = .001
          )
        ),
        priors = rw_priors,
        obs_interval = sim_rw_obs[[1]]$obs_interval
      ))
    },
    pattern = cross(map(rw_priors), sim_rw_obs),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient'
  ),
  
  # parameters for posterior when path lengths are exactly known
  tar_target(
    name = rw_post_exact,
    command = {
      N_exact = length(sim_rw$durations)
      list(list(
        parameters = c(shape = rw_priors$shape + sum(N_exact),
                       rate = rw_priors$rate + diff(range(sim_rw$times))),
        priors = rw_priors
      ))
    },
    pattern = map(rw_priors),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient'
  ),
  
  tar_target(
    name = rw_post_summaries_shortest, 
    command = {
      # extract parameters from model fits
      df = do.call(rbind, lapply(rw_post_shortest, function(res) {
        data.frame(post.shape = res$parameters['shape'], 
                   post.rate = res$parameters['rate'],
                   prior = res$priors$name, 
                   obs.interval = res$obs_interval)
      }))
      # compute hpds
      hpds = apply(df, 1, function(r) {
        HDInterval::hdi(qgamma, credMass = .95, 
                        shape = as.numeric(r['post.shape']), 
                        rate = as.numeric(r['post.rate']))
      }) 
      # enrich extracted information with posterior summaries
      df$post.mean = df$post.shape / df$post.rate
      df$hpd.lwr = hpds['lower',]
      df$hpd.upr = hpds['upper',]
      
      ggplot(df, aes(x = jitter(obs.interval), y = post.mean, ymin = hpd.lwr, 
                     ymax = hpd.upr, col = prior)) +
        geom_hline(yintercept = sim_rw_params$theta, lty = 3) + 
        geom_pointrange() +
        theme_few() + 
        theme(panel.border = element_blank())
      
    }
  ),
  
  tar_target(
    name = rw_post_summaries_samples, 
    command = {
      # extract posterior summaries from MC samples
      df = do.call(rbind, lapply(rw_post_samples, function(res) {
        hpd = HPDinterval(mcmc(res$post_samples))
        data.frame(post.mean = mean(res$post_samples), 
                   hpd.lwr = hpd[,'lower'],
                   hpd.upr = hpd[,'upper'],
                   prior = res$priors$name, 
                   obs.interval = res$obs_interval)
      }))
      
      ggplot(df, aes(x = jitter(obs.interval), y = post.mean, ymin = hpd.lwr, 
                     ymax = hpd.upr)) +
        geom_hline(yintercept = sim_rw_params$theta, lty = 3) + 
        geom_pointrange() +
        facet_wrap(~prior, scales = 'free') + 
        theme_few() + 
        theme(panel.border = element_blank())
      
    }
  ),
  
  tar_target(
    name = rw_post_summaries_hanks, 
    command = {
      # extract posterior summaries from MC samples
      df = do.call(rbind, lapply(rw_post_hanks, function(res) {
        hpd = HPDinterval(mcmc(res$post_samples$samples))
        data.frame(post.mean = mean(res$post_samples$samples), 
                   hpd.lwr = hpd[,'lower'],
                   hpd.upr = hpd[,'upper'],
                   prior = res$priors$name, 
                   obs.interval = res$obs_interval)
      }))
      
      ggplot(df, aes(x = jitter(obs.interval), y = post.mean, ymin = hpd.lwr, 
                     ymax = hpd.upr)) +
        geom_hline(yintercept = sim_rw_params$theta, lty = 3) + 
        geom_pointrange() +
        facet_wrap(~prior, scales = 'free') + 
        theme_few() + 
        theme(panel.border = element_blank())
      
    }
  ),
  
  tar_target(
    name = rw_post_summaries_hanks_uninformative, 
    command = {
      # extract posterior summaries from MC samples
      df = do.call(rbind, lapply(rw_post_hanks_uninformative, function(res) {
        hpd = HPDinterval(mcmc(res$post_samples$samples))
        data.frame(post.mean = mean(res$post_samples$samples), 
                   hpd.lwr = hpd[,'lower'],
                   hpd.upr = hpd[,'upper'],
                   prior = res$priors$name, 
                   obs.interval = res$obs_interval)
      }))
      
      ggplot(df, aes(x = jitter(obs.interval), y = post.mean, ymin = hpd.lwr, 
                     ymax = hpd.upr)) +
        geom_hline(yintercept = sim_rw_params$theta, lty = 3) + 
        geom_pointrange() +
        facet_wrap(~prior, scales = 'free') + 
        theme_few() + 
        theme(panel.border = element_blank())
      
    }
  ),
  
  tar_target(
    name = rw_post_summaries_combined, 
    command = {
      
      # load DTMC approximation results
      rw_post_dtmc_raw = rw_post_dtmc
      # rw_post_dtmc_raw = readRDS('rw_post_dtmc.rds')
      
      # extract parameters from model fits
      df = rbind(
        do.call(rbind, lapply(rw_post_mostlikely, function(res) {
          data.frame(post.shape = res$parameters['shape'], 
                     post.rate = res$parameters['rate'],
                     prior = res$priors$name, 
                     obs.interval = res$obs_interval, 
                     method = 'Most likely path')
        })),
        do.call(rbind, lapply(rw_post_shortest, function(res) {
          data.frame(post.shape = res$parameters['shape'], 
                     post.rate = res$parameters['rate'],
                     prior = res$priors$name, 
                     obs.interval = res$obs_interval, 
                     method = 'Shortest path')
        })),
        do.call(rbind, lapply(rw_post_exact, function(res) {
          data.frame(post.shape = res$parameters['shape'], 
                     post.rate = res$parameters['rate'],
                     prior = res$priors$name, 
                     obs.interval = 0, 
                     method = 'Truth known')
        }))
      )
      # compute hpds
      hpds = apply(df, 1, function(r) {
        HDInterval::hdi(qgamma, credMass = .95, 
                        shape = as.numeric(r['post.shape']), 
                        rate = as.numeric(r['post.rate']))
      }) 
      # enrich extracted information with posterior summaries
      df$post.mean = df$post.shape / df$post.rate
      df$hpd.lwr = hpds['lower',]
      df$hpd.upr = hpds['upper',]
      
      # remove terms that are not needed for full plotting
      df$post.shape = NULL
      df$post.rate = NULL
      
      df = rbind(
        df,
        do.call(rbind, lapply(rw_post_hanks, function(res) {
          hpd = HPDinterval(mcmc(res$post_samples$samples))
          data.frame(post.mean = mean(res$post_samples$samples), 
                   hpd.lwr = hpd[,'lower'],
                   hpd.upr = hpd[,'upper'],
                   prior = res$priors$name, 
                   obs.interval = res$obs_interval, 
                   method = 'Hanks')
        })),
        do.call(rbind, lapply(rw_post_hanks_uninformative, function(res) {
          hpd = HPDinterval(mcmc(res$post_samples$samples))
          data.frame(post.mean = mean(res$post_samples$samples), 
                     hpd.lwr = hpd[,'lower'],
                     hpd.upr = hpd[,'upper'],
                     prior = res$priors$name, 
                     obs.interval = res$obs_interval, 
                     method = 'Hanks-Uninformative')
        })),
        do.call(rbind, lapply(rw_post_samples, function(res) {
          hpd = HPDinterval(mcmc(res$post_samples))
          data.frame(post.mean = mean(res$post_samples), 
                     hpd.lwr = hpd[,'lower'],
                     hpd.upr = hpd[,'upper'],
                     prior = res$priors$name, 
                     obs.interval = res$obs_interval,
                     method = 'Model-based imputation')
        })),
        do.call(rbind, lapply(rw_post_dtmc_raw, function(res) {
          data.frame(
            prior = res$priors$name,
            obs.interval = res$obs_interval,
            method = 'DTMC approximation', 
            post.mean = res$posterior$post_mean,
            hpd.lwr = res$posterior$lower,
            hpd.upr = res$posterior$upper
          )
        }))
      )
      
      df = df %>% filter(
        method %in% c('DTMC approximation', 'Hanks', 'Model-based imputation',
                      'Truth known')
      )
      
      
      pl = ggplot(df %>% dplyr::filter(obs.interval < 3,
                                  method != 'Model-based imputation'), 
             aes(x = obs.interval, y = post.mean, ymin = hpd.lwr, 
                 ymax = hpd.upr, col = method)) +
        geom_hline(yintercept = sim_rw_params$theta, lty = 3) + 
        geom_line(alpha  = .4) + 
        geom_pointrange() +
        xlab('Time between observations') + 
        ylab(expression(theta)) + 
        scale_color_brewer('Imputation method', 
                           type = 'qual', palette = 'Dark2') + 
        facet_wrap(~prior, scales = 'free', ncol = 1) + 
        theme_few() + 
        ggtitle('Posterior means and HPDs by prior')
      
      pl2 = ggplot(df %>% dplyr::filter(obs.interval < 3), 
                  aes(x = obs.interval, y = post.mean, ymin = hpd.lwr, 
                      ymax = hpd.upr, col = method)) +
        geom_hline(yintercept = sim_rw_params$theta, lty = 3) + 
        geom_line(alpha  = .4) + 
        geom_pointrange() +
        xlab('Time between observations') + 
        ylab(expression(theta)) + 
        scale_color_brewer('Imputation method', 
                           type = 'qual', palette = 'Dark2') + 
        facet_wrap(~prior, scales = 'free', ncol = 1) + 
        theme_few() + 
        ggtitle('Posterior means and HPDs by prior')
      
      pl3 = ggplot(df, 
                   aes(x = obs.interval, y = post.mean, ymin = hpd.lwr, 
                       ymax = hpd.upr, col = method)) +
        geom_hline(yintercept = sim_rw_params$theta, lty = 3) + 
        geom_line(alpha  = .4) + 
        geom_pointrange() +
        xlab('Time between observations') + 
        ylab(expression(theta)) + 
        scale_color_brewer('Imputation method', 
                           type = 'qual', palette = 'Dark2') + 
        facet_wrap(~prior, scales = 'free', ncol = 1) + 
        theme_few() + 
        ggtitle('Posterior means and HPDs by prior')
    
      ggsave(pl, filename = 'posteriors_without_model.png',
             width = 12, height = 8)  
      ggsave(pl2, filename = 'posteriors_with_model.png',
             width = 12, height = 8)
      ggsave(pl3, filename = 'posteriors_with_model_longer_x.png',
             width = 12, height = 8)  
      
    }
  )
  
)
