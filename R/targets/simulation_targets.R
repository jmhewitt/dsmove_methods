simulation_targets = list(
  
  #
  # simple simulation
  #
  
  # prior distributions
  tar_target(
    name = sim_priors, 
    command = list(
      theta = c('shape' = 1, 'rate' = 1),
      betaAR = c('mean' = 0, 'sd' = 1e2)
    )
  ),
  
  # AR parameters to run in simulation
  tar_target(sim_betaAR, c(0,1)),
  
  # parameters to generate basic CTDS trajectories
  tar_target(
    name = sim_params, 
    command = list(
      t0 = 0,               # start time for simulation
      tf = 500,             # end time for simulation
      betaAR = sim_betaAR,  # directional persistence parameter
      beta = 0,             # log of transition rate
      dims = rep(1e3, 2),   # support for spatial domain
      x0 = rep(500, 2),     # location from which simulated trajectory begins
      x0.last = c(499, 500) # auto regressive loc. to set trajectory direction
    ),
    pattern = map(sim_betaAR)
  ),
  
  # realization of basic CTDS trajectory
  tar_target(
    name = sim_simple,
    command = {
      set.seed(2021)
      list(list(
        sim = ctds.quicksim(
          dims = sim_params$dims,
          beta_loc = sim_params$beta,
          beta_ar = sim_params$betaAR,
          v0 = sim_params$x0,
          t0 = sim_params$t0,
          tf = sim_params$tf,
          max.steps = 1e3,
          v0.last = sim_params$x0.last
        ),
        params = sim_params
      ))
    },
    pattern = map(sim_params)
  ),
  
  tar_target(obs_interval, c(2, 1, .5, .25)),
  
  # observe CTDS trajectory under increasing temporal resolution
  tar_target(
    name = sim_obs,
    command = list(list(
      obs = ctds.observe(
        states = sim_simple[[1]]$sim$states[-1,],
        times = sim_simple[[1]]$sim$times,
        t.obs = seq(
          from = sim_simple[[1]]$sim$times[1],
          to = tail(sim_simple[[1]]$sim$times, 1) + obs_interval,
          by = obs_interval
        )
      ),
      params = sim_simple[[1]]$params,
      obs_interval = obs_interval
    )),
    pattern = cross(obs_interval, sim_simple)
  ),

  tar_target(rep_batch, 1:10),
  
  tar_target(
    name = sim_fits_hanks,
    command = {

      post_samples = fit_hanks(
        params = list(beta_ar = 0, beta_loc = 0),
        niter = 1e4, priors = sim_priors,
        states = sim_obs[[1]]$obs$states,
        times = sim_obs[[1]]$obs$times,
        dims = sim_obs[[1]]$params$dims,
        reps = 100
      )

      r = list(list(
        samples = post_samples, 
        rep = rep_batch, 
        sim_obs = sim_obs[[1]]
      ))
      
      f = file.path('output', 'simulation')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '.rds', sep = ''))
      saveRDS(r, file = f)
      
      f
    },
    pattern = cross(sim_obs, rep_batch),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(
    name = sim_fits_hanks_univariate,
    command = {
      
      post_samples = fit_hanks(
        params = list(beta_ar = 0, beta_loc = 0),
        niter = 1e4, priors = sim_priors,
        states = sim_obs[[1]]$obs$states,
        times = sim_obs[[1]]$obs$times,
        dims = sim_obs[[1]]$params$dims,
        reps = 10,
        univariate = TRUE
      )
      
      r = list(list(
        samples = post_samples, 
        rep = rep_batch, 
        sim_obs = sim_obs[[1]]
      ))
      
      f = file.path('output', 'simulation')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '.rds', sep = ''))
      saveRDS(r, file = f)
      
      f
    },
    pattern = cross(sim_obs, rep_batch),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(
    name = sim_fit_dtmc_gapprox, 
    command = {
      
      message(paste(
        tar_name(), ': ', 
        ' betaAR: ', sim_obs[[1]]$params$betaAR,
        ' obs_ind:', sim_obs[[1]]$obs_interval,
        sep = ''
      ))
      
      gapprox = dtmc_ar_approximation(
        states = sim_obs[[1]]$obs$states, 
        times = sim_obs[[1]]$obs$times, 
        delta = .125/2,
        priors = sim_priors, 
        dims = sim_obs[[1]]$params$dims
      )
      
      list(list(
        gapprox = gapprox,
        sim_obs = sim_obs[[1]]
      ))
    },
    pattern = map(sim_obs),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(delta_seq, c(.125/2, .125, .25)),
  
  tar_target(
    name = sim_fit_dtmc_gapprox_delta_sensitivity, 
    command = {
      
      message(paste(
        tar_name(), ': ', 
        ' betaAR: ', sim_obs[[1]]$params$betaAR,
        ' obs_ind:', sim_obs[[1]]$obs_interval,
        sep = ''
      ))
      
      if(sim_obs[[1]]$obs_interval != 1) {
        return(NULL)
      }
      
      gapprox = dtmc_ar_approximation(
        states = sim_obs[[1]]$obs$states, 
        times = sim_obs[[1]]$obs$times, 
        delta = delta_seq,
        priors = sim_priors, 
        dims = sim_obs[[1]]$params$dims
      )
      
      list(list(
        gapprox = gapprox,
        sim_obs = sim_obs[[1]],
        delta = delta_seq
      ))
    },
    pattern = cross(delta_seq, map(sim_obs)),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(
    name = simulation_hanks_summaries,
    command = {
      
      sim_fits_hanks = readRDS('sim_fits_hanks.rds')

      burn = 1:5e3
      
      # extract identifiers for each batch output
      batch_ids = do.call(
        rbind,  lapply(1:length(sim_fits_hanks), function(ind) {
        res = sim_fits_hanks[[ind]]
        data.frame(batch = ind, 
          obs_interval = res$sim_obs$obs_interval,
          beta = res$sim_obs$params$beta,
          betaAR = res$sim_obs$params$betaAR)
      }))
      
      # determine unique groups
      groups = unique(batch_ids[, c('obs_interval', 'beta', 'betaAR')])
      groups = cbind(groups, gid = 1:nrow(groups))
      
      # associate group ids with the batch numbers
      batch_map = batch_ids %>% 
        dplyr::left_join(groups, by = c('obs_interval', 'beta', 'betaAR')) %>% 
        dplyr::select(batch, gid)
      
      # extract posterior summaries
      summaries = do.call(rbind, lapply(groups$gid, function(g) {
        # extract batches associated with group
        res_subset = sim_fits_hanks[
          batch_map %>% 
            dplyr::filter(gid == g) %>% 
            dplyr::select(batch) %>% 
            unlist()
        ]
        # aggregate posterior samples
        samples = do.call(rbind, lapply(res_subset, function(batches) {
          do.call(rbind, lapply(batches$samples$samples, function(res) {
            res$param_vec[-burn,]
          }))
        }))
        # transform parameter
        samples = cbind(samples, theta = exp(samples[,'log_theta']))
        # hpds 
        hpds = HPDinterval(mcmc(samples))
        # package results
        data.frame(
          param = colnames(samples),
          mean = colMeans(samples),
          lwr = hpds[,'lower'],
          upr = hpds[,'upper'],
          truth = c(unlist(res_subset[[1]]$sim_obs$params['betaAR']),
                    unlist(res_subset[[1]]$sim_obs$params['beta']),
                    exp(unlist(res_subset[[1]]$sim_obs$params['beta']))),
          method = 'Hanks',
          obs.interval = res_subset[[1]]$sim_obs$obs_interval,
          scenario = paste(
            'theta = ', exp(res_subset[[1]]$sim_obs$params$beta),
            ' betaAR = ', res_subset[[1]]$sim_obs$params$betaAR,
            sep = ''
          )
        )
      }))
      
      summaries
    }
  ),
  
  tar_target(
    name = simulation_hanks_univariate_summaries,
    command = {
      
      # AID MCMC sample files
      sample.files = dir(
        path = file.path('output', 'simulation'), 
        pattern = 'sim_fits_hanks_univariate_', 
        full.names = TRUE
      )
      
      samples.aggregate = do.call(rbind, lapply(sample.files, function(f) {
        x = readRDS(f)[[1]]
        burn = 1:5e3
        thin = 5
        samples = cbind(
          log_theta = x$samples$samples[[1]]$param_vec[-burn,'log_theta'],
          theta = exp(x$samples$samples[[1]]$param_vec[-burn,'log_theta'])
        )
        data.frame(
          log_theta = samples[seq(from = 1, to = nrow(samples), by = thin), 'log_theta'],
          theta = samples[seq(from = 1, to = nrow(samples), by = thin), 'theta'],
          rep = x$rep,
          obs_interval = x$sim_obs$obs_interval,
          beta = x$sim_obs$params$beta,
          betaAR = x$sim_obs$params$betaAR
        )
      }))
      
      # determine unique groups
      groups = unique(samples.aggregate[, c('obs_interval', 'beta', 'betaAR')])
      groups = cbind(groups, gid = 1:nrow(groups))

      # extract posterior summaries
      summaries = do.call(rbind, lapply(groups$gid, function(g) {
        # aggregate posterior samples associated with group
        res_subset = samples.aggregate %>% 
          filter(
            obs_interval == groups[g, 'obs_interval'],
            beta == groups[g, 'beta'],
            betaAR == groups[g, 'betaAR']
          )
        
        # extract parameter samples
        samples = res_subset[,c('log_theta','theta')]
        # hpds
        hpds = HPDinterval(mcmc(samples))
        # package results
        data.frame(
          param = c('log_theta', 'theta'),
          mean = colMeans(samples),
          lwr = hpds[,'lower'],
          upr = hpds[,'upper'],
          truth = c(
            unlist(res_subset$beta[1]),
            exp(unlist(res_subset$beta[1]))
          ),
          method = 'Hanks',
          obs.interval = res_subset$obs_interval[1],
          scenario = paste(
            'theta = ', exp(res_subset$beta[1]),
            ' betaAR = ', res_subset$betaAR[1],
            sep = ''
          )
        )
      })) %>% filter(scenario == 'theta = 1 betaAR = 0')
      
      summaries
    }
  ),
  
  tar_target(
    name = simulation_summaries_combined, 
    command = {
      
      sim_fit_dtmc_gapprox = readRDS('sim_fit_dtmc_gapprox.rds')
      
      # extract posterior summaries from model fits
      df = do.call(rbind, lapply(sim_fit_dtmc_gapprox, function(res) {
        rbind(
          data.frame(
            mean = res$gapprox$post_mean['theta'],
            lwr = res$gapprox$hpds['theta', 'lower'],
            upr = res$gapprox$hpds['theta', 'upper'],
            truth = exp(res$sim_obs$params$beta),
            obs.interval = res$sim_obs$obs_interval,
            scenario = paste(
              'theta = ', exp(res$sim_obs$params$beta),
              ' betaAR = ', res$sim_obs$params$betaAR,
              sep = ''
            ),
            method = 'DTMC',
            param = 'theta'
          ),
          data.frame(
            mean = res$gapprox$post_mean['betaAR'],
            lwr = res$gapprox$hpds['betaAR', 'lower'],
            upr = res$gapprox$hpds['betaAR', 'upper'],
            truth = res$sim_obs$params$betaAR,
            obs.interval = res$sim_obs$obs_interval,
            scenario = paste(
              'theta = ', exp(res$sim_obs$params$beta),
              ' betaAR = ', res$sim_obs$params$betaAR,
              sep = ''
            ),
            method = 'DTMC',
            param = 'betaAR'
          )
        )
      }))
      
      # append hanks results
      df = rbind(df, simulation_hanks_summaries)
      
      # 
      plot_eps = .02
      
      pl = ggplot(df, aes(x = obs.interval - (method == 'Hanks') * plot_eps, 
                          y = mean, ymin = lwr, ymax = upr, col = method)) + 
        facet_wrap(~scenario*param, nrow = 2, scales = 'free') + 
        geom_hline(mapping = aes(yintercept = truth), lty = 3) + 
        scale_x_continuous('Time between observations', 
                           breaks = sort(unique(df$obs.interval))) + 
        geom_pointrange() + 
        scale_color_brewer(type = 'qual', palette = 'Dark2') + 
        theme_few()
      
      pl
      
      browser()
      ggsave(pl, filename = 'posteriors_2param.pdf', width = 8, height = 5)
      
      
      df2 = df %>% 
        mutate(method = gsub('DTMC', 'Discretization', method),
               method = gsub('Hanks', 'Imputation', method),
               param = gsub('betaAR', 'Autocorrelation', param),
               param = gsub('theta', 'Speed', param)) %>% 
        filter(scenario == 'theta = 1 betaAR = 1')
      
      
      pl = ggplot(df2, aes(x = obs.interval - (method == 'Hanks') * plot_eps, 
                          y = mean, ymin = lwr, ymax = upr, col = method)) + 
        facet_wrap(~param, nrow = 2, scales = 'free') + 
        geom_hline(mapping = aes(yintercept = truth), lty = 3) + 
        scale_x_continuous('Time between observations', 
                           breaks = sort(unique(df$obs.interval))) + 
        geom_pointrange() + 
        scale_color_brewer('Approx. type', type = 'qual', palette = 'Dark2') + 
        theme_few() + 
        ylab('Est. Value')
      
      pl
      
      ggsave(pl, filename = 'posteriors_2param_simple.png', width = 8.2, 
             height = 6.5)
      
    }
  ),
  
  tar_target(
    name = sensitivity_summaries_combined, 
    command = {
      
      sim_fit_dtmc_gapprox_delta_sensitivity = lapply(
        grep(pattern = 'sim_fit_dtmc_gapprox_delta_sensitivity', 
             x = tar_objects(), 
             value = TRUE),
        function(x) { 
          eval(substitute(tar_read(x), env = list(x=x)))[[1]]
        }
      )
      
      sim_fit_dtmc_gapprox_delta_sensitivity = 
        sim_fit_dtmc_gapprox_delta_sensitivity[
          !sapply(sim_fit_dtmc_gapprox_delta_sensitivity, is.null)
        ]
      
      # extract posterior summaries from model fits
      df = do.call(rbind, lapply(sim_fit_dtmc_gapprox_delta_sensitivity, function(res) {
        rbind(
          data.frame(
            mean = res$gapprox$post_mean['theta'],
            lwr = res$gapprox$hpds['theta', 'lower'],
            upr = res$gapprox$hpds['theta', 'upper'],
            truth = exp(res$sim_obs$params$beta),
            obs.interval = res$sim_obs$obs_interval,
            scenario = paste(
              'theta = ', exp(res$sim_obs$params$beta),
              ' betaAR = ', res$sim_obs$params$betaAR,
              sep = ''
            ),
            delta = res$delta,
            method = 'DTMC',
            param = 'theta'
          ),
          data.frame(
            mean = res$gapprox$post_mean['betaAR'],
            lwr = res$gapprox$hpds['betaAR', 'lower'],
            upr = res$gapprox$hpds['betaAR', 'upper'],
            truth = res$sim_obs$params$betaAR,
            obs.interval = res$sim_obs$obs_interval,
            scenario = paste(
              'theta = ', exp(res$sim_obs$params$beta),
              ' betaAR = ', res$sim_obs$params$betaAR,
              sep = ''
            ),
            delta = res$delta,
            method = 'DTMC',
            param = 'betaAR'
          )
        )
      }))
      
      
      # 
      plot_eps = .02
      
      pl = ggplot(df, aes(x = delta, 
                          y = mean, ymin = lwr, ymax = upr, )) + 
        facet_wrap(~scenario*param, nrow = 2, scales = 'free') + 
        geom_hline(mapping = aes(yintercept = truth), lty = 3) + 
        scale_x_continuous(expression(delta)) +
        geom_pointrange() + 
        scale_color_brewer(type = 'qual', palette = 'Dark2') + 
        theme_few()
      
      pl
      ggsave(pl, filename = 'posteriors_2param.pdf', width = 8, height = 5)
    }
  ),
  
  tar_target(
    name = sim_marg_dtmc_approx, 
    command = {
      
      delta = .125/2
      
      gapprox = dtmc_ar_marginal_filtering(
        states = sim_obs[[1]]$obs$states, 
        times = sim_obs[[1]]$obs$times, 
        delta = delta,
        dims = sim_obs[[1]]$params$dims,
        pred_times = sim_obs[[1]]$obs$times[5] + delta * 0:3
      )
      
      list(list(
        gapprox = gapprox,
        sim_obs = sim_obs[[1]]
      ))
    },
    pattern = map(sim_obs),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(univariate_options, c(FALSE, TRUE)),
  
  tar_target(
    name = sim_dtmc_mcmc,
    command = {
      
      # unwrap input
      obs = sim_obs[[1]]
      
      # construct and run mcmc sampler
      res = fit_dtmc_mcmc(
        niter = 1e3,
        delta = .025,
        priors = sim_priors,
        t0 = obs$params$t0,
        tf = obs$params$tf,
        dims = obs$params$dims,
        states = obs$obs$states,
        times = obs$obs$times,
        univariate = univariate_options
      )
      
      res$obs = obs
      
      # save results to disk
      f = file.path('output', 'simulation')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '.rds', sep = ''))
      saveRDS(res, file = f)
      
      # return results
      f
    },
    pattern = cross(map(sim_obs), univariate_options),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(
    name = sim_dtmc_mcmc_catchup,
    command = {
      
      # unwrap input for this specific simulated dataset
      obs = sim_obs[[1]]
      
      # construct and run mcmc sampler
      res = fit_dtmc_mcmc(
        niter = 1e3,
        delta = .1,
        priors = sim_priors,
        t0 = obs$params$t0,
        tf = obs$params$tf,
        dims = obs$params$dims,
        states = obs$obs$states,
        times = obs$obs$times,
        univariate = univariate_options
      )
      
      res$obs = obs
      
      # save results to disk
      f = file.path('output', 'simulation')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '.rds', sep = ''))
      saveRDS(res, file = f)
      
      # return results
      f
    },
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  ),
  
  tar_target(
    name = sim_dtmc_mcmc_order,
    command = {
      
      # unwrap input
      obs = sim_obs[[1]]
      
      list(obs)
    },
    pattern = cross(map(sim_obs), univariate_options),
  ),
  
  tar_target(
    name = simulation_results_combined_univariate, 
    command = {
      
      # DTMC MCMC sample files
      sample.files = dir(
        path = file.path('output', 'simulation'), 
        pattern = 'sim_dtmc_mcmc_([A-z0-9]{8}|catchup)\\.rds', 
        full.names = TRUE
      )
      
      # compile results from DTMC MCMC samplers
      df = do.call(rbind, lapply(sample.files, function(f) {
        samples = readRDS(f)[[1]]
        if(!('obs' %in% names(samples))) {
          if(file.exists(gsub(pattern = '\\.rds', '_obs.rds', f))) {
            samples$obs = readRDS(gsub(pattern = '\\.rds', '_obs.rds', f))$obs
          } else {
            if(grepl('catchup', f)) {
              samples$obs = sim_obs[[1]]
            }
          }
        }
        samples$samples$param_vec = cbind(
          samples$samples$param_vec,
          'theta' = exp(
            samples$samples$param_vec[,'log_theta']
          )
        )
        hpd = HPDinterval(mcmc(samples$samples$param_vec))
        colnames(hpd) = c('lwr','upr')
        data.frame(
          param = colnames(samples$samples$param_vec),
          mean = colMeans(mcmc(samples$samples$param_vec)),
          betaARFixed = all(samples$samples$param_vec[,'betaAR'] == 0),
          hpd,
          truth = c(
            samples$obs$params$betaAR,
            samples$obs$params$beta,
            exp(samples$obs$params$beta)
          ),
          method = 'DTMC',
          obs.interval = samples$obs$obs_interval,
          scenario = paste(
            'theta = ', exp(samples$obs$params$beta),
            ' betaAR = ', samples$obs$params$betaAR,
            sep = ''
          )
        )
      }))
      
      # merge simulation results with hanks results
      df = rbind(
        df,
        cbind(simulation_hanks_summaries %>% mutate(method = 'AID'), 
              betaARFixed = FALSE),
        cbind(simulation_hanks_univariate_summaries %>% mutate(method = 'AID'),
              betaARFixed = TRUE)
      )
      
      # table version of results
      df %>% 
        filter(betaARFixed == TRUE,
               scenario == 'theta = 1 betaAR = 0',
               param == 'log_theta') %>%
        select(method, obs.interval, mean, lwr, upr) %>%
        mutate(mean = round(mean, 2),
               lwr = round(lwr, 2),
               upr = round(upr, 2)) %>% 
        arrange(method, -obs.interval)
      
      # visualize results!
      pl_univariate = ggplot(df %>% 
                    filter(scenario == 'theta = 1 betaAR = 0',
                           betaARFixed == TRUE) %>%
                    mutate(method = gsub('DTMC', 'State space', method),
                           # method = gsub('AID', 'Data augmentation', method),
                           param = gsub('betaAR', 'beta[1]', param),
                           param = gsub('^theta', 'beta[0]', param),
                           Approximation = method) %>%
                    filter(param == 'log_theta'), 
                  aes(x = obs.interval + 
                        ifelse(method =='State space', .01, -.01), 
                      y = mean, ymin = lwr, ymax = upr, col = Approximation)) + 
        # horizontal lines at truth 
        geom_hline(mapping = aes(yintercept = truth), lty = 3) + 
        # posterior means and HPD intervals
        geom_pointrange() + 
        # formatting
        xlab(expression('Time between observations'~(Delta))) + 
        ylab(expression(beta[0])) +
        scale_color_brewer(type = 'qual', palette = 'Dark2') + 
        theme_few() + 
        theme(axis.title.y = element_text(angle = 0, vjust = .5))
      
      f = file.path('output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '_univariate.png', sep=''))
      ggsave(pl_univariate, filename = f)
      
      # visualize results!
      pl_bivariate = ggplot(df %>% 
         filter(scenario == 'theta = 1 betaAR = 1',
                betaARFixed == FALSE) %>%
         mutate(method = gsub('DTMC', 'State space', method),
                method = gsub('AID', 'Data augmentation', method),
                param = gsub('betaAR', 'beta[1]', param),
                param = gsub('log_theta', 'beta[0]', param),
                Approximation = method) %>%
         filter(param %in% c('beta[0]', 'beta[1]')), 
       aes(x = obs.interval + 
             ifelse(method =='State space', .01, -.01), 
           y = mean, ymin = lwr, ymax = upr, col = Approximation)) + 
        # horizontal lines at truth 
        geom_hline(mapping = aes(yintercept = truth), lty = 3) + 
        # posterior means and HPD intervals
        geom_pointrange() + 
        # formatting
        xlab(expression('Time between observations'~(Delta))) + 
        ylab(expression(beta[0])) +
        scale_color_brewer(type = 'qual', palette = 'Dark2') + 
        theme_few() + 
        theme(axis.title.y = element_text(angle = 0, vjust = .5))
      
      pl_bivariate
      
      f = file.path('output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '.png', sep=''))
      ggsave(pl, filename = f)
      
      f
    }
  ),
  
  tar_target(
    name = simulation_results_combined_bivariate, 
    command = {
      
      # DTMC MCMC sample files
      sample.files = dir(
        path = file.path('output_bak_20211015', 'simulation'), 
        pattern = 'sim_dtmc_mcmc', 
        full.names = TRUE
      )
      
      # compile results from DTMC MCMC samplers
      df = do.call(rbind, lapply(sample.files, function(f) {
        samples = readRDS(f)[[1]]
        if(!('obs' %in% names(samples))) {
          if(file.exists(gsub(pattern = '\\.rds', '_obs.rds', f))) {
            samples$obs = readRDS(gsub(pattern = '\\.rds', '_obs.rds', f))$obs
          } else {
            if(grepl('catchup', f)) {
              samples$obs = sim_obs[[1]]
            }
          }
        }
        samples$samples$param_vec = cbind(
          samples$samples$param_vec,
          'theta' = exp(
            samples$samples$param_vec[,'log_theta']
          )
        )
        hpd = HPDinterval(mcmc(samples$samples$param_vec))
        colnames(hpd) = c('lwr','upr')
        data.frame(
          param = colnames(samples$samples$param_vec),
          mean = colMeans(mcmc(samples$samples$param_vec)),
          betaARFixed = all(samples$samples$param_vec[,'betaAR'] == 0),
          hpd,
          truth = c(
            samples$obs$params$betaAR,
            samples$obs$params$beta,
            exp(samples$obs$params$beta)
          ),
          method = 'DTMC',
          obs.interval = samples$obs$obs_interval,
          scenario = paste(
            'theta = ', exp(samples$obs$params$beta),
            ' betaAR = ', samples$obs$params$betaAR,
            sep = ''
          )
        )
      }))
      
      # merge simulation results with hanks results
      df = rbind(
        df,
        cbind(simulation_hanks_summaries %>% mutate(method = 'AID'), 
              betaARFixed = FALSE),
        cbind(simulation_hanks_univariate_summaries %>% mutate(method = 'AID'),
              betaARFixed = TRUE)
      )
      
      # table version of results
      print(
        df %>% 
          filter(betaARFixed == FALSE,
                 scenario == 'theta = 1 betaAR = 1',
                 method == 'DTMC',
                 param %in% c('log_theta')) %>%
          select(method, obs.interval, mean, lwr, upr, param) %>%
          mutate(mean = round(mean, 2),
                 lwr = round(lwr, 2),
                 upr = round(upr, 2)) %>% 
          arrange(param, method, -obs.interval)
      )
      
      # visualize results!
      pl_bivariate = ggplot(df %>% 
            filter(scenario == 'theta = 1 betaAR = 1',
                   betaARFixed == FALSE) %>%
            mutate(method = gsub('DTMC', 'State space', method),
                   param = gsub('betaAR', 'beta[1]', param),
                   param = gsub('log_theta', 'beta[0]', param),
                   Approximation = method) %>%
            filter(param %in% c('beta[0]', 'beta[1]')), 
          aes(x = obs.interval + 
                ifelse(method =='State space', .01, -.01), 
              y = mean, ymin = lwr, ymax = upr, col = Approximation)) + 
        # horizontal lines at truth 
        geom_hline(mapping = aes(yintercept = truth), lty = 3) + 
        # posterior means and HPD intervals
        geom_pointrange() + 
        # split plot
        facet_wrap(~param, nrow = 2, scales = 'free_y', strip.position = 'left', 
                   labeller = label_parsed) + 
        # formatting
        scale_x_continuous(breaks = c(1/4,1/2,1,2)) + 
        xlab(expression('Time between observations'~(Delta))) + 
        scale_color_brewer(type = 'qual', palette = 'Dark2') + 
        theme_few() + 
        theme(strip.text.y.left = element_text(angle = 0, vjust = .5),
              strip.placement = "outside",
              axis.title.y = element_blank())
      
      pl_bivariate
      
      f = file.path('output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      f = file.path(f, paste(tar_name(), '.png', sep=''))
      ggsave(pl_bivariate, filename = f)
      
      f
    }
  )
  
)
