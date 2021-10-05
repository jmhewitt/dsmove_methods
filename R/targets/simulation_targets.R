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
        reps = 10
      )

      list(list(
        samples = post_samples, 
        rep = rep_batch, 
        sim_obs = sim_obs[[1]]
      ))
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
        samples[,'log_theta'] = exp(samples[,'log_theta'])
        colnames(samples)[2] = 'theta'
        # hpds 
        hpds = HPDinterval(mcmc(samples))
        # package results
        data.frame(
          param = colnames(samples),
          mean = colMeans(samples),
          lwr = hpds[,'lower'],
          upr = hpds[,'upper'],
          truth = c(unlist(res_subset[[1]]$sim_obs$params['betaAR']),
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
      ggsave(pl, filename = 'posteriors_2param.pdf', width = 8, height = 5)
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
  
  tar_target(
    name = sim_dtmc_mcmc,
    command = {
      
      # number of MCMC samples to draw
      niter = 1e3
      
      # time discretization timestep
      delta = .05
      
      # unwrap input
      obs = sim_obs[[1]]
      
      # specify priors for analysis
      priors = sim_priors
      
      # discretize time
      tseq = seq(from = obs$params$t0, to = obs$params$tf, by = delta)
      
      # boilerplate dimension information for likelihood approximation fn.
      xcoords = 1:obs$params$dims[1]
      ycoords = 1:obs$params$dims[2]
      surface_heights = numeric(length = length(xcoords) * length(ycoords))
      
      # initialize flattened data structure in discrete time
      oseq = matrix(NA, nrow = length(tseq), ncol = 4)
      colnames(oseq) = c('x_ind', 'y_ind', 'x_coord', 'y_coord')
      
      # flatten data
      for(ind in 1:nrow(obs$obs$states)) {
        tgt_ind = which.min(abs(obs$obs$times[ind] - tseq))
        oseq[tgt_ind, c('x_coord', 'y_coord')] = obs$obs$states[ind,]
        oseq[tgt_ind, 'x_ind'] = which(oseq[tgt_ind,'x_coord'] == xcoords) - 1
        oseq[tgt_ind, 'y_ind'] = which(oseq[tgt_ind,'y_coord'] == ycoords) - 1
      }
      
      # gaussian approximation to the posterior 
      o = optim(par = c(0,.5), fn = function(theta) {
        dsmovetools2d:::ExactLocFilteredLL(
          init_dsts = matrix(oseq[1,c('x_ind','y_ind')], nrow = 4, ncol = 2, 
                             byrow = TRUE),
          init_srcs = matrix(
            c(oseq[1,c('x_ind','y_ind')] + c(0,1),
              oseq[1,c('x_ind','y_ind')] + c(0,-1),
              oseq[1,c('x_ind','y_ind')] + c(1,0),
              oseq[1,c('x_ind','y_ind')] + c(-1,0)),
            nrow = 4, ncol = 2, byrow = TRUE
          ), 
          init_log_probs = rep(log(1/4), 4),
          obs_x_coords = oseq[,'x_coord'],
          obs_y_coords = oseq[,'y_coord'],
          x_coords = xcoords, 
          y_coords = ycoords, 
          surface_heights = surface_heights, 
          log_self_tx = log(1-delta*exp(theta[2])), 
          betaAR = theta[1], 
          lptrunc = Inf
        )
      }, control = list(fnscale = -1), hessian = TRUE)
      
      # wrapper to evaluate log-posterior
      lpfn = function(x, x_ind, theta) {
        # Evaluate likelihood for model while holding all but one parameter fixed.
        # The parameter that is easiest to change is "x", which relates to theta 
        # via theta[x_ind] = x.
        # 
        # Parameters:
        #   x - univariate parameter for which to evaluate likelihood at
        #   x_ind - index in theta where "x" should be inserted
        #   theta - full parameter vector required to evaluate likelihood
        #   path - latent CTDS path along grid
        
        # merge parameter with theta
        theta[x_ind] = x
        
        # observation likelihood
        dsmovetools2d:::ExactLocFilteredLL(
          init_dsts = matrix(oseq[1,c('x_ind','y_ind')], nrow = 4, ncol = 2, 
                             byrow = TRUE),
          init_srcs = matrix(
            c(oseq[1,c('x_ind','y_ind')] + c(0,1),
              oseq[1,c('x_ind','y_ind')] + c(0,-1),
              oseq[1,c('x_ind','y_ind')] + c(1,0),
              oseq[1,c('x_ind','y_ind')] + c(-1,0)),
            nrow = 4, ncol = 2, byrow = TRUE
          ), 
          init_log_probs = rep(log(1/4), 4),
          obs_x_coords = oseq[,'x_coord'],
          obs_y_coords = oseq[,'y_coord'],
          x_coords = xcoords, 
          y_coords = ycoords, 
          surface_heights = surface_heights, 
          log_self_tx = log(1-delta*exp(theta[2])), 
          betaAR = theta[1], 
          lptrunc = Inf
        ) +
          # prior
          dnorm(x = theta[1], mean = priors$betaAR['mean'], 
                sd = priors$betaAR['sd'], log = TRUE) +
          dgamma(x = exp(theta[2]), shape = priors$theta['shape'], 
                 rate =  priors$theta['rate'], log = TRUE) + 
          # Jacobian to account for transformation in which sampling is done on 
          #   unconstrained space, but prior for theta[1] is defined on [0,\infty)
          jac.log(x = theta[2], log = TRUE)
      }
      
      #
      # gibbs sampler to explore posterior, initialized around mode
      #
      
      n_params = 2
      param_vec = o$par
      init_sd = sqrt(diag(solve(-o$hessian)))
      
      # construct MHRW samplers
      paramSamplers = lapply(1:n_params, function(ind) {
        dsmovetools::Mhrw1DAdaptive$new(
          x = param_vec[ind], sd = init_sd[ind],  lp = lpfn, C = .75, 
          alpha = .44, adaptive = TRUE
        )
      })
      
      # initialize output
      samples = list(
        param_vec = matrix(nrow = niter, ncol = length(param_vec)),
        lp = numeric(niter)
      )
      colnames(samples$param_vec) = c('betaAR', 'log_theta')
      
      # run sampler
      for(it in 1:niter) {
        
        message(it)
        
        # update model parameters
        for(i in 1:n_params) {
          update = paramSamplers[[i]]$sample(
            x_ind = i, theta = param_vec
          )
          param_vec[i] = update$x
        }
        
        # save samples
        samples$param_vec[it, ] = param_vec
        samples$lp[it] =  lpfn(x = param_vec, x_ind = 1:n_params, 
                               theta = param_vec)
        
      }
      
      
      # package results
      list(list(
        obs = obs,
        o = o,
        samples = samples
      ))
    },
    pattern = map(sim_obs),
    deployment = 'worker',
    storage = 'worker',
    memory = 'transient',
    error = 'continue'
  )
  
)
