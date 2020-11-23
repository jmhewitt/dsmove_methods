simulation_plan = drake_plan(

  # directory for simulation files
  sim_dir = target(
    file.path('output', 'simulation'),
    hpc = FALSE
  ),
  sim_plots = target(
    file.path('plots', 'simulation'),
    hpc = FALSE
  ),
  
  # build directories
  make_dirs = target({
    dir.create(sim_dir, recursive = TRUE)
    dir.create(sim_plots, recursive = TRUE)
  }, hpc = FALSE),
  
  # build domain and CTDS representation
  sim_domain = target(
    spatial_lattice_ctds(n_coord_dimensions = 2, 
                         cells_per_dimension = c(100, 100)),
    hpc = FALSE
  ),
  
  # simulation parameters
  sim_params = target(
    command = {
      list(
        beta_loc = matrix(-1, nrow = 1, ncol = 1),
        beta_dir = 0,
        beta_ar = 1,
        t0 = 0,
        tf = 200,
        max.steps = 1e3,
        weibull_shape = shape_param,
        seed = 2020
      )
    },
    transform = map(shape_param = 1),
    .tag_out = shape_param,
    hpc = FALSE
  ),
  
  prior_params = target(
    list(
      # maximum transition rate allowed for imputations. chosen such that 
      # max_lambda >> exp(max(sim_domain$Xloc %*% sim_params_1$beta_loc))
      max_lambda = 5,
      p_max = .99,
      beta_ar = list(mean = 0, sd = 1e2),
      beta_loc = list(mean = sim_params_1$beta_loc[1], sd = .1)
    ),
    hpc = FALSE
  ),
  
  # simulate a CTDS trajectory
  sim_trajectory = target(
    command = {
      res = list(
        sim = sim_ctds(sim_domain, sim_params$beta_loc, sim_params$beta_dir, 
                       sim_params$beta_ar, sim_params$t0, sim_params$tf, 
                       sim_params$max.steps, seed = sim_params$seed,
                       weibull_shape = sim_params$weibull_shape),
        params = sim_params
      )
      # save simulation object
      f = file.path(sim_dir, paste(id_chr(), '.rds', sep = ''))
      dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(res, file = f)
      f
    },
    transform = map(sim_params),
    format = 'file',
    hpc = FALSE
  ),

  # observe simulated trajectory
  sim_obs = target(
    command = {
      sim_pkg = readRDS(sim_trajectory)
      ctds_sim = sim_pkg$sim
      sim_params = sim_pkg$params
      # observe ctds trajectory at regular timepoints
      tseq = seq(from = ctds_sim$times[1], to = tail(ctds_sim$times, 1), 
                 by =  1 / obs_per_sec)
      n_obs = length(tseq)
      obs = ctds.observe(states = ctds_sim$states, times = ctds_sim$times, 
                         t.obs = tseq)
      # save object
      f = file.path(sim_dir, paste(id_chr(), '.rds', sep = ''))
      saveRDS(list(obs = obs, params = sim_pkg$params, nobs = n_obs), file = f)
      f
    },
    transform = cross(sim_trajectory, 
                      obs_per_sec = c(0.5, 1, 2, 4, 8)),
    format = 'file',
    hpc = FALSE
  ),
  
  # plot observations alongside true trajectory
  plot_sim_obs = target(
    command = {
      
      # load true trajectory and observations
      sim_pkg = readRDS(sim_trajectory)
      obs = readRDS(sim_obs)
      
      # extract spatial coordinates for trajectory and observations
      sim_coords = data.frame(
        sim_domain$coords[sim_pkg$sim$states,],
        time = sim_pkg$sim$times
      ) %>% 
        pivot_longer(cols = s1:s2, names_to = 'coord')
      obs_coords = data.frame(
        sim_domain$coords[obs$obs$states,],
        time = obs$obs$times
      ) %>% 
        pivot_longer(cols = s1:s2, names_to = 'coord')
      
      # build plot
      pl = ggplot(sim_coords, aes(x = time, y = value)) + 
        # complete trajectory
        geom_step() + 
        # observations
        geom_point(data = obs_coords, col = 'red', size = .1) + 
        # plot each spatial coordinate separately
        facet_wrap(~coord, nrow = 2, ncol = 1, scales = 'free',
                   strip.position = 'left') +
        theme_few() +  
        theme(strip.placement = 'outside', 
              strip.text.y = element_text(angle = 180)) + 
        xlab('Time') + 
        ylab('Coord.') + 
        ggtitle('CTDS trajectory (observations in red)')
        
      # save object
      f = file.path(sim_plots, paste(id_chr(), '.pdf', sep = ''))
      ggsave(pl, filename = f)
      f
    },
    transform = map(sim_obs),
    hpc = TRUE,
    format = 'file'
  ),
  
  # pre-compute family of path segments that could connect observations
  impute_segments = target(
    command = {
      # load observations
      sim = readRDS(sim_obs)
      obs = sim$obs
      # sample segments
      segments = sample_segments(states = obs$states, times = obs$times, 
                                 ctds_domain = sim_domain, 
                                 max_lambda = prior_params$max_lambda, 
                                 p_max = prior_params$p_max, nsamples = 100)
      # save object
      f = file.path(sim_dir, paste(id_chr(), '.rds', sep = ''))
      saveRDS(list(params = sim$params, segments = segments, obs = obs), 
              file = f)
      f
    },
    transform = map(sim_obs), 
    format = 'file',
    hpc = TRUE
  ),
  
  # path "seeds" used to explore the posterior space to initialize gibbs sampler
  useq = target({
    set.seed(2020)
    npaths = 10
    unique(pmin(c(0,
                  seq(from = 0, to = 1, length.out = npaths) +
                    runif(n = npaths, min = 0, max = 1/npaths),
                  1), 1
    ))
  }, hpc = FALSE),
  
  # use path seeds to sample initial paths and optimize associated model params.
  init_fits = target(
    command = {
      dat = readRDS(impute_segments)
      inits = init_integration(segments = dat$segments, obs = dat$obs, 
                                inits = list(beta_loc = 0, beta_ar = 0), 
                                priors = prior_params, niter = 50,  u = useq,
                                ctds_domain = sim_domain)
      list(list(inits = inits, file.segments = impute_segments))
    },
    transform = map(impute_segments),
    dynamic = map(useq),
    format = 'rds',
    hpc = TRUE
  ),
  
  # plot observations alongside true trajectory
  plot_init_fits = target(
    command = {
    
      ctds_paths = lapply(init_fits, function(init) {
        # extract and flatten optimized path
        path_components = init$inits$path_components
        epath = do.call(c, sapply(path_components, function(p) {
          unlist(p$epath)
        }))
        tpath = do.call(c, sapply(path_components, function(p) { p$tpath }))
        list(epath = epath, tpath = tpath, durations = diff(tpath))
      })
      
      # extract optimized initialization paths
      init_paths = lapply(init_fits, function(init) {
        # extract and flatten optimized path
        path_components = init$inits$path_components
        epath = do.call(c, sapply(path_components, function(p) {
          unlist(p$epath)
        }))
        tpath = do.call(c, sapply(path_components, function(p) { p$tpath }))
        spath = sim_domain$edge_df$to[epath]
        # convert optimized path to locations
        init_coords = data.frame(
          sim_domain$coords[spath,],
          time = tpath
        ) %>% 
          pivot_longer(cols = s1:s2, names_to = 'coord')
      })
      
      # load true trajectory and observations
      sim_pkg = readRDS(sim_trajectory)
      obs = readRDS(sim_obs)
      
      # extract spatial coordinates for trajectory and observations
      sim_coords = data.frame(
        sim_domain$coords[sim_pkg$sim$states,],
        time = sim_pkg$sim$times
      ) %>% 
        pivot_longer(cols = s1:s2, names_to = 'coord')
      obs_coords = data.frame(
        sim_domain$coords[obs$obs$states,],
        time = obs$obs$times
      ) %>% 
        pivot_longer(cols = s1:s2, names_to = 'coord')
      
      # build plots
      pl_inits = mapply(function(ipath, u) {
        df = rbind(
          cbind(sim_coords, type = 'Exact/Simulated'),
          cbind(ipath, type = 'Initial imputation')
        )
        ggplot(df, aes(x = time, y = value, col = type)) + 
          # complete trajectories
          geom_step() + 
          scale_color_brewer('Trajectory', type = 'qual', palette = 'Dark2') + 
          # observations
          geom_point(data = obs_coords, col = 'red', size = .1) + 
          # plot each spatial coordinate separately
          facet_wrap(~coord, nrow = 2, ncol = 1, scales = 'free',
                     strip.position = 'left') +
          theme_few() +  
          theme(strip.placement = 'outside', 
                strip.text.y = element_text(angle = 180)) + 
          xlab('Time') + 
          ylab('Coord.') + 
          ggtitle(paste('CTDS trajectory (observations in red; u=', u, ')', 
                        sep = ''))
      }, init_paths, useq, SIMPLIFY = FALSE)
      
      # save objects
      mapply(function(pl, u) {
        f = file.path(sim_plots, paste(id_chr(), '_u', u, '.pdf', sep = ''))
        ggsave(pl, filename = f)
        f
      }, pl_inits, useq)
      
    },
    transform = map(init_fits),
    hpc = TRUE,
    format = 'file'
  ),
  
  init_fits_summary = target(
    command = {
      
      # extract optimized parameters
      init_params = cbind(
        u = useq,
        do.call(rbind, lapply(init_fits, function(init) {
          data.frame(init$inits$params, ll = init$inits$ll)
        }))
      )
      
      # save as csv
      f = file.path(sim_plots, paste(id_chr(), '.csv', sep = ''))
      write.csv(init_params, file = f, row.names = FALSE)
      f
    },
    transform = map(init_fits),
    hpc = TRUE,
    format = 'file'
  ),
  
  # # approximate posterior for initial path with highest log-posterior
  # gibbs_fits = target(
  #   command = {
  #     
  #     # extract the best-fitting initial trajectory 
  #     selected_init = which.max(sapply(init_fits, function(init) {
  #       init$inits$ll
  #     }))
  #     init = init_fits[[selected_init]]
  #     
  #     # re-package the initial parameters
  #     init$inits$beta_ar = init$inits$params$beta_ar
  #     init$inits$beta_loc = init$inits$params$beta_loc
  #     init$inits$prop_cov = solve(-init$inits$hessian)
  #     
  #     # load family of path segments (and observations)
  #     segment.group = readRDS(init$file.segments)
  #     
  #     # define priors
  #     priors = list(
  #       beta_ar = list(mean = 0, sd = 1e2),
  #       beta_loc = list(mean = 0, sd = 1e2)
  #     )
  #     
  #     # set file output
  #     f = file.path(sim_dir, paste(id_chr(), '.rds', sep = ''))
  #     
  #     # sample save function
  #     checkpoint_fn = function(samples) {
  #       o = list(
  #         samples = samples,
  #         obs_per_sec = 1/unique(diff(segment.group$obs$times)),
  #         params = segment.group$params
  #       )
  #       saveRDS(o, file = f)
  #     }
  #     
  #     # # run gibbs sampler, receive file name with outputs
  #     # samples = fit_integration(
  #     #   segments = segment.group$segments, obs = segment.group$obs, niter = 1e4, 
  #     #   ncheckpoints = 1e2, inits = init$inits, priors = priors, 
  #     #   ctds_domain = sim_domain, checkpoint_fn = checkpoint_fn
  #     # )
  #     # 
  #     # # return 
  #     # checkpoint_fn(samples)
  #     
  #     # return file name
  #     f
  #   },
  #   transform = map(init_fits),
  #   format = 'file'
  # ),
  # 
  # # summarize approximate posterior
  # gibbs_summaries = target(
  #   command = {
  #     if(file.exists(gibbs_fits)) {
  #       # read mcmc output
  #       mcout = readRDS(gibbs_fits)
  #       # extract posterior samples for model parameters
  #       niter = sum(mcout$samples$ll != 0)
  #       nburn = 1e3
  #       m = mcmc(mcout$samples$param_vec[nburn:niter, , drop = FALSE])
  #       # aggregate and return posterior summaries
  #       hpd = HPDinterval(m)
  #       data.frame(
  #         parameter = c('beta_ar', 'beta_loc'),
  #         mean = colMeans(m),
  #         se = apply(m, 2, sd),
  #         lwr = hpd[, 1],
  #         upr = hpd[, 2],
  #         obs_per_second = mcout[['obs_per_sec']],
  #         truth = c(mcout$params[['beta_ar']], mcout$params[['beta_loc']])
  #       ) %>% 
  #         dplyr::mutate(covered = lwr <= truth & truth <= upr)
  #     } else {
  #       NULL
  #     }
  #   },
  #   transform = map(gibbs_fits),
  #   hpc = TRUE
  # ),
  
  # gibbs_summary_plots = target(
  #   command = {
  #     # merge simulation output
  #     df = dplyr::combine(gibbs_summaries)
  #     
  #     
  #   }, 
  #   transform = combine(gibbs_summaries),
  #   hpc = FALSE
  # )
  
  # fit model using crawl approximation
  sim_fit_hanks = target(
    command = {
      obs_pkg = readRDS(sim_obs)
      fit = fit_hanks(ctds_struct = sim_domain, ctds_obs = obs_pkg$obs)
      list(fit = fit, params = obs_pkg$params, nobs = obs_pkg$nobs)
    },
    transform = cross(sim_obs),
    hpc = TRUE
  ),
  
  # assemble crawl approximation results across timesteps
  sim_fit_hanks_summary = target(
    command = {
      summarize_fit_hanks(fit_output = sim_fit_hanks)
    },
    transform = map(sim_fit_hanks),
    hpc = TRUE
  ),
  
  aggregated_summaries = target(
    dplyr::bind_rows(sim_fit_hanks_summary),
    transform = combine(sim_fit_hanks_summary),
    hpc = TRUE
  ),

  plot_hanks_summary = target({
      pl = ggplot(aggregated_summaries,
           aes(x = tstep, y = Estimate, ymin = lwr, ymax = upr)) +
      geom_pointrange() +
      geom_hline(mapping = aes(yintercept = truth), lty = 3) +
      xlab('Time between observations') +
      ggtitle('Hanks recovery of true parameters (Truth at dotted line)') +
      facet_wrap(~param, scales = 'free') +
      theme_few()

    ggsave(pl, filename = file.path(sim_plots,
                                    paste(id_chr(), '.pdf', sep = '')))

  }, hpc = TRUE)
  
  # sim_mle_ests = target(
  #   sim_mle$est,
  #   transform = map(sim_mle)
  # ),
  
)
