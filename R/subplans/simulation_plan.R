simulation_plan = drake_plan(

  # directory for simulation files
  sim_dir = file.path('output', 'simulation'),
  sim_plots = file.path('plots', 'simulation'),
  
  # build directories
  make_dirs = {
    dir.create(sim_dir, recursive = TRUE)
    dir.create(sim_plots, recursive = TRUE)
  },
  
  # build domain and CTDS representation
  sim_domain = spatial_lattice_ctds(
    n_coord_dimensions = 2, 
    cells_per_dimension = c(100, 100)),
  
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
    transform = map(shape_param = c(1,5,10)),
    .tag_out = shape_param
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
    format = 'file'
  ),

  plot_durations = target(
    command = {
      sim_pkg = readRDS(sim_trajectory)
      ctds_sim = sim_pkg$sim
      pl = ggplot(data.frame(x = ctds_sim$durations), aes(x=x)) +
        stat_density(geom = 'line') +
        theme_few() +
        theme(panel.border = element_blank()) +
        xlab('Duration') +
        ylab('Density') + 
        ggtitle(paste('Weibull shape = ', sim_pkg$params$weibull_shape, 
                      sep = ''))
      ggsave(pl, filename = file.path(sim_plots,
                                      paste(id_chr(), '.png', sep = '')))
    },
    transform = map(sim_trajectory)
  ),

  # plot simulation
  sim_plot = target(
    command = {
      sim_pkg = readRDS(sim_trajectory)
      ctds_sim = sim_pkg$sim
      colnames(sim_domain$coords) = c('x', 'y')
      sim_domain$coords = data.frame(sim_domain$coords)
      pl = plot(x = ctds_sim, ctds_struct = sim_domain) +
        ggtitle(paste('Simulated trajectory (Weibull shape = ', 
                      sim_pkg$params$weibull_shape, ')',
                      sep = ''))
      ggsave(pl, filename = file.path(sim_plots,
                                      paste(id_chr(), '.png', sep = '')))
    },
    transform = map(sim_trajectory)
  ),
  
  # observe simulated trajectory
  sim_obs = target(
    command = {
      sim_pkg = readRDS(sim_trajectory)
      ctds_sim = sim_pkg$sim
      sim_params = sim_pkg$params
      n_obs = obs_per_sec * (sim_params$tf - sim_params$t0)
      obs = observe_ctds(ctds_sim, n_obs)
      # save object
      f = file.path(sim_dir, paste(id_chr(), '.rds', sep = ''))
      saveRDS(list(obs = obs, params = sim_pkg$params, nobs = n_obs), file = f)
      f
    },
    transform = cross(sim_trajectory, 
                      obs_per_sec = c(0.5, 1, 2, 4, 8)),
    format = 'file'
  ),
  
  # 
  # # plot likelihood surface
  # # sim_lik_surface = target(
  # #   plot_obs_lik(ctds_struct = sim_domain, obs = sim_obs, plot_dir = sim_plots,
  # #                beta_loc_seq = seq(from = -1.5, to = 2.5, length.out = 25),
  # #                beta_ar_seq = seq(from = -1, to = 5, length.out = 25),
  # #                # beta_loc_seq = seq(from = 2, to = 8, length.out = 25),
  # #                # beta_ar_seq = seq(from = -1, to = 1, length.out = 25),
  # #                beta_loc, beta_ar),
  # #   dynamic = map(sim_obs), 
  # #   max_expand = 1
  # # ),
  # 
  
  # fit model using crawl approximation
  sim_fit_hanks = target(
    command = {
      obs_pkg = readRDS(sim_obs)
      fit = fit_hanks(ctds_struct = sim_domain, ctds_obs = obs_pkg$obs, 
                weibull_est = weibull_est)
      list(fit = fit, weibull_est = weibull_est, 
           params = obs_pkg$params, nobs = obs_pkg$nobs)
    },
    transform = cross(sim_obs, weibull_est = c(TRUE, FALSE))
  ),
  
  # assemble crawl approximation results across timesteps
  sim_fit_hanks_summary = target(
    command = {
      summarize_fit_hanks(fit_output = sim_fit_hanks)
    },
    transform = map(sim_fit_hanks)
  ),
  
  aggregated_summaries = target(
    dplyr::bind_rows(sim_fit_hanks_summary),
    transform = combine(sim_fit_hanks_summary, by = weibull_est)
  ),

  plot_hanks_summary = {
    pl = ggplot(aggregated_summaries,
           aes(x = tstep, y = Estimate, ymin = lwr, ymax = upr,
               col = weibull_est)) +
      geom_pointrange() +
      geom_hline(mapping = aes(yintercept = truth), lty = 3) +
      xlab('Time between observations') +
      ggtitle('Hanks recovery of true parameters (Truth at dotted line)') +
      facet_grid(param~weibull_truth, scales = 'free') +
      theme_few()

    ggsave(pl, filename = file.path(sim_plots,
                                    paste(id_chr(), '.pdf', sep = '')))

  }

  # 
  # # compare imputed trajectories to truth
  # # imputed_plots = target(
  # #   plot_imputations(ctds_struct = sim_domain, trajectory = sim_trajectory,
  # #                    imputations = sim_fit_hanks$ctmc.list,
  # #                    raster.coords = sim_fit_hanks$raster.coords,
  # #                    tstep = sim_fit_hanks$tstep,
  # #                    plot_dir = sim_plots),
  # #   dynamic = map(sim_fit_hanks)
  # # ),
  # 
  # # compare shortest path imputation to truth
  # # shortest_path_plots = target(
  # #   plot_shortest_path(ctds_struct = sim_domain, trajectory = sim_trajectory, 
  # #                      obs = sim_obs, plot_dir = sim_plots, 
  # #                      beta_loc = beta_loc),
  # #   transform = map(sim_obs)
  # # ),
  # 
  # # fit model to completely observed trajectory
  # # sim_fit_exact = target(
  # #   fit_exact(ctds_struct = sim_domain, trajectory = sim_trajectory, 
  # #             beta_loc = beta_loc, beta_dir = beta_dir, beta_ar = beta_ar)
  # # ),
  # 
  # # assemble exact fit results
  # # sim_fit_exact_summary = 
  # #   data.frame(
  # #     param = c('beta_loc', 'beta_ar'), 
  # #     truth = c(beta_loc, beta_ar),
  # #     est = sim_fit_exact$par[-2],
  # #     se = sqrt(diag(solve(-sim_fit_exact$hessian[-2,-2])))
  # #   ) %>% 
  # #   dplyr::mutate(
  # #     lwr = est - 1.96 * se,
  # #     upr = est + 1.96 * se,
  # #     covered = (lwr <= truth) & (truth <= upr)
  # #   )
  # 
)
