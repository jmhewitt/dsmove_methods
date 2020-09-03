simulation_plan = drake_plan(

  # directory for simulation files
  sim_dir = file_out(!!file.path('output', 'simulation')),
  sim_plots = file_out(!!file.path('plots', 'simulation')),
  
  # simulation domain dimensions
  n_coord_dimensions = 2,
  cells_per_dimension = c(100, 100),
  
  # build domain and CTDS representation
  sim_domain = spatial_lattice_ctds(n_coord_dimensions, cells_per_dimension),
  
  # simulation parameters
  beta_loc = matrix(-1, nrow = 1, ncol = 1),
  beta_dir = 0,
  beta_ar = 1,
  t0 = 0,
  tf = 200,
  max.steps = 1e3,
  obs_per_sec = c(0.5, 1, 2, 4, 8),
  n_obs = obs_per_sec * (tf - t0),
  
  # simulate a CTDS trajectory
  sim_trajectory = target(
    sim_ctds(sim_domain, beta_loc, beta_dir, beta_ar, t0, tf, 
             max.steps, sim_dir),
    format = 'file'
  ),
  
  # look at average duration in simulated trajectory
  sim_avg_duration = {
    ctds_sim = readRDS(sim_trajectory)
    mean(ctds_sim$durations)
  },
  
  # plot simulation
  sim_plot = {
    ctds_sim = readRDS(sim_trajectory)
    colnames(sim_domain$coords) = c('x', 'y')
    sim_domain$coords = data.frame(sim_domain$coords)
    pl = plot(x = ctds_sim, ctds_struct = sim_domain) + 
      ggtitle('Simulated trajectory')
    ggsave(pl, filename = file.path(sim_plots, 
                                    paste(id_chr(), '.png', sep = '')))
  },
  
  # observe simulated trajectory
  sim_obs = target(
    observe_ctds(sim_trajectory, n_obs, sim_dir),
    dynamic = map(n_obs),
    format = 'file'
  ),
  
  # plot likelihood surface
  sim_lik_surface = target(
    plot_obs_lik(ctds_struct = sim_domain, obs = sim_obs, plot_dir = sim_plots,
                 beta_loc_seq = seq(from = -1.5, to = 2.5, length.out = 25),
                 beta_ar_seq = seq(from = -1, to = 5, length.out = 25),
                 # beta_loc_seq = seq(from = 2, to = 8, length.out = 25),
                 # beta_ar_seq = seq(from = -1, to = 1, length.out = 25),
                 beta_loc, beta_ar),
    dynamic = map(sim_obs), 
    max_expand = 1
  ),
  
  # fit model using crawl approximation
  sim_fit_hanks = target(
    fit_hanks(ctds_struct = sim_domain, obs = sim_obs, plot_dir = sim_plots),
    dynamic = map(sim_obs)
  ),
  
  # assemble crawl approximation results across timesteps
  sim_fit_hanks_summary = target(
    summarize_fit_hanks(fit = sim_fit_hanks, beta_loc = beta_loc, 
                        beta_ar = beta_ar),
    dynamic = map(sim_fit_hanks)
  ),
  
  plot_hanks_summary = {
    pl = ggplot(sim_fit_hanks_summary, 
           aes(x = tstep, y = Estimate, ymin = lwr, ymax = upr)) + 
      geom_pointrange() + 
      geom_hline(mapping = aes(yintercept = Estimate), 
                 data = data.frame(
                   param = c('(Intercept)', 'crw'),
                   Estimate = c(beta_loc, beta_ar)), 
                 lty = 3) + 
      xlab('Time between observations') +
      ggtitle('Hanks recovery of true parameters (Truth at dotted line)') + 
      facet_wrap(~param) +
      theme_few()
    
    ggsave(pl, filename = file.path(sim_plots, 
                                    paste(id_chr(), '.pdf', sep = '')))
      
  },
  
  # compare imputed trajectories to truth
  imputed_plots = target(
    plot_imputations(ctds_struct = sim_domain, trajectory = sim_trajectory,
                     imputations = sim_fit_hanks$ctmc.list,
                     raster.coords = sim_fit_hanks$raster.coords,
                     tstep = sim_fit_hanks$tstep,
                     plot_dir = sim_plots),
    dynamic = map(sim_fit_hanks)
  ),
  
  # fit model to completely observed trajectory
  sim_fit_exact = target(
    fit_exact(ctds_struct = sim_domain, trajectory = sim_trajectory, 
              beta_loc = beta_loc, beta_dir = beta_dir, beta_ar = beta_ar)
  ),
  
  # assemble exact fit results
  sim_fit_exact_summary = 
    data.frame(
      param = c('beta_loc', 'beta_ar'), 
      truth = c(beta_loc, beta_ar),
      est = sim_fit_exact$par[-2],
      se = sqrt(diag(solve(-sim_fit_exact$hessian[-2,-2])))
    ) %>% 
    dplyr::mutate(
      lwr = est - 1.96 * se,
      upr = est + 1.96 * se,
      covered = (lwr <= truth) & (truth <= upr)
    )
  
)
