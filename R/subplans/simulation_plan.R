simulation_plan = drake_plan(

  # directory for simulation files
  sim_dir = file_out(!!file.path('output', 'simulation')),
  sim_plots = file_out(!!file.path('plots', 'simulation')),
  
  # simulation domain dimensions
  n_coord_dimensions = 2,
  cells_per_dimension = c(5, 5),
  
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
  
  # observe simulated trajectory
  sim_obs = target(
    observe_ctds(sim_trajectory, n_obs, sim_dir),
    dynamic = map(n_obs),
    format = 'file'
  ),
  
  # plot likelihood surface
  sim_lik_surface = target(
    plot_obs_lik(ctds_struct = sim_domain, obs = sim_obs, plot_dir = sim_plots,
                 # beta_loc_seq = seq(from = -1, to = 5, length.out = 25),
                 # beta_ar_seq = seq(from = -3, to = 5, length.out = 25),
                 beta_loc_seq = seq(from = -1.5, to = 10, length.out = 25),
                 beta_ar_seq = seq(from = -10, to = 3, length.out = 25),
                 beta_loc, beta_ar),
    dynamic = map(sim_obs)
  )
  
)
