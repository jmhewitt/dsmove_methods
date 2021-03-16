simulation_targets = list(
  
  #
  # simple simulation
  #
  
  # parameters to generate basic CTDS trajectories
  tar_target(
    name = sim_params, 
    command = list(
      t0 = 0,             # start time for simulation
      tf = 100,           # end time for simulation
      betaAR = 1,         # directional persistence parameter
      beta = 0,           # log of transition rate
      dims = rep(100, 2), # support for spatial domain
      x0 = rep(50, 2),    # location from which simulated trajectory begins
      x0.last = c(49, 50) # auto regressive loc. to set direction for trajectory
  )),
  
  # realization of basic CTDS trajectory
  tar_target(
    name = sim_simple, 
    command = ctds.quicksim(
      dims = sim_params$dims, 
      beta_loc = sim_params$beta, 
      beta_ar = sim_params$betaAR, 
      v0 = sim_params$x0, 
      t0 = sim_params$t0, 
      tf = sim_params$tf, 
      max.steps = 1e3, 
      v0.last = sim_params$x0.last
  )),
    
  # observe CTDS trajectory under increasing temporal resolution
  tar_target(
    name = sim_obs,
    command = list(
      ctds.observe(
        states = sim_simple$states[-1,],
        times = sim_simple$times,
        t.obs = seq(
          from = sim_simple$times[1],
          to = tail(sim_simple$times, 1) + obs_interval,
          by = obs_interval
        )
    )),
    pattern = map(obs_interval)
  ),
  
  tar_target(obs_interval, c(1, .5, .25))
  
  
)
