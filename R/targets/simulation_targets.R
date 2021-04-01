simulation_targets = list(
  
  #
  # simple simulation
  #
  
  # parameters to generate basic CTDS trajectories
  tar_target(
    name = sim_params, 
    command = list(
      t0 = 0,               # start time for simulation
      tf = 500,             # end time for simulation
      betaAR = 1,           # directional persistence parameter
      beta = 0,             # log of transition rate
      dims = rep(1e3, 2),   # support for spatial domain
      x0 = rep(500, 2),     # location from which simulated trajectory begins
      x0.last = c(499, 500) # auto regressive loc. to set trajectory direction
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
  
  tar_target(obs_interval, c(5, 1, .5, .25)),
  
  tar_target(
    name = sim_path_segments,  
    # cue = tar_cue(mode = 'always'),
    command = {
      # remove outer layer of wrapping
      sim_obs = sim_obs[[1]]
      # sample path segments
      list(lapply(1:(nrow(sim_obs$states) - 1), function(obs_ind) {
        list(sample_path_segments(
          x0 = c(sim_obs$states[obs_ind,], 1), 
          xf = c(sim_obs$states[obs_ind + 1,], 1), 
          t0 = sim_obs$times[obs_ind],
          tf = sim_obs$times[obs_ind + 1],
          dims = c(sim_params$dims, 1), 
          max_tx_rate = 2,
          high_quantile = .99,
          computational_max_steps = 1e3, 
          surface_heights = rep(0, prod(sim_params$dims)),
          domain_heights = 1,
          nsegments = 100
        ))
      }))
    }, 
    pattern = map(sim_obs)
  ),
  
  tar_target(
    name = sim_fits, 
    command = {

      priors = list(
        beta_ar = list(mean = 0, sd = 1e2),
        beta = list(mean = 0, sd = 1e2)
      )
      
      post_samples = fit_single(
        niter = 1e4, priors = priors, 
        path_fam = sim_path_segments[[1]], dims = sim_params$dims, 
        seg_start_times = sim_obs[[1]]$times, 
        params = list(beta_ar = 0, beta_loc = 0),
        init_sd = rep(.05, 2)
      )
      
      list(list(samples = post_samples, rep = reps, sim_obs = sim_obs,
                obs_interval = obs_interval))
    }, 
    pattern = cross(map(sim_path_segments, sim_obs, obs_interval), reps), 
    deployment = 'worker'
  ),
  
  tar_target(reps, 1:100)

  
)
