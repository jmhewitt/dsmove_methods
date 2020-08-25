observe_ctds = function(sim_trajectory, n_obs, sim_dir) {
  
  # load simulated trajectory
  ctds_sim = readRDS(sim_trajectory)
  
  # observe trajectory
  ctds_obs = ctds.observe(
    states = ctds_sim$states, times = ctds_sim$times, 
    t.obs = seq(from = ctds_sim$times[1], 
                to = ctds_sim$times[length(ctds_sim$times)], 
                length.out = n_obs))
  
  # save object
  f = file.path(sim_dir, paste(id_chr(), '.rds', sep = ''))
  dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(ctds_obs, file = f)
  
  f
}


