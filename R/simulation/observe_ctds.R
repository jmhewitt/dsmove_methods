observe_ctds = function(sim_trajectory, n_obs, sim_dir) {
  
  # load simulated trajectory
  ctds_sim = readRDS(sim_trajectory)
  
  # tseq = seq(from = ctds_sim$times[1], 
  #            to = ctds_sim$times[length(ctds_sim$times)], 
  #            length.out = n_obs)
  
  tseq = sort(runif(n = n_obs, min = ctds_sim$times[1], 
               max = ctds_sim$times[length(ctds_sim$times)]))
  
  # observe trajectory
  ctds_obs = ctds.observe(
    states = ctds_sim$states, times = ctds_sim$times, 
    t.obs = tseq)
  
  # save object
  f = file.path(sim_dir, paste(id_chr(), '.rds', sep = ''))
  dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(ctds_obs, file = f)
  
  f
}


