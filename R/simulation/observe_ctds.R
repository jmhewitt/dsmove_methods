observe_ctds = function(ctds_sim, n_obs) {
  
  tseq = sort(runif(n = n_obs, min = ctds_sim$times[1], 
               max = ctds_sim$times[length(ctds_sim$times)]))
  
  # observe trajectory
  ctds_obs = ctds.observe(
    states = ctds_sim$states, times = ctds_sim$times, 
    t.obs = tseq)
  
  ctds_obs
}


