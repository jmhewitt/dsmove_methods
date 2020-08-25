sim_ctds = function(sim_domain, beta_loc, beta_dir, beta_ar, t0, tf, 
                    max.steps, sim_dir) {
  
  # define central coordinate as average coordinate value in all dimensions
  coord.center = colMeans(apply(sim_domain$coords, 2, range))
  
  # identify index of location closest to coord.center
  center.ind = which.min(apply(
    sweep(x = sim_domain$coords, MARGIN = 2, STATS = coord.center, FUN = '-'),
    1, function(r) sum(r^2)
  ))
  
  # simulate trajectory
  ctds_sim = ctds.fwdsim(ctds_struct = sim_domain, beta_loc = beta_loc, 
                         beta_dir = beta_dir, v0 = center.ind, t0 = t0, tf = tf, 
                         max.steps = max.steps, beta_ar = beta_ar, 
                         v0.last = NULL)
  
  # save simulation object
  f = file.path(sim_dir, 'ctds_sim.rds')
  dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(ctds_sim, file = f)
  
  f
}


