library(nimble)

# central location coordinates
box.center = c(s1 = 50, s2 = 50)

# central location index
box.ind = which(
  coords[,'s1'] == box.center['s1'] & 
    coords[,'s2'] == box.center['s2']
)

ctds_struct$coords = data.frame(ctds_struct$coords)
ctds_struct$coords$x = ctds_struct$coords$s1
ctds_struct$coords$y = ctds_struct$coords$s2

ctds_sim = ctds.fwdsim(ctds_struct = ctds_struct, 
                       beta_loc = matrix(1, nrow = 1, ncol = 1), 
                       beta_dir = 0, v0 = box.ind, t0 = 0, tf = 100, 
                       max.steps = 1e3, beta_ar = 0, v0.last = NULL)


document('packages/dsmovetools/')

ctds_obs = ctds.observe(states = ctds_sim$states, times = ctds_sim$times, 
                        t.obs = seq(from = ctds_sim$times[1], 
                                    to = ctds_sim$times[length(ctds_sim$times)], 
                                    length.out = 1800))



# impute path from observation
imputed = ctds.shortest_impute(states = ctds_obs$states, times = ctds_obs$times, 
                               ctds_struct = ctds_struct)

plot.ctds_realization(x = ctds_sim, ctds_struct = ctds_struct,
                      ctds_obs = ctds_obs)

plot.ctds_realization(x = imputed, ctds_struct = ctds_struct,
                      ctds_obs = ctds_obs)

# document('packages/dsmovetools/')


# cctds_nbhd_ll = compileNimble(ctds_nbhd_ll)



if(any(is.na(unlist(imputed)))) {
  stop('Failed imputation')
} else {
  
  o = optim(par = c(0), fn = function(theta) {
    cctds_nbhd_ll(x = imputed$states, durations = imputed$durations, 
                  N = length(imputed$states), 
                  inedges_by_loc = do.call(c, ctds_struct$in_edges_inds), 
                  inloc_start = c(1, 1 + cumsum(ctds_struct$in_degree)), 
                  tolocs_by_edge = ctds_struct$edge_df$to, 
                  fromlocs_by_edge = ctds_struct$edge_df$from, 
                  outedges_by_loc = do.call(c, ctds_struct$out_edges_inds), 
                  loc_start = c(1, 1 + cumsum(ctds_struct$out_degree)), 
                  Xloc = ctds_struct$Xloc, 
                  betaLoc = matrix(theta[1], nrow = 1, ncol = 1), 
                  Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df), ncol = 1), 
                  betaDir = 0, W = ctds_struct$w_ij, betaAR = 0, 
                  nbrlocs_by_loc = do.call(c, nbs.local), 
                  nbrlocs_start = c(1, 1 + cumsum(sapply(nbs.local, length))), 
                  log = TRUE)
  }, method = 'BFGS', control = list(fnscale = -1), hessian = TRUE)
  
  print(o$par)
  
  (solve(-o$hessian))
  chol(solve(-o$hessian))
  
  
}

