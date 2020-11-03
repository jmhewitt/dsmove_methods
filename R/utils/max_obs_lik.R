max_obs_lik = function(ctds_struct, obs) {
  
  # load observations
  dat = readRDS(obs)  
  ctds_obs = dat$obs
  
  # extract timestep
  tstep = round(mean(diff(ctds_obs$times)), 2)
  
  #
  # approximate likelihood
  #

  # impute path from observation
  imputed = ctds.shortest_impute(states = ctds_obs$states, 
                                 times = ctds_obs$times, 
                                 ctds_struct = ctds_struct)

  # find MLE
  o = optim(par = c(0,0), function(theta) {
    cctds_nbhd_ll(x = imputed$states, durations = imputed$durations, 
                  N = length(imputed$states), 
                  inedges_by_loc = do.call(c, ctds_struct$in_edges_inds), 
                  inloc_start = c(1, 1 + cumsum(ctds_struct$in_degree)), 
                  tolocs_by_edge = ctds_struct$edge_df$to, 
                  fromlocs_by_edge = ctds_struct$edge_df$from, 
                  outedges_by_loc = do.call(c, ctds_struct$out_edges_inds), 
                  loc_start = c(1, 1 + cumsum(ctds_struct$out_degree)), 
                  Xloc = ctds_struct$Xloc, 
                  betaLoc = matrix(theta[-1], nrow = 1, ncol = 1), 
                  Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df), ncol = 1), 
                  betaDir = 0, W = ctds_struct$w_ij, betaAR = theta[1], 
                  nbrlocs_by_loc = do.call(c, ctds_struct$nbs.local), 
                  nbrlocs_start = c(1, 1 + cumsum(sapply(ctds_struct$nbs.local, 
                                                         length))), 
                  log = TRUE)
  }, method = 'BFGS', control = list(fnscale = -1, maxit = 500), hessian = TRUE)
  
  # extract uncertainty estimate
  cov = solve(-o$hessian)
  
  # summarize fit
  est = data.frame(est = o$par, se = sqrt(diag(cov)), 
                   truth = unlist(dat$params[c('beta_ar', 'beta_loc')]),
                   param = c('beta_ar', 'beta_loc'),
                   weibull_shape = dat$params$weibull_shape,
                   tstep = tstep) %>%
    dplyr::mutate(lwr = est - 1.96 * se, upr = est + 1.96 * se,
                  covered = (lwr <= truth) & (truth <= upr))
  
  # package results
  list(
    o = o,
    est = est
  )
}