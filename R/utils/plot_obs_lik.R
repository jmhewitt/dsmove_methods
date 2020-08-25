plot_obs_lik = function(ctds_struct, obs, plot_dir, beta_loc_seq, beta_ar_seq,
                        beta_loc, beta_ar) {
  
  # load observations
  ctds_obs = readRDS(obs)
  
  # extract timestep
  tstep = round(diff(ctds_obs$times[1:2]), 2)

  
  #
  # approximate likelihood
  #
  
  # impute path from observation
  imputed = ctds.shortest_impute(states = ctds_obs$states, 
                                 times = ctds_obs$times, 
                                 ctds_struct = ctds_struct)

  
  #
  # evaluate likelihood at grid of locations
  #
  
  gr = expand.grid(beta_loc = beta_loc_seq, beta_ar = beta_ar_seq)
  
  ll = apply(gr, 1, function(r) {
    cctds_nbhd_ll(x = imputed$states, durations = imputed$durations, 
                  N = length(imputed$states), 
                  inedges_by_loc = do.call(c, ctds_struct$in_edges_inds), 
                  inloc_start = c(1, 1 + cumsum(ctds_struct$in_degree)), 
                  tolocs_by_edge = ctds_struct$edge_df$to, 
                  fromlocs_by_edge = ctds_struct$edge_df$from, 
                  outedges_by_loc = do.call(c, ctds_struct$out_edges_inds), 
                  loc_start = c(1, 1 + cumsum(ctds_struct$out_degree)), 
                  Xloc = ctds_struct$Xloc, 
                  betaLoc = matrix(r['beta_loc'], nrow = 1, ncol = 1), 
                  Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df), ncol = 1), 
                  betaDir = 0, W = ctds_struct$w_ij, betaAR = r['beta_ar'], 
                  nbrlocs_by_loc = do.call(c, ctds_struct$nbs.local), 
                  nbrlocs_start = c(1, 1 + cumsum(sapply(ctds_struct$nbs.local, 
                                                         length))), 
                  log = TRUE)
  })
  
  # likelihood surface plot
  pl = ggplot(data.frame(gr, ll = ll), 
              aes(x = beta_loc, y = beta_ar, z = (ll), fill = (ll))) + 
    # likelihood raster, with contours
    geom_raster() + 
    geom_contour(col = 'grey60') + 
    geom_text_contour(col = 'white') +
    # true parameters
    geom_point(data = data.frame(beta_loc = beta_loc[1], beta_ar = beta_ar,
                                 ll = 0), 
               pch = 4, col = 'white') + 
    # formatting and labels
    scale_fill_viridis('', direction = -1) + 
    xlab(expression(beta[loc])) + 
    ylab(expression(beta[AR])) + 
    ggtitle('Log-likelihood surface (Truth at x)') +
    theme_few() + 
    theme(panel.border = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5))
  
  # save plot
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(pl, filename = file.path(plot_dir, paste('sim_lik_tstep_', tstep, 
                                                  '.pdf', sep = '')))
 
  pl
}