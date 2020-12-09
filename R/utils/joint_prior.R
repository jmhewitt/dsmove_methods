joint_prior = function(beta_loc, beta_ar, ctds_struct, penalty_rate = 1/.5,
                       log = FALSE) {
  #
  # Parameters:
  #  penalty_rate - 1/expected number of reversals per unit time to anticipate
  
  # precompute common infinitesimal generator elements
  outedges_by_loc = do.call(c, ctds_struct$out_edges_inds)
  loc_start = c(1, 1 + cumsum(ctds_struct$out_degree))
  
  # index of central-most coordinate
  ind.center = which.min(
    apply(sweep(x = ctds_struct$coords, MARGIN = 2, 
                STATS = colMeans(ctds_struct$coords), FUN = '-'),
          1, function(r) sum(r^2))
  )
  
  # arbitrarily select an edge, which implies a previous location
  e_cur = ctds_struct$in_edges_inds[[ind.center]][1]
  # get outbound edges from the central-most coordinate
  edges = ctds_struct$out_edges_inds[[ind.center]]
  
  # identify current location
  v_cur = ctds_struct$edge_df$to[e_cur]
  # identify locations associated with outbound edges
  locs = ctds_struct$edge_df$to[edges]
  
  # identify the edge associated with a reversal of movement
  reversal_ind = which(
    ctds_struct$edge_df$from[e_cur] == ctds_struct$edge_df$to[edges]
  )
  
  # basis for joint prior
  expected_reversals = function(beta_ar, beta_loc) {
    
    mapply(function(beta_ar, beta_loc) {
      # get local transition parameters via infinitesimal generator extract
      A_cols = c(e_cur, edges)
      A = local_generator(locs = c(v_cur, locs), row_edges = e_cur,
                          col_edges = A_cols,
                          tolocs_by_edge = ctds_struct$edge_df$to,
                          fromlocs_by_edge = ctds_struct$edge_df$from,
                          outedges_by_loc = outedges_by_loc,
                          loc_start = loc_start, Xloc = ctds_struct$Xloc,
                          betaLoc = matrix(beta_loc),
                          Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df),
                                        ncol = 1),
                          betaDir = matrix(0), W = ctds_struct$w_ij,
                          betaAR = beta_ar)
      
      # expected number of transitions per unit time
      -A[1] *
      # scaled by expected proportion or reversals
      A[-1][reversal_ind] / sum(A[-1])
      
    }, beta_ar, beta_loc)
  }
  
  dexp(x = expected_reversals(beta_ar = beta_ar, beta_loc = beta_loc), 
       log = log, rate = penalty_rate)
}