ll_exact = function(epath, durations, beta_loc, beta_dir, beta_ar, 
                    ctds_struct) {
  
  # precompute common infinitesimal generator elements
  outedges_by_loc = do.call(c, ctds_struct$out_edges_inds)
  loc_start = c(1, 1 + cumsum(ctds_struct$out_degree))

  # initialize output
  nll = 0
  
  # reformat/extract AR component
  beta_ar_computational = ifelse(is.null(beta_ar), 0, beta_ar)

  # extract number of edges visited by trajectory
  npath = length(epath)
  
  # build likelihood by looping over trajectory
  for(i in 1:npath) { 
    
    # extract current edge
    e_cur = epath[i]
    
    # location in spatial domain associated with current edge
    v_cur = ctds_struct$edge_df$to[e_cur]
    
    # edges and destination locations that could have been transitioned to
    edges = ctds_struct$out_edges_inds[[v_cur]]
    locs = ctds_struct$edge_df$to[edges]
    
    # get local transition parameters via infinitesimal generator extract
    A_cols = c(e_cur, edges)
    A = local_generator(locs = c(v_cur, locs), row_edges = e_cur,
                        col_edges = A_cols,
                        tolocs_by_edge = ctds_struct$edge_df$to,
                        fromlocs_by_edge = ctds_struct$edge_df$from,
                        outedges_by_loc = outedges_by_loc,
                        loc_start = loc_start, Xloc = ctds_struct$Xloc,
                        betaLoc = beta_loc,
                        Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df),
                                      ncol = 1),
                        betaDir = beta_dir, W = ctds_struct$w_ij,
                        betaAR = beta_ar_computational)
    
    # aggregate transition likelihood
    if(i+1 <= npath) {
      # identify column that denotes the state transitioned to
      tx_ind = which(epath[i+1] == A_cols)
      # aggregate likelihood
      nll = nll + log(A[tx_ind]) - log(sum(A[-1]))
    }
    
    # aggregate holding time likelihood
    if(is.finite(durations[i])) {
      nll = nll + dexp(x = durations[i], rate = -A[1], log = TRUE)
    }
    
  }

  nll
}
