fit_exact = function(ctds_struct, trajectory, beta_loc, beta_dir, beta_ar) {
  
  # load observations
  ctds_trajectory = readRDS(trajectory)
  
  # precompute common infinitesimal generator elements
  outedges_by_loc = do.call(c, ctds_struct$out_edges_inds)
  loc_start = c(1, 1 + cumsum(ctds_struct$out_degree))
  
  # log-likelihood for completely observed trajectory
  ll = function(beta_loc, beta_dir, beta_ar) {
    
    # initialize output
    nll = 0
    
    # reformat/extract AR component
    beta_ar_computational = ifelse(is.null(beta_ar), 0, beta_ar)

    # extract number of edges visited by trajectory
    npath = length(ctds_trajectory$edges)

    # build likelihood by looping over trajectory
    for(i in 2:npath) {

      # edge from which trajectory moved
      e_prev = ctds_trajectory$edges[i-1]
      
      # location in spatial domain associated with previous edge
      v_prev = ctds_struct$edge_df$to[e_prev]
      
      # edges that could have been transitioned to
      edges = ctds_struct$out_edges_inds[[v_prev]]
      
      # destination locations associated with transition edges
      locs = ctds_struct$edge_df$to[edges]
      
      # get local transition parameters via infinitesimal generator extract
      A_cols = c(e_prev, edges)
      A = local_generator(locs = c(v_prev, locs), row_edges = e_prev,
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
      
      # identify column that denotes the state transitioned to 
      tx_ind = which(ctds_trajectory$edges[i] == A_cols)
      
      # aggregate likelihood
      nll = nll + 
        # holding time
        dexp(x = ctds_trajectory$durations[i-1], rate = -A[1], log = TRUE) + 
        # transition probability
        log(A[tx_ind]) - log(sum(A[-1]))
    }

    nll
  }

  # parameter MLE
  o = optim(par = c(beta_loc, beta_dir, beta_ar), function(theta) {
    ll(beta_loc = matrix(theta[1:length(beta_loc)], ncol = 1), 
       beta_dir = matrix(theta[length(beta_loc) + 1:length(beta_dir)], ncol = 1), 
       beta_ar =  theta[length(theta)])
  }, method = 'BFGS', control = list(fnscale = -1), hessian = TRUE)
  
  o
}
