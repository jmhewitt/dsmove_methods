sample_path_segments = function(x0, xf, t0, tf, max_tx_rate, high_quantile, 
                                dims, tgt_max_steps, computational_max_steps, 
                                surface_heights, domain_heights, nsegments) {
  # Sample bridged random walks on discrete spaces when given initial and final
  # conditions, and information to restrict the maximum number of transitions.
  
  # get maximum number of steps to take subject to reachability constraints
  pp_lambda = (tf - t0) * max_tx_rate
  tgt_max_steps = qpois(p = high_quantile, lambda = pp_lambda)
  
  # sample paths
  paths = dsmovetools:::SampleConstrainedBridgedRWPathFamily(
    a0coords = x0 - 1, dstcoords = xf - 1, dims = dims, steps = tgt_max_steps, 
    max_steps = computational_max_steps, 
    surface_heights = as.numeric(surface_heights), 
    domain_heights = domain_heights,
    n = nsegments
  )
  
  # unwrap 0-based indexing from c++
  paths$path = lapply(paths$path, function(path) { path + 1 })
  
  # log-probability of generating path
  paths$lp = paths$log_weights
  path_len = sapply(paths$path, function(path) nrow(path) - 1)
  paths$lp_length = dpois(x = path_len, lambda = pp_lambda, log = TRUE)
  
  # initialize path sampling weight for paths; sampling weight will be inversely
  # proportional to number of times path appears in set
  paths$log_weights = rep(0, length(paths$log_weights))
  
  #
  # aggregate duplicate paths to reduce storage overhead
  #
  
  # primary loop to aggregate duplicate path segments
  path_ind = 1
  while(path_ind <= length(paths$path)) {
    # secondary loop to look for matches
    path_ind2 = path_ind + 1
    while(path_ind2 <= length(paths$path)) {
      # process match
      if(identical(paths$path[[path_ind]], paths$path[[path_ind2]])) {
        # aggregate log mass
        paths$log_weights[path_ind] = log_add(
          log_a = paths$log_weights[path_ind], 
          log_b = paths$log_weights[path_ind2]
        )
        # remove duplicate path
        paths$path = paths$path[-path_ind2]
        paths$log_weights = paths$log_weights[-path_ind2]
        paths$lp = paths$lp[-path_ind2]
        paths$lp_length = paths$lp_length[-path_ind2]
        # realign path_ind2 counter with shift
        path_ind2 = path_ind2 - 1
      }
      # increment secondary loop counter
      path_ind2 = path_ind2 + 1
    }
    # increment primary loop counter
    path_ind = path_ind + 1
  }
  
  # standardize log weights
  paths$log_weights = paths$log_weights - log_sum(paths$log_weights)
  
  paths
}