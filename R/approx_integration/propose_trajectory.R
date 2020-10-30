propose_trajectory = function(states, times, ctds_domain, segments, beta_loc) {
  # Use segments to propose a trajectory consistent with observations
  #
  # Parameters:
  #  max_lambda - prior information: upper bound on maximum transition rate
  #  p_max - upper bound on Poisson transition tail
  #  nsamples - number of segments to simulate between each timepoint
  #    (duplicates will be upweighted)
  
  # randomly sample an initial edge
  in_edges = ctds_domain$in_edges_inds[[states[1]]]
  epath = ifelse(length(in_edges) > 1, sample(x = in_edges, size = 1), 
                 in_edges)
  
  # set initial time
  tpath = times[1]
  
  # build onto the path in edge-space
  for(i in 1:(length(states)-1)) {
    
    # extract sampling weights and indexes for segment
    pathwts = do.call(rbind, lapply(1:length(segments[[i]]), function(k) {
      s = segments[[i]][[k]]
      npaths = ifelse(is.null(s$paths), 1, nrow(s$paths))
      data.frame(lengthind = k, pathind = 1:npaths, w = s$weights)
    }))
    
    # sample a path segment
    pathrow = ifelse(nrow(pathwts) > 1, 
                     sample(x = nrow(pathwts), size = 1, prob = pathwts$w),
                     1)
    
    # extract path (in edge space)
    r = pathwts[pathrow, ]
    p = unlist( segments[[i]][[r$lengthind]]$paths[r$pathind, ] )
    
    # extend path
    if(!is.null(p)) {
      epath = c(epath, p)
      # TODO: update sampling s.t. it is from a non-homogeneous PP
      tpath = c(tpath, 
                sort(runif(n = length(p), min = times[i], max = times[i+1])))
    }
    
  }
  
  # convert path to state-space
  spath = ctds_domain$edge_df$to[epath]
  
  # TODO: get a version of this where we a) only return the segment id's, and 
  #       b) we only sample individual segments
  
  # # TODO: eventually remove this debugging check
  # # verify that the imputed path is compatible with the observations
  # obs.imputed = ctds.observe(states = spath, times = tpath, t.obs = times)
  # stopifnot(
  #   all(identical(obs.imputed$times, times),
  #       identical(obs.imputed$states, states))
  # )
  
  # package and return trajectory
  list(states = spath, times = tpath)
}