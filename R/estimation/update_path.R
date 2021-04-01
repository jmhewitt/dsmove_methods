update_path = function(init_prev_loc, path_fam, path_inds, seg_start_times,
                       seg_times, seg_ind_wts, seg_time_wts, dims, betaAR, 
                       beta) {
  # Metropolis-Hastings update of a segment of the latent trajectory
  #
  # Parameters:
  #  init_prev_loc - last unique location visited before path_segment
  #  path_fam - nested list of path segments that connect to each other. outer
  #    list has family path segments for each timepoint, innermost data has 
  #    the specific paths stored as a matrix, each row of which is an adjacent 
  #    location.
  #  path_inds - vector of path segments in path_fam that define the path
  #  seg_start_times - start time for each path segment
  #  seg_times - list in which each element is a vector of times at which 
  #    transitions occur
  #  seg_wts - vector of sampling weights for each selected segment
  #  dims - dimensions of spatial domain
  #  betaAR - directional persistence parameter
  #  beta - (log) speed parameter
  
  n_segs = length(path_fam)
  
  # must always be aware of last location transitioned from
  prev_loc = init_prev_loc
  
  # update path segments incrementally
  for(seg_ind in 1:n_segs) {
    
    # family of segments and time range for current period
    pf = path_fam[[seg_ind]][[1]]
    t0 = seg_start_times[seg_ind]
    tf = seg_start_times[seg_ind + 1]
    
    # extract path segment
    prop_seg = path_inds[seg_ind]
    
    # number of timepoints to sample for proposal
    n_sampled_times = nrow(pf$path[[prop_seg]]) - 1
    # random times
    prop_times = c(t0, sort(runif(n = n_sampled_times, min = t0, max = tf)))
    
    # sampling likelihood for times
    prop_time_wts = lfactorial(n_sampled_times) - n_sampled_times * log(tf - t0)
    
    # update path (no accept/reject needed since proposal is full-conditional)
    seg_times[[seg_ind]] = prop_times
    seg_time_wts[seg_ind] = prop_time_wts
    
    # update prev_loc if segment changed locations
    seg_len = nrow(pf$path[[path_inds[seg_ind]]])
    if(seg_len > 1) {
      seg = path_fam[[seg_ind]][[1]]$path[[path_inds[seg_ind]]]
      row_nums = 1:seg_len
      prev_loc = seg[seg_len - 1, , drop = FALSE]
    }
    
  }
  
  list(
    path_inds = path_inds,
    seg_times = seg_times, 
    seg_ind_wts = seg_ind_wts,
    seg_time_wts = seg_time_wts,
    init_prev_loc = init_prev_loc
  )
}