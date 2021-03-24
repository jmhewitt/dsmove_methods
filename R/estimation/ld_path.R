ld_path = function(path_fam, path_inds, init_prev_loc, dims, betaAR) {
  # utility function: log-likelihood for an entire path sequence
  #
  # Parameters:
  #  path_fam - nested list of path segments that connect to each other. outer
  #    list has family path segments for each timepoint, innermost data has 
  #    the specific paths stored as a matrix, each row of which is an adjacent 
  #    location.
  #  path_inds - vector specifying which path segments from path_fam are used to 
  #    define the path
  #  init_prev_loc - location immediately before first location in path, for 
  #    AR movement modeling
  #  dims - dimensions of spatial domain
  #  betaAR - directional persistence parameter
  
  # number of segments in path
  n_segs = length(path_inds)
  
  # must always be aware of last location transitioned from
  prev_loc = init_prev_loc
  
  # initialize log-likelihood result
  ld = 0
  
  # aggregate likelihood over path segments
  for(seg_ind in 1:n_segs) {
    # extract path segment
    seg = path_fam[[seg_ind]][[1]]$path[[path_inds[seg_ind]]]
    seg_len = nrow(seg)
    # aggregate segment likelihood
    ld = ld + ld_path_segment(path_segment = seg, prev_loc = prev_loc, 
                              dims = dims, betaAR = betaAR)
    # update prev_loc if segment changed locations
    if(seg_len > 1) {
      row_nums = 1:seg_len
      prev_loc = seg[seg_len - 1, , drop = FALSE]
    }
  }
  
  ld
}
