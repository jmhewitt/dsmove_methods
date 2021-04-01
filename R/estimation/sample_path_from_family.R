sample_path_from_family = function(path_fam, dims, seg_start_times, 
                                   unif_segment_lengths = NULL) {
  # Sample a path from a list of bridged random walk segments.  Optionally, 
  # ensure that the sampled path has a similar length across all segments.
  # 
  # Parameters:
  #  path_fam - nested list of path segments that connect to each other. outer
  #    list has family path segments for each timepoint, innermost data has 
  #    the specific paths stored as a matrix, each row of which is an adjacent 
  #    location.
  #  dims - dimensions of spatial domain
  #  seg_start_times - start time for each path segment
  #  unif_segment_lengths - if NULL, then the length of each path segment is 
  #    sampled independently; otherwise, a quantile between [0,1) used to choose 
  #    each path segment's length wrt. the distribution of the path segment's 
  #    possible lengths, which ensures that the path lengths are similar across 
  #    segments.
  
  # number of path segments to sample
  n_segs = length(path_fam)
  
  #
  # sample path segments
  #
  
  # path segments sampled independently
  seg_inds = sapply(path_fam, function(paths) {
    sample.gumbeltrick(log.p = paths[[1]]$lp)
  })
  
  # sampling probability for each segment
  seg_ind_wts = sapply(1:n_segs, function(i) {
    # probability of generating sample 
    path_fam[[i]][[1]]$lp[seg_inds[i]]
  })
  
  # sample (initial) arrival times for all locations in each path segment
  seg_times = lapply(1:n_segs, function(seg_ind) {
    # segment length
    seg_len = nrow(path_fam[[seg_ind]][[1]]$path[[seg_inds[seg_ind]]])
    # segment time range
    t0 = seg_start_times[seg_ind]
    tf = seg_start_times[seg_ind + 1]
    # random times
    c(t0, sort(runif(n = seg_len - 1, min = t0, max = tf)))
  })
  
  # likelihood of sampling each segment under independence proposal scheme
  seg_time_wts = sapply(1:n_segs, function(seg_ind) {
    # number of timepoints sampled
    n_sampled_times = nrow(
      path_fam[[seg_ind]][[1]]$path[[seg_inds[seg_ind]]]
    ) - 1
    # duration of segment
    dt = diff(seg_start_times[seg_ind + 0:1])
    # sampling likelihood for times
    lfactorial(n_sampled_times) - n_sampled_times * log(dt)
  })
  
  # randomly initialize a last location for first observation
  init_nbhd = dsmovetools:::TestRookNeighborhood(
    dims = dims, x = path_fam[[1]][[1]]$path[[1]][1, , drop = FALSE] - 1
  ) + 1
  init_prev_loc = init_nbhd[sample(nrow(init_nbhd), size = 1), , drop = FALSE]
  
  # package sample
  list(
    path_inds = seg_inds,
    seg_times = seg_times,
    seg_ind_wts = seg_ind_wts,
    seg_time_wts = seg_time_wts,
    init_prev_loc = init_prev_loc
  )
}
