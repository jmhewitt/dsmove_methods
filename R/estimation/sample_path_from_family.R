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
  
  # sample path segments
  if(is.null(unif_segment_lengths)) {
    # path segments with independent segment lengths
    seg_inds = sapply(path_fam, function(paths) {
      sample.gumbeltrick(log.p = paths[[1]]$log_weights) 
    })
  } else {
    # path segments s.t. segment lengths are dependent
    seg_inds = sapply(path_fam, function(paths, u) {
      # Sample path s.t. the path lengths are random, but sampled using the same 
      # random variate for each path segment; the specific path segments are 
      # sampled randomly, conditionally on the segment's sampled path length
      # 
      # Parameters:
      #  paths - family of path segments that connect two observations
      #  u - common variate used to draw path lengths across all path segments
      
      # lengths of paths in family
      segment_lengths = sapply(paths[[1]]$path, nrow)
      # aggregate log-mass across paths with identical length
      length_log_prob = aggregate(
        paths[[1]]$log_weights, by = list(len = segment_lengths), FUN = log_sum
      )
      # normalize probabilities, back-transform to linear scale
      length_log_prob$p = exp(length_log_prob$x - log_sum(length_log_prob$x))
      length_log_prob$cdf = cumsum(length_log_prob$p)
      # sample a path length using common random variate across all segments
      seg_len = length_log_prob$len[min(which(length_log_prob$cdf >= u))]
      # identify indices of paths with sampled length
      valid_paths = which(segment_lengths == seg_len)
      # sample a path with given length; return the path's index in paths[[1]]
      valid_paths[
        sample.gumbeltrick(log.p = paths[[1]]$log_weights[valid_paths])
      ]
    }, u = unif_segment_lengths)
  }
  
  # sample arrival times for all locations in each path segment
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
  seg_wts = sapply(1:n_segs, function(seg_ind) {
    # number of timepoints sampled
    n_sampled_times = nrow(
      path_fam[[seg_ind]][[1]]$path[[seg_inds[seg_ind]]]
    ) - 1
    # duration of segment
    dt = diff(seg_start_times[seg_ind + 0:1])
    # sampling weight for path segment
    path_fam[[seg_ind]][[1]]$log_weights[seg_inds[seg_ind]] + 
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
    seg_wts = seg_wts,
    init_prev_loc = init_prev_loc
  )
}
