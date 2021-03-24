ld_path_segment = function(path_segment, prev_loc, dims, betaAR) {
  # utility function: log-likelihood for sequence of locations in path segment
  #
  # Parameters:
  #  path_segment - matrix of locations visited in path; each row is a loc.
  #  prev_loc - last unique location visited before path_segment
  #  dims - dimensions of spatial domain
  #  betaAR - directional persistence parameter
  
  path_len = nrow(path_segment)
  row_nums = 1:path_len
  
  # ll for paths with no transitions
  if(path_len == 1) {
    return(0)
  }
  
  # log-likelihood for initial transition
  ll = dsmovetools:::TxModelLd(
    prev_loc = prev_loc - 1,
    cur_loc = subset(path_segment, row_nums == 1, drop = FALSE) - 1, 
    dst_loc = subset(path_segment, row_nums == 2, drop = FALSE) - 1,
    dims = dims, betaAR = betaAR
  )
  
  # aggregate mass for remaining transitions
  if(path_len > 2) {
    for(i in 2:(path_len-1)) {
      ll = ll + dsmovetools:::TxModelLd(
        prev_loc = subset(path_segment, row_nums == (i - 1), drop = FALSE) - 1, 
        cur_loc = subset(path_segment, row_nums == i, drop = FALSE) - 1, 
        dst_loc = subset(path_segment, row_nums == (i + 1), drop = FALSE) - 1,
        dims = dims, betaAR = betaAR
      )
    }
  }
  
  ll
}