std_is_wts = function(w, log = TRUE) {
  # standardize importance sampling weights
  #
  # Parameters:
  #  w - importance weights
  #  log - TRUE if w is actually the log of the importance weights
  #
  # Return:
  #  if log==TRUE, the logs of the standardized importance weights
   
  # standardize and return importance weights on linear scale
  if(!log) {
    return(w / sum(w))
  }
  
  # initialize log of sum of importance weights
  log_wts_sum = -Inf
  
  # aggregate importance sampling mass
  for(lw in w) {
    log_wts_sum = log_add(log_a = log_wts_sum, log_b = lw)
  }
  
  w - log_wts_sum
}