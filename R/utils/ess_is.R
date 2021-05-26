ess_is = function(w) {
  # effective sample size for importance samples, given standardized weights.
  # Follows Givens and Hoeting eq's. 6.29 and 6.31.
  #
  # Parameters:
  #  w - standardized importance weights
  
  n = length(w)
  cvsq = n * sum(w^2) - 1
  n / (1 + cvsq)
}
