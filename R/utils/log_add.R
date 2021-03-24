log_add = function(log_a, log_b) {
  x = log_a - log_b
  exp_x = exp(x)
  if(is.infinite(exp_x)) {
    return(log_a)
  }
  if(exp_x == 0) {
    return(log_b)
  }
  return(log_b + log(1 + exp_x))
}
