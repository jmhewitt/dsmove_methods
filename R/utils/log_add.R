log_add = function(log_a, log_b) {
  # evaluate log(c) = log(a + b) when given log(a) and log(b)
  if(all(is.infinite(log_a), is.infinite(log_b))) {
    return(-Inf)
  }
  x = log_a - log_b
  exp_x = exp(x)
  if(is.infinite(exp_x)) {
    warning('exp(log(a) - log(b)) is infinite')
    return(log_a)
  }
  if(exp_x == 0) {
    warning('exp(log(a) - log(b)) is 0')
    return(log_b)
  }
  return(log_b + log(1 + exp_x))
}

log_sum = function(x) {
  # evaluate log(c) = log(x1 + ... + xn) when given log(x1), ..., log(xn)
  
  log_res = x[1]
  len = length(x)
  if(len > 1) {
    for(i in 2:len) {
      log_res = log_add(log_res, x[i])
    }
  }
  log_res
}