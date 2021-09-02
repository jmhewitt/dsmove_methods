# tools for working with log-posteriors evaluated on gridded domains

log_C_gridded = function(x, y, log_fxy, dx, dy) {
  # return log(C) where C is the normalization constant for a density f(x,y),
  # and log(f(x,y)) is provided on a regular grid with spacing dx, dy.
  # 
  # Parameters:
  #  x - vector of x coordinates associated with log(f(x,y)) values
  #  y - vector of y coordinates associated with log(f(x,y)) values
  #  log_fxy - vector of log(f(x,y)) values
  #  dx - distance between grid cells in x
  #  dy - distance between grid cells in y
  
  log_sum(log_fxy + log(dx) + log(dy))
}

log_density_gridded = function(x, y, ld, dx, dy) {
  # return the normalized log-density values for a log-density surface 
  # log(f(x,y)) evaluated on a grid when provided information about the grid 
  # and prior
  # 
  # Parameters: 
  #  x - vector of x coordinates associated with log-density values
  #  y - vector of y coordinates associated with log-density values
  #  ll - vector of log-density values
  #  dx - distance between grid cells in x
  #  dy - distance between grid cells in y

  ld - log_C_gridded(x = x, y = y, log_fxy = ld, dx = dx, dy = dy)
}

log_marginal_x = function(x, y, ld, dy) { 
  # evaluate the marginal log(f(x)) = log(\int f(x,y) dy ) when provided the 
  # joint density log(f(x,y)) evaluated on a grid
  #
  # Parameters:
  #  x - vector of x coordinates associated with log-density values
  #  y - vector of y coordinates associated with log-density values
  #  ld - vector of log-density values
  #  dy - distance between grid cells in y
  
  data.frame(x = x, y = y, ld = ld) %>% 
    group_by(x) %>% 
    summarise(
      ld = log_sum(ld + log(dy))
    ) %>% 
    ungroup()
}

mean.marginal = function(x, p) {
  # Mean when given marginal parameter values and probabilities
  # 
  # Parameters:
  #  x - grid of values
  #  p - probability of value
  
  sum(x * p)
}

var.marginal = function(x, p) {
  # Variance when given marginal parameter values and probabilities
  # 
  # Parameters:
  #  x - grid of values
  #  p - probability of value
  
  mu = mean.marginal(x = x, p = p)
  
  sum((x-mu)^2 * p)
}

hpd.marginal = function(x, p, level = .95) {
  # HPD interval when given marginal parameter values and probabilities
  # 
  # Parameters:
  #  x - grid of values
  #  p - probability of value
  
  # ensure input is ordered
  o = order(x)
  x = x[o]
  p = p[o]
  
  cdf = cumsum(p)
  
  lwr.ind = 1
  upr.ind = min(which(cdf - cdf[lwr.ind] >= level))
  len = x[upr.ind] - x[lwr.ind]
  for(i in 2:length(x)) {
    valid = cdf - cdf[i] >= level
    if(any(valid)) {
      upr.ind.prop = min(which(valid))
      len.prop = x[upr.ind.prop] - x[i]
      if(len.prop < len) {
        lwr.ind = i
        upr.ind = upr.ind.prop
      }
    }
  }
  
  x[c(lwr.ind, upr.ind)]
}