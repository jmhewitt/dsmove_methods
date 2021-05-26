hpd_is = function(x, w, prob = 0.95) {
  # HPD intervals based on importance samples
  # 
  # Parameters:
  #  x - importance samples
  #  w - standardized importance weights
  #  prob - probability level

  if(sum(w) != 1) {
    stop('Importance weights are not standardized.')
  }  
  
  # sort samples and weights
  o = order(x)
  x.sorted = x[o]
  cdf.sorted = cumsum(w[o])
  
  # initialize an interval with appropriate mass
  lwr = x.sorted[1]
  upr = x.sorted[min(which((cdf.sorted - cdf.sorted[1]) >= prob))]
  
  # compare against all other intervals with appropriate mass
  for(i in 2:length(x.sorted)) {
    # search for another interval with appropriate mass
    lwr.prop = x.sorted[i]
    inds = which((cdf.sorted - cdf.sorted[i]) >= prob)
    # replace current interval with new interval if shorter
    if(length(inds)>0) {
      upr.prop = x.sorted[min(inds)]
      if((upr.prop - lwr.prop) < (upr - lwr)) {
        lwr = lwr.prop
        upr = upr.prop
      }
    }
  }
  
  # return HPD interval
  c(lower = lwr, upper = upr)
}