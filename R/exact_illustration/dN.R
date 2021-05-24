dN = function(N, s0, sf, t0, tf, shape, rate, log = FALSE) {
  # Parameters:
  #  N - vector of step counts to evaluate density for
  #  s0 - location at which process starts
  #  sf - location at which process ends
  #  t0 - time at which process starts
  #  tf - time at which process ends
  #  shape - shape parameter for prior Gamma distribution on tx. rate
  #  rate - rate parameter for prior Gamma distribution on tx. rate
  #  log - return log density if TRUE
  
  ds = sf - s0
  dt = tf - t0
  
  ld = sapply(N, function(N) {
    p1 = (N + ds[1] + ds[2])/2
    p2 = (N + ds[1] - ds[2])/2
    if(any(p1 != round(p1), p2 != round(p2))) {
      return(-Inf)
    }
    # bridged random walk component of density
    lchoose(n = N, k = p1) + 
    lchoose(n = N, k = p2) -
    N * log(4) + 
    # number of transitions component of density
    lgamma(shape + N) - N * log(1 + rate / dt) -
    lfactorial(N)
  })
  
  if(log) { ld } else { exp(ld) }
}