dtmc_approximation = function(states, times, delta, priors, niter) {
  # Parameters:
  #  states - matrix of observation states/coordinates
  #  times - times at which observations are made
  #  delta - the discretization timestep to use
  #  niter - number of MCMC samples to draw

  # number of discretization timesteps between observations
  nsteps = diff(times[1:2]) / delta
  
  # build discrete space dimensions to approximate movement on infinite lattice
  dims = c(5e3,5e3,1)
  # height of domain surface
  zsurf = matrix(0, nrow = dims[1], ncol = dims[2])
  # height of vertical layer
  zval = 1
  
  # translate coordinates to discrete approximation domain
  states.shifted = 500 - states
  
  # add a dummy height coordinate to shifted states
  states.shifted = cbind(states.shifted, 0)
  
  # log-likelihood for transitions
  ll = function(self_tx_prob) {
    res = 0
    for(ind in 1:(nrow(states.shifted)-1)) {
      # forward-filter distribution
      af = dsmovetools:::FFRWLightLogConstrainedSelfTx(
        a0coords = states.shifted[ind,, drop = FALSE], 
        log_a0values = 0, dims = dims, steps = nsteps, 
        surface_heights = zsurf, domain_heights = zval, 
        log_self_tx = log(self_tx_prob)
      )
      # extract log-probability at destination; aggregate likelihood
      res = res + af[
        rowSums(sapply(1:3, function(coord) {
          af[,coord] == states.shifted[ind+1, coord]
        })) == 3,
        4
      ]
    }
    res
  }
  
  # develop gaussian approximation to posterior, at least for the moment
  o = optim(par = 1/delta/2, fn = function(theta) {
    # likelihood
    ll(self_tx_prob = 1 - theta * delta) + 
    # prior
    dgamma(x = theta, shape = priors$shape, rate = priors$rate, log = TRUE)
  }, method = 'Brent', control = list(fnscale = -1), lower = 0, upper = 1/delta, 
  hessian = TRUE)
  
  # extract initial gaussian approximation
  gapprox = list(
    post_mean = o$par,
    post_var = -1/o$hessian
  )
  
  # enrich with HPDs 
  gapprox$lower = qnorm(
    p = .025, 
    mean = gapprox$post_mean, 
    sd = sqrt(gapprox$post_var),
    lower.tail = TRUE
  )
  gapprox$upper = qnorm(
    p = .025, 
    mean = gapprox$post_mean, 
    sd = sqrt(gapprox$post_var),
    lower.tail = FALSE
  )
  
  # return results
  gapprox
}
