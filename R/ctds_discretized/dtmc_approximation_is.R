dtmc_approximation_is = function(states, times, delta, priors, niter,
                                 gapprox, var.inflate = 1.2) {
  # Importance sampler using proposals from the Gaussian approximation 
  # to the posterior.
  # 
  # Parameters:
  #  states - matrix of observation states/coordinates
  #  times - times at which observations are made
  #  delta - the discretization timestep to use
  #  niter - number of MCMC samples to draw
  #  gapprox - List containing mean and variance of Gaussian approximation
  #  var.inflate - ammount by which to increase the posterior approximation 
  #    variance in order to develop an importance sampler that better covers 
  #    the tails

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
  
  # log-posterior for model
  lp = function(theta) {
    # enforce parameter support
    if(theta > 1/delta) {
      return(-Inf)
    }
    if(theta < 0) {
      return(-Inf)
    }
    # likelihood
    ll(self_tx_prob = 1 - theta * delta) + 
    # prior
    dgamma(x = theta, shape = priors$shape, rate = priors$rate, log = TRUE)
  }
  
  # initialize posterior samples and sampler state
  prop.sd = sqrt(gapprox$post_var * var.inflate)
  samples = rep(NA, niter)
  weights = rep(NA, niter)
  samples_ll = rep(NA, niter)

  # importance sample
  for(it in 1:niter) {
    
    # proposal and likelihood
    prop = rnorm(n = 1, mean = gapprox$post_mean, sd = prop.sd)
    lprop = lp(prop)

    # IS weight
    w = lprop - dnorm(x = prop, mean = gapprox$post_mean, sd = prop.sd, 
                      log = TRUE)
      
    # save samples and weights
    samples[it] = prop
    weights[it] = w
    samples_ll[it] = lprop
  }
  
  # return results
  list(
    samples = samples,
    weights = weights,
    ll = samples_ll
  )
}
