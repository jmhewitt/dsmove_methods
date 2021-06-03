dtmc_ar_approximation = function(states, times, delta, priors, niter = 1e3, 
                                 dims) {
  # Compute Gaussian approximation to the posterior
  #
  # Parameters:
  #  states - matrix of observation states/coordinates
  #  times - times at which observations are made
  #  delta - the discretization timestep to use
  #  niter - number of optimization steps
  #  dims - state space
  
  #
  # construct dummy height variables and information, to align with code
  #
  
  # append dummy variable and dimension information
  dims = c(dims, 1)
  states = cbind(states, 0)
  # height of domain surface
  zsurf = matrix(0, nrow = dims[1], ncol = dims[2])
  # height of vertical layer
  zval = 1
  
  # discretized timesteps
  tseq = seq(from = times[1], to = tail(times,1), by = delta)
  
  # map observations to discretized timesteps
  obs_coords = matrix(NA, nrow = length(tseq), ncol = 3)
  for(ind in 1:length(times)) {
    obs_coords[which(tseq == times[ind]),] = states[ind,]
  }
  
  # uniform mass over neighborhood structure for initial state
  a0_prev_coords = dsmovetools:::TestZConstrainedRookNeighborhood(
    dims = dims, 
    x = states[1,], 
    zfield = zsurf, 
    zvals = zval
  )
  log_a0val = rep(1, nrow(a0_prev_coords))
  log_a0val = log(log_a0val / sum(log_a0val))
  
  # log-likelihood for observations
  ll = function(theta) {
    dsmovetools:::ARFilteredLL(
      a0 = matrix(states[1,], 
                  nrow = nrow(a0_prev_coords), 
                  ncol = ncol(states), byrow = TRUE),
      a0_prev_coords = a0_prev_coords, 
      obs_coords = obs_coords, 
      log_a0val = log_a0val, 
      dims = dims, 
      surface_heights = zsurf, domain_heights = zval, 
      log_self_tx = log(1 - theta[1] * delta), betaAR = theta[2]
    ) + 
      dgamma(x = theta[1], 
             shape = priors$theta['shape'], 
             rate = priors$theta['rate'], 
             log = TRUE) + 
      dnorm(x = theta[2], 
            mean = priors$betaAR['mean'], 
            sd = priors$betaAR['sd'], 
            log = TRUE)
  }
  
  # develop gaussian approximation to posterior
  o = optim(
    par = c(1/delta/2, 0), 
    fn = ll, 
    control = list(fnscale = -1, maxit = niter), 
    method = 'L-BFGS-B', 
    lower = c(0 + 1e-1,-10), 
    upper = c(1/delta - 1e-1, 10), 
    hessian = TRUE
  )
  
  # extract initial gaussian approximation
  gapprox = list(
    post_mean = o$par,
    post_cov = solve(-o$hessian)
  )
  gapprox$post_sd = sqrt(diag(gapprox$post_cov))
  names(gapprox$post_mean) = c('theta', 'betaAR')
  names(gapprox$post_sd) = names(gapprox$post_mean)
  
  # enrich with HPDs 
  gapprox$hpds = cbind(
    lower = qnorm(
      p = .025, 
      mean = gapprox$post_mean, 
      sd = gapprox$post_sd,
      lower.tail = TRUE
    ),
    upper = qnorm(
      p = .025, 
      mean = gapprox$post_mean, 
      sd = gapprox$post_sd,
      lower.tail = FALSE
    )
  )
  rownames(gapprox$hpds) = names(gapprox$post_mean)
  
  # return results
  gapprox
}
