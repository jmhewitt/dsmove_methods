fit_dtmc_mcmc = function(niter, delta, priors, t0, tf, dims, states, times,
                         univariate = FALSE) {
  # MCMC approximatino to the posterior, using a DTMC likelihood approximation
  #
  # Parameters:
  #  niter - number of MCMC samples to draw
  #  delta - time discretization timestep
  #  priors - list of priors for model
  #  t0 - initial time of simulation
  #  tf - final time of simulation
  #  dims - number of coordinates for discrete space
  #  states - observed locations
  #  times - times of observed locations
  #  univariate - FALSE to estimate speed parameter only (betaAR assumed 0)

  # discretize time
  tseq = seq(from = t0, to = tf, by = delta)
  
  # boilerplate dimension information for likelihood approximation fn.
  xcoords = 1:dims[1]
  ycoords = 1:dims[2]
  surface_heights = numeric(length = length(xcoords) * length(ycoords))
  
  # initialize flattened data structure in discrete time
  oseq = matrix(NA, nrow = length(tseq), ncol = 4)
  colnames(oseq) = c('x_ind', 'y_ind', 'x_coord', 'y_coord')
  
  # flatten data
  for(ind in 1:nrow(states)) {
    tgt_ind = which.min(abs(times[ind] - tseq))
    oseq[tgt_ind, c('x_coord', 'y_coord')] = states[ind,]
    oseq[tgt_ind, 'x_ind'] = which(oseq[tgt_ind,'x_coord'] == xcoords) - 1
    oseq[tgt_ind, 'y_ind'] = which(oseq[tgt_ind,'y_coord'] == ycoords) - 1
  }
  
  # set initial parameters
  param_vec = c(0, .5)
  
  # set optimization specs
  if(univariate) {
    est_param = c(FALSE, TRUE)
  } else {
    est_param = c(TRUE, TRUE)
  }
  n_params = sum(est_param)
  
  # wrapper to evaluate log-posterior
  lpfn = function(x, x_ind, theta) {
    # Evaluate likelihood for model while holding all but one parameter fixed.
    # The parameter that is easiest to change is "x", which relates to theta 
    # via theta[x_ind] = x.
    # 
    # Parameters:
    #   x - univariate parameter for which to evaluate likelihood at
    #   x_ind - index in theta where "x" should be inserted
    #   theta - full parameter vector required to evaluate likelihood
    #   path - latent CTDS path along grid
    
    # merge parameter with theta
    theta[x_ind] = x
    
    # observation likelihood
    dsmovetools2d:::ExactLocFilteredLL(
      init_dsts = matrix(oseq[1,c('x_ind','y_ind')], nrow = 4, ncol = 2, 
                         byrow = TRUE),
      init_srcs = matrix(
        c(oseq[1,c('x_ind','y_ind')] + c(0,1),
          oseq[1,c('x_ind','y_ind')] + c(0,-1),
          oseq[1,c('x_ind','y_ind')] + c(1,0),
          oseq[1,c('x_ind','y_ind')] + c(-1,0)),
        nrow = 4, ncol = 2, byrow = TRUE
      ), 
      init_log_probs = rep(log(1/4), 4),
      obs_x_coords = oseq[,'x_coord'],
      obs_y_coords = oseq[,'y_coord'],
      x_coords = xcoords, 
      y_coords = ycoords, 
      surface_heights = surface_heights, 
      log_self_tx = log(1-delta*exp(theta[2])), 
      betaAR = theta[1], 
      lptrunc = Inf
    ) +
      # prior
      dnorm(x = theta[1], mean = priors$betaAR['mean'], 
            sd = priors$betaAR['sd'], log = TRUE) +
      dgamma(x = exp(theta[2]), shape = priors$theta['shape'], 
             rate =  priors$theta['rate'], log = TRUE) + 
      # Jacobian to account for transformation in which sampling is done on 
      #   unconstrained space, but prior for theta[1] is defined on [0,\infty)
      jac.log(x = theta[2], log = TRUE)
  }
  
  # gaussian approximation to the posterior 
  o = optim(par = param_vec[est_param], fn = function(theta) {
    lpfn(x = theta, x_ind = which(est_param), theta = param_vec)
  }, control = list(fnscale = -1), hessian = TRUE)
  
  #
  # gibbs sampler to explore posterior, initialized around mode
  #
  
  param_vec[est_param] = o$par
  init_sd = sqrt(diag(solve(-o$hessian)))
  
  # construct MHRW samplers
  paramSamplers = lapply(1:length(which(est_param)), function(ind) {
    Mhrw1DAdaptive$new(x = param_vec[which(est_param)[ind]], 
                       sd = init_sd[ind], 
                       lp = lpfn, C = .75, alpha = .44, adaptive = TRUE)
  })
  
  # initialize output
  samples = list(
    param_vec = matrix(nrow = niter, ncol = length(param_vec)),
    lp = numeric(niter)
  )
  colnames(samples$param_vec) = c('betaAR', 'log_theta')
  
  # run sampler
  for(it in 1:niter) {
    
    message(it)
  
    # update model parameters
    for(i in 1:length(which(est_param))) {
      update = paramSamplers[[i]]$sample(
        x_ind = which(est_param)[i], theta = param_vec
      )
      param_vec[which(est_param)[i]] = update$x
    }
    
    # save samples
    samples$param_vec[it, ] = param_vec
    samples$lp[it] =  lpfn(x = param_vec, x_ind = 1:length(param_vec), 
                           theta = param_vec)
    
  }
  
  # package results
  res = list(list(
    o = o,
    samples = samples
  ))
  
}
