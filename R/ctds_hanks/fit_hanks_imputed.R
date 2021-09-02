fit_hanks_imputed = function(
  imputed_trajectories, niter, cell_size, prior
) {
  # Approximate posterior distribution for the 2 parameter CTDS model using 
  # multiply-imputed complete trajectories as the posterior approx. strategy.
  #
  # Parameters:
  #  imputed_trajectories - list of multiply-imputed, completely observed
  #     trajectories. each should be a data.frame containing one row for each 
  #    transition.  each row should have minimal columns "lon_ind", "lat_ind",
  #    "duration", and "movement_dir"
  #  niter - number of MCMC samples to draw
  #  prior - function(betaAR, speed) to evaluate log-prior at parameter values
  
  # log-posterior
  lpfn = function(x, x_ind, theta, trajectory) {
    # Evaluate log-posterior while optionally parameters fixed.
    # The parameter that is easiest to change is "x", which relates to theta 
    # via theta[x_ind] = x.
    
    # merge parameter with theta
    theta[x_ind] = x
    
    # extract parameters
    betaAR = theta[1]
    speed = exp(theta[2])
    
    # likelihood
    dsmovetools2d:::completeObsLL(
      betaAR = betaAR, 
      speed = speed, 
      cell_size = cell_size, 
      durations = trajectory$duration, 
      direction_of_movement = as.integer(
        as.integer(trajectory$direction_of_movement) - 1
      )
    ) + 
    # prior
    prior(betaAR, speed) + 
    # Jacobian to account for transformation in which sampling is done on 
    #   unconstrained space, but prior for theta[2] is defined on [0,\infty)
    jac.log(x = theta[2], log = TRUE)
  }
  
  
  # 
  # analyze each of the imputed trajectories
  #
  
  samples_out = list()
  
  reps = length(imputed_trajectories)
  
  for(rep in 1:reps) {
    
    n_params = 2
    
    # optimize initial parameters
    o = optim(par = c(0,0), fn = function(theta) {
      lpfn(x = theta, x_ind = 1:n_params, theta = theta, 
           trajectory = imputed_trajectories[[rep]])
    }, control = list(fnscale = - 1), hessian = TRUE)
    
    param_vec = o$par
    init_sd = sqrt(diag(solve(-o$hessian)))
    
    # construct MHRW samplers
    paramSamplers = lapply(1:n_params, function(ind) {
      dsmovetools::Mhrw1DAdaptive$new(x = param_vec[ind], sd = init_sd[ind], 
                         lp = lpfn, C = .75, alpha = .44, adaptive = TRUE)
    })
    
    # initialize output
    samples = list(
      param_vec = matrix(nrow = niter, ncol = length(param_vec)),
      lp = numeric(niter)
    )
    colnames(samples$param_vec) = c('betaAR', 'log_speed')
    
    # run sampler
    for(it in 1:niter) {
      
      # update model parameters
      for(i in 1:n_params) {
        update = paramSamplers[[i]]$sample(
          x_ind = i, theta = param_vec, trajectory = imputed_trajectories[[rep]]
        )
        param_vec[i] = update$x
      }
      
      # save samples
      samples$param_vec[it, ] = param_vec
      samples$lp[it] =  lpfn(x = param_vec, x_ind = 1:n_params, 
                             theta = param_vec, 
                             trajectory = imputed_trajectories[[rep]])
      
    }
    
    samples_out[[rep]] = samples
    
  } 
  
  samples_out
}
