fit_hanks = function(params, niter, priors, states, times, dims, reps) {
  # Approximate posterior distribution for the 2 parameter CTDS model
  #
  # Parameters:
  #  params - list containing initial model parameters
  #  niter - number of MCMC samples to draw
  #  priors - list containing specification of model priors
  #  states - matrix of observation states/coordinates
  #  times - times at which observations are made
  #  dims - dimensions of spatial domain
  
  # following package use from help('ctmcmove')
  
  # extract coordinate/time triples
  xyt = cbind(states, times)
  colnames(xyt) = c('x', 'y', 't')
  xy = xyt[, c('x', 'y')]
  x = xyt[, 'x']
  y = xyt[, 'y']
  t = xyt[, 't']
  
  
  #
  # fit functional movement model to telemetry data
  #
  
  ## Define the knots of the spline expansion.
  ##
  ## Problems with fitting the functional movement model can often be fixed by
  ## varying the spacing of the knots.
  knots = seq(min(t),max(t),by=1)
  ## create B-spline basis vectors used to approximate the path
  b=create.bspline.basis(c(min(t),max(t)),breaks=knots,norder=3)
  ## define the sequence of times on which to sample the imputed path
  tpred=seq(min(t),max(t),by=1/4)
  
  ## Fit latent Gaussian model using MCMC
  out=mcmc.fmove(xy,t,b,tpred,QQ="CAR",n.mcmc=400,a=1,r=1,num.paths.save=reps)
  
  
  #
  # create rasters
  #
  
  X = raster(
    matrix(1, nrow = dims[1], ncol = dims[2]), 
    xmn = .5,
    xmx = dims[1],
    ymn = .5,
    ymx = dims[2]
  )
  
  crs(X) = "+proj=utm +zone=32"
  
  coords = round(coordinates(X))

  # log-density for complete trajectory
  ld_path = function(theta, states, durations) {
    # Parameters
    #   theta - model parameters on the transformed scale, so all elements of 
    #     theta are supported on \mathbb{R}
    
    # initial state of path
    prev_loc = states[1,]
    
    # initialize log-density
    ld = 0
    
    # aggregate log-density across transitions
    for(ind in 2:(nrow(states)-1)) {
      # evaluate grid transition density
      dtx = dsmovetools:::TxModelLd(
        cur_loc = states[ind,], prev_loc = prev_loc, dims = dims, 
        betaAR = theta[1], dst_loc = states[ind+1,]
      )
      # skip transitions that are invalid due to the linear interpolation
      if(!is.finite(dtx)) {
        dtx = 0
      }
      # aggregate grid transition density
      ld = ld + dtx + 
        # and duration density
        dexp(x = durations[ind], rate = exp(theta[-1]), log = TRUE)
      
      # update current state
      prev_loc = states[ind,]
    }
    
    ld
  }
  
  # wrapper to evaluate log-posterior
  lpfn = function(x, x_ind, theta, states, durations) {
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
    
    # path likelihood
    ld_path(theta = theta, states = states, durations = durations) + 
    # prior
    dnorm(x = theta[1], mean = priors$betaAR['mean'], sd = priors$betaAR['sd'],
          log = TRUE) +
    dgamma(x = exp(theta[-1]), shape = priors$theta['shape'], 
           rate =  priors$theta['rate'], log = TRUE) + 
    # Jacobian to account for transformation in which sampling is done on 
    #   unconstrained space, but prior for theta[-1] is defined on [0,\infty)
    jac.log(x = theta[-1], log = TRUE)
  }
  
  
  # 
  # analyze each of the imputed trajectories
  #
  
  samples_out = list()
  
  for(rep in 1:reps) {
    
    # discretize and format path
    path=out$pathlist[[rep]]
    ctmc=path2ctmc(path$xy,path$t, X, method="LinearInterp")
    ctmc$states = coordinates(coords[ctmc$ec,])
    
    # read in initial parameters
    param_vec = c(params$beta_ar, params$beta_loc)
    n_params = length(param_vec)
    
    # optimize initial parameters
    o = optim(par = param_vec, fn = function(theta) {
      lpfn(x = theta, x_ind = 1:n_params, theta = theta, states = ctmc$states, 
           durations = ctmc$rt)
    }, control = list(fnscale = - 1), hessian = TRUE)
    if(o$convergence == 0) {
      param_vec = o$par
      init_sd = sqrt(diag(solve(-o$hessian)))
    }
    
    # construct MHRW samplers
    paramSamplers = lapply(1:n_params, function(ind) {
      Mhrw1DAdaptive$new(x = param_vec[ind], sd = init_sd[ind], 
                         lp = lpfn, C = .75, alpha = .44, adaptive = TRUE)
    })
    
    # initialize output
    samples = list(
      param_vec = matrix(nrow = niter, ncol = length(param_vec)),
      lp = numeric(niter),
      ctmc = ctmc
    )
    colnames(samples$param_vec) = c('betaAR', 'log_theta')
    
    # run sampler
    for(it in 1:niter) {
      
      # update model parameters
      for(i in 1:n_params) {
        update = paramSamplers[[i]]$sample(
          x_ind = i, theta = param_vec, states = ctmc$states, 
          durations = ctmc$rt
        )
        param_vec[i] = update$x
      }
      
      # save samples
      samples$param_vec[it, ] = param_vec
      samples$lp[it] =  lpfn(x = param_vec, x_ind = 1:n_params, 
                             theta = param_vec, states = ctmc$states, 
                             durations = ctmc$rt)
      
    }
    
    samples_out[[rep]] = samples
    
  } 
  
  
  list(
    samples = samples_out,
    tstep = mean(diff(times))
  )
}
