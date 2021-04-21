fit_hanks = function(params, niter, priors, states, times, dims, reps) {
  # Parameters:
  #  params - list containing initial model parameters
  #  niter - number of MCMC samples to draw
  #  priors - list containing specification of model priors
  #  states - matrix of observation states/coordinates
  #  times - times at which observations are made
  #  dims - dimensions of spatial domain
  
  # following package use from help('ctmcmove')
  
  # targets::tar_load(sim_obs)
  # targets::tar_load(sim_params)
  # library(ctmcmove)
  # 
  # states = sim_obs[[1]]$states
  # times = sim_obs[[1]]$times
  # dims = sim_params$dims
  
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

  # log-density for trajectory
  ld_path = function(theta, states, durations) {
    
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
    dnorm(x = theta[1], mean = priors$beta_ar$mean, sd = priors$beta_ar$sd,
          log = TRUE) +
    sum(dnorm(x = theta[-1], mean = priors$beta$mean, sd = priors$beta$sd,
              log = TRUE))
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
    colnames(samples$param_vec) = c('beta_ar', 'beta_loc')
    
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
      # samples$path[[it]] = path
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
