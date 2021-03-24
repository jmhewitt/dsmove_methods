fit = function(path, init_prev_loc, params, niter, priors, path_fam, dims, 
               init_prop_cov, seg_start_times) {
  # Parameters:
  #  path - list containing initial path segments
  #  init_prev_loc - location immediately before first location in path, for 
  #    AR movement modeling
  #  params - list containing initial model parameters
  #  niter - number of MCMC samples to draw
  #  priors - list containing specification of model priors
  #  path_fam - nested list of path segments that connect to each other. outer
  #    list has family path segments for each timepoint, innermost data has 
  #    the specific paths stored as a matrix, each row of which is an adjacent 
  #    location.
  #  dims - dimensions of spatial domain
  #  init_prop_cov - initial proposal cov. matrix for model parameter RW sampler
  #  seg_start_times - start time for each path segment

  # read in initial parameters
  param_vec = c(params$beta_ar, params$beta_loc)
  
  # initialize output
  samples = list(
    path = vector('list', niter),
    param_vec = matrix(nrow = niter, ncol = length(param_vec)),
    lp = numeric(niter)
  )
  
  #
  # initialize MH-RW sampler
  #  
  
  # wrapper to evaluate log-posterior
  lpfn = function(theta, path) {
    # path likelihood
    ld_path(path_fam = path_fam, path_inds = path$path_inds, 
            init_prev_loc = init_prev_loc, dims = dims, betaAR = theta[1]) + 
    # transition time likelihood
    sum(dexp(x = diff(unlist(path$seg_times)), rate = exp(theta[-1]), 
             log = TRUE)) + 
    # prior
    dnorm(x = theta[1], mean = priors$beta_ar$mean, sd = priors$beta_ar$sd, 
          log = TRUE) + 
    sum(dnorm(x = theta[-1], mean = priors$beta$mean, sd = priors$beta$sd, 
              log = TRUE))
  }
  
  # construct MHRW sampler
  paramSampler = MhrwAdaptive$new(
    x = param_vec, mu = param_vec, Sigma = init_prop_cov, 
    lambda = c(1, 1), lp = lpfn, C = .1, alpha = 1, alpha_star = .44, 
    adaptive = TRUE, adaptation_frequency = Inf
  )
  
  for(it in 1:niter) {
    
    message(paste('Step: ', it, '\n'))
    
    # update model parameters
    update = paramSampler$sample(path = path)
    param_vec = update$x
    
    # update path
    path = update_path(
      init_prev_loc = init_prev_loc, path_fam = path_fam, 
      path_inds = path$path_inds, 
      seg_start_times = seg_start_times, 
      seg_times = path$seg_times, seg_wts = path$seg_wts, 
      dims = dims, betaAR = param_vec[1], beta = param_vec[-1]
    )
    
    # save samples
    samples$path[[it]] = path
    samples$param_vec[it, ] = param_vec
    samples$lp[it] =  lpfn(param_vec, path)
    
    print(c(param_vec, samples$lp[it]))
    
    print(paramSampler$print())
  }
  
  # return final samples
  samples
}