fit_init = function(init_prev_loc, params, niter, priors, path_fam, dims, 
                    seg_start_times, unif_segment_lengths, 
                    optim.method = 'BFGS') {
  # Parameters:
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
  #  seg_start_times - start time for each path segment
  #  unif_segment_lengths - if NULL, then the length of each path segment is 
  #    sampled independently; otherwise, a quantile between [0,1) used to choose 
  #    each path segment's length wrt. the distribution of the path segment's 
  #    possible lengths, which ensures that the path lengths are similar across 
  #    segments.

  # sample initial path
  path = sample_path_from_family(
    path_fam = path_fam, dims = dims, seg_start_times = seg_start_times, 
    unif_segment_lengths = unif_segment_lengths
  )
  
  # extract initial location
  init_prev_loc = path$init_prev_loc
  
  # read in initial parameters
  param_vec = c(params$beta_ar, params$beta_loc)
  
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

  # optimize inital path and parameters  
  for(it in 1:niter) {
    
    message(paste('Step: ', it, '\n'))
    
    # optimize model parameters
    o = optim(par = param_vec, fn = lpfn, path = path, 
              control = list(fnscale = -1), method = optim.method)
    if(o$convergence == 0) {
      param_vec = o$par
    } else {
      warning(paste('Parameter optimization failed with error code', 
                    o$convergence, sep = ' '))
    }
    
    # update path
    path = update_path(
      init_prev_loc = init_prev_loc, path_fam = path_fam, 
      path_inds = path$path_inds, 
      seg_start_times = seg_start_times, 
      seg_times = path$seg_times, seg_wts = path$seg_wts, 
      dims = dims, betaAR = param_vec[1], beta = param_vec[-1]
    )
    
    print(c(param_vec, lpfn(param_vec, path)))

  }
  
  # return optimized parameters
  samples = list(
    path = path,
    param_vec = param_vec,
    lp =  lpfn(param_vec, path)
  )
}