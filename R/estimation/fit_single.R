fit_single = function(init_prev_loc = NULL, params, niter, priors, path_fam, dims, 
                    unif_segment_lengths = NULL, init_sd, seg_start_times) {
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
  #  init_prop_cov - initial proposal cov. matrix for model parameter RW sampler
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
  n_params = length(param_vec)
  
  # initialize output
  samples = list(
    path = path,
    param_vec = matrix(nrow = niter, ncol = length(param_vec)),
    lp = numeric(niter)
  )
  colnames(samples$param_vec) = c('beta_ar', 'beta_loc')
  
    
  #
  # initialize MH-RW sampler
  #  
  
  lp_filtered = function(theta, path_inds, seg_times, path_fam, init_prev_loc,
                         dims) {
    
    # TODO: re-implement this in c++
    # TODO: allow versions that only return path or timing likelihood components
    
    # initial state of path
    prev_loc = init_prev_loc
    cur_loc = path_fam[[1]][[1]]$path[[path_inds[1]]][1,]
    cur_time = seg_times[[1]][1]
    
    # initialize log-density
    ld = 0
    
    # aggregate log-density over path segments
    for(ind in 1:length(path_inds)) {
      # extract path segment from path family
      path_seg = path_fam[[ind]][[1]]$path[[path_inds[ind]]]
      # process each step in path segment
      for(seg_ind in 1:nrow(path_seg)) {
        # extract current location in path segment
        dst_loc = path_seg[seg_ind,]
        # process cell transitions within path
        if(!identical(cur_loc, dst_loc)) {
          # aggregate grid transition density
          ld = ld + dsmovetools:::TxModelLd(
            cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
            betaAR = theta[1], dst_loc = dst_loc
          ) + 
            # and duration density
            dexp(x = seg_times[[ind]][seg_ind] - cur_time, 
                 rate = exp(theta[-1]), log = TRUE)
          # update current state
          prev_loc = cur_loc
          cur_loc = dst_loc
          cur_time = seg_times[[ind]][seg_ind]
        }
      }
    }
    
    ld
  }
  
  # wrapper to evaluate log-posterior
  lpfn = function(x, x_ind, theta, path) {
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
    lp_filtered(theta = theta, path_inds = path$path_inds, 
                seg_times = path$seg_times, path_fam = path_fam, 
                init_prev_loc = init_prev_loc, dims = dims) + 
    # prior
    dnorm(x = theta[1], mean = priors$beta_ar$mean, sd = priors$beta_ar$sd,
          log = TRUE) +
    sum(dnorm(x = theta[-1], mean = priors$beta$mean, sd = priors$beta$sd,
              log = TRUE))
  }
  
  # optimize initial parameters
  o = optim(par = param_vec, fn = function(theta) {
    lpfn(x = theta, x_ind = 1:n_params, theta = theta, path = path)
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
  
  for(it in 1:niter) {
    
    # update model parameters
    for(i in 1:n_params) {
      update = paramSamplers[[i]]$sample(
        x_ind = i, theta = param_vec, path = path
      )
      param_vec[i] = update$x
    }
    
    # update path
    path = update_path(
      init_prev_loc = init_prev_loc, path_fam = path_fam,
      path_inds = path$path_inds,
      seg_start_times = seg_start_times,
      seg_times = path$seg_times, seg_ind_wts = path$seg_ind_wts,
      seg_time_wts = path$seg_time_wts,
      dims = dims, betaAR = param_vec[1], beta = param_vec[-1]
    )
    
    # save samples
    # samples$path[[it]] = path
    samples$param_vec[it, ] = param_vec
    samples$lp[it] =  lpfn(x = param_vec, x_ind = 1:n_params, 
                           theta = param_vec, path = path)
    
  }
  
  # return final samples
  samples
}