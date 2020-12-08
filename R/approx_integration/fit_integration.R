fit_integration = function(segments, obs, inits, priors, niter, ctds_domain,
                           ncheckpoints, checkpoint_fn) {
  
  # extract sampling weights and indexes for all segments
  pathwts = lapply(segments, function(sfam) {
    do.call(rbind, lapply(1:length(sfam), function(lengthind) {
      s = sfam[[lengthind]]
      npaths = ifelse(is.null(s$paths), 1, nrow(s$paths))
      data.frame(lengthind = lengthind, pathind = 1:npaths, w = s$weights)
    }))
  })
  
  #
  # read in initial data augmentation path
  #
  
  path_components = inits$path_components
  
  
  #
  # sampler
  #
  
  # read in initial parameters
  param_vec = c(inits$beta_ar, inits$beta_loc)
  
  # initialize output
  samples = list(
    path = vector('list', niter),
    param_vec = matrix(nrow = niter, ncol = length(param_vec)),
    ll = numeric(niter),
    mem = vector('list', niter)
  )
  
  # flatten path
  epath = do.call(c, sapply(path_components, function(p) {
    unlist(p$epath)
  }))
  tpath = do.call(c, sapply(path_components, function(p) { p$tpath }))

  #
  # initialize MH-RW sampler
  #
  
  # wrapper to evaluate log-posterior 
  lpfn = function(theta, epath, durations) {
    # likelihood
    ll_exact(epath = epath, durations = durations,
             beta_loc = matrix(theta[-1], ncol = 1), beta_dir = 0,
             beta_ar = theta[1], ctds_struct = ctds_domain) + 
    # prior
    dnorm(x = theta[1], mean = priors$beta_ar$mean, 
          sd = priors$beta_ar$sd, log = TRUE) + 
      dnorm(x = theta[-1], mean = priors$beta_loc$mean, 
            sd = priors$beta_loc$sd, log = TRUE)
  }
  
  # construct MHRW sampler
  paramSampler = MhrwAdaptive$new(
    x = param_vec, mu = param_vec, Sigma = inits$prop_cov, 
    lambda = c(1, 1), lp = lpfn, C = .1, alpha = 1, alpha_star = .44, 
    adaptive = TRUE, adaptation_frequency = 1
  )
  
  # initialize data likelihood
  ll = ll_exact(epath = epath, durations = diff(tpath),
                beta_loc = matrix(param_vec[-1], ncol = 1),
                beta_dir = 0, beta_ar = param_vec[1],
                ctds_struct = ctds_domain)
  
  # set checkpoint indices
  checkpoint.inds = seq(from = 1, to = niter, length.out = ncheckpoints)[-1]
    
  for(it in 1:niter) {
    
    # update model parameters
    update = paramSampler$sample(epath =  epath, durations = diff(tpath))
    param_vec = update$x
    
    #
    # update path
    #
    
    message(paste('Step: ', it, '\n'))
    
    # sample path segments
    for(i in 2:length(path_components)) {

      # propose new path segment
      p = pathwts[[i-1]]
      seg_ind = ifelse(nrow(p) > 1, sample(x = nrow(p), size = 1, prob = p$w), 1)
      prop_components = list(path_ind = seg_ind)

      # extract path
      r = p[seg_ind, ]
      prop_components$epath = segments[[i-1]][[r$lengthind]]$paths[r$pathind,]
      # sample times
      prop_components$tpath = sort(runif(
        n = length(prop_components$epath),
        min = obs$times[i-1],
        max = obs$times[i]
      ))

      # flatten proposed path
      prop_components_all = path_components
      prop_components_all[[i]] = prop_components
      epath_prop = do.call(c, sapply(prop_components_all, function(p) {
        unlist(p$epath)
      }))
      tpath_prop = do.call(c, sapply(prop_components_all, function(p) {
        p$tpath
      }))
      
      if(!identical(epath_prop, epath)) {

        # extract original and proposed path segments
        x0 = extract_segment(epath = epath, tpath = tpath, 
                             tmin = obs$times[i-1], tmax = obs$times[i])
        x_prop = extract_segment(epath = epath_prop, tpath = tpath_prop, 
                                 tmin = obs$times[i-1], tmax = obs$times[i])
        
        # compute likelihood on original and changed path segments
        ll0 = lpfn(theta = param_vec, epath = x0$epath, 
                   durations = diff(x0$tpath))
        ll_prop = lpfn(theta = param_vec, epath = x_prop$epath, 
                       durations = diff(x_prop$tpath))
        
        ltrange = log(obs$times[i] - obs$times[i-1])
        ntimesprop = length(prop_components$tpath)
        ntimespath = length(path_components[[i]]$tpath)
        
        # accept/reject 
        logR = ll_prop - ll0 + 
          # path proposal
          log(pathwts[[i-1]]$w[path_components[[i]]$path_ind]) -
          log(pathwts[[i-1]]$w[prop_components$path_ind]) + 
          # times proposal
          (lfactorial(ntimespath) - ntimespath * ltrange ) -
          (lfactorial(ntimesprop) - ntimesprop * ltrange )
        accept = log(runif(n = 1)) < logR
        
        if(accept) {
          path_components = prop_components_all
          ll = ll_exact(epath = epath_prop, durations = diff(tpath_prop),
                        beta_loc = matrix(param_vec[-1], ncol = 1),
                        beta_dir = 0, beta_ar = param_vec[1],
                        ctds_struct = ctds_domain)
          epath = epath_prop
          tpath = tpath_prop
        }
        
      }

    }
    
    # save samples
    samples$path[[it]] = cbind(epath = epath, tpath = tpath)
    samples$param_vec[it, ] = param_vec
    samples$ll[it] =  ll
    samples$mem[[it]] = sapply(ls(), object_size)
    
    # checkpoint operations
    if(it %in% checkpoint.inds) {
      checkpoint_fn(samples)
    }
    
  }
  
  # return final samples
  samples
}