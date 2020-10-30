fit_integration = function(segments, obs, inits, priors, niter, ctds_domain) {
  
  # extract sampling weights and indexes for all segments
  pathwts = lapply(segments, function(sfam) {
    do.call(rbind, lapply(1:length(sfam), function(lengthind) {
      s = sfam[[lengthind]]
      npaths = ifelse(is.null(s$paths), 1, nrow(s$paths))
      data.frame(lengthind = lengthind, pathind = 1:npaths, w = s$weights)
    }))
  })
  
  #
  # initialize data augmentation path
  #
  
  path_components = vector('list', length(segments) + 1)
  
  # randomly sample an initial edge, and set initial time
  in_edges = ctds_domain$in_edges_inds[[obs$states[1]]]
  path_components[[1]]$path_ind = ifelse(
    length(in_edges) > 1, sample(x = 1:length(in_edges), size = 1), 1
  )
  path_components[[1]]$epath = in_edges[path_components[[1]]$path_ind]
  path_components[[1]]$tpath = obs$times[1]
  
  # sample path segments
  for(i in 2:length(path_components)) {
    p = pathwts[[i-1]]
    # sample path segment
    # seg_ind = ifelse(nrow(p) > 1, sample(x = nrow(p), size = 1, prob = p$w), 1)
    seg_ind = 1
    path_components[[i]]$path_ind = seg_ind
    # extract path
    r = p[seg_ind, ]
    path_components[[i]]$epath = segments[[i-1]][[r$lengthind]]$paths[r$pathind,]
    # sample times
    path_components[[i]]$tpath = sort(runif(
      n = length(path_components[[i]]$epath), 
      min = obs$times[i-1], 
      max = obs$times[i]
    ))
  }
  
  
  #
  # sampler
  #
  
  # initialize state
  param_vec = c(inits$beta_ar, inits$beta_loc)
  
  # initialize output
  samples = list(
    path = vector('list', niter),
    param_vec = matrix(nrow = niter, ncol = length(param_vec)),
    ll = numeric(niter)
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
    ll_exact(epath = epath, durations = durations,
             beta_loc = matrix(theta[-1], ncol = 1), beta_dir = 0,
             beta_ar = theta[1], ctds_struct = ctds_domain)
  }
  
  # posterior mode for initial conditions
  o = optim(par = param_vec, function(theta) {
    lpfn(theta = theta, epath = epath, durations = diff(tpath))
  }, control = list(fnscale = -1), method = 'BFGS', hessian = TRUE)
  
  # update initial parameters
  param_vec = o$par
  
  # construct MHRW sampler
  paramSampler = MhrwAdaptive$new(
    x = param_vec, mu = o$par, Sigma = solve(-o$hessian), lambda = c(1, 1), 
    lp = lpfn, C = .1, alpha = 1, alpha_star = .44, adaptive = TRUE, 
    adaptation_frequency = 1
  )
  
  # compute initial likelihood
  ll = lpfn(theta = param_vec, epath = epath, durations = diff(tpath))
  
  for(it in 1:niter) {
    
    # update model parameters
    update = paramSampler$sample(epath =  epath, durations = diff(tpath))
    param_vec = update$x
    ll = update$lp
    
    #
    # update path
    #
    
    message(paste('Step: ', it, '\n'))
    
    # sample path segments
    for(i in 2:length(path_components)) {

      # # flatten current path, and get ll for segment
      # epath_seg = do.call(c, sapply(path_components[(i-1):(i+1)], function(p) {
      #   unlist(p$epath)
      # }))
      # tpath_seg = do.call(c, sapply(path_components[(i-1):(i+1)],
      #                               function(p) { p$tpath }))
      # ll_seg = ll_exact(epath = epath_seg, durations = diff(tpath_seg),
      #                   beta_loc = matrix(param_vec[-1], ncol = 1), beta_dir = 0,
      #                   beta_ar = param_vec[1], ctds_struct = ctds_domain)

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
        
        if(any(length(x0$tpath) == 1,
               length(x_prop$tpath) == 1)) {
          browser()
        }
        
        # compute likelihood on original and changed path segments
        ll0 = ll_exact(epath = x0$epath, durations = diff(x0$tpath),
                       beta_loc = matrix(param_vec[-1], ncol = 1),
                       beta_dir = 0, beta_ar = param_vec[1],
                       ctds_struct = ctds_domain)
        ll_prop = ll_exact(epath = x_prop$epath, durations = diff(x_prop$tpath),
                           beta_loc = matrix(param_vec[-1], ncol = 1),
                           beta_dir = 0, beta_ar = param_vec[1],
                           ctds_struct = ctds_domain)
        
        # # compute likelihood only on changed portion of path
        # ll_prop = ll_exact(epath = epath_prop, durations = diff(tpath_prop),
        #                    beta_loc = matrix(param_vec[-1], ncol = 1),
        #                    beta_dir = 0, beta_ar = param_vec[1],
        #                    ctds_struct = ctds_domain)
        # 
        # 
        # 
        # # accept/reject
        # logR = ll_prop - ll +
        #   log(pathwts[[i-1]]$w[path_components[[i]]$path_ind]) -
        #   log(pathwts[[i-1]]$w[prop_components$path_ind])
        
        # accept/reject 
        logR = ll_prop - ll0 + 
          log(pathwts[[i-1]]$w[path_components[[i]]$path_ind]) -
          log(pathwts[[i-1]]$w[prop_components$path_ind])
        
        
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
    
    # # maximum possible length of trajectories
    # sum(unlist(sapply(segments, function(s) {
    #   ncol(s[[length(s)]]$paths)
    # })))
    # 
    # # length of sampled trajectories
    # plot(sapply(samples$path, function(p) {
    #   nrow(p)
    # })[-(1:60)], type = 'l', ylab = 'path len')
    # 
    
    message(ll)
    print(param_vec)
  }
  
  browser()
  
  samples
}