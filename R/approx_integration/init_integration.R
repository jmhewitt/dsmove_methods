init_integration = function(segments, obs, inits, priors, niter, ctds_domain,
                            u = runif(n = 1)) {
  # Parameters:
  #  u - common variate used to draw path lengths across all segments for this 
  #      path
  
  # extract sampling weights and indexes for all segments
  pathwts = lapply(1:length(segments), function(segind) {
    sfam = segments[[segind]]
    do.call(rbind, lapply(1:length(sfam), function(lengthind) {
      s = sfam[[lengthind]]
      if(is.null(s$paths)) {
        npaths = 1
        pathlen = 0
        s$ll = 0
      } else {
        npaths = nrow(s$paths)
        pathlen = ncol(s$paths)
      }
      data.frame(w = s$weights, pathlen = pathlen, pathind = 1:npaths, 
                 lengthind = lengthind, ll = s$ll)
    }))
  })
  
  # aggregate path weights by path length
  pathwts.aggregated = lapply(pathwts, function(p) {
    p %>% 
      dplyr::group_by(pathlen) %>% 
      dplyr::summarise(w = sum(w)) %>% 
      dplyr::ungroup()
  })
  
  
  #
  # initialize data augmentation path
  #
  
  # initialize path segment container
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
    
    segind = i-1
    p = pathwts[[segind]]
    
    # extract length of sampled segment
    plen = pathwts.aggregated[[segind]] %>% 
      dplyr::mutate(cdf = cumsum(w), exceeds = cdf >= u) %>% 
      dplyr::filter(exceeds == TRUE) %>% 
      dplyr::slice(1) %>% 
      dplyr::select(pathlen) %>% 
      unlist()
  
    # sample segment conditional on segment length
    seg_ind = p %>% 
      dplyr::mutate(rowind = 1:dplyr::n()) %>%
      dplyr::filter(pathlen == plen) %>% 
      dplyr::sample_n(size = 1, weight = w) %>% 
      dplyr::select(rowind) %>%
      unlist()
  
    path_components[[i]]$path_ind = seg_ind
    
    # extract path
    r = p[seg_ind, ]
    path_components[[i]]$epath = segments[[segind]][[r$lengthind]]$paths[r$pathind,]
    
    # sample times
    path_components[[i]]$tpath = sort(runif(
      n = length(path_components[[i]]$epath), 
      min = obs$times[i-1], 
      max = obs$times[i]
    ))
    
  }
   
   
  #
  # optimize path and associated model params. via sim.-annealing-like scheme
  #
    
  # initialize state
  param_vec = c(inits$beta_ar, inits$beta_loc)
  
  # flatten path
  epath = do.call(c, sapply(path_components, function(p) {
    unlist(p$epath)
  }))
  tpath = do.call(c, sapply(path_components, function(p) { p$tpath }))
  
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
  
  # posterior mode for initial conditions
  o = optim(par = param_vec, function(theta) {
    lpfn(theta = theta, epath = epath, durations = diff(tpath))
  }, control = list(fnscale = -1), method = 'BFGS', hessian = TRUE)
  
  # update initial parameters and extract initial likelihood
  param_vec = o$par
  ll = o$value
  hessian = o$hessian
  
  # TODO: add in convergence criteria
  
  for(it in 1:niter) {
    
    #
    # update path
    #
    
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
          # # exact path proposal
          # pathwts[[i-1]]$ll[path_components[[i]]$path_ind] -
          # pathwts[[i-1]]$ll[prop_components$path_ind] + 
          # times proposal
          (lfactorial(ntimespath) - ntimespath * ltrange ) -
          (lfactorial(ntimesprop) - ntimesprop * ltrange )
          
        accept = logR > 0
        
        if(accept) {
          path_components = prop_components_all
          epath = epath_prop
          tpath = tpath_prop
        }
        
      }

    }
    
    # update model parameters
    o = optim(par = param_vec, function(theta) {
      lpfn(theta = theta, epath = epath, durations = diff(tpath))
    }, control = list(fnscale = -1), method = 'BFGS', hessian = TRUE)
    
    # update parameters and likelihood
    param_vec = o$par
    ll = o$value
    hessian = o$hessian
    
  }
  
  list(
    path_components = path_components,
    params = list(
      beta_ar = param_vec[1], 
      beta_loc = param_vec[-1]
    ),
    ll = ll,
    hessian = hessian
  )
}
