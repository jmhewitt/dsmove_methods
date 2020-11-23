fit_hanks = function(ctds_struct, ctds_obs, weibull_est) {
  
  # following package use from help('ctmcmove')
  
  # extract coordinate/time triples
  xyt = cbind(ctds_struct$coords[ctds_obs$states,], ctds_obs$times)
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
  knots = seq(min(t),max(t),by=1/4)
  ## create B-spline basis vectors used to approximate the path
  b=create.bspline.basis(c(min(t),max(t)),breaks=knots,norder=3)
  ## define the sequence of times on which to sample the imputed path
  tpred=seq(min(t),max(t),by=1/24/60)
  
  ## Fit latent Gaussian model using MCMC
  out=mcmc.fmove(xy,t,b,tpred,QQ="CAR",n.mcmc=400,a=1,r=1,num.paths.save=30)
  
  
  #
  # creating rasters
  #
  
  X = raster(
    matrix(ctds_struct$Xloc, nrow = length(unique(ctds_struct$coords[, 's1']))), 
    xmn = min(ctds_struct$coords[, 's1']),
    xmx = max(ctds_struct$coords[, 's1']), 
    ymn = min(ctds_struct$coords[, 's2']), 
    ymx = max(ctds_struct$coords[, 's2'])
  )
  
  crs(X) = "+proj=utm +zone=32"
  
  
  #
  # Turn CTMC discrete path into latent Poisson GLM data
  #
  
  glm.list=list()
  ctmc.list = list()
  path.list = list()
  
  P = 5
  
  # impute and format P paths
  for(i in 1:P) {
    
    cat(i," ")
    
    # discretize and format path
    path=out$pathlist[[i]]
    ctmc=path2ctmc(path$xy,path$t, X, method="ShortestPath")
    
    glm.list[[i]]=ctmc2glm(ctmc, X, X)
    
    # save path 
    ctmc.list[[i]] = ctmc
    path.list[[i]] = path
    
    ## remove transitions that are nearly instantaneous
    ##  (These are essentially outliers in the following regression analyses)
    idx.0=which(glm.list[[i]]$tau<10^-5)
    if(length(idx.0)>0){
      glm.list[[i]]=glm.list[[i]][-idx.0,]
    }
    glm.list[[i]]$t=glm.list[[i]]$t-min(glm.list[[i]]$t)
    
  }
  
  

  
  ##
  ## Stack the P imputations together
  ##
  
  glm.data=glm.list[[1]]
  for(i in 2:P){
    glm.data=rbind(glm.data, glm.list[[i]])
  }
  
  ##########################################################################
  ##
  ## 6. Fit Poisson GLM
  ##    (here we are fitting all "M" paths simultaneously,
  ##     giving each one a weight of "1/M")
  ##
  ##########################################################################
  
  X = cbind(1, glm.data$crw)
  offsets = log(glm.data$tau)
  fit = optim(par = c(0,0), fn = function(theta) {
    ll = 0
    log_lambdas = X %*% theta[1:2]
    for(i in 1:(nrow(glm.data)/4)) {
      inds = (i-1)*4 + 1:4
      lambdas = exp(log_lambdas[inds])
      lambdas_total = sum(lambdas)
      tau = glm.data$tau[inds[1]]
      ll = ll +
        # cell transition
        sum(log_lambdas[inds] * glm.data$z[inds]) - log(lambdas_total) +
        # exponential cell duration
        dexp(x = tau, rate = exp(theta[1]), log = TRUE)
    }
    ll/P
  }, control = list(fnscale = -1), hessian = TRUE)
  
  list(
    fit = fit,
    tstep = mean(diff(ctds_obs$times)),
    ctmc.list = ctmc.list,
    # glm.list = glm.list,
    # path.list = path.list,
    raster.coords = coordinates(X)
  )
}
