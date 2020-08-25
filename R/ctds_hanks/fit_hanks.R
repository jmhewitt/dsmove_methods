fit_hanks = function(ctds_struct, obs) {
  
  # load observations
  ctds_obs = readRDS(obs)
  
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
    ctds_struct$Xloc, 
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
  
  P = 5
  
  # impute and format P paths
  for(i in 1:P) {
    
    cat(i," ")
    
    # discretize and format path
    path=out$pathlist[[i]]
    ctmc=path2ctmc(path$xy,path$t, X, method="LinearInterp")
    glm.list[[i]]=ctmc2glm(ctmc, X, X)
    
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
  
  fit.SWL=glm(z ~ crw, weights = rep(1/P, nrow(glm.data)), 
              family = "poisson", offset = log(tau), data = glm.data)
  
  fit.SWL
}
