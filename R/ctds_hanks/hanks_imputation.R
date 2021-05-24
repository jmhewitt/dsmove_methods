hanks_imputation = function(states, times, reps, samples_per_rep, priors) {
  # Parameters:
  #  states - matrix of observation states/coordinates
  #  times - times at which observations are made
  #  reps - number of imputed samples to draw
  #  samples_per_rep - number of MC samples to draw for \theta per imputation

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
  knots = seq(min(t),max(t),by=1/2)
  ## create B-spline basis vectors used to approximate the path
  b=create.bspline.basis(c(min(t),max(t)),breaks=knots,norder=3)
  ## define the sequence of times on which to sample the imputed path
  tpred=seq(min(t),max(t),by=.01)

  ## Fit latent Gaussian model using MCMC
  out=mcmc.fmove(xy,t,b,tpred,QQ="CAR",n.mcmc=1000,a=1,r=1,num.paths.save=reps)
  
  #
  # analyze each of the imputed trajectories
  #

  samples_out = c()

  for(rep in 1:reps) {

    # discretize path
    path=out$pathlist[[rep]]
    path.lattice = round(path$xy)
    
    # get length of discretized path
    n_imputed = sum(
      diff(path.lattice[,1]) | diff(path.lattice[,2])
    )
    
    # conditional posterior for theta, given path
    post_params = c(shape = priors$shape + n_imputed, 
                    rate = priors$rate + diff(range(path$t)))
    
    # draw posterior samples
    samples_out = c(
      samples_out, rgamma(n = samples_per_rep, shape = post_params['shape'], 
                          rate = post_params['rate'])
    )

  }

  list(
    samples = samples_out
  )
}
