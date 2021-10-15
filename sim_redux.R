library(dsmovetools2d)
library(bisque)

targets::tar_load(sim_obs)

priors = targets::tar_read(sim_priors)

obs = sim_obs$sim_obs_457d4bc5

delta = .05
tseq = seq(from = obs$params$t0, to = obs$params$tf, by = delta)

oseq = matrix(NA, nrow = length(tseq), ncol = 4)
colnames(oseq) = c('x_ind', 'y_ind', 'x_coord', 'y_coord')

xcoords = 1:1e3
ycoords = 1:1e3
surface_heights = numeric(length = length(xcoords) * length(ycoords))

for(ind in 1:nrow(obs$obs$states)) {
  tgt_ind = which.min(abs(obs$obs$times[ind] - tseq))
  oseq[tgt_ind, c('x_coord', 'y_coord')] = obs$obs$states[ind,]
  oseq[tgt_ind, 'x_ind'] = which(oseq[tgt_ind,'x_coord'] == xcoords) - 1
  oseq[tgt_ind, 'y_ind'] = which(oseq[tgt_ind,'y_coord'] == ycoords) - 1
}

# devtools::document('packages/dsmovetools2d/')

dsmovetools2d:::ExactLocFilteredLL(
  init_dsts = matrix(oseq[1,c('x_ind','y_ind')], nrow = 4, ncol = 2, 
                     byrow = TRUE),
  init_srcs = matrix(
    c(oseq[1,c('x_ind','y_ind')] + c(0,1),
      oseq[1,c('x_ind','y_ind')] + c(0,-1),
      oseq[1,c('x_ind','y_ind')] + c(1,0),
      oseq[1,c('x_ind','y_ind')] + c(-1,0)),
    nrow = 4, ncol = 2, byrow = TRUE
  ), 
  init_log_probs = rep(log(1/4), 4),
  obs_x_coords = oseq[,'x_coord'],
  obs_y_coords = oseq[,'y_coord'],
  x_coords = xcoords, 
  y_coords = ycoords, 
  surface_heights = surface_heights, 
  log_self_tx = log(.9), 
  betaAR = 0, 
  lptrunc = Inf
)

o = optim(par = c(0,.5), fn = function(theta) {
  print(theta)
  dsmovetools2d:::ExactLocFilteredLL(
    init_dsts = matrix(oseq[1,c('x_ind','y_ind')], nrow = 4, ncol = 2, 
                       byrow = TRUE) ,
    init_srcs = matrix(
      c(oseq[1,c('x_ind','y_ind')] + c(0,1),
        oseq[1,c('x_ind','y_ind')] + c(0,-1),
        oseq[1,c('x_ind','y_ind')] + c(1,0),
        oseq[1,c('x_ind','y_ind')] + c(-1,0)),
      nrow = 4, ncol = 2, byrow = TRUE
    ) - 1, 
    init_log_probs = rep(log(1/4), 4),
    obs_x_coords = oseq[,'x_coord'],
    obs_y_coords = oseq[,'y_coord'],
    x_coords = xcoords, 
    y_coords = ycoords, 
    surface_heights = surface_heights, 
    log_self_tx = log(1-delta*exp(theta[2])), 
    betaAR = theta[1], 
    lptrunc = Inf
  )
}, control = list(fnscale = -1), hessian = TRUE)


# wrapper to evaluate log-posterior
lpfn = function(x, x_ind, theta) {
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
  
  # observation likelihood
  dsmovetools2d:::ExactLocFilteredLL(
    init_dsts = matrix(oseq[1,c('x_ind','y_ind')], nrow = 4, ncol = 2, 
                       byrow = TRUE) ,
    init_srcs = matrix(
      c(oseq[1,c('x_ind','y_ind')] + c(0,1),
        oseq[1,c('x_ind','y_ind')] + c(0,-1),
        oseq[1,c('x_ind','y_ind')] + c(1,0),
        oseq[1,c('x_ind','y_ind')] + c(-1,0)),
      nrow = 4, ncol = 2, byrow = TRUE
    ) - 1, 
    init_log_probs = rep(log(1/4), 4),
    obs_x_coords = oseq[,'x_coord'],
    obs_y_coords = oseq[,'y_coord'],
    x_coords = xcoords, 
    y_coords = ycoords, 
    surface_heights = surface_heights, 
    log_self_tx = log(1-delta*exp(theta[2])), 
    betaAR = theta[1], 
    lptrunc = Inf
  ) +
  # prior
  dnorm(x = theta[1], mean = priors$betaAR['mean'], sd = priors$betaAR['sd'],
        log = TRUE) +
  dgamma(x = exp(theta[2]), shape = priors$theta['shape'], 
         rate =  priors$theta['rate'], log = TRUE) + 
  # Jacobian to account for transformation in which sampling is done on 
  #   unconstrained space, but prior for theta[1] is defined on [0,\infty)
  jac.log(x = theta[2], log = TRUE)
}

n_params = 2
param_vec = o$par
init_sd = sqrt(diag(solve(-o$hessian)))

# construct MHRW samplers
paramSamplers = lapply(1:n_params, function(ind) {
  dsmovetools::Mhrw1DAdaptive$new(x = param_vec[ind], sd = init_sd[ind], 
                     lp = lpfn, C = .75, alpha = .44, adaptive = TRUE)
})

niter = 1e2

# initialize output
samples = list(
  param_vec = matrix(nrow = niter, ncol = length(param_vec)),
  lp = numeric(niter)
)
colnames(samples$param_vec) = c('betaAR', 'log_theta')

# run sampler
for(it in 1:niter) {
  
  message(it)
  
  # update model parameters
  for(i in 1:n_params) {
    update = paramSamplers[[i]]$sample(
      x_ind = i, theta = param_vec
    )
    param_vec[i] = update$x
  }
  
  # save samples
  samples$param_vec[it, ] = param_vec
  samples$lp[it] =  lpfn(x = param_vec, x_ind = 1:n_params, 
                         theta = param_vec)
  
}

#
# posterior analysis
#

library(coda)

plot(samples$lp, type = 'l')

m = mcmc(samples$param_vec)
m[,2] = exp(m[,2])
colnames(m)[2] = 'theta'

effectiveSize(m)

plot(m[,1])
plot(m[,2])

HPDinterval(m)

unlist(obs$params[c('betaAR','beta')])
