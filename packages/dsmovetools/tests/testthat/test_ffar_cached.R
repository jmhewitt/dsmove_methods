context('ffar_cached.cpp C++ methods')

test_that('Validating random walk self-transition filtering on a grid', {

  # number of dimensions
  ndim = 3
  # number of coordinates in each dimension
  dims = c(100,100,1)
  # coordinates
  coords = expand.grid(x = 0:(dims[1]-1), y = 0:(dims[2]-1), z = 0)
  # self-transition probability
  self_tx_prob = .7
  # directional persistence parameter
  betaAR = .5
  # initial probability mass
  p0 = 1

  # height of vertical layer
  zval = 1

  # height of domain surface
  zsurf = matrix(0, nrow = dims[1], ncol = dims[2])
  
  
  #
  # 7-step diffusion
  #

  x0 = matrix(c(4,1,0), nrow = 1)
  x0_prev = matrix(c(5,1,0), nrow = 1)
  nsteps = 7

  # forward-filter via C++
  af = FFRWLightLogConstrainedSelfTxAR(
    a0 = x0, a0_prev_coords = x0_prev, log_a0val = log(p0), dims = dims, 
    steps = nsteps, surface_heights = zsurf, domain_heights = zval, 
    log_self_tx = log(self_tx_prob), betaAR = betaAR
  )
  
  # forward-filter via cached C++
  af.alt = FFRWLightLogConstrainedSelfTxARCached(
    a0 = x0, a0_prev_coords = x0_prev, log_a0val = log(p0), dims = dims, 
    steps = nsteps, surface_heights = zsurf, domain_heights = zval, 
    log_self_tx = log(self_tx_prob), betaAR = betaAR
  )
  
  library(microbenchmark)
  
  microbenchmark(
    FFRWLightLogConstrainedSelfTxAR(
      a0 = x0, a0_prev_coords = x0_prev, log_a0val = log(p0), dims = dims, 
      steps = nsteps, surface_heights = zsurf, domain_heights = zval, 
      log_self_tx = log(self_tx_prob), betaAR = betaAR
    )
  )
  
  microbenchmark(
    FFRWLightLogConstrainedSelfTxARCached(
      a0 = x0, a0_prev_coords = x0_prev, log_a0val = log(p0), dims = dims, 
      steps = nsteps, surface_heights = zsurf, domain_heights = zval, 
      log_self_tx = log(self_tx_prob), betaAR = betaAR
    )
  )
  
  # check equivalence  
  expect_equal(af, af.alt)
})
