context('ffar.cpp C++ methods')

test_that('Validating random walk self-transition filtering on a grid', {

  # number of dimensions
  ndim = 3
  # number of coordinates in each dimension
  dims = c(10,10,1)
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
  
  # iterative probability diffusion, implemented in R
  pstep = function(x0, x0_prev, log_a0val, dims, steps, log_self_tx, betaAR) {
    
    # log of transition probability
    log_tx = log(1 - exp(log_self_tx))
      
    # initialize diffused probability vector with initial conditions
    res = matrix(c(x0, x0_prev, log_a0val), nrow = 1)
    
    # run diffusions
    for(i in 1:steps) {
      # unaggregated diffusion probabilities within step
      tmp = matrix(nrow = 0, ncol = ncol(res))
      # forward-diffuse mass from each atom in the current diffused prob. vector
      for(j in 1:nrow(res)) {
        # reachable destinations and forward probabilities
        p = TxModelParams(cur_loc = res[j, 1:length(x0)], 
                          prev_loc = res[j, length(x0) + 1:length(x0)], 
                          dims = dims, betaAR = betaAR)
        # reweight forward probabilities by source and self-tx. probs
        lp = p[, ncol(p)] + log_tx + res[j, ncol(res)]
        # format output
        tmp = rbind(
          # append with mass diffused from other atoms
          tmp,
          # diffused mass from leaving current atom
          cbind(p[,-ncol(p)], 
                matrix(data = res[j,1:length(x0)], 
                       nrow = nrow(p), 
                       ncol = length(x0), byrow = TRUE), 
                lp),
          # diffused mass from staying in current atom (self-transition)
          matrix(c(res[j,-ncol(res)], log_self_tx + res[j, ncol(res)]), 
                 nrow = 1)
        )
      }
      # aggregate transition probabilities
      res = data.frame(tmp) %>% 
        group_by(V1, V2, V3, V4, V5, V6) %>% 
        summarise(lp = log(sum(exp(lp)))) %>% 
        ungroup() %>% as.matrix()
    }
    
    # package results
    dimnames(res) = NULL
    res
  }
  
  
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
  
  # forward-filter via R
  af.alt = pstep(
    x0 = x0, x0_prev = x0_prev, log_a0val = log(p0), dims = dims, 
    steps = nsteps, log_self_tx = log(self_tx_prob), betaAR = betaAR
  )
  
  # check equivalence  
  expect_equal(af, af.alt)
})
