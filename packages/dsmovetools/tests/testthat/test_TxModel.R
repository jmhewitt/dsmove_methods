context('TxModel.cpp C++ methods')

test_that('Validating directional persistence gives expected results', {

  cur_loc = c(50, 50)
  prev_loc = c(49, 50)
  dims = c(100, 100)
  
  betaAR = 1
  betaAR.equal = 0
  betaAR.large = 1e2
  
  p = TxModelParams(cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
                    betaAR = betaAR)
  
  p.equal = TxModelParams(cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
                          betaAR = betaAR.equal)
  
  p.large = TxModelParams(cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
                          betaAR = betaAR.large)
  
  # extremely strong preference to continue motion in current direction
  nzi = (p.large[,1] == 51) & (p.large[,2] == 50)
  expect_identical(
    exp(p.large[,3][nzi]),
    1
  )
  
  # no directional preferences
  expect_identical(
    rep(1/nrow(p.equal), nrow(p.equal)),
    exp(p.equal[,3])
  )
  
  # un-shifted probabilities, to recover betaAR
  #  NOTE: This test can break if the neighborhood enumeration changes
  expect_identical(
    as.numeric(scale(p[,3], center = TRUE, scale = FALSE)),
    c(-1, 0, 1, 0)
  )
  
  # probabilities are properly normalized
  expect_identical(
    sum(exp(p[,3])), 
    1
  )
  
})

test_that('Validating cell transition sampler', {
  
  cur_loc = c(50, 50)
  prev_loc = c(49, 50)
  dims = c(100, 100)
  
  betaAR = 1
  betaAR.equal = 0
  
  # exact transition parameters
  p = TxModelParams(cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
                    betaAR = betaAR)
  p.equal = TxModelParams(cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
                          betaAR = betaAR.equal)
  
  
  # crude/quick monte carlo test sampler implementation
  nrep = 1e4
  set.seed(2021)
  
  p.samples = replicate(n = nrep,
    TxModelSample(cur_loc = cur_loc, prev_loc = prev_loc,dims = dims, 
                  betaAR = betaAR)
  )
  p.samples.equal = replicate(n = nrep,
    TxModelSample(cur_loc = cur_loc, prev_loc = prev_loc,dims = dims, 
                  betaAR = betaAR.equal)
  )
  
  # monte carlo probabilities
  p.empirical = table(p.samples[1,], p.samples[2,]) / nrep
  p.empirical.equal = table(p.samples.equal[1,], p.samples.equal[2,]) / nrep
  

  # quick-estimate probabilities are accurate to two decimal places
  expect_equal(
    exp(p[,3]),
    apply(p, 1, function(r) {
      p.empirical[as.character(r[1]), as.character(r[2])]
    }), 
    tolerance = 1e-2
  )
  
  # ...on different distributions
  expect_equal(
    exp(p.equal[,3]),
    apply(p.equal, 1, function(r) {
      p.empirical.equal[as.character(r[1]), as.character(r[2])]
    }), 
    tolerance = 1e-2
  )
  
})

test_that('Validating cell transition density function', {
  
  cur_loc = c(50, 50)
  prev_loc = c(49, 50)
  dims = c(100, 100)
  
  betaAR = 1
  betaAR.equal = 0
  
  # exact transition parameters
  p = TxModelParams(cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
                    betaAR = betaAR)
  p.equal = TxModelParams(cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
                          betaAR = betaAR.equal)
  
  # recover log densities by the output location
  
  expect_identical(
    apply(p[(nrow(p):1),], 1, function(r) { 
      TxModelLd(cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
                betaAR = betaAR, dst_loc = r[1:length(dims)])
    }),
    rev(p[,3])
  )
  
  expect_identical(
    apply(p.equal, 1, function(r) { 
      TxModelLd(cur_loc = cur_loc, prev_loc = prev_loc, dims = dims, 
                betaAR = betaAR.equal, dst_loc = r[1:length(dims)])
    }),
    p.equal[,3]
  )
  
})
