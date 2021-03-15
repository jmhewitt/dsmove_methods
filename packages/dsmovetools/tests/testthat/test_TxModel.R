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
