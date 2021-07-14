context('CTDS2DProbs.h C++ methods')


test_that('Validating directional persistence gives expected results', {
  
  # define grid
  lons = seq(from = 0, to = 40, length.out = 100)
  lats = seq(from = 50, to = 70, length.out = 100)
  
  # arbitrary height-field
  surface_heights = runif(length(lons) * length(lats))
  
  # coordinates (lon,lat) format
  init_dsts = rbind(
    c(50,50)
  )
  init_srcs = rbind(
    c(49,50)
  )
  
  betaAR = 1
  betaAR.equal = 0
  betaAR.large = 1e2
  
  p = LogTxProbs(
    lons = lons, lats = lats, surface_heights = surface_heights, 
    lon_from_ind = init_srcs[1], lat_from_ind = init_srcs[2], 
    lon_to_ind = init_dsts[1], lat_to_ind = init_dsts[2], 
    betaAR = betaAR
  )
  
  p.equal = LogTxProbs(
    lons = lons, lats = lats, surface_heights = surface_heights, 
    lon_from_ind = init_srcs[1], lat_from_ind = init_srcs[2], 
    lon_to_ind = init_dsts[1], lat_to_ind = init_dsts[2], 
    betaAR = betaAR.equal
  )
  
  p.large = LogTxProbs(
    lons = lons, lats = lats, surface_heights = surface_heights, 
    lon_from_ind = init_srcs[1], lat_from_ind = init_srcs[2], 
    lon_to_ind = init_dsts[1], lat_to_ind = init_dsts[2], 
    betaAR = betaAR.large
  )
  
  # extremely strong preference to continue motion in current direction
  nzi = (p.large[,'lon_to_ind'] == 51) & (p.large[,'lat_to_ind'] == 50)
  expect_identical(
    exp(p.large[,'log_prob'][nzi]),
    1
  )
  
  # no directional preferences
  expect_identical(
    rep(1/nrow(p.equal), nrow(p.equal)),
    exp(p.equal[,'log_prob'])
  )
  
  # un-shifted probabilities, to recover betaAR
  #  NOTE: This test can break if the neighborhood enumeration changes
  expect_equal(
    as.numeric(scale(p[,'log_prob'], center = TRUE, scale = FALSE)),
    c(-1,0,0,1)
  )
  
  # probabilities are properly normalized
  expect_equal(
    sum(exp(p[,'log_prob'])), 
    1
  )
  
})

