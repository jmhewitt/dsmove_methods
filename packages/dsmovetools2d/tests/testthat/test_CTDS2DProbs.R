context('CTDS2DProbs.h C++ methods')


test_that('Validating directional persistence gives expected results', {
  
  # define grid
  lons = seq(from = 0, to = 40, length.out = 87)
  lats = seq(from = 50, to = 70, length.out = 89)
  
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
  expect_equal(
    exp(p.large[,'log_prob'][nzi]),
    1
  )
  
  # no directional preferences
  expect_equal(
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


test_that('Validating directional persistence respects boundaries', {
  
  # define grid
  lons = seq(from = 0, to = 40, length.out = 103)
  lats = seq(from = 50, to = 70, length.out = 100)
  
  # arbitrary height-field
  surface_heights = runif(length(lons) * length(lats))
  
  # coordinates (lon,lat) format
  init_dsts = rbind(
    c(0,0)
  )
  init_srcs = rbind(
    c(0,1)
  )
  
  betaAR = 1
  
  p = LogTxProbs(
    lons = lons, lats = lats, surface_heights = surface_heights, 
    lon_from_ind = init_srcs[1], lat_from_ind = init_srcs[2], 
    lon_to_ind = init_dsts[1], lat_to_ind = init_dsts[2], 
    betaAR = betaAR
  )
  
  # expected neighborhood structure
  expect_identical(
    structure(
      c(0, 0, 0, 0, 0, 1, 1, 0), 
      .Dim = c(2L, 4L), 
      .Dimnames = list(
        NULL, c("lon_from_ind", "lat_from_ind", "lon_to_ind", "lat_to_ind")
      )
    ), 
    p[,1:4]
  )
  
  # more likely to move in new direction rather than reverse direction
  expect_lt(p[1,'log_prob'], p[2,'log_prob'])
})

test_that('Validating directional persistence respects height boundaries', {
  
  set.seed(2019)
  
  # define grid
  lons = seq(from = 0, to = 40, length.out = 101)
  lats = seq(from = 50, to = 70, length.out = 100)
  
  # arbitrary height-field
  surface_heights = matrix(runif(length(lons) * length(lats)), 
                           nrow = length(lats))
  
  # height above which locations are undefined
  max_height = .8
  
  # find a tall location to transition to/from
  tall_src = which(surface_heights > max_height, arr.ind = TRUE)[1e3,]
  
  
  # surface_heights column-major list of heights at each grid point,
  # where latitudes form the rows, and longitudes form the columns
  
  # coordinates (lon,lat) format, shifted for 0-indexing
  init_dsts = c(
    'lon_ind' = as.numeric(tall_src['col']) - 1, 
    'lat_ind' = as.numeric(tall_src['row']) - 1
  )
  
  init_srcs = init_dsts + c(0,1)
  
  # verify the init_dst location is not a valid location
  expect_gt(
    surface_heights[init_dsts['lat_ind']+1, init_dsts['lon_ind']+1], 
    max_height
  )
  
  # force the init_src location to be valid
  surface_heights[
    init_srcs['lat_ind']+1, init_srcs['lon_ind']+1
  ] = max_height - .1
  
  betaAR = 1
  
  p = LogTxProbsElevation(
    lons = lons, lats = lats, surface_heights = surface_heights, 
    lon_from_ind = init_srcs['lon_ind'], lat_from_ind = init_srcs['lat_ind'], 
    lon_to_ind = init_dsts['lon_ind'], lat_to_ind = init_dsts['lat_ind'], 
    betaAR = betaAR, min_elevation = 0, max_elevation = .8
  )
  
  # verify no transitions are possible from an invalid location
  expect_equal(
    nrow(p), 0
  )
  
  # transition probs. when state direction is reversed, starting from a good loc.
  p.rev = LogTxProbsElevation(
    lons = lons, lats = lats, surface_heights = surface_heights, 
    lon_from_ind = init_dsts['lon_ind'], lat_from_ind = init_dsts['lat_ind'], 
    lon_to_ind = init_srcs['lon_ind'], lat_to_ind = init_srcs['lat_ind'], 
    betaAR = betaAR, min_elevation = 0, max_elevation = .8
  )
  
  # verify no transitions are made to the bad location
  expect_false(
    any(
      (p.rev[,'lon_to_ind'] == init_dsts['lon_ind']) & 
      (p.rev[,'lat_to_ind'] == init_dsts['lat_ind'])
    )
  )
  
  # verify total prob mass is 1
  expect_equal(
    sum(exp(p.rev[,'log_prob'])), 1
  )
  
  # verify all transitions are made to good locations
  expect_true(
    all(
      apply(p.rev, 1, function(r) {
        surface_heights[r['lat_to_ind']+1, r['lon_to_ind']+1]
      }) < max_height
    )
  )
  
})
