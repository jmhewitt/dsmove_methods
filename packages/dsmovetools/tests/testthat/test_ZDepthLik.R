context('DepthLik.h C++ methods')

test_that('Validating basic input and output for depth likelihoods', {
  
  dims = c(13, 13, 16)
  
  zvals = 1:16
  
  zfield = matrix(0, nrow = dims[1], ncol = dims[2])

  # raise the 3D space's vertical surface by placing "mountains" where the
  # complex 2D spatial domain is not valid, here for loc's near the edge and
  # at the average zval height
  coords = expand.grid(1:dims[1], 1:dims[2])
  for(i in 1:nrow(coords)) {
    if(sum(coords[i,]) < 7) {
      zfield[i] = mean(zvals) + 1e-10
    }
  }
  
  obs_depths = c(10,NA,10,0,0)
  
  # coordinate has reasonable depth wrt. observation
  expect_equal(
    DepthLikEval(dims = dims, coords = c(0,0,0), obs_depths = obs_depths,
                 zfield = as.numeric(zfield), zvals = zvals, ind = 0),
    0
  )
  
  # coordinate does not have reasonable depth wrt. observation
  expect_equal(
    DepthLikEval(dims = dims, coords = c(0,0,0), obs_depths = obs_depths,
                 zfield = as.numeric(zfield), zvals = zvals, ind = 4),
    -Inf
  )
  
  # coordinate has reasonable depth wrt. observation
  expect_equal(
    DepthLikEval(dims = dims, coords = c(10,10,0), obs_depths = obs_depths,
                 zfield = as.numeric(zfield), zvals = zvals, ind = 4),
    0
  )  
  
  # na observations handled correctly
  expect_equal(
    DepthLikEval(dims = dims, coords = c(0,0,0), obs_depths = obs_depths,
                 zfield = as.numeric(zfield), zvals = zvals, ind = 1),
    -Inf
  )
  
})