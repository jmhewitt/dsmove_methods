context('ZConstrainedRookNeighborhood.h C++ methods')

test_that('Validating constrained 2D rook neighborhoods on a grid', {
  
  # we represent a complex 2D spatial domain via a 3D space that only has one 
  # vertical height layer.  ultimately, we place "impossibly high" mountains at 
  # locations in the 3D space where the 2D spatial domain is not defined.
  dims = c(16, 16, 1)
  
  # define the numerical height of the 3D space's vertical layer
  zvals = 0
  
  # valid coordinates in the 3D space have a vertical height that meets or 
  # exceeds the 3D space's vertical surface.  start by allowing all locations.
  zfield = matrix(zvals, nrow = dims[1], ncol = dims[2])
  
  # raise the 3D space's vertical surface by placing "mountains" where the 
  # complex 2D spatial domain is not valid, here for loc's near the edge
  coords = expand.grid(1:dims[1], 1:dims[2])
  for(i in 1:nrow(coords)) {
    if(sum(coords[i,]) < 7) {
      zfield[i] = zvals + 1e-10
    }
  }
  
  # invalid location
  expect_identical(
    TestZConstrainedRookNeighborhood(dims = dims, x = c(1,1,0), 
                                     zfield = as.numeric(zfield), 
                                     zvals = zvals),
    matrix(0, ncol = 3, nrow = 1)[-1, , drop = FALSE]
  )
  
  # point along complex boundary
  expect_identical(
    TestZConstrainedRookNeighborhood(dims = dims, x = c(2,3,0), 
                                     zfield = as.numeric(zfield), 
                                     zvals = zvals),
    matrix(c(3,3,0,2,4,0), ncol = 3, byrow = TRUE)
  )
  
  # unconstrained points
  expect_identical(
    TestZConstrainedRookNeighborhood(dims = dims, x = c(7,7,0), 
                                     zfield = as.numeric(zfield), 
                                     zvals = zvals),
    TestRookNeighborhood(dims = dims, x = c(7,7,0))
  )
  
})

test_that('Validating constrained 3D rook neighborhoods on a grid', {
  
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
  
  # invalid location, height is too low
  expect_identical(
    TestZConstrainedRookNeighborhood(dims = dims, x = c(1,1,7), 
                                     zfield = as.numeric(zfield), 
                                     zvals = zvals),
    matrix(0, ncol = 3, nrow = 1)[-1, , drop = FALSE]
  )
  
  # all valid locations, neighborhood entirely above height constraints
  expect_identical(
    TestZConstrainedRookNeighborhood(dims = dims, x = c(1,1,9), 
                                     zfield = as.numeric(zfield), 
                                     zvals = zvals),
    TestRookNeighborhood(dims = dims, x = c(1,1,9))
  )
  
  # point along height boundary
  expect_identical(
    TestZConstrainedRookNeighborhood(dims = dims, x = c(1,1,8), 
                                     zfield = as.numeric(zfield), 
                                     zvals = zvals),
    TestRookNeighborhood(dims = dims, x = c(1,1,8))[-3,]
  )
  
  # unconstrained points
  expect_identical(
    TestZConstrainedRookNeighborhood(dims = dims, x = c(7,7,0), 
                                     zfield = as.numeric(zfield), 
                                     zvals = zvals),
    TestRookNeighborhood(dims = dims, x = c(7,7,0))
  )
  
})