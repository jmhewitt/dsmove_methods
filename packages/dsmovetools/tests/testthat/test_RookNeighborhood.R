context('RookNeighborhood.h C++ methods')

test_that('Validating rook neighborhoods on a grid', {

  #
  # 1D neighborhood tests
  #
  
  # within-boundary point
  expect_identical(
    TestRookNeighborhood(dims = 16, x = 1),
    matrix(c(0, 2), ncol = 1)
  )
  
  # lower-boundary point
  expect_identical(
    TestRookNeighborhood(dims = 16, x = 0),
    matrix(c(1), ncol = 1)
  )
  
  # upper-boundary point
  expect_identical(
    TestRookNeighborhood(dims = 16, x = 15),
    matrix(14, ncol = 1)
  )
  
  # point outside boundary
  expect_identical(
    TestRookNeighborhood(dims = 16, x = 20),
    matrix(0, ncol = 1, nrow = 1)[-1, , drop = FALSE]
  )
  
  #
  # 2D neighborhood tests
  #
  
  dims = c(100,100)

  # points outside boundaries
  expect_identical(
    rbind(
      TestRookNeighborhood(dims = dims, x = c(0,-1)),
      TestRookNeighborhood(dims = dims, x = c(300,0)),
      TestRookNeighborhood(dims = dims, x = c(300,-4))
    ),
    matrix(0, ncol = 2, nrow = 1)[-1, , drop = FALSE]
  )
  
  # lower corner
  expect_identical(
    TestRookNeighborhood(dims = dims, x = c(0,0)),
    matrix(c(1,0,0,1), ncol = 2, byrow = TRUE)
  )
  
  # upper corner
  expect_identical(
    TestRookNeighborhood(dims = dims, x = c(99,99)),
    matrix(c(98,99, 99,98), ncol = 2, byrow = TRUE)
  )
  
  # inner cell
  expect_identical(
    TestRookNeighborhood(dims = dims, x = c(12,20)),
    matrix(c(11,20, 12,19, 13,20, 12,21), ncol = 2, byrow = TRUE)
  )
  
})
