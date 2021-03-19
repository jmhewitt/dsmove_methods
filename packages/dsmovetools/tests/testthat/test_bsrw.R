context('bsrw.cpp C++ methods')

test_that('Validate bridged random walk sampling behaves as expected', {
  
  # size of discrete spatial domain
  dims = c(100,100, 1)
  
  # coordinates in 0-indexed format, for c++
  x0 = c(50, 50, 1) - 1
  xf = c(70, 70, 1) - 1
  # x0 = c(1, 1, 1) - 1
  # xf = x0
  
  # minimum number of steps to take
  steps = 100
  
  # spatial domain is unconstrained
  zvals = 1
  zfield = matrix(0, nrow = dims[1], ncol = dims[2])
  
  # sample path
  path = SampleConstrainedBridgedRWPath(
    a0coords = x0, dstcoords = xf, dims = dims, steps = steps, max_steps = 1e3, 
    surface_heights = as.numeric(zfield), domain_heights = zvals
  )
  
  # path satisfies length constraint
  expect_gte(
    nrow(path),
    steps
  )
  
  # path follows rook neighborhood geometry, with exactly one change per step
  expect_identical(
    rep(1, nrow(path) - 1),
    rowSums(abs(apply(path, 2, diff)))
  )
  
  # path starts at x0
  expect_identical(
    x0,
    path[1,]
  )
  
  # path ends at xf
  expect_identical(
    xf,
    path[nrow(path),]
  )
  
})
