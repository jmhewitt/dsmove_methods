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
  sampled_path = SampleConstrainedBridgedRWPath(
    a0coords = x0, dstcoords = xf, dims = dims, steps = steps, max_steps = 1e3, 
    surface_heights = as.numeric(zfield), domain_heights = zvals
  )
  
  # path satisfies length constraint
  expect_gte(
    nrow(sampled_path),
    steps
  )
  
  # path follows rook neighborhood geometry, with exactly one change per step
  expect_identical(
    rep(1, nrow(sampled_path) - 1),
    rowSums(abs(apply(sampled_path, 2, diff)))
  )
  
  # path starts at x0
  expect_identical(
    x0,
    sampled_path[1,]
  )
  
  # path ends at xf
  expect_identical(
    xf,
    sampled_path[nrow(sampled_path),]
  )
  
})


test_that('Validate sampling families of bridged random walks', {
  
  set.seed(2021)
  
  # number of paths to sample
  n = 10
  
  # size of discrete spatial domain
  dims = c(100,100, 1)
  
  # coordinates in 0-indexed format, for c++
  x0 = c(50, 50, 1) - 1
  xf = c(70, 70, 1) - 1
  
  # minimum number of steps to take
  steps = 100
  
  # spatial domain is unconstrained
  zvals = 1
  zfield = matrix(0, nrow = dims[1], ncol = dims[2])
  
  # sample family of bridged paths
  path_fam = SampleConstrainedBridgedRWPathFamily(
    a0coords = x0, dstcoords = xf, dims = dims, steps = steps, max_steps = 1e3, 
    surface_heights = as.numeric(zfield), domain_heights = zvals, n = n
  )
  
  # sampled paths should have different lengths
  expect_gt(
    length(unique(sapply(path_fam$path, function(sampled_path) {
      nrow(sampled_path)
    }))),
    1
  )
  
  # validate each path
  invisible(lapply(path_fam$path, function(sampled_path) {
    # path follows rook neighborhood geometry, with exactly one change per step
    expect_identical(
      rep(1, nrow(sampled_path) - 1),
      rowSums(abs(apply(sampled_path, 2, diff)))
    )
    
    # path starts at x0
    expect_identical(
      x0,
      sampled_path[1,]
    )
    
    # path ends at xf
    expect_identical(
      xf,
      sampled_path[nrow(sampled_path),]
    )
  }))
  
  
  #
  # validate that paths of equal length are uniformly sampled
  #
  
  # number of times each path length appears
  path_lengths = sapply(path_fam$path, nrow)
  path_length_counts = table(path_lengths)
  
  # test requires that some path lengths are seen multiple times 
  expect_true(
    any(path_length_counts > 1)
  )
  
  # extract duplicate path lengths
  duplicate_lengths = as.numeric(
    names(path_length_counts)[which(path_length_counts > 1)]
  )
  
  # ensure paths of equal length have same sampling weight
  for(len in duplicate_lengths) {
    # ensure paths of equal length are not all the same
    expect_false(
      identical(
        path_fam$path[which(path_lengths == len)][[1]],
        path_fam$path[which(path_lengths == len)][[2]]
      )
    )
    # sampling weights for equal-length paths
    len_wts = path_fam$log_weights[which(path_lengths == len)]
    # crudely sampling weights are identical
    expect_equal(0, sd(len_wts))
  }
  
})
