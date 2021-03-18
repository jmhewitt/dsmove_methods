context('ffrw.cpp C++ methods')

test_that('Validating random walk forward filtering on a grid: 1D simple', {

  # number of dimensions
  ndim = 1
  # number of coordinates in each dimension
  dims = 20
  
  # initial probability mass
  coords = matrix(c(0,1), ncol = ndim)
  p0 = rep(.5,2)
  
  # 1-step diffusion
  expect_identical(
    TestFFRW(a0coords = coords, a0values = p0, dims = dims, steps = 1),
    matrix(c(0,.25, 1,.5, 2,.25), ncol = 2, byrow = TRUE)
  )
  
})

test_that('Validating random walk forward filtering on a grid: 1D, log-scale', {
  
  # number of dimensions
  ndim = 1
  # number of coordinates in each dimension
  dims = 20
  # number of steps to take
  nsteps = 100
  
  # initial probability mass
  coords = matrix(c(0,1), ncol = ndim)
  p0 = rep(.5,2)
  lp0 = log(p0)
  
  # standard evaluation of diffused mass
  af_std = TestFFRWLight(a0coords = coords, a0values = p0, dims = dims, 
                         steps = nsteps)
  
  # diffused probability mass accumulated on log scale
  af_alt =TestFFRWLightLog(a0coords = coords, log_a0values = lp0, dims = dims, 
                           steps = nsteps)
  
  # back-transform probabilities to (natural) linear scale
  af_alt[,ncol(af_alt)] = exp(af_alt[,ncol(af_alt)])
  
  # n-step diffusion
  expect_equal(af_std, af_alt)
})

test_that('Validating random walk forward filtering on a grid: 2D simple', {

  # number of dimensions
  ndim = 2
  # number of coordinates in each dimension
  dims = c(10,10)
  # coordinates
  coords = expand.grid(x = 0:(dims[1]-1), y = 0:(dims[2]-1))
  
  # template for neighbors
  nbr.template = matrix(c(-1,0, 1,0, 0,-1, 0,1), ncol = 2, byrow = TRUE)
  
  # 1-step transition matrix
  P = matrix(0, nrow = nrow(coords), ncol = nrow(coords))
  for(i in 1:nrow(coords)) {
    for(j in 1:nrow(nbr.template)) {
      nbr = coords[i,] + nbr.template[j,]
      nbr.col = which(
        (unlist(nbr[1]) == coords[,1]) & (unlist(nbr[2]) == coords[,2])
      )
      P[i,nbr.col] = 1
    }
  }
  P = sweep(x = P, MARGIN = 1, STATS = rowSums(P), FUN = '/')
  
  # lexicographic ordering fn. from statnet.common package
  # author: Pavel N. Krivitsky
  order.matrix<-function(..., na.last = TRUE, decreasing=FALSE){
    x <- list(...)[[1L]]
    do.call(base::order,c(lapply(seq_len(ncol(x)), function(i) x[,i]), 
                          na.last=na.last, 
                          decreasing=decreasing)
    )
  }
  
  # extract probabilities
  pstep = function(steps, loc) {
    # diffused transition matrix
    Pm = expm::`%^%`(P, steps)
    # row in which to find output
    loc.row = which(
      (unlist(loc[1]) == coords[,1]) & (unlist(loc[2]) == coords[,2])
    )
    # output coordinates
    nnz = which(Pm[loc.row,] > 0)
    res = cbind(coords[nnz,], Pm[loc.row, nnz])
    o = order.matrix(res)
    res[o,]
  }
  
  # initial probability mass
  x0 = matrix(c(4,8), ncol = ndim)
  p0 = 1
  
  # 1-step diffusion
  nsteps = 1
  expect_equivalent(
    TestFFRW(a0coords = x0, a0values = p0, dims = dims, steps = nsteps),
    as.matrix(pstep(steps = nsteps, loc = x0))
  )
  
  # 7-step diffusion
  nsteps = 7
  expect_equivalent(
    TestFFRW(a0coords = x0, a0values = p0, dims = dims, steps = nsteps),
    as.matrix(pstep(steps = nsteps, loc = x0))
  )
  
})

test_that('Validating random walk forward filtering on a grid: 2D multiple', {
  
  # NOTE: This test is slow because the P matrix benchmark is computationally
  #  expensive both to construct and use in dense matrix multiplications.
  #  Test is made faster by reducing dims from c(100,100) to c(25,25)
  
  # number of dimensions
  ndim = 2
  # number of coordinates in each dimension
  dims = c(25,25)
  # coordinates
  coords = expand.grid(x = 0:(dims[1]-1), y = 0:(dims[2]-1))
  
  # template for neighbors
  nbr.template = matrix(c(-1,0, 1,0, 0,-1, 0,1), ncol = 2, byrow = TRUE)
  
  # 1-step transition matrix
  P = matrix(0, nrow = nrow(coords), ncol = nrow(coords))
  for(i in 1:nrow(coords)) {
    for(j in 1:nrow(nbr.template)) {
      nbr = coords[i,] + nbr.template[j,]
      nbr.col = which(
        (unlist(nbr[1]) == coords[,1]) & (unlist(nbr[2]) == coords[,2])
      )
      P[i,nbr.col] = 1
    }
  }
  P = sweep(x = P, MARGIN = 1, STATS = rowSums(P), FUN = '/')
  
  # lexicographic ordering fn. from statnet.common package
  # author: Pavel N. Krivitsky
  order.matrix<-function(..., na.last = TRUE, decreasing=FALSE){
    x <- list(...)[[1L]]
    do.call(base::order,c(lapply(seq_len(ncol(x)), function(i) x[,i]), 
                          na.last=na.last, 
                          decreasing=decreasing)
    )
  }
  
  # extract probabilities
  pstep = function(steps, loc, p0) {
    # initial mass vector, in state-space format
    x0 = numeric(nrow(P))
    loc.inds = apply(loc, 1, function(r) {
      which(
        (unlist(r[1]) == coords[,1]) & (unlist(r[2]) == coords[,2])
      )
    })
    x0[loc.inds] = p0
    # diffused transition matrix pre-scaled by initial mass
    Pm = x0 %*% expm::`%^%`(P, steps)
    # aggregate output
    nnz = which(Pm > 0)
    res = cbind(coords[nnz,], Pm[nnz])
    # sort and return
    o = order.matrix(res)
    res[o,]
  }
  
  # initial probability mass
  x0 = matrix(c(4,8, 1,0), ncol = ndim, byrow = TRUE)
  p0 = runif(2)
  p0 = p0 / sum(p0)
  
  # 1-step diffusion
  nsteps = 1
  expect_equivalent(
    TestFFRW(a0coords = x0, a0values = p0, dims = dims, steps = nsteps),
    as.matrix(pstep(steps = nsteps, loc = x0, p0 = p0))
  )
  
  # 700-step diffusion
  nsteps = 700
  expect_equivalent(
    TestFFRW(a0coords = x0, a0values = p0, dims = dims, steps = nsteps),
    as.matrix(pstep(steps = nsteps, loc = x0, p0 = p0))
  )
})

test_that('Validating random walk forward filtering on a grid: 2D two-step', {
  
  # number of dimensions
  ndim = 2
  # number of coordinates in each dimension
  dims = c(10,10)
  # coordinates
  coords = expand.grid(x = 0:(dims[1]-1), y = 0:(dims[2]-1))
  
  # template for neighbors
  nbr.template = matrix(c(-1,0, 1,0, 0,-1, 0,1), ncol = 2, byrow = TRUE)
  
  # 1-step transition matrix
  P = matrix(0, nrow = nrow(coords), ncol = nrow(coords))
  for(i in 1:nrow(coords)) {
    for(j in 1:nrow(nbr.template)) {
      nbr = coords[i,] + nbr.template[j,]
      nbr.col = which(
        (unlist(nbr[1]) == coords[,1]) & (unlist(nbr[2]) == coords[,2])
      )
      P[i,nbr.col] = 1
    }
  }
  P = sweep(x = P, MARGIN = 1, STATS = rowSums(P), FUN = '/')
  
  # lexicographic ordering fn. from statnet.common package
  # author: Pavel N. Krivitsky
  order.matrix<-function(..., na.last = TRUE, decreasing=FALSE){
    x <- list(...)[[1L]]
    do.call(base::order,c(lapply(seq_len(ncol(x)), function(i) x[,i]), 
                          na.last=na.last, 
                          decreasing=decreasing)
    )
  }
  
  # extract probabilities
  pstep = function(steps, loc, p0) {
    # initial mass vector, in state-space format
    x0 = numeric(nrow(P))
    loc.inds = apply(loc, 1, function(r) {
      which(
        (unlist(r[1]) == coords[,1]) & (unlist(r[2]) == coords[,2])
      )
    })
    x0[loc.inds] = p0
    # diffused transition matrix pre-scaled by initial mass
    Pm = x0 %*% expm::`%^%`(P, steps)
    # aggregate output
    nnz = which(Pm > 0)
    res = cbind(coords[nnz,], Pm[nnz])
    # sort and return
    o = order.matrix(res)
    res[o,]
  }
  
  # initial probability mass
  x0 = matrix(c(4,8, 1,0), ncol = ndim, byrow = TRUE)
  p0 = runif(2)
  p0 = p0 / sum(p0)
  
  
  # 34-step diffusion
  nsteps = 38
  
  # test c++ implementation using two calls to filtering methods
  x1 = TestFFRW(a0coords = x0, a0values = p0, dims = dims, steps = nsteps/2)
  xf = TestFFRW(a0coords = x1[,1:ndim], a0values = x1[,-(1:ndim)], dims = dims, 
                steps = nsteps/2)
  
  expect_equivalent(
    xf, 
    as.matrix(pstep(steps = nsteps, loc = x0, p0 = p0))
  )
})

test_that('Validating feasibility on a large, 3D grid', {
  
  # number of dimensions
  ndim = 3
  # number of coordinates in each dimension
  dims = c(1128, 1287, 300)
  
  # initial probability mass
  x0 = matrix(c(4,8,0, 1,0,150), ncol = ndim, byrow = TRUE)
  p0 = runif(2)
  p0 = p0 / sum(p0)
  
  nsteps = 100
  xf = TestFFRW(a0coords = x0, a0values = p0, dims = dims, steps = nsteps)
  
  expect_equal(sum(xf[,4]), 1)
})

test_that('Validating feasibility on a large, 3D grid with alternate methods', {

  # This test is slow since the "Light" algorithm does lots of memory swapping,
  # but the memory overhead is mostly fixed.
  
  # number of dimensions
  ndim = 3
  # number of coordinates in each dimension
  dims = c(1128, 1287, 300)

  # initial probability mass
  x0 = matrix(c(4,8,0, 1,0,150), ncol = ndim, byrow = TRUE)
  p0 = runif(2)
  p0 = p0 / sum(p0)

  # use smaller number of steps for faster package testing
  # nsteps = 200
  nsteps = 100
  xf = TestFFRWLightLog(a0coords = x0, log_a0values = log(p0), dims = dims, 
                        steps = nsteps)

  expect_equal(sum(exp(xf[,4])), 1)
})

test_that('Validating lightweight implementation of forward filtering', {
  
  # number of dimensions
  ndim = 2
  # number of coordinates in each dimension
  dims = c(100, 100)
  
  # initial probability mass
  x0 = matrix(c(4,8, 1,0), ncol = ndim, byrow = TRUE)
  p0 = runif(2)
  p0 = p0 / sum(p0)
  
  nsteps = 50
  
  expect_equal(
    TestFFRWLight(a0coords = x0, a0values = p0, dims = dims, steps = nsteps),
    TestFFRW(a0coords = x0, a0values = p0, dims = dims, steps = nsteps)
  )
  
})

test_that('Validating forward filtering to destination is successful', {
  
  # number of dimensions
  ndim = 3
  # numer of coordinates in each dimension
  dims = c(100,100,1)
  
  # initial and final locations
  x0 = matrix(c(5,5,0), nrow = 1)
  xf = matrix(c(7,7,0), nrow = 1)
  xf_far = matrix(c(25,25,0), nrow = 1)
  xf_close = matrix(c(6,5,0), nrow = 1)
  
  # minimum number of steps required for diffusion
  nsteps  = 10
  
  # initial probability mass 
  p0_log = 0
  
  # height of vertical layer
  zval = 1
  
  # height of domain surface
  zsurf = matrix(0, nrow = dims[1], ncol = dims[2])
  
  # diffuse from source to near destination
  af = FFRWLogConstrainedDst(
    a0coords = x0, dstcoords = xf, log_a0values = p0_log, dims = dims, 
    steps = nsteps, max_steps = 1e2, surface_heights = zsurf, 
    domain_heights = zval
  )
  
  # diffuse from source to near destination, maintaining reachability info. 
  af_reachable = FFRWLogConstrainedDstReachable(
    a0coords = x0, dstcoords = xf, log_a0values = p0_log, dims = dims, 
    steps = nsteps, max_steps = 1e2, surface_heights = zsurf, 
    domain_heights = zval
  )
  
  # diffuse from source to far destination
  af_far = FFRWLogConstrainedDst(
    a0coords = x0, dstcoords = xf_far, log_a0values = p0_log, dims = dims, 
    steps = nsteps, max_steps = 1e2, surface_heights = zsurf, 
    domain_heights = zval
  )
  
  # diffuse from source to near destination
  af_close = FFRWLogConstrainedDst(
    a0coords = x0, dstcoords = xf_close, log_a0values = p0_log, dims = dims, 
    steps = nsteps, max_steps = 1e2, surface_heights = zsurf, 
    domain_heights = zval
  )
  
  
  dst_attainable = function(pvecs, dst) {
    # TRUE when dst location has nonzero mass at a diffusion
    sapply(pvecs, function(p) {
      # see if any of the locations in diffusion match dst
      any(apply(p, 1, function(r) {
        all(r[1:length(dims)] == dst)
      }))
    })
  }
  
  # all diffusions have non-zero probability of reaching dst location
  expect_true(
    all(
      any(dst_attainable(af, xf)),
      any(dst_attainable(af_close, xf_close)),
      any(dst_attainable(af_far, xf_far))
    )
  )
  
  # the diffusions that reach destination are properly recorded
  expect_identical(
    which(dst_attainable(af_reachable[[1]], xf)) - 1,
    af_reachable[[2]]
  )
  
  # all diffusions have required minimum length
  expect_true(
    all(
      length(af) > nsteps,
      length(af_close) > nsteps,
      length(af_far) > nsteps
    )
  )
  
  # close diffusion does not exceed minimum length
  expect_true(
    length(af_close) == (nsteps + 1)
  )
  
})
