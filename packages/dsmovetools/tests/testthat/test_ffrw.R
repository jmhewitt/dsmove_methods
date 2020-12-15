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