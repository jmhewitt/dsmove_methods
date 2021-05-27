context('SparseNDimArray.h C++ methods')

test_that('Validating basic input and output', {

  # number of dimensions
  dim = 3
  # number of coordinates
  n = 10
  
  # random coordinates and values
  coords = matrix(runif(n*dim), nrow = n)
  x = runif(n)
  
  # lexicographic ordering fn. from statnet.common package
  # author: Pavel N. Krivitsky
  order.matrix<-function(..., na.last = TRUE, decreasing=FALSE){
    x <- list(...)[[1L]]
    do.call(base::order,c(lapply(seq_len(ncol(x)), function(i) x[,i]), 
                          na.last=na.last, 
                          decreasing=decreasing)
    )
  }
  
  # lexicographic ordering for coordinates
  o = order.matrix(coords)
  
  # verify r/w data to the SparseNDimArray class in a predictable manner
  expect_identical(
    cbind(coords[o,,drop = FALSE], x[o]),
    TestSparseNdimArrayReadWrite(coords = coords, values = x)
  )
  
})

test_that('Validating input and output of coordinate pairs', {
  
  # number of dimensions
  dim = 2
  # number of coordinates
  n = 10
  
  # random coordinates and values
  coords1 = matrix(1:(dim), nrow = n, ncol = dim)
  coords2 = cbind(
    c(2,2,4,4,6,6,7,7,5,5),
    c(1,2,1,2,1,2,1,2,1,2)
  )
  x = runif(n)
  
  # lexicographic ordering fn. from statnet.common package
  # author: Pavel N. Krivitsky
  order.matrix<-function(..., na.last = TRUE, decreasing=FALSE){
    x <- list(...)[[1L]]
    do.call(base::order,c(lapply(seq_len(ncol(x)), function(i) x[,i]), 
                          na.last=na.last, 
                          decreasing=decreasing)
    )
  }
  
  # lexicographic ordering for coordinates
  o = order.matrix(cbind(coords1, coords2))
  reordered = cbind(coords1,coords2,x)[o,]
  dimnames(reordered) = NULL
  
  # verify r/w data to the SparseNDimArray class in a predictable manner
  expect_identical(
    TestBivariateSparseNdimArrayReadWrite(
      coords1 = coords1, coords2 = coords2, values = x
    ),
    reordered
  )
  
})
