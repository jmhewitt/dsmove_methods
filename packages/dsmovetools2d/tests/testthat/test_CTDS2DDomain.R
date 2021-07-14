context('CTDS2DDomain.h C++ methods')

test_that('Validating basic input and output', {
  
  # lexicographic ordering fn. from statnet.common package
  # author: Pavel N. Krivitsky
  order.matrix<-function(..., na.last = TRUE, decreasing=FALSE){
    x <- list(...)[[1L]]
    do.call(base::order,c(lapply(seq_len(ncol(x)), function(i) x[,i]), 
                          na.last=na.last, 
                          decreasing=decreasing)
    )
  }
  
  # define grid
  lons = seq(from = 0, to = 40, length.out = 100)
  lats = seq(from = 50, to = 70, length.out = 100)
  
  # arbitrary height-field
  surface_heights = runif(length(lons) * length(lats))
  
  # coordinates
  init_dsts = rbind(
    c(30,30),
    c(25,19)
  )
  init_srcs = rbind(
    c(29,30),
    c(25,20)
  )
  
  # random probabilities
  log_probs = log(runif(nrow(init_dsts)))
  
  # output from c++
  cpp_out = TestCTDS2DDomainIO(
    lons = lons, lats = lats, surface_heights = surface_heights,
    init_dsts = init_dsts, init_srcs = init_srcs, log_probs = log_probs
  )
  
  # format for verification data
  r_check = cbind(init_srcs, init_dsts, log_probs)
  
  expect_equivalent(
    cpp_out[order.matrix(cpp_out),],
    r_check[order.matrix(r_check),]
  )
  
})

