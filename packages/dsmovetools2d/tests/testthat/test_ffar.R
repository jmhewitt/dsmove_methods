context('ffar.h C++ methods')


test_that('Validating self-transition filtering on a grid', {
  
  # define grid
  lons = seq(from = 0, to = 40, length.out = 100)
  lats = seq(from = 50, to = 70, length.out = 100)
  
  # arbitrary height-field
  surface_heights = runif(length(lons) * length(lats))
  
  # coordinates (lon,lat) format
  init_dsts = rbind(
    c(50,50)
  )
  init_srcs = rbind(
    c(49,50)
  )
  
  # initial log probability
  init_log_probs = rep(log(1/nrow(init_dsts)), nrow(init_dsts))
  
  # model parameters
  betaAR = 1
  log_self_tx = log(.7)
  
  # diffusion count
  nsteps = 7
  
  cpp_out = FF_DTMC(
    lons = lons, lats = lats, surface_heights = surface_heights,
    init_dsts = init_dsts, init_srcs = init_srcs,
    init_log_probs = init_log_probs, steps = 7, log_self_tx = log_self_tx,
    betaAR = betaAR
  )
  
  
  #
  # verification
  #
  
  # iterative probability diffusion, implemented in R
  pstep = function(x0, x0_prev, log_a0val, dims, steps, log_self_tx, betaAR) {
    
    # log of transition probability
    log_tx = log(1 - exp(log_self_tx))
    
    # initialize diffused probability vector with initial conditions
    res = matrix(c(x0, x0_prev, log_a0val), nrow = 1)
    
    # run diffusions
    for(i in 1:steps) {
      # unaggregated diffusion probabilities within step
      tmp = matrix(nrow = 0, ncol = ncol(res))
      # forward-diffuse mass from each atom in the current diffused prob. vector
      for(j in 1:nrow(res)) {
        # reachable destinations and forward probabilities
        p = dsmovetools:::TxModelParams(
          cur_loc = res[j, 1:length(x0)], 
          prev_loc = res[j, length(x0) + 1:length(x0)], 
          dims = dims, betaAR = betaAR
        )
        # reweight forward probabilities by source and self-tx. probs
        lp = p[, ncol(p)] + log_tx + res[j, ncol(res)]
        # format output
        tmp = rbind(
          # append with mass diffused from other atoms
          tmp,
          # diffused mass from leaving current atom
          cbind(p[,-ncol(p)], 
                matrix(data = res[j,1:length(x0)], 
                       nrow = nrow(p), 
                       ncol = length(x0), byrow = TRUE), 
                lp),
          # diffused mass from staying in current atom (self-transition)
          matrix(c(res[j,-ncol(res)], log_self_tx + res[j, ncol(res)]), 
                 nrow = 1)
        )
      }
      # aggregate transition probabilities
      res = data.frame(tmp) %>% 
        dplyr::group_by(V1, V2, V3, V4, V5, V6) %>% 
        dplyr::summarise(lp = log(sum(exp(lp)))) %>% 
        dplyr::ungroup() %>% 
        as.matrix()
    }
    
    # package results
    dimnames(res) = NULL
    res
  }
  
  # forward-filter via R, rearranging output to match column order
  af.alt = pstep(
    x0 = cbind(init_dsts, 0), x0_prev = cbind(init_srcs, 0), 
    dims = c(length(lons), length(lats), 1), log_a0val = init_log_probs, 
    steps = nsteps, log_self_tx = log_self_tx, betaAR = betaAR
  )[,c(4,5,1,2,7)]
  
  # lexicographic ordering fn. from statnet.common package
  # author: Pavel N. Krivitsky
  order.matrix<-function(..., na.last = TRUE, decreasing=FALSE){
    x <- list(...)[[1L]]
    do.call(base::order,c(lapply(seq_len(ncol(x)), function(i) x[,i]), 
                          na.last=na.last, 
                          decreasing=decreasing)
    )
  }
  
  expect_equivalent(
    af.alt[order.matrix(af.alt),],
    cpp_out[order.matrix(cpp_out),]
  )
  
})