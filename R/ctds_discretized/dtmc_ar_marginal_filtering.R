dtmc_ar_marginal_filtering = function(states, times, delta, dims, pred_times) {
  # Compute marginal filtering distributions
  #
  # Parameters:
  #  states - matrix of observation states/coordinates
  #  times - times at which observations are made
  #  delta - the discretization timestep to use
  #  dims - state space
  #  pred_times - times at which marginal filtering distributions should be 
  #   computed; must be aligned with the discretized times implied by delta
  
  
  # TODO: make this more compatible with the coefficient structure we will end 
  # up using
  theta = c(1,0)
  
  #
  # construct dummy height variables and information, to align with code
  #
  
  # append dummy variable and dimension information
  dims = c(dims, 1)
  states = cbind(states, 0)
  # height of domain surface
  zsurf = matrix(0, nrow = dims[1], ncol = dims[2])
  # height of vertical layer
  zval = 1
  
  # discretized timesteps
  tseq = seq(from = times[1], to = tail(times,1), by = delta)
  
  # verify prediction times are estimable wrt time discretization
  if(!all(pred_times %in% tseq)) {
    stop('All pred_times must be compatible with time discretization.')
  }
  
  # prediction times must be increasing to work with the c++ code
  pred_times = sort(pred_times)
  
  # map prediction times to discretized timesteps
  pred_steps = sapply(pred_times, function(t) {
    which(t == tseq)
  })
  
  # map observations to discretized timesteps
  obs_coords = matrix(NA, nrow = length(tseq), ncol = 3)
  for(ind in 1:length(times)) {
    obs_coords[which(tseq == times[ind]),] = states[ind,]
  }
  
  # uniform mass over neighborhood structure for initial state
  a0_prev_coords = dsmovetools:::TestZConstrainedRookNeighborhood(
    dims = dims, 
    x = states[1,], 
    zfield = zsurf, 
    zvals = zval
  )
  log_a0val = rep(1, nrow(a0_prev_coords))
  log_a0val = log(log_a0val / sum(log_a0val))
  
  # uniform mass over neighborhood structure for final state
  af_prev_coords = dsmovetools:::TestZConstrainedRookNeighborhood(
    dims = dims, 
    x = states[nrow(states),], 
    zfield = zsurf, 
    zvals = zval
  )
  log_afval = rep(1, nrow(af_prev_coords))
  log_afval = log(log_afval / sum(log_afval))
  
  # forward predictive distribution
  x_pred = dsmovetools:::ARPredDist(
    a0 = matrix(states[1,],
                nrow = nrow(a0_prev_coords),
                ncol = ncol(states), byrow = TRUE),
    a0_prev_coords = a0_prev_coords,
    pred_steps = pred_steps,
    obs_coords = obs_coords,
    log_a0val = log_a0val,
    dims = dims,
    surface_heights = zsurf, domain_heights = zval,
    log_self_tx = log(1 - theta[1] * delta), betaAR = theta[2]
  )
  
  # backward information filters
  x_bif = dsmovetools:::ARBackInfoFilteringDist(
    a0 = matrix(states[nrow(states),], 
                nrow = nrow(af_prev_coords), 
                ncol = ncol(states), byrow = TRUE),
    a0_prev_coords = af_prev_coords, 
    pred_steps = pred_steps,
    obs_coords = obs_coords, 
    log_a0val = log_afval, 
    dims = dims, 
    surface_heights = zsurf, domain_heights = zval, 
    log_self_tx = log(1 - theta[1] * delta), betaAR = theta[2]
  )
  
  # merge distributions, yielding marginal (assumes coords are in 3D)
  res = mapply(function(pred, bif) {
    inner_join(
      x = data.frame(pred),
      y = data.frame(bif),   
      by = paste('X', 1:6, sep = '')
    ) %>% 
      mutate(
        lp = X7.x + X7.y,
        lp = lp - log_sum(lp)
      ) %>% 
      select(X1_cur = X1, X2_cur = X2, X1_prev = X4, X2_prev = X5, lp)
  }, x_pred, x_bif, SIMPLIFY = FALSE)
  
  # collapse distributions according to the observable location
  res2 = lapply(res, function(d) {
    d %>% 
      group_by(X1 = X1_cur, X2 = X2_cur) %>% 
      summarise(lp = log_sum(lp), .groups = 'keep') %>% 
      ungroup()
  })
  
  # label with timepoints
  names(res) = paste('Prediction time: ', pred_times, sep = '')
  names(res2) = names(res)
  
  list(
    marginal_directional = res,
    marginal_location = res2
  )
}
