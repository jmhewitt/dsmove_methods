joint_prior = function(beta_loc, beta_ar, dims, penalty_rate = 1/.5,
                       log = FALSE) {
  # Prior that penalizes the expected number of reversals per unit time in 
  # paths.  The function can accept a vector of beta_loc (intercept terms) 
  # values and a vector of beta_ar values, so that the joint prior can be 
  # studied more easily.
  # 
  # Parameters:
  #  penalty_rate - 1/expected number of reversals per unit time to anticipate
  #  dims - dimensions of spatial domain
  #  beta_loc - covariates driving transition rates; assumed that intercept is 
  #    stored in beta_loc[1] and defines a meaningful "base" movement speed

  # index of central-most coordinate, 0-referenced for c++
  ind.center = floor(dims/2)
  
  # first neighbor of central-most coordinate, 0-referenced from c++
  nbr = dsmovetools:::TestRookNeighborhood(dims = dims, x = ind.center)[
    1, , drop = FALSE
  ]

  # principal feature of joint process used to define joint prior
  expected_reversals = exp(mapply(function(beta_loc, beta_ar) {
    # (log) expected num. transitions per unit time (mean for exp dist'n.)
  - beta_loc[1] +
    # (log) scaled by expected proportion of reversals (i.e., via reversal prob)
      dsmovetools:::TxModelLd(
        cur_loc = ind.center, prev_loc = nbr, dims = dims, betaAR = beta_ar, 
        dst_loc = nbr
      )
  }, beta_loc, beta_ar))
  
  # penalize expected number of reversals per unit time
  dexp(x = expected_reversals, log = log, rate = penalty_rate)
}