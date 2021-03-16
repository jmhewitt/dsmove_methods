#' Forward simulation of a simple CTDS model
#' 
#' Forward-samples from the generating distributions for the CTDS model until 
#' either the object is observable at time \code{tf}, or a maximum number of 
#' transitions is realized (\code{max.steps}).
#'   
#' @param ctds_struct Output of \code{build_ctds}; a representation of a CTDS 
#'   domain formatted for computing
#' @param beta_loc coefficients for location-based drivers of movement
#' @param beta_dir coefficients for direction-based drivers of movement
#' @param beta_ar strength of directional-persistence parameter
#' @param v0 the node from which sampling should start
#' @param v0.last the last-exited node before \code{v0}, if using directional
#'   persistence---e.g., if \code{beta_ar != 0}.
#' @param t0 the time at which sampling begins
#' @param tf the time after which sampling should end
#' @param max.steps the maximum number of transitions to simulate
#' @param weibull.shape the shape parameter used to generate holding times; 
#'   if \code{weibull.shape==1}, then the distribution is exactly the 
#'   exponential distribution
#'   
#' 
#' @importFrom stats rexp
#' 
#' @example examples/ctds.quicksim.R
#' 
#' @export
#' 
ctds.quicksim = function(dims, beta_loc, beta_ar, v0, t0, tf, 
                         max.steps = Inf, v0.last = NULL) {
  
  if(tf < t0) {
    stop('tf occurs before t0; no sampling is required!')
  }
  
  if(is.null(v0.last)) {
    # get neighborhood of initial position
    nbhd = TestRookNeighborhood(dims = dims, x = v0)
    # randomly sample a starting location from neighbors
    v0.last = nbhd[sample(x = nrow(nbhd), size = 1),]
  }
  
  # initialize output
  vpath = matrix(c(v0.last, v0), nrow = 2, byrow = TRUE)
  times = t0
  
  # initialize state tracking
  v = v0
  vprev = v0.last
  tcurrent = t0
  
  # extract ar parameter
  beta_ar_computational = ifelse(is.null(beta_ar), 0, beta_ar)
  
  # convert beta_loc to rate
  tx_rate = exp(beta_loc)

  # forward sample
  it = 1
  while(it < max.steps) {
    
    it = it + 1
    
    # sample and update transition times
    dur = rexp(n = 1, rate = tx_rate)
    tcurrent = tcurrent + dur
    times = c(times, tcurrent)
    
    # sample and update state
    vnew = TxModelSample(cur_loc = v, prev_loc = vprev, dims = dims, 
                         betaAR = beta_ar_computational)
    vpath = rbind(vpath, vnew)
    vprev = v
    v = vnew
    
    # stop sampling if time is exceeded 
    if(tcurrent >= tf) {
      break
    }
    
  }
  
  #
  # package results
  #
  
  rownames(vpath) = NULL
  
  res = list(states = vpath, times = times, durations = diff(times))
  
  class(res) = 'ctds_realization'
  
  res
}
