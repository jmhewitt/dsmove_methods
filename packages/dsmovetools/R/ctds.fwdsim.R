#' Forward simulation of a CTDS object
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
#' @export
#' 
ctds.fwdsim = function(ctds_struct, beta_loc, beta_dir, v0, t0, tf, 
                       max.steps = Inf, beta_ar = NULL, v0.last = NULL,
                       weibull.shape = 1) {
  
  if(tf < t0) {
    stop('tf occurs before t0; no sampling is required!')
  }
  
  # sampling is conducted entirely in edge-space
  
  if(is.null(v0.last)) {
    # randomly sample a starting edge
    e0 = sample(x = ctds_struct$in_edges_inds[[v0]], size = 1)
  } else {
    stop('Non-NULL v0.last argument not implemented yet.')
  }
  
  # initialize output
  edge_path = e0
  times = t0
  
  # initialize state tracking
  e = e0
  tcurrent = t0
  
  # precompute common simulation elements
  outedges_by_loc = do.call(c, ctds_struct$out_edges_inds)
  loc_start = c(1, 1 + cumsum(ctds_struct$out_degree))
  
  # extract ar parameter
  beta_ar_computational = ifelse(is.null(beta_ar), 0, beta_ar)

  # forward sample
  it = 1
  while(it < max.steps) {
    
    it = it + 1
    
    # extract current location in spatial domain
    v = ctds_struct$edge_df$to[e]

    # edges that can be transitioned to
    edges = ctds_struct$out_edges_inds[[v]]

    # destination locations associated with transition edges
    locs = ctds_struct$edge_df$to[edges]

    # get local transition parameters
    A = local_generator(locs = c(v, locs), row_edges = e, 
                        col_edges = c(e, edges), 
                        tolocs_by_edge = ctds_struct$edge_df$to,
                        fromlocs_by_edge = ctds_struct$edge_df$from, 
                        outedges_by_loc = outedges_by_loc, 
                        loc_start = loc_start, Xloc = ctds_struct$Xloc, 
                        betaLoc = beta_loc, 
                        Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df), 
                                      ncol = 1), 
                        betaDir = beta_dir, W = ctds_struct$w_ij, 
                        betaAR = beta_ar_computational)
    
    # sample and update transition times
    dur = rweibull(n = 1, scale = -1/A[1], shape = weibull.shape)
    tcurrent = tcurrent + dur
    times = c(times, tcurrent)
    
    # sample and update state
    if(length(edges) == 1) {
      e = edges
    } else {
      e = sample(x = edges, size = 1, prob = A[-1])
    }
    edge_path = c(edge_path, e)
    
    # stop sampling if time is exceeded 
    if(tcurrent >= tf) {
      break
    }
    
  }
  
  # extract state path
  states = ctds_struct$edge_df$to[edge_path]
  
  #
  # package results
  #
  
  res = list(states = states, times = times, durations = diff(times),
             edges = edge_path)
  
  class(res) = 'ctds_realization'
  
  res
}
