#' Shortest-path CTDS imputation
#'
#' This method returns a CTDS imputation that uses the shortest path to connect 
#' observations.  Since the exact in-edge is unknown, this will be randomly 
#' selected.  The method can be used to initialize latent state for Gibbs 
#' sampling.  
#'   
#' @param ctds_struct Output of \code{build_ctds}; a representation of a 
#'   CTDS domain formatted for computing
#' @param states Sequence of nodes at which a CTDS is observed
#' @param times Times at which observations were recorded
#' @param output.einds \code{TRUE} to include the edge indices associated with 
#'   each transition
#' 
#' @importFrom igraph shortest_paths
#' 
#' @export
#' 
ctds.shortest_impute = function(states, times, ctds_struct, 
                                output.einds = FALSE) {
 
  # number of observations
  n = length(states)
  
  # build imputation by segments
  imputed = do.call(rbind, lapply(1:(n-1), function(i) {
    if(states[i] == states[i+1]) {
      # no transition was observed; placeholder to duplicate previous state
      states.imputed = NA
    } else  {
      # states associated with one shortest directed path between observations
      states.imputed = as.numeric(shortest_paths(graph = ctds_struct$graph, 
                                                 from = states[i], 
                                                 to = states[i+1], mode = 'out', 
                                                 output = 'epath')$epath[[1]])
    }
    # placeholder for fully equally-spaced transition times
    times.imputed = rep(NA, length(states.imputed))
    times.imputed[length(times.imputed)] = times[i+1]
    # package segments
    data.frame(states = states.imputed, times = times.imputed)
  }))
  
  
  #
  # infill state at observation times where no transitions were observed
  #
  
  # indices which need infill
  repeated.states = which(is.na(imputed$states))
  if(length(repeated.states) > 0) {
    # indices which can provide infill
    valid.states = setdiff(1:length(imputed$states), repeated.states)
    # indices to be used for infill
    infill.states = sapply(repeated.states, function(r) {
      valid.states[max(which(valid.states < r))]
    })
    # infill missing state
    imputed$states[repeated.states] = imputed$states[infill.states]
  }
  
  
  #
  # finish imputation
  #
  
  # sample a random starting state
  init_inds = ctds_struct$in_edges_inds[[states[1]]]
  if(length(init_inds) > 1) {
    init_state = sample(x = init_inds, size = 1)
  } else {
    init_state = init_inds
  }
  
  # prepend starting state and time
  imputed = list(states = c(init_state, imputed$states), 
                 times = c(times[1], imputed$times))
  
  # linearly interpolate missing times
  time.inds = which(!is.na(imputed$times))
  for(i in 1:(length(time.inds)-1)) {
    inds = time.inds[i]:time.inds[i+1]
    n.inds = length(inds)
    imputed$times[inds] = seq(from = imputed$times[inds[1]], 
                              to = imputed$times[inds[n.inds]],
                              length.out = n.inds)
  }
  
  # remove self-transitions from imputation
  true.tx = c(TRUE, diff(imputed$states) != 0)
  imputed$states = imputed$states[true.tx]
  imputed$times = imputed$times[true.tx]
  
  if(output.einds) {
    imputed$einds = imputed$states
  }
  
  # back-transform edges to nodes, for plotting, etc.
  imputed$states = ctds_struct$edge_df$to[imputed$states]
  
  imputed$durations = diff(imputed$times)
  
  class(imputed) = 'ctds_realization'
  
  imputed
}
