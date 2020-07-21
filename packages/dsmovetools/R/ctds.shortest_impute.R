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
#' 
#' @importFrom igraph shortest_paths
#' 
#' @export
#' 
ctds.shortest_impute = function(states, times, ctds_struct) {
 
  # number of observations
  n = length(states)
  
  # edge-states associated with shortest directed paths between observations
  imputed = lapply(1:(n-1), function(i) {
    as.numeric(shortest_paths(graph = ctds_struct$graph, from = states[i], 
                              to = states[i+1], mode = 'out', 
                              output = 'epath')$epath[[1]])
  })
  
  # imputed sequence of nodes visited
  vpath = c(states[1], ctds_struct$edge_df$to[do.call(c, imputed)])
  
  # imputed transition times
  times.imputed = do.call(c, sapply(1:(n-1), function(i) {
    ntx = length(imputed[[i]])
    if(ntx > 0) {
      seq(from = times[i], to = times[i+1], length.out = ntx)
    } else if(i==1) {
      times[1]
    } else {
      NULL
    }
  }))
  
  # package results
  imputed = list(
    states = vpath,
    times = times.imputed,
    durations = diff(times.imputed)
  )
  
  class(imputed) = 'ctds_realization'
  
  imputed
}
