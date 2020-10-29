#' Forward-filter a CTDS state distribution without observations or params.
#' 
#' Forward-filtering follows via a random walk along the CTDS graph structure
#' 
#' @param ctds_struct Output of \code{build_ctds}; a representation of a CTDS 
#'   domain formatted for computing
#' @param a0 Initial state distribution vector (sparse)
#' @param steps Number of transitions to forward filter for
#' 
#' @return A list \code{steps + 1} entries.  The initial state distribution 
#'   \code{a0} is contained in the first entry, and each successive entry has 
#'   the state distribution after one additional diffusion according to a 
#'   random walk along the spatial domain's graph, specified in 
#'   \code{ctds_struct}.
#' 
#' @export
#' 
ctds.ffrw = function(ctds_struct, a0, steps) {
  
  # initialize output
  a = vector(mode = 'list', length = steps + 1)
  a[[1]] = as.numeric(a0)
  
  # forward-diffuse initial mass
  if(length(a) > 1) {
    for(i in 2:length(a)) {
      # initialize diffused mass vector
      b = numeric(ctds_struct$nedges)
      # forward-diffuse mass from non-zero elements of last probability vector
      nzinds = which(a[[i-1]] > 0)
      for(j in nzinds) {
        # extract edge and probability associated with sparse vector index
        edge = j
        prob = a[[i-1]][j]
        # spatial location associated with edge
        loc = ctds_struct$edge_df$to[edge]
        # edges associated with neighboring spatial location
        out_edges = ctds_struct$out_edges_inds[[loc]]
        nout = length(out_edges)
        # diffuse and aggregate mass
        # NOTE: sparse vector math is not used b/c sparse vector addition is a 
        #  computational bottleneck from the Matrix package
        b[out_edges] = b[out_edges] + prob * rep(1/nout, nout)
      }
      a[[i]] = b
    }
  }
  
  a
}
