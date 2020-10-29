#' Backward-sample a CTDS path without observations or params.
#' 
#' Backward-sampling follows via a random walk along the CTDS graph structure
#' 
#' @param n number of paths to sample
#' @param ctds_struct Output of \code{build_ctds}; a representation of a CTDS 
#'   domain formatted for computing
#' @param af Final state distribution vector (sparse)
#' @param ff Forward-filtering distributions, which induces a path length.
#' 
#' @return A vector of edges, representing a path
#' 
#' @export
#' 
ctds.bsrw = function(n, ctds_struct, af, ff) {
  
  # path length is one less than ff length since ff[[1]] is initial mass dist'n.
  len = length(ff) - 1
  
  if(len < 1) {
    stop("Haven't safety-checked code for use with path lengths of 0.")
  }
  
  # prepare to backward sample: extract diffused mass at af
  b = ff[[len+1]] * as.numeric(af)
  
  # backward-sample n paths
  paths = replicate(n, {
    # initialize path
    p = numeric(length = len)
    nzinds = which(b > 0)
    p[len] = ifelse(length(nzinds) > 1, 
                    sample(x = nzinds, size = 1, prob = b[nzinds]),
                    nzinds)
    # backward-sample
    if(len > 1) {
      for(k in (len-1):1) {
        # last spatial movement associated with previous backwards-sampled edge
        loc = ctds_struct$edge_df$from[p[k+1]]
        # edges consistent with last spatial movement
        in_edges = ctds_struct$in_edges_inds[[loc]]
        nin = length(in_edges)
        # condition ff dist'n. on edges that yield contiguous path.
        # no addt'l. weighting b/c all transitions are RW transitions
        bs = numeric(ctds_struct$nedges)
        bs[in_edges] = ff[[k+1]][in_edges]
        # sample next edge
        nzinds = which(bs > 0)
        p[k] = ifelse(length(nzinds) > 1, 
                      sample(x = nzinds, size = 1, prob = bs[nzinds]),
                      nzinds)
      }
    }
    # return path
    p
  })
  
  paths
}
