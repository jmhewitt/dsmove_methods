#' Extract directed CTDS edges that begin and end at locations in locs
#' 
#' @name filter_edges
#'
#' @param locs vector of locations for which edges should be filtered
#' @param inedges_by_loc edges sorted/grouped by their endpoint location
#' @param loc_start index of first in-edge for a location in 
#'   \code{inedges_by_loc}
#' @param fromlocs_by_edge location from which each edge begins
#' 
#' @return vector of edges that begin and end at locations in \code{locs}
#' 
#' @import nimble
#' 

NULL

#' 
#' @export
#' @rdname filter_edges
#' 
filter_edges = nimbleFunction(
  run = function(locs = integer(1), inedges_by_loc = integer(1), 
                 loc_start = integer(1), fromlocs_by_edge = integer(1)) {
    
    # return a vector of edges
    returnType(integer(1))
    
    nlocs <- length(locs)
    
    # maximum number of edges that may be returned
    nedges_max <- 0
    for(i in 1:nlocs) {
      loc <- locs[i]
      # add number of in-edges associated with "loc"
      nedges_max <- nedges_max + loc_start[loc+1] - loc_start[loc]
    }
    
    # initialize output and count of filtered edges
    edges <- integer(length = nedges_max, init = FALSE)
    nedges <- 0
    
    # loop over locations
    for(i in 1:nlocs) {
      loc <- locs[i]
      # loop over edges that end at loc
      for(j in loc_start[loc]:(loc_start[loc+1]-1)) {
        in_edge <- inedges_by_loc[j]
        # check if in_edge begins at a location in locs
        if(intContains(x = fromlocs_by_edge[in_edge], vec = locs)) {
          # add to collection of filtered edges
          nedges <- nedges + 1
          edges[nedges] <- in_edge
        }
      }
    }
    
    if(nedges > 0) {
      return(edges[1:nedges])
    } else {
      stop("No filtered edges to return.")
    }
  }
)