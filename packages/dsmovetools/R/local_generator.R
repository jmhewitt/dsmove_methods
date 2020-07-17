#' Construct a locally-constrained infinitesimal generator matrix
#' 
#' @name local_generator
#'
#' @param locs vector of locations associated with \code{edges}
#' @param row_edges vector of edges to build generator matrix for
#' @param col_edges vector of edges to build generator matrix for
#' @param tolocs_by_edge location at which each edge in domain ends
#' @param fromlocs_by_edge location from which each edge in domain begins
#' @param outedges_by_loc edges in domain sorted by their starting location
#' @param loc_start index of first out-edge for a location in 
#'   \code{outedges_by_loc}
#' @param Xloc covariates for location-based drivers of movement across domain
#' @param Xdir (edge-based) covariates for directional drivers of movement 
#'   across domain
#' @param W unit-vector matrix for whole domain; rows define each edge's 
#'   direction of movement
#' @param betaLoc coefficients for location-based drivers of movement
#' @param betaDir coefficients for directional drivers of movement
#' @param betaAR scalar coefficient for directional persistence in movement
#' 
#' 
#' @import nimble
#' 

NULL

#' 
#' @export
#' @rdname local_generator
#' 
local_generator = nimbleFunction(
  run = function(locs = integer(1), row_edges = integer(1), 
                 col_edges = integer(1), tolocs_by_edge = integer(1), 
                 fromlocs_by_edge = integer(1), outedges_by_loc = integer(1), 
                 loc_start = integer(1), Xloc = double(2), betaLoc = double(1), 
                 Xdir = double(2), betaDir = double(1), W = double(2), 
                 betaAR = double(0)) {
    
    # return an infinitesimal generator matrix
    returnType(double(2))
    
    nlocs <- length(locs)
    nrows <- length(row_edges)
    ncols <- length(col_edges)
    
    # zero-initialize local infinitesimal generator matrix
    A <- matrix(0, nrow = nrows, ncol = ncols)
    
    # pre-compute common, location-based drivers of movement
    lambda_loc <- exp((Xloc[locs,] %*% betaLoc)[,1])
    
    # loop over matrix rows
    for(i in 1:nrows) {
      edge_i <- row_edges[i]
      
      # get location at which edge_i ends
      to_loc <- tolocs_by_edge[edge_i]
      
      # row_edges and col_edges may have different ordering; determine where 
      # edge_i exists in col_edges, so we can store the infinitesimal 
      # transition rate in the correct location
      identity_column <- intWhich(x = edge_i, vec = col_edges)
      
      # set infinitesimal rate at which edge_i is being left
      A[i, identity_column] <- - lambda_loc[intWhich(x = to_loc, vec = locs)]
      
      # initialize off-diagonal mass counter
      offdiagMass <- 0
      
      # extract direction of movement associated with edge_i (for AR component)
      w_i <- W[edge_i,]
      
      # loop over all out-edges from to_loc
      for(j in loc_start[to_loc]:(loc_start[to_loc+1]-1)) {
        edge_j <- outedges_by_loc[j]
        
        # restrict to "off-diagonal" entries for A
        if(edge_i != edge_j) {
          # process out-edge edge_j if it has a column in our matrix
          if(intContains(x = edge_j, vec = col_edges)) {
            dst_loc <- tolocs_by_edge[edge_j]
            
            # unstandardized infinitesimal rate for transition to_loc -> dst_loc
            lambda <- exp(
              # directional drivers of movement
              (Xdir[edge_j,] %*% betaDir)[1,1] + 
                # directional persistence
                (t(w_i) %*% W[edge_j,])[1,1] * betaAR
            )
            
            # set unstandardized rate at which edge_i is transitioning to edge_j
            A[i,intWhich(x = edge_j, vec = col_edges)] <- lambda
            
            # aggregate unstandardized off-diagonal mass
            offdiagMass <- offdiagMass + lambda
          }
        }
      }
      
      # standardize off-diagonal entries so that the rowsum is 0
      for(j in 1:ncols) {
        if(j!=identity_column) {
          if(A[i,j] != 0) {
            A[i,j] <- - A[i,j] / offdiagMass * A[i, identity_column]
          }
        }
      }
      
    }
    
    return(A)
  }
)