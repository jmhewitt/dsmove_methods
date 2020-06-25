source('sketches/2020_06_24_initial_fns/domain.R')

library(nimble)

#
# determine support of local transition matrix
#

# central location coordinates
box.center = c(s1 = 50, s2 = 50)

# central location index
box.ind = which(
  coords[,'s1'] == box.center['s1'] & 
    coords[,'s2'] == box.center['s2']
)

# location indices to include in local matrix exponential
nbs.inds = nbs.local[[box.ind]]

# Some naming conventions for sparse vector structures: INDEXTYPE_ORDERING


# extract directed CTDS edges that begin and end at locations in locs
local_edges = filter_edges(locs = nbs.inds, 
                           inedges_by_loc = 
                             do.call(c, ctds_struct$in_edges_inds), 
                           loc_start = c(1, 1 + cumsum(ctds_struct$in_degree)), 
                           fromlocs_by_edge = ctds_struct$edge_df$from)

# returns the first index ind such that vec[ind] == x, or 0 if x is not in vec
intWhich = nimbleFunction(
  run = function(x = integer(0), vec = integer(1)) {
    
    returnType(integer(0))
    
    res <- 0
    
    N <- length(vec)
    
    i <- 1
    searching <- TRUE
    while(searching) {
      
      if(vec[i]==x) {
        res <- i
        searching <- FALSE
      }
      
      i <- i + 1
      if(i > N) { searching <- FALSE }
    }
    
    return(res)
  }
)


# workflow to figure out which edges edge_start can transition to:
#   1) look up to_loc associated with edge_start
#   2) return outedges associated with to_loc
# 
# ultimately also need to index the outedges into appropriate columns in matrix
# it would be interesting if we could actually just figure things out by working 
# with the edge_df matrix.  the cleanest thing might simply be to write a 
# "which" indexing function in nimble.  it will be a crude search, but it will 
# be simple to write and require few assumptions about the inputs.

local_generator = nimbleFunction(
  run = function(locs = integer(1), edges = integer(1), 
                 tolocs_by_edge = integer(1), fromlocs_by_edge = integer(1),
                 outedges_by_loc = integer(1), loc_start = integer(1), 
                 Xloc = double(2), betaLoc = double(1), Xdir = double(2), 
                 betaDir = double(1), W = double(2), betaAR = double(0)) {
    # Parameters:
    #  locs - vector of locations associated with edges
    #  edges - vector of edges to build generator matrix for
    #  tolocs_by_edge - location at which each edge ends
    #  fromlocs_by_edge - location from which each edge begins
    #  outedges_by_loc - edges sorted/grouped by their starting location
    #  loc_start - index of first out-edge for a location in outedges_by_loc
    #  Xloc - covariates for location-based drivers of movement
    #  Xdir - (edge-based) covariates for directional drivers of movement
    #  W - unit-vector matrix; rows define each edge's direction of movement
    #  betaLoc - coefficients for location-based drivers of movement
    #  betaDir - coefficients for directional drivers of movement
    #  betaAR - scalar coefficient for directional persistence in movement
    
    # return an infinitesimal generator matrix
    returnType(double(2))
    
    nlocs <- length(locs)
    nedges <- length(edges)
    
    # zero-initialize local infinitesimal generator matrix
    A <- matrix(0, nrow = nedges, ncol = nedges)
    
    # pre-compute common, location-based drivers of movement
    lambda_loc <- exp((Xloc[locs,] %*% betaLoc)[,1])
    
    # loop over matrix rows
    for(i in 1:nedges) {
      edge_i <- edges[i]
      
      # get location at which edge_i ends
      to_loc <- tolocs_by_edge[edge_i]
      
      # set infinitesimal rate at which edge_i is being left
      A[i,i] <- - lambda_loc[intWhich(x = to_loc, vec = locs)]
      
      # initialize off-diagonal mass counter
      offdiagMass <- 0
      
      # extract direction of movement associated with edge_i (for AR component)
      w_i <- W[edge_i,]
        
      # loop over all out-edges from to_loc
      for(j in loc_start[to_loc]:(loc_start[to_loc+1]-1)) {
        edge_j <- outedges_by_loc[j]
        
        # restrict to off-diagonal entries for A
        if(edge_i != edge_j) {
          # process out-edge edge_j if it has a column in our matrix
          if(intContains(x = edge_j, vec = edges)) {
            dst_loc <- tolocs_by_edge[edge_j]
            
            # unstandardized infinitesimal rate for transition to_loc -> dst_loc
            lambda <- exp(
              # directional drivers of movement
              (Xdir[edge_j,] %*% betaDir)[1,1] + 
              # directional persistence
              (t(w_i) %*% W[edge_j,])[1,1] * betaAR
            )
            
            # set unstandardized rate at which edge_i is transitioning to edge_j
            A[i,intWhich(x = edge_j, vec = edges)] <- lambda
             
            # aggregate unstandardized off-diagonal mass
            offdiagMass <- offdiagMass + lambda
          }
        }
      }
      
      # standardize off-diagonal entries so that the rowsum is 0
      for(j in 1:nedges) {
        if(i!=j) {
          if(A[i,j] != 0) {
            A[i,j] <- - A[i,j] / offdiagMass * A[i,i]
          }
        }
      }
      
    }
    
    return(A)
  }
)

clocal_generator = compileNimble(local_generator)

A = clocal_generator(edges = local_edges, 
                    tolocs_by_edge = ctds_struct$edge_df$to,
                fromlocs_by_edge = ctds_struct$edge_df$from,
                outedges_by_loc = do.call(c, ctds_struct$out_edges_inds),
                loc_start = c(1, 1 + cumsum(ctds_struct$out_degree)), 
                locs = nbs.inds, Xloc = ctds_struct$Xloc, 
                betaLoc = matrix(1, nrow = 1, ncol = 1), 
                Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df), ncol = 1), 
                betaDir = 0, W = ctds_struct$w_ij, betaAR = .5)


#
# compute likelihood!  ...it's a little scary because with directional 
# persistence we will need to compute O(n) matrix exponentials per likelihood 
# evaluation, where n is the number of observations.  but this might not be 
# so bad since it is perhaps similar to the amount of work required to evaluate 
# likelihoods in the basic deep dive model
#
# we will also need to do some extra work to make sure that we transfer the 
# overlapping probability mass from one transition matrix to the next

# impute path from observation
imputed = ctds.shortest_impute(states = ctds_obs$states, times = ctds_obs$times, 
                               ctds_struct = ctds_struct)
