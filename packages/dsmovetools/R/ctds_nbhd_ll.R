#' Extract directed CTDS edges that begin and end at locations in locs
#' 
#' @name ctds_nbhd_ll
#'
#' @param x vector of states visited
#' 
#' 
#' @import nimble
#' 

NULL

#' 
#' @export
#' @rdname ctds_nbhd_ll
#' 
ctds_nbhd_ll = nimbleFunction(
  run = function(x = integer(1), durations = double(1), N = integer(0),
                 inedges_by_loc = integer(1), inloc_start = integer(1),
                 tolocs_by_edge = integer(1), fromlocs_by_edge = integer(1),
                 outedges_by_loc = integer(1), loc_start = integer(1), 
                 Xloc = double(2), betaLoc = double(1), Xdir = double(2), 
                 betaDir = double(1), W = double(2), betaAR = double(0),
                 nbrlocs_by_loc = integer(1), nbrlocs_start = integer(1),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    
    ll <- 0
    
    # use observation likelihood to initialize probabilities for latent state
    edges <- inedges_by_loc[inloc_start[x[1]]:(inloc_start[x[1]+1] - 1)]
    n_edges <- length(edges)
    p <- rep(1/n_edges, n_edges)
    
    # loop over observations
    for(i in 1:N) {

      # observation likelihood
      edges_new <- inedges_by_loc[inloc_start[x[i]]:(inloc_start[x[i]+1] - 1)]
      n_edges_new <- length(edges_new)
      p_new <- rep(1/n_edges_new, n_edges_new)

      # local infinitesimal generator for transition being processed
      locs <- nbrlocs_by_loc[nbrlocs_start[x[i]]:(nbrlocs_start[x[i]+1] - 1)]
      local_edges <- filter_edges(locs = locs, inedges_by_loc = inedges_by_loc,
                                  loc_start = loc_start,
                                  fromlocs_by_edge = fromlocs_by_edge)
      A <- local_generator(
        locs = locs, row_edges = local_edges, col_edges = local_edges,
        tolocs_by_edge = tolocs_by_edge, fromlocs_by_edge = fromlocs_by_edge,
        outedges_by_loc = outedges_by_loc, loc_start = loc_start, Xloc = Xloc,
        betaLoc = betaLoc, Xdir = Xdir, betaDir = betaDir, W = W,
        betaAR = betaAR
      )

      # map current state probabilities to local edge list
      n_local_edges <- length(local_edges)
      p_local <- rep(0, n_local_edges)
      for(j in 1:n_edges) {
        ind <- intWhich(x = edges[j], vec = local_edges)
        if(ind > 0) {
          p_local[ind] <- p[j]
        }
      }

      # update current state probabilities with observation likelihood
      for(j in 1:n_edges_new) {
        ind <- intWhich(x = edges_new[j], vec = local_edges)
        if(ind >0) {
          p_local[ind] <- p_local[ind] * p_new[j]
        }
      }

      # accumulate likelihood
      mass <- sum(p_local)
      ll <- ll + log(mass)
      p_local <- p_local / mass

      # diffuse and update state probabilities
      if(i < N) {

        # flatten infinitesimal generator
        ne <- n_local_edges^2
        Aentries <- numeric(length = ne, init = FALSE)
        ind <- 1
        for(k in 1:n_local_edges) {
          for(j in 1:n_local_edges) {
            Aentries[ind] <- A[j,k]
            ind <- ind + 1
          }
        }

        # matrix exponential
        expA <- expocall_gpadm(H = Aentries[1:ne], t = durations[i],
                               nrows = n_local_edges, ncols = n_local_edges)

        # update probability state storage
        n_edges <- n_local_edges
        edges <- local_edges
        p <- numeric(length = n_local_edges, init = FALSE)
        p[1:n_local_edges] <- p_local[1:n_local_edges] %*%
          expA[1:n_local_edges, 1:n_local_edges]
      }

    }
    
    if(log) { return(ll) } 
    else { return(exp(ll)) }
})