#' Build structures for manipulating CTDS models
#'
#' @param g Either an \code{igraph} or an adjacency matrix that defines the 
#'   neighborhood structure of a spatial domain with \eqn{n} locations.
#' @param Xloc \eqn{n x p} Matrix of location-based drivers of movement.  
#'   Each row defines covariates for a single location.  If \code{NULL}, then an 
#'   intercept-only matrix will be automatically generated instead.
#' @param coords \eqn{n x d} Matrix associating coordinates with each location.  
#'   Used for determining covariates for directional drivers of motion.  
#'   If \code{NULL}, then the returned \code{ctds} object will not be able to 
#'   fit models with directional drivers of motion.
#' @param Vdir Either one \eqn{n x d} matrix, or a list of such matrices, 
#'   each of which defines one directional driver of movement.  For example, 
#'   the \code{j}th row in the \code{i}th \code{Vdir} matrix specifies a vector.  
#'   The vector defines the strength and magnitude of the \code{i}th directional 
#'   driver of movement at location \code{j}.
#'
#' @return A \code{ctds_struct} object containing information needed to compute 
#'   on the CTDS model.  One of the most important components of
#'   \code{ctds_struct} is the data.frame \code{edge_df}.  The data.frame is 
#'   structured so that all of the out-edges from a node appear as a consecutive 
#'   block of rows.
#' 
#' @example examples/build_ctds.R
#' 
#' @import igraph
#' 
#' @export
#' 
build_ctds = function(g, Xloc = NULL, coords = NULL, Vdir = NULL) {
  
  # directed igraph representation of spatial domain
  if(!inherits(g, 'igraph')) {
    g = graph_from_adjacency_matrix(adjmatrix = g, mode = 'directed')
  } else {
    if(!is_directed(g)) {
      g = as.directed(g, mode = 'mutual')
    }
  }
  
  # number of nodes and edges in graph
  n = vcount(g)
  nedges = ecount(g)
  
  # edge-list enumeration of edges
  edge_df = as_data_frame(g, 'edges')
  
  # in/out-degree of all nodes
  out_degree = degree(g, mode = 'out')
  in_degree = degree(g, mode = 'in')

  # vector with index of each node's first out-edge; 
  # has n+1 entries, and last entry counts number of edges in graph + 1
  out_edges_start = c(0, cumsum(out_degree)) + 1
  
  # indices for outbound edges associated with each node
  out_edges_inds = lapply(1:n, function(v) {
    seq(from = out_edges_start[v], to = out_edges_start[v+1] - 1)
  })
  
  # indices for edges that lead in to each node
  in_edges_inds = lapply(1:n, function(v) { 
    which(edge_df$to==v)
  })
  
  # intercept-only location-based drivers of motion
  if(is.null(Xloc)) {
    Xloc = matrix(1, nrow = vcount(g), ncol = 1)
  }
  
  # compute unit-vectors that point in the direction of neighboring nodes
  if(is.null(coords)) {
    w_ij = NULL 
  } else { 
    # (euclidean) displacements from source nodes to destination nodes
    dv = coords[edge_df$to,] - coords[edge_df$from,]
    # (euclidean) directions from source nodes to destination nodes
    w_ij = as.matrix(sweep(dv, 1, sqrt(rowSums(dv^2)), '/'))
    # trim space because we don't need dimnames
    dimnames(w_ij) = NULL
  }
  
  # compute directional drivers of movement
  if(is.null(Vdir)) {
    q_lij = NULL
  } else {
    # basic input validation and munging
    if(is.null(coords)) {
      stop('Must provide coords to compute directional drivers of movement.')
    }
    if(!inherits(Vdir, 'list')) {
      Vdir = list(Vdir)
    }
    # larger drivers of movement where alignment is stronger between the 
    # covariate Vl and direction to neighbor w_ij
    q_lij = do.call(cbind, lapply(Vdir, function(Vl) {
      rowSums(Vl[edge_df$from,] * w_ij)
    }))
  }
  
  
  #
  # package results
  #
  
  res = list(
    graph = g,
    coords = coords,
    nedges = nedges,
    edge_df = edge_df,
    in_degree = in_degree,
    out_degree = out_degree,
    out_edges_inds = out_edges_inds,
    in_edges_inds = in_edges_inds,
    Xloc = Xloc,
    w_ij = w_ij,
    q_lij = q_lij
  )
  
  class(res) = 'ctds_struct'
  
  res
}