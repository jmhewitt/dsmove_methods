# define spatial domain, and basic ctds quantities

library(devtools)
document('packages/dsmovetools/')


#
# construct spatial lattice
#

library(sp)

# number of coordinate dimensions
dim = 2

# number of cells per dimension
n.dim = 5

# construct lattice
gr = GridTopology(cellcentre.offset = rep(0,dim), cellsize = rep(1,dim), 
                  cells.dim = rep(n.dim,dim))

# extract coordinates
coords = coordinates(gr)


#
# define transition neighborhoods
#

library(spdep)

# single-step neighborhood for CTDS model
nbs.tx = grid2nb(grid = gr, queen = FALSE, nb = TRUE, self = FALSE)

# neighborhood for local matrix exponentials
nbs.local = grid2nb(grid = gr, queen = TRUE, nb = TRUE, self = TRUE)


#
# graph-structure to define directed edges for CTDS model
#

library(igraph)

g = make_graph(
  edges = do.call(c, sapply(1:length(nbs.tx), function(from_loc) { 
    do.call(as.numeric, strsplit(
      # flatten edges into from_loc -> to_loc format for igraph::make_graph
      x = paste(from_loc, nbs.tx[[from_loc]], collapse = ' '), 
      split = ' '
    ))
  }))
)

nbs.extended = neighborhood(graph = g, order = 10, mode = 'out')

#
# precompute quantities used to determine CTDS transition rates
#

ctds_struct = build_ctds(g = g, coords = coords)

