spatial_lattice_ctds = function(n_coord_dimensions, cells_per_dimension) {
  
  #
  # define spatial domain
  #
  
  # construct lattice
  gr = GridTopology(cellcentre.offset = rep(0, n_coord_dimensions), 
                    cellsize = rep(1, n_coord_dimensions), 
                    cells.dim = cells_per_dimension)
  
  # extract coordinates
  coords = coordinates(gr)
  
  
  #
  # define transition neighborhoods
  #
  
  # single-step neighborhood for CTDS model
  nbs.tx = grid2nb(grid = gr, queen = FALSE, nb = TRUE, self = FALSE)
  
  # neighborhood for local matrix exponentials
  nbs.local = grid2nb(grid = gr, queen = TRUE, nb = TRUE, self = TRUE)
  
  
  #
  # graph-structure to define directed edges for CTDS model
  #
  
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
  ctds_struct$nbs.local = nbs.local
  
  ctds_struct
}
