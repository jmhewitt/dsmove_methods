library(igraph)

# number of nodes (n) along each dimension (d)
n = 100
d = 2

# build 2D lattice in igraph with directed edges
g = make_lattice(dimvector = rep(n,d), directed = TRUE, mutual = TRUE)

# associate euclidean coordinates with lattice points
coords = expand.grid(x = 1:n, y = 1:n)

# directional driver of motion toward center of lattice
Vdir = t(apply(coords, 1, function(x) { x - rep(n/2, d) }))

# initialize structures for working with ctds models
ctds = build_ctds(g = g, Xloc = NULL, coords = coords, Vdir = Vdir)
