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


# workflow to figure out which edges edge_start can transition to:
#   1) look up to_loc associated with edge_start
#   2) return outedges associated with to_loc
# 
# ultimately also need to index the outedges into appropriate columns in matrix
# it would be interesting if we could actually just figure things out by working 
# with the edge_df matrix.  the cleanest thing might simply be to write a 
# "which" indexing function in nimble.  it will be a crude search, but it will 
# be simple to write and require few assumptions about the inputs.

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
