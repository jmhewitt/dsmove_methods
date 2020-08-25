library(nimble)

#
# determine support of local transition matrix
#

# central location coordinates
box.center = c(s1 = 5, s2 = 5)

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

A = local_generator(row_edges = local_edges, col_edges = local_edges,
                    tolocs_by_edge = ctds_struct$edge_df$to,
                    fromlocs_by_edge = ctds_struct$edge_df$from,
                    outedges_by_loc = do.call(c, ctds_struct$out_edges_inds),
                    loc_start = c(1, 1 + cumsum(ctds_struct$out_degree)), 
                    locs = nbs.inds, Xloc = ctds_struct$Xloc, 
                    betaLoc = matrix(1, nrow = 1, ncol = 1), 
                    Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df), 
                                  ncol = 1), 
                    betaDir = 0, W = ctds_struct$w_ij, betaAR = .5)

ctds_struct$coords = data.frame(ctds_struct$coords)
ctds_struct$coords$x = ctds_struct$coords$s1
ctds_struct$coords$y = ctds_struct$coords$s2

# ctds_sim = ctds.fwdsim(ctds_struct = ctds_struct, 
#                        beta_loc = matrix(.01, nrow = 1, ncol = 1), 
#                        beta_dir = 0, v0 = box.ind, t0 = 0, tf = 100, 
#                        max.steps = 200, beta_ar = 1, v0.last = NULL)


ctds_sim = ctds.fwdsim(ctds_struct = ctds_struct, 
                       beta_loc = matrix(-1, nrow = 1, ncol = 1), 
                       beta_dir = 0, v0 = box.ind, t0 = 0, tf = 200, 
                       max.steps = 1e3, beta_ar = 1, v0.last = NULL)


ctds_obs = ctds.observe(states = ctds_sim$states, times = ctds_sim$times, 
                        t.obs = seq(from = ctds_sim$times[1], 
                                    to = ctds_sim$times[length(ctds_sim$times)], 
                                    length.out = 100))



#
# compute likelihood!  ...it's a little scary because with directional 
# persistence we will need to compute O(n) matrix exponentials per likelihood 
# evaluation, where n is the number of observations.  but this might not be 
# so bad since it is perhaps similar to the amount of work required to evaluate 
# likelihoods in the basic deep dive model
#
# we will also need to do some extra work to make sure that we transfer the 
# overlapping probability mass from one transition matrix to the next



# ctds_obs = ctds.observe(states = ctds_sim$states, times = ctds_sim$times, 
#                         t.obs = seq(from = ctds_sim$times[1], 
#                                     to = ctds_sim$times[length(ctds_sim$times)], 
#                                     length.out = 50))



# impute path from observation
imputed = ctds.shortest_impute(states = ctds_obs$states, times = ctds_obs$times, 
                               ctds_struct = ctds_struct)

plot.ctds_realization(x = ctds_sim, ctds_struct = ctds_struct, 
                      ctds_obs = ctds_obs)

plot.ctds_realization(x = imputed, ctds_struct = ctds_struct, 
                      ctds_obs = ctds_obs)

# document('packages/dsmovetools/')



# next steps (for likelihood):
#  1) determine obervation likelihood for latent ctds states from imputed
#  2) HMM-like forward/backward algorithm to diffuse probability mass, yielding 
#      likelihood... this will be slightly complicated by the need to transfer 
#      only the overlapping components of the densities for the 3x3 cells.

# observation likelihood is uniformly distributed over these edges
ctds_struct$in_edges_inds[imputed$states[1]]

# location indices to include in local matrix exponential
nbs.inds = nbs.local[[imputed$states[1]]]

local_edges = filter_edges(locs = nbs.inds, 
                           inedges_by_loc = 
                             do.call(c, ctds_struct$in_edges_inds), 
                           loc_start = c(1, 1 + cumsum(ctds_struct$in_degree)), 
                           fromlocs_by_edge = ctds_struct$edge_df$from)

A = local_generator(row_edges = local_edges, col_edges = local_edges,
                    tolocs_by_edge = ctds_struct$edge_df$to,
                    fromlocs_by_edge = ctds_struct$edge_df$from,
                    outedges_by_loc = do.call(c, ctds_struct$out_edges_inds),
                    loc_start = c(1, 1 + cumsum(ctds_struct$out_degree)), 
                    locs = nbs.inds, Xloc = ctds_struct$Xloc, 
                    betaLoc = matrix(1, nrow = 1, ncol = 1), 
                    Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df), 
                                  ncol = 1), 
                    betaDir = 0, W = ctds_struct$w_ij, betaAR = .5)

library(Matrix)
pmat = expm(A * diff(imputed$times[1:2]))


# likelihood steps
#  1) evaluate observation likelihood for current state
#          - uniform mass over all the edges that end at observation location

#  2) build neighborhood-based local generator (in edge-space) for current 
#     state.  by design, this should be a superset of the observation likelihood
#          - we have demo code for this already

#  3) map current state probabilities to local generator.  some mass will be 
#     lost, but this reflects an observation likelihood effect rather than an 
#     issue with the movement approximation.

#  4) map observation likelihood to local generator

#  5) update mapped state probabilities with observation likelihood; update 
#     overall likelihood contributions

#  6) diffuse updated state probabilities via matrix exponential




document('packages/dsmovetools/')


cctds_nbhd_ll = compileNimble(ctds_nbhd_ll)



ll = cctds_nbhd_ll(x = imputed$states, durations = imputed$durations, 
                  N = length(imputed$states), 
                  inedges_by_loc = do.call(c, ctds_struct$in_edges_inds), 
                  inloc_start = c(1, 1 + cumsum(ctds_struct$in_degree)), 
                  tolocs_by_edge = ctds_struct$edge_df$to, 
                  fromlocs_by_edge = ctds_struct$edge_df$from, 
                  outedges_by_loc = do.call(c, ctds_struct$out_edges_inds), 
                  loc_start = c(1, 1 + cumsum(ctds_struct$out_degree)), 
                  Xloc = ctds_struct$Xloc, 
                  betaLoc = matrix(2, nrow = 1, ncol = 1), 
                  Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df), ncol = 1), 
                  betaDir = 0, W = ctds_struct$w_ij, betaAR = 0, 
                  nbrlocs_by_loc = do.call(c, nbs.local), 
                  nbrlocs_start = c(1, 1 + cumsum(sapply(nbs.local, length))), 
                  log = TRUE)
