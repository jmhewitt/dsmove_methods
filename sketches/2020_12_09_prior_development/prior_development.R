library(ggplot2)
library(metR)
library(mvtnorm)
library(dplyr)
library(ggthemes)
library(viridis)
library(coda)



#
# load ctds structures
#

source("./packages.R")
invisible(lapply(list.files("./R", full.names = TRUE, recursive = TRUE), source))
source('R/plan.R')

ctds_struct = readd('sim_domain')

# precompute common infinitesimal generator elements
outedges_by_loc = do.call(c, ctds_struct$out_edges_inds)
loc_start = c(1, 1 + cumsum(ctds_struct$out_degree))


# index of central-most coordinate
ind.center = which.min(
  apply(sweep(x = ctds_struct$coords, MARGIN = 2, 
              STATS = colMeans(ctds_struct$coords), FUN = '-'),
        1, function(r) sum(r^2))
)

# arbitrarily select an edge, which implies a previous location
e_cur = ctds_struct$in_edges_inds[[ind.center]][1]
# get outbound edges from the central-most coordinate
edges = ctds_struct$out_edges_inds[[ind.center]]

# identify current location
v_cur = ctds_struct$edge_df$to[e_cur]
# identify locations associated with outbound edges
locs = ctds_struct$edge_df$to[edges]
  
# identify the edge associated with a reversal of movement
reversal_ind = which(
  ctds_struct$edge_df$from[e_cur] == ctds_struct$edge_df$to[edges]
)


#
# standard, joint normal prior
#


xseq = seq(from = -3, to = 3, length.out = 20)
yseq = seq(from = -3, to = 3, length.out = 20)

# marginal means and variances
muvec = c(0, 0)
sdvec = c(1e2, 1e2)
# bivariate correlation
rho = -.9999
# rho = 0

# associated covariance matrix
Sigma = matrix(c(sdvec[1]^2, rho * prod(sdvec), rho * prod(sdvec), sdvec[2]^2),
               nrow = 2)

# evaluate likelihood across grid
df = expand.grid(beta_loc = xseq, beta_ar = yseq)
df$ll = dmvnorm(x = df, mean = muvec, sigma = Sigma, log = TRUE)

contour_col = 'grey50'
contour_text_col = 'grey80'

ggplot(df, aes(x = beta_loc, y = beta_ar, fill = ll, z = ll)) + 
  # log-prior surface
  geom_raster() + 
  # contours and labels
  geom_contour2(col = contour_col) + 
  geom_text_contour(col = contour_text_col) + 
  # "undesired" points
  geom_point(x = 2, y = -2, col = 'white', pch = 4) + 
  # formatting
  scale_fill_viridis(direction = -1) + 
  theme_few() + 
  theme(panel.border = element_blank()) + 
  coord_equal()

# the end of the story is that with non-informative marginals, there needs to 
# be extremely strong prior correlation in order to impart meaningful 
# penalization onto the joint distribution


#
# alternate prior: only penalize expected number of reversals per unit time
#  by placing an exponential distribution on the expected number of reversals...
#  of course, this will induce a prior on the marginal values, so we might need 
#  to consider formal copula methods to adjust
#

xseq = seq(from = -10, to = 10, length.out = 20)
yseq = seq(from = -10, to = 10, length.out = 20)

xseq = seq(from = -3, to = 3, length.out = 20)
yseq = seq(from = -3, to = 3, length.out = 20)

# marginal means and variances
muvec = c(0, 0)
sdvec = c(1e2, 1e2)
# penalty strength
lambda = 0

# basis for joint prior
expected_reversals = function(beta_ar, beta_loc) {
  
  mapply(function(beta_ar, beta_loc) {
    # get local transition parameters via infinitesimal generator extract
    A_cols = c(e_cur, edges)
    A = local_generator(locs = c(v_cur, locs), row_edges = e_cur,
                        col_edges = A_cols,
                        tolocs_by_edge = ctds_struct$edge_df$to,
                        fromlocs_by_edge = ctds_struct$edge_df$from,
                        outedges_by_loc = outedges_by_loc,
                        loc_start = loc_start, Xloc = ctds_struct$Xloc,
                        betaLoc = matrix(beta_loc),
                        Xdir = matrix(0, nrow = nrow(ctds_struct$edge_df),
                                      ncol = 1),
                        betaDir = matrix(0), W = ctds_struct$w_ij,
                        betaAR = beta_ar)
    
    # expected number of transitions per unit time
    -A[1] *
    # scaled by expected proportion or reversals
    A[-1][reversal_ind] / sum(A[-1])
    
  }, beta_ar, beta_loc)
}

# set penalty strength (i.e., 1/expected # flips per time we want to tolerate)
# penalty_rate = 1/.5
penalty_rate = 1/.125

# evaluate likelihood across grid
df = expand.grid(beta_loc = xseq, beta_ar = yseq)
df$ll = dexp(x = expected_reversals(beta_ar = df$beta_ar, beta_loc = df$beta_loc), 
             log = TRUE, rate = penalty_rate)

contour_col = 'grey50'
contour_text_col = 'grey80'

ggplot(df, aes(x = beta_loc, y = beta_ar, fill = ll, z = ll)) + 
  # log-prior surface
  geom_raster() + 
  # contours and labels
  geom_contour2(col = contour_col) + 
  geom_text_contour(col = contour_text_col) + 
  # "undesired" points
  geom_point(x = 2, y = -2, col = 'white', pch = 4) + 
  # axes
  geom_hline(yintercept = 0, col = 'white', lty = 3) + 
  geom_vline(xintercept = 0, col = 'white', lty = 3) + 
  # formatting
  scale_fill_viridis(direction = -1) + 
  theme_few() + 
  theme(panel.border = element_blank()) + 
  coord_equal()



#
# rejection sample the proposed prior distribution to learn about its properties
#

# specify proposal distribution range
proposal_sd = c(3,3)

# compute rejection sampling ratio
M = dexp(x = 0, rate = penalty_rate) / 
  prod(dnorm(x = muvec, mean = muvec, sd = proposal_sd))
  
# draw from prior distribution
rejection.sample = t(replicate(n = 1e4, expr = {
  
  accept = FALSE
  
  while(!accept) {
    # propose coordinate
    x = rnorm(n = 2, mean = muvec, sd = proposal_sd)
    
    # evaluate log-proposal density
    g = sum(dnorm(x = x, mean = muvec, sd = proposal_sd, log = TRUE))
    
    # evaluate log-target density
    f = dexp(x = expected_reversals(beta_ar = x[2], beta_loc = x[1]), 
             log = TRUE, rate = penalty_rate)
    
    # rejection step
    accept = log(runif(1)) < f - g - log(M)
    if(accept) {
      break
    }
  }
  
  x
}))
colnames(rejection.sample) = c('beta_loc', 'beta_ar')

# verify that the rejection sample resembles the density contours
ggplot(df, aes(x = beta_loc, y = beta_ar, fill = exp(ll), z = exp(ll))) + 
  # log-prior surface
  geom_raster() + 
  # contours and labels
  geom_contour2(col = contour_col) + 
  geom_text_contour(col = contour_text_col) + 
  # "undesired" points
  geom_point(x = 2, y = -2, col = 'white', pch = 4) + 
  geom_point(x = 1.5, y = -2.6, col = 'white', pch = 4) + 
  # axes
  geom_hline(yintercept = 0, col = 'white', lty = 3) + 
  geom_vline(xintercept = 0, col = 'white', lty = 3) + 
  # rejection-sampled estmiate
  stat_density_2d(data = data.frame(beta_loc = rejection.sample[,'beta_loc'],
                               beta_ar = rejection.sample[,'beta_ar'],
                               ll = 0)) + 
  # plot limits
  xlim(-3,3) + 
  ylim(-3,3) + 
  # formatting
  scale_fill_viridis(direction = -1) + 
  theme_few() + 
  theme(panel.border = element_blank()) + 
  coord_equal()

# look at marginals and summaries
summary(mcmc(rejection.sample))
HPDinterval(mcmc(rejection.sample))


# we see that the prior induces otherwise informative marginals, but the range 
# is probably fairly reasonable.
ggarrange(
  ggplot(data.frame(beta_loc = rejection.sample[,'beta_loc']), aes(x=beta_loc)) +
    stat_density(geom = 'line') + 
    theme_few() + 
    theme(panel.border = element_blank()),
  ggplot(data.frame(beta_ar = rejection.sample[,'beta_ar']), aes(x=beta_ar)) +
    stat_density(geom = 'line') + 
    theme_few() + 
    theme(panel.border = element_blank()),
  ncol = 2, nrow = 1
)
