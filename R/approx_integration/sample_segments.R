sample_segments = function(states, times, ctds_domain, max_lambda, p_max,
                           nsamples) {
  # Sample path segments, where nsamples paths are returned and the path lengths 
  # are uniformly distributed.  Note that this is not the same as randomly 
  # sampling from among all paths that connect observations.
  #
  # Parameters:
  #  max_lambda - prior information: upper bound on maximum transition rate
  #  p_max - upper bound on Poisson transition tail
  #  nsamples - number of segments to simulate between each timepoint
  #    (duplicates will be upweighted)
  
  # extract CTDS domain's graph structure
  g = ctds_domain$graph
  V = V(g)
  E = E(g)
  
  # time between observations
  dt = diff(times)
  
  # determine range of path lengths between observations
  dist_ranges = do.call(rbind, lapply(1:(length(states)-1), function(i) { c(
    # lengths of shortest paths (i.e., num. tx's.) between observed states
    distances(graph = g, v = V[states[i]], to = V[states[i+1]], mode = 'out'),
    # prior bound on lengths of longest paths between observed states
    qpois(p = p_max, lambda = max_lambda * dt[i])
  )}))
  colnames(dist_ranges) = c('min', 'max')
  
  # ensure that the sampling range for path lengths is valid
  unordered_transitions = dist_ranges[, 'min'] > dist_ranges[, 'max']
  if(any(unordered_transitions)) {
    dist_ranges[unordered_transitions, 'max'] = 
      dist_ranges[unordered_transitions, 'min']
    warning('Unordered transitions in segment samples.')
  }

  # sample path segments between observed states
  segments = lapply(1:(length(states)-1), function(i) {
    # distribute mass uniformly among edges associated with observation
    edges = ctds_domain$in_edges_inds[[states[i]]]
    nedges = length(edges)
    a0 = sparseVector(i = edges, x = rep(1/nedges, nedges), length = length(E))
    # distribute mass uniformly among edges associated with next observation
    edges = ctds_domain$in_edges_inds[[states[i+1]]]
    nedges = length(edges)
    af = sparseVector(i = edges, x = rep(1/nedges, nedges), length = length(E))
    # get forward-filtering RW distributions for path lengths to sample over
    pathlens = seq(from = dist_ranges[i, 'min'], to = dist_ranges[i, 'max'])
    ff = ctds.ffrw(ctds_struct = ctds_domain, a0 = a0, steps = tail(pathlens, 1))
    # remove path lengths that are not feasible for connecting a0 to af
    valid_lengths = which(
      sapply(ff, function(v) { sum(v * as.numeric(af)) > 0 })
    ) - 1
    pathlens = pathlens[pathlens %in% valid_lengths]
    # uniformly distribute segments to sample over path lengths
    npaths = rmultinom(n = 1, size = nsamples, prob = rep(1, length(pathlens)))
    # remove path lengths that will not be sampled
    nzpaths = npaths > 0
    npaths = npaths[nzpaths]
    pathlens = pathlens[nzpaths]
    # sample path segments
    lapply(1:length(pathlens), function(j) {
      if(pathlens[j] == 0) {
        list(paths = NULL, weights = npaths[j] / nsamples)
      } else {
        # backward-sample paths
        paths = data.frame(matrix(
          ctds.bsrw(n = npaths[j], ctds_struct = ctds_domain, af = af,
                    ff = ff[1:(pathlens[j]+1)]), 
          nrow = npaths[j], byrow = TRUE
        ))
        colnames(paths) = paste('e', 1:pathlens[j], sep = '')
        # count and consolidate duplicate paths
        paths.aggregated = aggregate(cbind(paths[0], numdup = 1), paths, length)
        # format return
        list(paths = paths.aggregated[, 1:pathlens[j], drop = FALSE],
             weights = paths.aggregated$numdup / nsamples)
      }
    })
  })
  
  segments
}