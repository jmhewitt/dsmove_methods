# developing a path proposal distribution that is roughly uniform over total 
# path length

segments = readRDS(readd('impute_segments_sim_obs_0.5_sim_trajectory_sim_params_1'))$segments

ctds_domain = readd('sim_domain')

# extract sampling weights and indexes for all segments
pathwts = lapply(1:length(segments), function(segind) {
  sfam = segments[[segind]]
  do.call(rbind, lapply(1:length(sfam), function(lengthind) {
    s = sfam[[lengthind]]
    if(is.null(s$paths)) {
      npaths = 1
      pathlen = 0
    } else {
      npaths = nrow(s$paths)
      pathlen = ncol(s$paths)
    }
    data.frame(w = s$weights, pathlen = pathlen, pathind = 1:npaths, 
               lengthind = lengthind)
  }))
})

pathwts.aggregated = lapply(pathwts, function(p) {
  p %>% 
    dplyr::group_by(pathlen) %>% 
    dplyr::summarise(w = sum(w)) %>% 
    dplyr::ungroup()
})

pathlen.support = c(
  # shortest possible path
  min = sum(sapply(pathwts, function(d) d$pathlen[1] )),
  # longest possible path
  max = sum(sapply(pathwts, function(d) tail(d$pathlen, 1) ))
)


# draw random paths
paths = replicate(n = 10, expr = {
  
  # common variate used to draw path lengths across all segments
  u = runif(n = 1)
  
  # randomly draw segments
  path.segments = lapply(1:length(pathwts), function(segind) {
    # extract length of sampled segment
    rowind = min(which(cumsum(pathwts.aggregated[[segind]]$w) >= u))
    len = pathwts.aggregated[[segind]]$pathlen[rowind]
    # sample segment conditional on segment length
    d = pathwts[[segind]] %>% dplyr::filter(pathlen == len)
    rowind = sample(x = 1:nrow(d), size = 1, prob = d$w)
    
    unlist(segments[[segind]][[d$lengthind[rowind]]]$paths[d$pathind[rowind],])
  })
  
  
})



pathlens.alt = replicate(n = 1e4, expr = {
  
  u = runif(n = 1)
  
  sum(sapply(pathwts.aggregated, function(d) {
    d$pathlen[min(which(cumsum(d$w) >= u))]
  }))
  
})

hist(pathlens.alt)
summary(pathlens.alt)
pathlen.support
