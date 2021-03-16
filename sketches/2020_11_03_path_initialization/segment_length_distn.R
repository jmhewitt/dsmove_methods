segments = readRDS(readd('impute_segments_sim_obs_0.5_sim_trajectory_sim_params_1'))$segments

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
    data.frame(w = s$weights, pathlen = pathlen, segment = segind)
  }))
})

# Monte Carlo sample path lengths
pathlens = replicate(n = 1e4, expr = {
  
  sum(sapply(pathwts, function(d) {
    sample(x = d$pathlen, size = 1, prob = d$w)
  }))
  
})

pathlen.support = c(
  # shortest possible path
  min = sum(sapply(pathwts, function(d) d$pathlen[1] )),
  # longest possible path
  max = sum(sapply(pathwts, function(d) tail(d$pathlen, 1) ))
)

# of course... although path lengths are uniformly distributed over each 
# segment, the total path length is the sum of the segment lengths, and we see 
# a CLT-like effect since we are looking at the sum of many independent samples
hist(pathlens)

qqnorm(pathlens, pch = '.')
qqline(pathlens)

# as expected, however, we see that the mean path length is about equal to the 
# mean of the segments.  the issue, of course, is that we don't like the 
# current distribution of sampled path lengths
round(c(mean(pathlens), mean(pathlen.support)))


# TODO: See if we can get a more uniform sampling of the path lengths
#  Try... using a single x ~ U(0,1) to seed inverse CDF sampling at each 
#  segment, or something like this.  The idea is that a single x ~ U should 
#  a) yield a more predictable path length, that is furthermore b) s.t. the 
#  path segments are all uniformly "longer" or "shorter" possible connecting 
#  segments.  Such a strategy will not allow us to uniformly sample from paths 
#  of a fixed length (which can be comprised of some shortest-possible and 
#  longest-possible segment lengths), but should still afford a reasonable 
#  approach to exploring the space of possible path lengths.  we might even 
#  prefer this method since it *should* yield a path with relatively 
#  paths 

pathwts.aggregated = lapply(pathwts, function(p) {
  p %>% 
    dplyr::group_by(pathlen) %>% 
    dplyr::summarise(w = sum(w)) %>% 
    dplyr::ungroup()
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
