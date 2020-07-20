#' Simulate observing a complete CTDS trajectory at specified timepoints
#' 
#' Given the components of a completely observed CTDS trajectory, extract the 
#' location of a trajectory at the target times \code{t.obs}.
#' 
#' @param states Complete record of states visited
#' @param times Times at which each of \code{states} was visited
#' @param t.obs Times at which the trajectory should be observed
#' 
#' 
#' @export
#' 
ctds.observe = function(states, times, t.obs) {
  tind = findInterval(t.obs, times)
  
  res = list(
    states = states[tind],
    times = t.obs
  )
  
  class(res) = 'ctds_observations'
  
  res
}