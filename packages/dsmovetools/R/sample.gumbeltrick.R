#' Sample a categorical variable when given log-probabilities
#' 
#' @param log.p log of categorical probabilitites; the probabilities do not 
#'   need to be normalized
#' 
#' @references https://timvieira.github.io/blog/post/2014/07/31/gumbel-max-trick/  
#' 
#' @export
#' 
sample.gumbeltrick = function(log.p) {
  # sample gumbel variates
  g = -log(-log(runif(n = length(log.p))))
  # return index of gumbel-max
  which.max(log.p + g)
}
