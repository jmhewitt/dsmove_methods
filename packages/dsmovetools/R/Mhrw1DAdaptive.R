#' Adaptive Random Walk Metropolis-Hastings (A-RWMH) sampler class
#' 
#' A-RWMH sampling for univariate densities, following 
#' Andrieu and Thoms (2008).  A-RWMH samplers need to maintain information 
#' about the target distribution being explored between sampler calls.  When 
#' A-RWMH samplers are used as a Gibbs step, required state is a nuisance to 
#' track and increases the likelihood for implementation errors.  This class is 
#' designed to make such samplers accessible by automatically/internally 
#' managing the state and required adaptation.
#' 
#' @references Andrieu, Christophe, and Johannes Thoms. "A tutorial on adaptive 
#'   MCMC." Statistics and computing 18.4 (2008): 343-373.
#' 
#' @export
#' 
#' @importFrom R6 R6Class
#' 
Mhrw1DAdaptive = R6Class('Mhrw1DAdaptive', 
                       
  private = list(
    # current location of sampler
    .x = NULL,
    # function to compute log of target density
    .lp = NULL,
    # adaptation and proposal components
    .n = 0,              # number of samples drawn
    .sd = NULL,          # proposal s.d.
    .C = NULL,
    .alpha = 0,          # current acceptance rate
    .alpha_star = NULL,  # target acceptance rate
    .adaptive = NULL     # TRUE to make the sampler adapt the proposal s.d.
  ),
  
  public = list(
    
    #' @description
    #' Output the current values of the sampler's state, and adaptation
    #'
    print = function() {
      c(
        sd = private$.sd,
        alpha = private$.alpha,
        alpha_star = private$.alpha_star
      )
    },
    
    #' @description
    #' Initialize a new Adaptive Random Walk Metropolis-Hastings sampler
    #' 
    #' @param x Initial location for sampler
    #' @param sd Initial proposal standard deviation
    #' @param lp Function to evaluate log of target density, given \code{x}
    #' @param C Base step-size scale.  Note: Setting C=1 can lead to adaptation 
    #'   that occurs to quickly, yielding degenerate proposal covariance 
    #'   matrices.
    #' @param alpha Target univariate acceptance rate for adaptation
    #' @param adaptive \code{TRUE} to adapt proposal covariance \code{Sigma}, 
    #'   \code{FALSE} otherwise.
    #' 
    initialize = function(x, sd, lp, C = 1, alpha = .44, adaptive = TRUE) {

      # copy initialization arguments
      private$.x = x
      private$.lp = lp
      private$.sd = sd
      private$.C = C
      private$.alpha_star = alpha
      private$.adaptive = adaptive
      
    },

    
    #' @description
    #' Use RWMH to propose, then accept/reject, a new value from the target 
    #' distribution.  If the sampler was configured with \code{adaptive=TRUE}, 
    #' then the sampler will automatically adapt as well.
    #' 
    #' @param ... Additional arguments to pass to the function that evaluates
    #'   the log-posterior
    #' 
    sample = function(...) {
      
      # extract current value
      x0 = private$.x

      # random variates for proposal
      z = rnorm(n = 1, sd = private$.sd)

      # assemble proposal
      x = x0 + z

      # metropolis ratio
      #   Recall: proposal is symmetric, so proposal density is not needed
      lp.x0 = private$.lp(x0, ...)
      lp.x = private$.lp(x, ...)
      lR = lp.x  - lp.x0

      # accept/reject
      accept = log(runif(1)) <= lR
      if(accept) {
        private$.x = x
      }
      
      # update acceptance rate
      private$.n = private$.n + 1
      private$.alpha = private$.alpha + (accept - private$.alpha) / private$.n
      
      #
      # adaptation
      #

      if(private$.adaptive) {
      
        adaptScale = private$.C / sqrt(private$.n)
        adaptShift = private$.alpha - private$.alpha_star
        private$.sd = private$.sd * exp( adaptScale * adaptShift )

      }
      
      # package results
      return(list(x = private$.x, accepted = accept, 
                  lp = ifelse(accept, lp.x, lp.x0)))
    }
    
  )
)