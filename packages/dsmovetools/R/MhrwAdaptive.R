#' Adaptive Random Walk Metropolis-Hastings (A-RWMH) sampler class
#' 
#' A-RWMH sampling for multivariate densities, following Algorithm 6 in 
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
MhrwAdaptive = R6Class('MhrwAdaptive', 
                       
  private = list(
    # sampler dimension
    .n = NULL,
    # current location of sampler
    .x = NULL,
    # estimate of target distribution's mean and covariance matrix
    .mu = NULL,
    .Sigma = NULL,
    .Sigma_chol = NULL,
    # function to compute log of target density
    .lp = NULL,
    # adaptation scale
    .lambda = NULL,
    # adaptation step sizes (p. 356, "Automatic choice of the stepsizes")
    .adaptation_count = 0,
    .adaptation_frequency = NULL,
    .C = NULL,
    .alpha = NULL,
    .alpha_star = NULL,
    .adaptive = NULL,
    .calls = 0
  ),
  
  public = list(
    
    #' @description
    #' Output the current values of the sampler's state, and adaptation
    #'
    print = function() {
      list(
        mu = private$.mu,
        Sigma = private$.Sigma,
        lambda = private$.lambda
      )
    },
    
    #' @description
    #' Initialize a new Adaptive Random Walk Metropolis-Hastings sampler
    #' 
    #' @param x Initial location for sampler
    #' @param mu Initial guess of target distribution's mean
    #' @param Sigma Initial guess of target distribution's covariance matrix
    #' @param lambda Initial scaling factors, for adaptation
    #' @param lp Function to evaluate log of target density, given \code{x}
    #' @param C Base step-size scale.  Note: Setting C=1 can lead to adaptation 
    #'   that occurs to quickly, yielding degenerate proposal covariance 
    #'   matrices.
    #' @param alpha Power of the diminishing adaptation step size's decay.
    #'   alpha must be in the range (1/(1+lambda), 1]
    #' @param alpha_star Target univariate acceptance rate for adaptation
    #' @param adaptive \code{TRUE} to adapt proposal covariance \code{Sigma}, 
    #'   \code{FALSE} otherwise.
    #' @param adaptation_frequency How often to adapt sampler parameters
    #' 
    initialize = function(x, mu, Sigma, lambda, lp, C = .75, alpha = 1, 
                          alpha_star = .44, adaptive = TRUE,
                          adaptation_frequency = 1) {
      
      # copy initialization arguments
      private$.x = x
      private$.mu = mu
      private$.Sigma = Sigma
      private$.lambda = lambda
      private$.lp = lp
      private$.C = C
      private$.alpha = alpha
      private$.alpha_star = alpha_star
      private$.adaptive = adaptive
      private$.adaptation_frequency = adaptation_frequency
      
      # initialize sampler internals
      private$.n = length(x)
      private$.Sigma_chol = chol(private$.Sigma)
      
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
      z = sqrt(private$.lambda) *
        (private$.Sigma_chol %*% rnorm(n = private$.n))

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


      #
      # componentwise adaptation
      #

      if(private$.adaptive) {
        
        private$.calls = private$.calls + 1
        
        if(private$.calls %% private$.adaptation_frequency == 0 ) {
          
          # adaptation step size
          private$.adaptation_count = private$.adaptation_count + 1
          gamma = private$.C / private$.adaptation_count^private$.alpha
          
          for(i in 1:private$.n) {
            
            # componentwise accept/reject rate
            x = x0
            x[i] = x[i] + z[i]
            ar.local = min(1, exp(private$.lp(x, ...) - lp.x0))
            
            # update componentwise scale
            log_lambda = log(private$.lambda[i])
            log_lambda = log_lambda + gamma * (ar.local - private$.alpha_star)
            private$.lambda[i] = exp(log_lambda)
            
          }
          
          # update estimate of target mean
          z = private$.x - private$.mu
          private$.mu = private$.mu + gamma * z
          
          # update estimate of target covariance
          private$.Sigma = private$.Sigma + gamma * (z %*% t(z) - private$.Sigma)
          private$.Sigma_chol = chol(private$.Sigma)
          
        }
        
      }
      
      # package results
      return(list(x = private$.x, accepted = accept, 
                  lp = ifelse(accept, lp.x, lp.x0)))
    }
    
  )
)