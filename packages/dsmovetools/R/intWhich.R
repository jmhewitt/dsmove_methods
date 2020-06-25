#' Find the index at which \code{x} first appears in \code{vec}
#' 
#' @name intWhich
#'
#' @param x value to look for in \code{vec}
#' @param vec integer vector to search over
#' 
#' @return the first index \code{ind} such that \code{vec[ind] == x}, or 0 if 
#'   \code{x} is not in \code{vec}
#' 
#' @import nimble
#' 

NULL

#' 
#' @export
#' @rdname intWhich
#' 
intWhich = nimbleFunction(
  run = function(x = integer(0), vec = integer(1)) {
    
    returnType(integer(0))
    
    res <- 0
    
    N <- length(vec)
    
    i <- 1
    searching <- TRUE
    while(searching) {
      
      if(vec[i]==x) {
        res <- i
        searching <- FALSE
      }
      
      i <- i + 1
      if(i > N) { searching <- FALSE }
    }
    
    return(res)
  }
)