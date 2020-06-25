#' Nimble version of a basic contains function
#' 
#' Function returns TRUE if the integer vector \code{vec} contains the element 
#' \code{x}
#' 
#' @name intContains
#'
#' @param x value to look for in \code{vec}
#' @param vec integer vector to search over
#' 
#' @import nimble
#' 

NULL

#' 
#' @export
#' @rdname intContains
#' 
# 
intContains = nimbleFunction(
  run = function(x = integer(0), vec = integer(1)) {
    
    returnType(logical(0))
    
    contains <- FALSE
    
    N <- length(vec)
    
    i <- 1
    searching <- TRUE
    while(searching) {
      
      if(vec[i] == x) {
        contains <- TRUE
        searching <- FALSE
      }
      
      i <- i + 1
      if(i > N) { searching <- FALSE }
    }
    
    return(contains)
  }
)
