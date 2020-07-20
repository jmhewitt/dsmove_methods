#' Extract directed CTDS edges that begin and end at locations in locs
#' 
#' @name expocall_gpadm
#'
#' @param x vector of states visited
#' 
#' 
#' @import nimble
#' 

NULL

#' 
#' @export
#' @rdname expocall_gpadm
#' 
expocall_gpadm = nimbleFunction(
  run = function(H = double(1), t = double(0), nrows = integer(0), 
                 ncols = integer(0)) {
    
    t_wrap = numeric(1, init = FALSE)
    t_wrap[1] <- t
    
    returnType(double(2))
    
    ideg <- integer(6, length = 1)
    m <- integer(nrows, length = 1)
    ldh <- integer(ncols, length = 1)
    lwsp <- integer(4*m*m+ideg+1 + 1, length = 1)
    wsp <- numeric(0, length = lwsp[1])
    ipiv <- integer(0, length = m[1])
    iexph <- integer(1)
    ns <- integer(1)
    iflag <- integer(1)
    
    expokit_gpadm(ideg, m, t_wrap, H, ldh, wsp, lwsp, ipiv, iexph, ns, iflag)
    
    ind_start <- iexph[1]
    ind_end <- ind_start + nrows * nrows - 1
    
    expm <- matrix(wsp[ind_start:ind_end], nrow = nrows, ncol = ncols)
    
    return(expm)
  }
)