#' Extract directed CTDS edges that begin and end at locations in locs
#' 
#' @name expokit_gpadm
#'
#' @param x vector of states visited
#' 
#' 
#' @import nimble
#' 

NULL

#' 
#' @export
#' @rdname expokit_gpadm
#' 
expokit_gpadm = nimbleExternalCall(
  prototype = function(ideg = integer(1), m = integer(1), t = double(1),
                       H = double(1), ldh = integer(1), wsp = double(1),
                       lwsp = integer(1), ipiv = integer(1),
                       iexph = integer(1), ns = integer(1),
                       iflag = integer(1)){},
  returnType = void(),
  Cfun = 'wrapdgpadm_',
  headerFile = file.path(find.package('dsmovetools'), 'inst', 'expokit.h'),
  oFile = file.path(find.package('dsmovetools'), 'inst',
                    'expokit_flattened.o')
)