#' Plot observations of a CTDS object
#' 
#' Method allows for underlaying complete trajectories
#' 
#'   
#' @param x \code{ctds_observation} object for plotting
#' @param ctds_struct \code{ctds_struct} object with domain information to 
#'   assist with plotting
#' @param pt_size When plotting a spatial domain defined via point locations, 
#'   \code{pt_size} controls the size of the points in the grid
#' @param pt_col When plotting a spatial domain defined via point locations, 
#'   \code{pt_size} controls the color of the points in the grid
#' @param jitter.factor CTDS realizations may be plotted with their coordinates 
#'   jittered slightly to aid in readability.  This argument is a number that
#'   specifies the \code{factor} argument to \code{base::jitter}.
#' @param ctds_realization \code{ctds_realization} object containing 
#'   observations of the \code{ctds_realization} to plot, or a list of such 
#'   objects to plot.
#' @param realization_alpha amount of transparency with which to underlay the 
#'   \code{ctds_realization} object
#'   
#' @param ... Additional plotting arguments 
#' 
#' 
#' @import ggplot2 ggthemes
#' 
#' @export
#'
plot.ctds_observations = function(
  x, ctds_struct, pt_size = .1, pt_col = 'grey60', jitter.factor = 0, 
  ctds_realization = NULL, realization_alpha = .55, ...) {
  
  # associate (jittered) coordinates with every observation in x
  obs.df = data.frame(
    x = jitter(ctds_struct$coords$x[x$states], factor = jitter.factor),
    y = jitter(ctds_struct$coords$y[x$states], factor = jitter.factor),
    ind = 1:length(x$states)
  )
  
  # build plot
  pl = ggplot(obs.df, aes(x = x, y = y, col = ind)) + 
    # plot spatial domain as points
    geom_point(data = ctds_struct$coords, mapping = aes(x = x, y = y), 
               size = pt_size, col = pt_col, inherit.aes = FALSE) + 
    # formatting
    theme_few() + 
    scale_color_distiller(type = 'seq', palette = 'Reds', direction = 1) +
    guides(col = 'none')
  
  # add realization layer
  if(!is.null(ctds_realization)) {
    
    if(inherits(ctds_realization, 'ctds_realization')) {
      ctds_realization.list = list(ctds_realization)
    } else {
      ctds_realization.list = ctds_realization
    }
    
    for(i in 1:length(ctds_realization.list)) {
      
      ctds_realization = ctds_realization.list[[i]]
      
      path.df = data.frame(
        x = jitter(ctds_struct$coords$x[ctds_realization$states], 
                   factor = jitter.factor),
        y = jitter(ctds_struct$coords$y[ctds_realization$states], 
                   factor = jitter.factor),
        ind = 1:length(ctds_realization$states)
      )
      
      # add coordinate of next location in tidy format
      path.df$xend = c(path.df$x[-1], NA)
      path.df$yend = c(path.df$y[-1], NA)
      
      pl = pl + 
        # underlay travel of path
        geom_segment(mapping = aes(x = x, y = y, xend = xend, yend = yend), 
                     data = path.df, linejoin = 'bevel', lineend = 'round', 
                     inherit.aes = FALSE, alpha = realization_alpha) +
        # overlay points at state changes
        geom_point(mapping = aes(x = x, y = y), data = path.df, 
                   inherit.aes = FALSE, alpha = realization_alpha)
      
    }
    
  }
  
  # add observation layer
  pl = pl + geom_point()
  
  pl
}