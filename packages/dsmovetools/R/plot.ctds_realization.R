#' Plot a realization of a CTDS object
#' 
#' Method allows for overlaying observations
#'   
#' @param x \code{ctds_realization} object for plotting
#' @param ctds_struct \code{ctds_struct} object with domain information to 
#'   assist with plotting
#' @param pt_size When plotting a spatial domain defined via point locations, 
#'   \code{pt_size} controls the size of the points in the grid
#' @param pt_col When plotting a spatial domain defined via point locations, 
#'   \code{pt_size} controls the color of the points in the grid
#' @param jitter.factor CTDS realizations may be plotted with their coordinates 
#'   jittered slightly to aid in readability.  This argument is a number that
#'   specifies the \code{factor} argument to \code{base::jitter}.
#' @param ctds_obs \code{ctds_observations} object containing observations of 
#'   the \code{ctds_realization} to plot.
#' @param ... Additional plotting arguments 
#' 
#' 
#' @import ggplot2 ggthemes
#' 
#' @export
#'
plot.ctds_realization = function(
  x, ctds_struct, pt_size = .1, pt_col = 'grey60', jitter.factor = 0, 
  ctds_obs = NULL, ...) {
  
  # associate (jittered) coordinates with every observation in x
  path.df = data.frame(
    x = jitter(ctds_struct$coords$x[x$states], factor = jitter.factor),
    y = jitter(ctds_struct$coords$y[x$states], factor = jitter.factor),
    ind = 1:length(x$states)
  )
  
  # add coordinate of next location in tidy format
  path.df$xend = c(path.df$x[-1], NA)
  path.df$yend = c(path.df$y[-1], NA)
  
  # build plot
  pl = ggplot(path.df, aes(x = x, y = y, xend = xend, yend = yend)) + 
    # plot spatial domain as points
    geom_point(data = ctds_struct$coords, mapping = aes(x = x, y = y), 
               size = pt_size, col = pt_col, inherit.aes = FALSE) + 
    # underlay path of arrows to show travel of path
    geom_segment(linejoin = 'bevel', lineend = 'round') +
    # overlay points at state changes
    geom_point(mapping = aes(col = ind)) +
    # formatting
    theme_few() + 
    scale_color_distiller(type = 'seq', palette = 'Greys', direction = 1) +
    guides(col = 'none')
  
  # add observation layer
  if(!is.null(ctds_obs)) {
    
    obs.df = data.frame(
      x = jitter(ctds_struct$coords$x[ctds_obs$states], factor = jitter.factor),
      y = jitter(ctds_struct$coords$y[ctds_obs$states], factor = jitter.factor)
    )
    
    pl = pl + geom_point(data = obs.df, mapping = aes(x = x, y = y), 
                         col = 'red', inherit.aes = FALSE)
  }
  
  pl
}