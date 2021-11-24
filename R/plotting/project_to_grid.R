project_to_grid = function(
  data, 
  coord = c('lon', 'lat'), 
  data_crs = CRS('+proj=longlat +datum=WGS84 +no_defs'), 
  rasterObj, 
  plot_crs = crs(rasterObj)
) {
  # prepare spatial data for plotting with ggplot2.  create a data.frame with 
  # lon/lat coordinates projected onto a gridded product, noting that the 
  # gridded product itself may need to be projected first in order to prepare 
  # data for plotting.  
  #
  # if the gridded product needs to be projected first, then the grid will be 
  # interpolated to ensure it still forms a regular grid in projected space.
  #
  # Parameters:
  #   data - tidy data.frame with data to plot
  #   coord - name of columns in data containing lon and lat coords (or x, y)
  #   data_crs - crs of source data
  #   rasterObj - gridded object to project coordinates onto
  #   plot_crs - plotting crs

  plot_crs = crs(plot_crs)
  
  # wrap data in a spatial format to allow projection
  data_spdf = SpatialPoints(
    coords = data[,coord],
    proj4string = data_crs
  )
  
  # project data
  data_projected = spTransform(data_spdf, plot_crs)
  
  # project the raster if needed
  if(!identical(crs(rasterObj), plot_crs)) {
    message('Projecting plotting domain in rasterObj...')
    rasterObj = projectRaster(from = rasterObj, crs = plot_crs)
    message('Projection complete.')
  }
  
  # map the projected data coordinates to the raster grid
  proj_coords = coordinates(rasterObj)
  data_coords = coordinates(data_projected)
  for(ind in 1:2) {
    gridvals = unique(proj_coords[,ind])
    data[, paste(coord[ind], '.projected', sep = '')] = sapply(
      X = data_coords[,ind], 
      FUN = function(x) {
        gridvals[which.min(abs(x - gridvals))]
      }
    )
  }
  
  # return enriched data.frame
  data
}
