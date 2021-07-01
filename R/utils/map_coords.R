map_coords = function(lon, lat, lon_grid, lat_grid, coord_grid) {
  # Returns the index of the closest gridded longitude and latitude
  # 
  # Parameters:
  #  lon - vector of target longitudes
  #  lat - vector of target latitudes
  #  lon_grid - vector of unique longitudes in grid
  #  lat_grid - vector of unique latitudes in grid
  
  do.call(rbind, lapply(1:length(lon), function(ind) {
    coord_ind = which.min(
      sp::spDists(
        x = matrix(c(lon[ind], lat[ind]), nrow = 1), 
        y = coord_grid, 
        longlat = TRUE
      )
    )
    data.frame(
      lon_ind = which(coord_grid[coord_ind,1] == lon_grid),
      lat_ind = which(coord_grid[coord_ind,2] == lat_grid)
    )
  }))
}
