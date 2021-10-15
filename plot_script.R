library(raster)
library(sp)
library(metR)
library(ggplot2)
library(ggthemes)
library(ggnewscale)
library(lubridate)
library(dplyr)

# load depth data
template_bins = read.csv('data/BeakedWhale/template_bins.csv')
depths_raw = read.csv('data/BeakedWhale/ZcTag095_DUML_series_20200703.csv')
depths_raw$Date = as.POSIXct(x = depths_raw$Date, tz = 'UTC', 
                             origin = '1970-01-01 00:00.00 UTC')

depths_raw$depth_mapped = template_bins$center[sapply(depths_raw$Depth, function(d) {
  which.min(abs(d - template_bins$center))
})]

post_marginal = readRDS('post_marginal_map.rds')
post_crawl = readRDS('crawl_discretized_post_locations.rds')


targets::tar_load(cee_start)

# times at which posterior for marginal surface should be made
pred_times = strptime(
  x = gsub(pattern = 'Prediction time: ', replacement = '', 
           x = names(post_marginal$marginal_location)),
  format = '%Y-%m-%d %H:%M:%S', 
  tz = 'GMT'
)

# target prediction time
ptime = 4

targets::tar_load(whale_domain)

# CTDS-based posterior location
x = post_marginal$marginal_location[[ptime]]

# crawl-based posterior location
x_crawl = do.call(rbind, lapply(post_crawl, function(imputation) {
  ind = which(imputation$time == pred_times[ptime])
  imputation[ind, c('mapped_lon', 'mapped_lat')]
}))


#
# map data to grid
#

# map the observed lat-lons to the bathymetry lon/lat grid
# (needed to set initial filtering locations)
coords = coordinates(whale_domain)
lons = unique(coords[,1])
lats = unique(coords[,2])

# extract coordinates from CTDS via time discretization
x$lon = lons[x$lon_to_ind + 1]
x$lat = lats[x$lat_to_ind + 1]

# # extract coordinates from CTDS via time discretization
# x$lon = lons[x$lon_to_ind]
# x$lat = lats[x$lat_to_ind]

# x$lon_to_ind = NULL
# x$lat_to_ind = NULL
# 
# saveRDS(x, file = 'zc93_post_location_2019-08-19_19-15-00.rds')
 
# plotting projection
plot_crs = '+proj=aea +lat_1=34.5 +lat_2=36 +lon_0=-75'


locs.unprojected = SpatialPointsDataFrame(
  coords = x[,c('lon','lat')],
  data = data.frame(lp = x$lp),
  proj4string = crs(whale_domain)
)

locs.projected = spTransform(
  locs.unprojected, CRS(plot_crs)
)

locs.crawl.projected = spTransform(
  x = SpatialPoints(
    coords = x_crawl, proj4string = crs(whale_domain)
  ), 
  CRSobj = CRS(plot_crs)
)



# subsample data, for faster plotting
bathy.downsampled = sampleRegular(x = whale_domain, size = 1e6, asRaster = TRUE)

# projected bathymetry coordinates, for plotting
bathy.projected = projectRaster(from = bathy.downsampled, crs = plot_crs)

# separate land from bathymetry
land.mask = bathy.projected
land.mask[land.mask <= 0] = NA
land.mask[land.mask > 0] = 1
bathy.projected[bathy.projected > 0] = NA
bathy.projected = bathy.projected * -1


coords.projected = coordinates(bathy.projected)
lons.projected = unique(coords.projected[,1])
lats.projected = unique(coords.projected[,2])

df.locs = data.frame(locs.projected)
df.locs$lon.gridded = sapply(df.locs$lon, function(lon) {
  lons.projected[which.min(abs(lon - lons.projected))]
})
df.locs$lat.gridded = sapply(df.locs$lat, function(lat) {
  lats.projected[which.min(abs(lat - lats.projected))]
})

df.locs.crawl = data.frame(locs.crawl.projected)
df.locs.crawl$lon.gridded = sapply(df.locs.crawl$mapped_lon, function(lon) {
  lons.projected[which.min(abs(lon - lons.projected))]
})
df.locs.crawl$lat.gridded = sapply(df.locs.crawl$mapped_lat, function(lat) {
  lats.projected[which.min(abs(lat - lats.projected))]
})


# munge raster to tidy format
df = data.frame(coordinates(bathy.projected), 
                depth = raster::values(bathy.projected))
df.land = data.frame(coordinates(land.mask), 
                     mask = raster::values(land.mask))

df$depth[df$depth > 2000] = NA

# plot gps track against bathymetry and coastline
pl = ggplot(mapping = aes(x = x, y = y)) + 
  # bathymetry
  geom_raster(mapping = aes(fill = depth), data = df, col = 1) +
  scale_fill_distiller('Depth (m)', palette = 'Blues', direction = -1,
                       na.value = 'transparent', trans = 'reverse') +
  # land mask
  new_scale_fill() + 
  guides(fill = 'none') + 
  geom_raster(mapping = aes(fill = mask), data = df.land) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  # marginal posterior support
  geom_point(mapping = aes(x = lon.gridded, y = lat.gridded), 
             data = df.locs, 
             inherit.aes = FALSE,
             pch = 1,
             size = 2) +
  # marginal posterior density
  geom_point(mapping = aes(x = lon.gridded, y = lat.gridded, col = exp(lp)), 
             data = df.locs %>% filter(lp >= log(.05)), 
             inherit.aes = FALSE,
             size = 2) +
  scale_color_viridis_c('Prob.', direction = -1) + 
  xlim(15000,28000) + 
  ylim(3762500,3775000 + 1e3) + 
  # crawl posterior sample for location
  geom_point(mapping = aes(x = lon.gridded, y = lat.gridded),
             data = df.locs.crawl,
             inherit.aes = FALSE,
             size = 1, col = 1, alpha = .8) +
  # bathymetry contours
  geom_contour2(mapping = aes(z = depth), data = df, col = 'grey50',
                breaks = seq(500,2000,by=100)) + 
  geom_text_contour(mapping = aes(z = depth), data = df, 
                    breaks = seq(500,2000,by=100)) + 

# 1 0.3694757 24508.85 3765940     TRUE
# 2 0.5087359 25262.96 3765942     TRUE
  
  # # expanded neighborhood points
  # geom_point(
  #   data = data.frame(coordinates(expanded_nbhd_coords.projected)), 
  #   mapping = aes(x = lon, y = lat), inherit.aes = FALSE, col = 2,
  #   size = .5
  # ) + 
  # observations
  # geom_path(mapping = aes(x = Longitude, y = Latitude), 
  #           inherit.aes = FALSE, alpha = .3, lwd = .5,
  #           data = data.frame(animalobs.projected)) +
  # geom_point(mapping = aes(x = Longitude, y = Latitude), 
  #            inherit.aes = FALSE,
  #            data = data.frame(animalobs.projected), size = .5) +
  # map extent labels
  # scale_x_continuous(
  #   breaks = range(coordinates(obserrors.projected)[,1]) + 50e3 * c(-1,1),
  #   limits = range(coordinates(obserrors.projected)[,1]) + 50e3 * c(-1,1),
  #   oob = scales::oob_keep,
  #   labels = telefit:::lon_trans()$format(
  #     extent_vec[c('min_lon', 'max_lon')]
  #   )) +
  # scale_y_continuous(
  #   breaks = range(coordinates(obserrors.projected)[,2]) + 1e3 * c(-1,1),
  #   limits = range(coordinates(obserrors.projected)[,2]) + 1e3 * c(-1,1),
  #   oob = scales::oob_keep,
  #   labels = telefit:::lat_trans()$format(
  #     extent_vec[c('min_lat', 'max_lat')]
  #   )) +
  # formatting
  coord_equal() + 
  theme_few() + 
  theme(panel.border = element_blank()) + 
  xlab('Easting (m)') + 
  ylab('Northing (m)') +
  ggtitle(paste('Zc093 posterior location at ', 
                strftime(pred_times[ptime], '%H:%M:%S', tz = 'UTC'), sep = ''), 
          subtitle = '(Posterior support for MAP parameters)')

# depths_raw %>% filter(Date == pred_times[ptime])

pl

ggsave(pl, filename = 'zc093_post_support_map.png')
