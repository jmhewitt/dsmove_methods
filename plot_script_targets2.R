library(raster)
library(sp)
library(metR)
library(ggplot2)
library(ggthemes)
library(ggnewscale)
library(lubridate)
library(dplyr)

# load data
targets::tar_load(whale_pkg)
pkg = whale_pkg$whale_pkg_5bc05e7c

# load posteriors
post_marginal = readRDS('marginal_locations.rds')
# post_crawl = readRDS('crawl_discretized_post_locations.rds')
post_crawl = readRDS('crawl_discretized_post_locations100.rds')

# extract prediction times
pred_times = as.POSIXct(
  x = sapply(post_marginal, function(r) r$meta$time), 
  tz = 'UTC', origin = '1970-01-01 00:00.00 UTC'
)

targets::tar_load(cee_start)

# # manifest of marginal distributions
# do.call(rbind, lapply(post_marginal, function(r) r$meta))

# target prediction time
ptime_without_depth = 18
ptime_with_depth = 6

targets::tar_load(whale_domain)

# CTDS-based posterior location
x_with_depth = post_marginal[[ptime_with_depth]]$posterior
x_without_depth = post_marginal[[ptime_without_depth]]$posterior

# normalize posterior
x_with_depth$lp = x_with_depth$lp - dsmovetools2d:::log_sum_c(x_with_depth$lp)
x_without_depth$lp = x_without_depth$lp - 
  dsmovetools2d:::log_sum_c(x_without_depth$lp)

# crawl-based posterior location
x_crawl = do.call(rbind, lapply(post_crawl, function(imputation) {
  ind = which.min(abs(imputation$time - pred_times[ptime_without_depth]))
  imputation[ind, c('mapped_lon', 'mapped_lat')]
}))

# condense to empirical probability distribution
x_crawl = x_crawl %>% 
  mutate(lp = log(1 /n())) %>%
  group_by(mapped_lon, mapped_lat) %>%
  summarise(lp = dsmovetools2d:::log_sum_c(lp)) %>% 
  ungroup()

#
# map data to grid
#

# map the observed lat-lons to the bathymetry lon/lat grid
# (needed to set initial filtering locations)
coords = pkg$grid$coords
lons = pkg$grid$lons
lats = pkg$grid$lats

# extract coordinates from CTDS via time discretization
x_with_depth$lon = lons[x_with_depth$lon_to_ind]
x_with_depth$lat = lats[x_with_depth$lat_to_ind]
x_without_depth$lon = lons[x_without_depth$lon_to_ind]
x_without_depth$lat = lats[x_without_depth$lat_to_ind]

# plotting projection
plot_crs = '+proj=aea +lat_1=34.5 +lat_2=36 +lon_0=-75'


locs.unprojected.depth = SpatialPointsDataFrame(
  coords = x_with_depth[,c('lon','lat')],
  data = data.frame(lp = x_with_depth$lp),
  proj4string = crs(whale_domain)
)

locs.projected.depth = spTransform(
  locs.unprojected.depth, CRS(plot_crs)
)

locs.unprojected.without_depth = SpatialPointsDataFrame(
  coords = x_without_depth[,c('lon','lat')],
  data = data.frame(lp = x_without_depth$lp),
  proj4string = crs(whale_domain)
)

locs.projected.without_depth = spTransform(
  locs.unprojected.without_depth, CRS(plot_crs)
)

locs.crawl.projected = spTransform(
  x = SpatialPointsDataFrame(
    coords = x_crawl, data = data.frame(lp = x_crawl$lp),
    proj4string = crs(whale_domain)
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

df.locs.depth = data.frame(locs.projected.depth)
df.locs.depth$lon.gridded = sapply(df.locs.depth$lon, function(lon) {
  lons.projected[which.min(abs(lon - lons.projected))]
})
df.locs.depth$lat.gridded = sapply(df.locs.depth$lat, function(lat) {
  lats.projected[which.min(abs(lat - lats.projected))]
})

df.locs.without_depth = data.frame(locs.projected.without_depth)
df.locs.without_depth$lon.gridded = sapply(df.locs.without_depth$lon, function(lon) {
  lons.projected[which.min(abs(lon - lons.projected))]
})
df.locs.without_depth$lat.gridded = sapply(df.locs.without_depth$lat, function(lat) {
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

df.contour = df

# df$depth[
#   df$depth <
#     abs(pkg$mapped_data[which(pkg$mapped_times == pred_times[ptime]), 'depth'])
# ] = NA

df.locs.combined = rbind(
  cbind(df.locs.depth %>% 
          select(lon.gridded, lat.gridded, lp) %>% 
          filter(lp >= log(.005)), 
        Model = 'Time discretization (Surface + Depth)'),
  cbind(df.locs.without_depth %>% 
          select(lon.gridded, lat.gridded, lp) %>% 
          filter(lp >= log(.005)), 
        Model = 'Time discretization (Surface only)'),
  cbind(df.locs.crawl %>% select(lon.gridded, lat.gridded, lp), 
        Model = 'CTCRW as AID')
)
df.locs.combined$Model = factor(
  df.locs.combined$Model, 
  levels = sort(unique(df.locs.combined$Model))[c(1,3,2)]
)



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
  # marginal posterior density
  geom_point(mapping = aes(x = lon.gridded, y = lat.gridded, col = exp(lp)), 
             data = df.locs.combined, 
             inherit.aes = FALSE,
             size = 2) +
  # formatting
  facet_wrap(~Model) + 
  scale_color_viridis_c('Prob.', direction = -1) + 
  xlim(15000,28000) + 
  ylim(3765000 - 1e3,3775000 + 2e3) + 
  # bathymetry contours
  geom_contour2(mapping = aes(z = depth), data = df.contour, 
                col = 'grey50', breaks = seq(500,2000,by=100)) + 
  geom_text_contour(mapping = aes(z = depth), data = df.contour, 
                    breaks = seq(500,2000,by=100)) + 
  # formatting
  coord_equal() + 
  theme_few() + 
  # theme(panel.border = element_blank()) + 
  xlab('Easting (m)') + 
  ylab('Northing (m)') +
  ggtitle(paste(pkg$tagID, ' posterior location at ', 
                strftime(pred_times[ptime_with_depth], '%H:%M:%S', tz = 'UTC'), 
                ', observed at ', 
                format(
                  pkg$mapped_data[
                    which(pkg$mapped_times == pred_times[ptime_with_depth]), 'depth'
                  ], big.mark = ','
                ), 'm deep',
                sep = ''))

pl

ggsave(pl, filename = 'zc095_post_loc_comparison.png', 
       width = 18, height = 9, dpi = 'print')


#
# plot the raw data 
#

locs.unprojected.obs = SpatialPointsDataFrame(
    coords = pkg$data$loc[,c('lon','lat')],
    data = data.frame(time = strftime(x = pkg$data$loc$time, 
                                      format = '%H:%M')),
    proj4string = crs(whale_domain)
)

locs.projected.obs = spTransform(
  locs.unprojected.obs, CRS(plot_crs)
)

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
  # observations
  geom_path(mapping = aes(x = lon, y = lat), 
             data = data.frame(locs.projected.obs), 
             inherit.aes = FALSE) +
  geom_point(mapping = aes(x = lon, y = lat), 
             data = data.frame(locs.projected.obs), 
             inherit.aes = FALSE,
             size = 2) +
  # formatting
  xlim(15000,28000) + 
  ylim(3765000 - 1e3,3775000 + 2e3) + 
  # bathymetry contours
  geom_contour2(mapping = aes(z = depth), data = df.contour, 
                col = 'grey50', breaks = seq(500,2000,by=100)) + 
  geom_text_contour(mapping = aes(z = depth), data = df.contour, 
                    breaks = seq(500,2000,by=100)) + 
  # observation times
  ggrepel::geom_label_repel(mapping = aes(x = lon, y = lat, label = time), 
                            data = data.frame(locs.projected.obs), 
                            inherit.aes = FALSE, fill = '#3182bd33') +
  # formatting
  coord_equal() + 
  theme_few() + 
  # theme(panel.border = element_blank()) + 
  xlab('Easting (m)') + 
  ylab('Northing (m)')


pl

ggsave(pl, filename = 'zc095_obs_data.png', 
       width = 6, height = 4, dpi = 'print')
