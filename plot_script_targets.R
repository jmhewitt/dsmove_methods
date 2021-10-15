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
post_marginal = readRDS('post_marginal_map_updated.rds')
post_crawl = readRDS('crawl_discretized_post_locations.rds')

# extract prediction times
pred_times = post_marginal$times
tz(pred_times) = 'UTC'

targets::tar_load(cee_start)

# target prediction time
ptime = 6

targets::tar_load(whale_domain)

# CTDS-based posterior location
x = post_marginal$marginal_location[[ptime]]

# crawl-based posterior location
x_crawl = do.call(rbind, lapply(post_crawl, function(imputation) {
  ind = which.min(abs(imputation$time - pred_times[ptime]))
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
x$lon = lons[x$lon_to_ind]
x$lat = lats[x$lat_to_ind]

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

df.contour = df

df$depth[
  df$depth <
    pkg$mapped_data[which(pkg$mapped_times == pred_times[ptime]), 'depth']
] = NA

df.locs.combined = rbind(
  cbind(df.locs %>% 
          select(lon.gridded, lat.gridded, lp) %>% 
          filter(lp >= log(.005)), Model = 'Time discretization'),
  cbind(df.locs.crawl %>% select(lon.gridded, lat.gridded, lp), 
        Model = 'CRAWL as AID')
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
  # # marginal posterior support
  # geom_point(mapping = aes(x = lon.gridded, y = lat.gridded),
  #            data = df.locs,
  #            inherit.aes = FALSE,
  #            pch = 1,
  #            size = 2) +
  # marginal posterior density
  # geom_point(mapping = aes(x = lon.gridded, y = lat.gridded, col = exp(lp)), 
  #            data = df.locs %>% filter(lp >= log(.005)), 
  #            inherit.aes = FALSE,
  #            size = 2) +
  geom_point(mapping = aes(x = lon.gridded, y = lat.gridded, col = exp(lp)), 
             data = df.locs.combined, 
             inherit.aes = FALSE,
             size = 2) +
  facet_wrap(~Model) + 
  scale_color_viridis_c('Prob.', direction = -1) + 
  xlim(15000,28000) + 
  ylim(3762500,3775000 + 1e3) + 
  # crawl posterior sample for location
  # geom_point(mapping = aes(x = lon.gridded, y = lat.gridded, col = exp(lp)), 
  #            data = df.locs.crawl, 
  #            inherit.aes = FALSE,
  #            size = 2, shape = 18) +
  # bathymetry contours
  geom_contour2(mapping = aes(z = depth), data = df.contour, 
                col = 'grey50', breaks = seq(500,2000,by=100)) + 
  geom_text_contour(mapping = aes(z = depth), data = df.contour, 
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
  # theme(panel.border = element_blank()) + 
  xlab('Easting (m)') + 
  ylab('Northing (m)') +
  ggtitle(paste(pkg$tagID, ' posterior location at ', 
                strftime(pred_times[ptime], '%H:%M:%S', tz = 'UTC'), 
                ', observed at ', 
                format(
                  pkg$mapped_data[
                    which(pkg$mapped_times == pred_times[ptime]), 'depth'
                  ], big.mark = ','
                ), 'm deep',
                sep = ''), 
          subtitle = '(Posterior near MAP parameters; feasible depths in blue)')



pl

ggsave(pl, filename = 'zc093_post_loc_comparison.pdf', 
       width = 18, height = 9)
