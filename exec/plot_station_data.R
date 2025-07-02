library(data.table)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(scico)
library(downscaleToPoint)
library(ncdf4)

# Define all necessary paths
# ------------------------------------------------------------------------------
data_dir = file.path(here::here(), "raw_data")
image_dir = file.path(data_dir, "images")
if (!dir.exists(image_dir)) dir.create(image_dir)

meta_path = file.path(data_dir, "meta.rds")

# Load the meta data
# ------------------------------------------------------------------------------
station_meta = readRDS(file.path(data_dir, "meta.rds"))

# Find out which stations that have good precip data, good tmean data or both
station_meta[, let(
  has_tmean = n_tmean > 200,
  has_precip = n_precip > 200 & n_good_precip_flag > 200 & n_unique_precip > 40,
  tag = NA_character_
)]
station_meta[has_tmean & has_precip, let(tag = "TP")]
station_meta[has_tmean & !has_precip, let(tag = "T")]
station_meta[has_precip & !has_tmean, let(tag = "P")]

# Load the CPRCM data (this must be run on hauktjern)
# ------------------------------------------------------------------------------

cprcm_dir = "/bigdisk3/pro/cprcm/"
cprcm_raw_files = list.files(file.path(cprcm_dir, "raw"), full.names = TRUE)

cprcm_raw_files = grep("ERAINT", cprcm_raw_files, value = TRUE)
cprcm_raw_files = grep("v1_day_\\d", cprcm_raw_files, value = TRUE)

nc = nc_open(cprcm_raw_files[1])
nc_lon = ncvar_get(nc, "lon")
nc_lat = ncvar_get(nc, "lat")
nc_close(nc)

cprcm_poly = cbind(
  c(
    nc_lon[nrow(nc_lon), 1:ncol(nc_lon)],
    nc_lon[nrow(nc_lon):1, ncol(nc_lon)],
    nc_lon[1, ncol(nc_lon):1],
    nc_lon[1:nrow(nc_lon), 1]
  ),
  c(
    nc_lat[nrow(nc_lon), 1:ncol(nc_lon)],
    nc_lat[nrow(nc_lon):1, ncol(nc_lon)],
    nc_lat[1, ncol(nc_lon):1],
    nc_lat[1:nrow(nc_lon), 1]
  )
)
cprcm_poly = st_polygon(list(cprcm_poly))
cprcm_poly = st_sf(geometry = st_sfc(cprcm_poly), crs = 4326)

# Plot a map with counts of the number of weather stations inside different hexagons
# ------------------------------------------------------------------------------

plot_data = rbind(
  station_meta[tag %in% c("T", "TP")][, let(tag = "T")],
  station_meta[tag %in% c("P", "TP")][, let(tag = "P")]
)
plot_data[, let(
  tag = factor(
    tag,
    levels = c("T", "P"),
    labels = c("Temperature", "Precipitation")
  )
)]
plot_data = sf::st_as_sf(
  plot_data,
  coords = c("lon", "lat"),
  crs = sf::st_crs(4326)
)
plot_data$lon = sf::st_coordinates(plot_data)[, 1]
plot_data$lat = sf::st_coordinates(plot_data)[, 2]

map = rnaturalearth::ne_countries(returnclass = "sf", scale = 110)

plot = ggplot() +
  geom_sf(data = map) +
  geom_hex(data = plot_data, aes(x = lon, y = lat), bins = 15) +
  geom_sf(data = map, fill = NA, color = "white") +
  geom_sf(data = cprcm_poly, fill = NA, color = "white", linewidth = 1, linetype = "dashed") +
  scale_fill_viridis_c() +
  facet_wrap(~tag) +
  coord_sf(
    xlim = sf::st_bbox(plot_data)[c(1, 3)] + c(-.5, 1.5),
    ylim = sf::st_bbox(plot_data)[c(2, 4)] + c(-1, 1),
    crs = sf::st_crs(plot_data)
  ) +
  guides(col = "none") +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black"),
    text = element_text(size = 15),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")
  ) +
  labs(col = "", x = "", y = "", fill = "Count") +
  scale_x_continuous(
    breaks = seq(-20, 40, by = 20),
    labels = paste0("$", c(20, 0, 20, 40), "^\\circ$", c("W", "", "E", "E")),
    minor_breaks = seq(-10, 30, by = 10)
  ) +
  scale_y_continuous(
    breaks = seq(40, 70, by = 10),
    labels = paste0("$", seq(40, 70, by = 10), "^\\circ$N"),
    minor_breaks = seq(35, 65, by = 10)
  )

plot_tikz(
  file = file.path(image_dir, "map.pdf"),
  plot = plot,
  width = 10,
  height = 4.5
)
