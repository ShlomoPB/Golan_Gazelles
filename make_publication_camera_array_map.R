#!/usr/bin/env Rscript

# Map the camera array and sampling strata from project-owned coordinates.
# This deliberately avoids adding an unverified external boundary layer.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(readxl)
  library(sf)
})

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  coordinate_file = "Camera Coordinates.xlsx",
  site_cov_file = "SiteCovs_060324.csv",
  output_prefix = "publication_camera_array_map"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

coordinates <- read_excel(arg_map$coordinate_file)
names(coordinates)[1:3] <- c("camera_label", "latitude", "longitude")
coordinates <- coordinates %>%
  transmute(
    Camera.ID = as.integer(gsub("[^0-9]", "", as.character(camera_label))),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  ) %>%
  filter(!is.na(Camera.ID), is.finite(latitude), is.finite(longitude))

site_covariates <- read_csv(arg_map$site_cov_file, show_col_types = FALSE) %>%
  transmute(
    Camera.ID = as.integer(Camera.ID),
    sampling_stratum = Polygon.Name,
    stratum_type = Polygon.Type
  )

sites <- site_covariates %>%
  left_join(coordinates, by = "Camera.ID")

if (nrow(sites) != 60 || n_distinct(sites$Camera.ID) != 60) {
  stop("Expected 60 unique cameras after joining site covariates and coordinates.")
}
if (any(!is.finite(sites$latitude)) || any(!is.finite(sites$longitude))) {
  stop("At least one camera lacks a valid coordinate.")
}

sites_sf <- st_as_sf(
  sites,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(32636)

strata_sf <- sites_sf %>%
  group_by(sampling_stratum, stratum_type) %>%
  summarise(
    n_cameras = n(),
    geometry = st_convex_hull(st_union(geometry)),
    .groups = "drop"
  )

label_points <- suppressWarnings(st_point_on_surface(strata_sf))
label_xy <- st_coordinates(label_points)
label_data <- st_drop_geometry(label_points) %>%
  mutate(x = label_xy[, 1], y = label_xy[, 2])

bbox <- st_bbox(sites_sf)
x_range <- unname(bbox["xmax"] - bbox["xmin"])
y_range <- unname(bbox["ymax"] - bbox["ymin"])
x_pad <- 0.08 * x_range
y_pad <- 0.06 * y_range
x_min <- unname(bbox["xmin"] - x_pad)
x_max <- unname(bbox["xmax"] + x_pad)
y_min <- unname(bbox["ymin"] - y_pad)
y_max <- unname(bbox["ymax"] + y_pad)

scale_length_m <- 10000
scale_x0 <- x_min + 0.06 * (x_max - x_min)
scale_y0 <- y_min + 0.045 * (y_max - y_min)
north_x <- x_max - 0.07 * (x_max - x_min)
north_y0 <- y_min + 0.06 * (y_max - y_min)
north_y1 <- north_y0 + 0.10 * (y_max - y_min)

map_plot <- ggplot() +
  geom_sf(
    data = strata_sf,
    fill = "#DDE4D5",
    color = "#778473",
    linewidth = 0.55,
    alpha = 0.72
  ) +
  geom_sf(
    data = sites_sf,
    shape = 21,
    fill = "#173F66",
    color = "white",
    stroke = 0.55,
    size = 2.7
  ) +
  geom_text_repel(
    data = label_data,
    aes(x, y, label = sampling_stratum),
    size = 2.8,
    fontface = "bold",
    color = "#263337",
    box.padding = 0.38,
    point.padding = 0.20,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 20260818,
    segment.color = "#87918A",
    segment.size = 0.35
  ) +
  annotate(
    "segment",
    x = scale_x0,
    xend = scale_x0 + scale_length_m,
    y = scale_y0,
    yend = scale_y0,
    linewidth = 1.15,
    color = "#1F292D"
  ) +
  annotate(
    "segment",
    x = c(scale_x0, scale_x0 + scale_length_m),
    xend = c(scale_x0, scale_x0 + scale_length_m),
    y = scale_y0 - 420,
    yend = scale_y0 + 420,
    linewidth = 0.7,
    color = "#1F292D"
  ) +
  annotate(
    "text",
    x = scale_x0 + scale_length_m / 2,
    y = scale_y0 + 1100,
    label = "10 km",
    size = 3.2,
    color = "#1F292D"
  ) +
  annotate(
    "segment",
    x = north_x,
    xend = north_x,
    y = north_y0,
    yend = north_y1,
    linewidth = 0.8,
    arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed")
  ) +
  annotate(
    "text",
    x = north_x,
    y = north_y1 + 1100,
    label = "N",
    fontface = "bold",
    size = 4
  ) +
  coord_sf(
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max),
    expand = FALSE,
    datum = st_crs(32636)
  ) +
  scale_x_continuous(labels = function(x) format(x / 1000, trim = TRUE)) +
  scale_y_continuous(labels = function(x) format(x / 1000, trim = TRUE)) +
  labs(
    x = "Easting (km)",
    y = "Northing (km)",
    title = "Camera-trap sampling array in the Golan Heights",
    subtitle = "Sixty cameras in 12 five-camera sampling strata",
    caption = paste0(
      "Points show camera locations; shaded hulls and labels show sampling strata.\n",
      "WGS 84 / UTM Zone 36N. No regional or political boundary is inferred."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E8E6", linewidth = 0.35),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4B5A5E"),
    plot.caption = element_text(hjust = 0, color = "#5B686C", lineheight = 1.05)
  )

write_csv(st_drop_geometry(sites_sf), paste0(arg_map$output_prefix, "_sites.csv"))
write_csv(st_drop_geometry(strata_sf), paste0(arg_map$output_prefix, "_sampling_strata.csv"))

ggsave(
  paste0(arg_map$output_prefix, ".png"),
  map_plot,
  width = 7.4,
  height = 8.5,
  dpi = 320
)
ggsave(
  paste0(arg_map$output_prefix, ".pdf"),
  map_plot,
  width = 7.4,
  height = 8.5
)

metadata <- c(
  paste0("Generated: ", format(Sys.time(), tz = "UTC"), " UTC"),
  paste0("Coordinate source: ", arg_map$coordinate_file),
  paste0("Sampling-stratum source: ", arg_map$site_cov_file),
  "Input CRS: EPSG:4326",
  "Plot CRS: EPSG:32636 (WGS 84 / UTM zone 36N)",
  "Scope: camera-array and sampling-stratum layout only; no boundary layer used."
)
writeLines(metadata, paste0(arg_map$output_prefix, "_metadata.txt"))

cat("Saved provisional publication camera-array map with prefix:", arg_map$output_prefix, "\n")
