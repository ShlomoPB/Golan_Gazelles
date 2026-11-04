library(readxl)
library(sf)

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  species = "fawn",
  count_file = "gazelle_juvenile_abundance_60x7.csv",
  output_file = "gazelle_fawn_spNMix_data.rds",
  site_cov_file = "SiteCovs_060324.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  coords_file = "Camera Coordinates.xlsx"
)

for (arg in args) {
  pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(pieces) == 2) {
    key <- sub("^--", "", pieces[1])
    arg_map[[key]] <- pieces[2]
  }
}

species_label <- arg_map$species
count_file <- arg_map$count_file
site_cov_file <- arg_map$site_cov_file
obs_cov_file <- arg_map$obs_cov_file
coords_file <- arg_map$coords_file
output_file <- arg_map$output_file

read_counts <- function(path) {
  x <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  stopifnot("camera_name" %in% names(x))
  y <- as.matrix(x[paste0("X", 1:7)])
  storage.mode(y) <- "integer"
  rownames(y) <- x$camera_name
  y
}

z_scale <- function(x) {
  sd_x <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sd_x) || sd_x == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / sd_x
}

to_site_by_round <- function(df, column, site_ids) {
  matrix(
    df[[column]],
    nrow = length(site_ids),
    ncol = 7,
    byrow = TRUE,
    dimnames = list(site = site_ids, replicate = paste0("X", 1:7))
  )
}

y <- read_counts(count_file)
site_ids <- rownames(y)

site_covs_raw <- read.csv(site_cov_file, na.strings = c("NA", ""), check.names = FALSE)
site_covs_raw <- site_covs_raw[match(as.integer(site_ids), site_covs_raw$Camera.ID), ]
if (any(is.na(site_covs_raw$Camera.ID))) {
  stop("Site covariates are missing one or more camera IDs found in the count file.")
}

abund_keep <- c(
  "Nature.Reserves",
  "Open.Spaces",
  "Built.and.Roads",
  "null.jackal.abundance",
  "null.boar.abundance",
  "Null.Wolf.Abundance",
  "Wolf.pack.abundance",
  "Minefields",
  "Firing.Zones",
  "Proximity.to.Water",
  "human.abundance"
)

abund_covs <- site_covs_raw[abund_keep]
abund_covs[] <- lapply(abund_covs, as.numeric)
abund_covs_scaled <- as.data.frame(lapply(abund_covs, z_scale))
names(abund_covs_scaled) <- paste0(make.names(names(abund_covs)), ".z")
abund_covs <- cbind(abund_covs, abund_covs_scaled)
rownames(abund_covs) <- site_ids

obs_covs_raw <- read.csv(obs_cov_file, na.strings = c("NA", ""), check.names = FALSE)
obs_covs_raw <- obs_covs_raw[order(obs_covs_raw$`Camera Number`, obs_covs_raw$`Round Number`), ]

expected <- expand.grid(`Camera Number` = as.integer(site_ids), `Round Number` = 1:7)
expected <- expected[order(expected$`Camera Number`, expected$`Round Number`), ]
if (!identical(obs_covs_raw$`Camera Number`, expected$`Camera Number`) ||
    !identical(obs_covs_raw$`Round Number`, expected$`Round Number`)) {
  stop("Observation covariates are not exactly one row per camera x round in site order.")
}

days_operational <- to_site_by_round(obs_covs_raw, "Days.Operational", site_ids)
cattle_sightings <- to_site_by_round(obs_covs_raw, "Cattle.Sightings", site_ids)

det_covs <- list(
  days_operational = days_operational,
  cattle_sightings = cattle_sightings,
  log_days_operational.z = matrix(
    z_scale(log(days_operational)),
    nrow = nrow(days_operational),
    ncol = ncol(days_operational),
    dimnames = dimnames(days_operational)
  ),
  cattle_sightings.z = matrix(
    z_scale(cattle_sightings),
    nrow = nrow(cattle_sightings),
    ncol = ncol(cattle_sightings),
    dimnames = dimnames(cattle_sightings)
  ),
  path_type = to_site_by_round(obs_covs_raw, "Path.Type", site_ids),
  camera_model = to_site_by_round(obs_covs_raw, "Camera.Model", site_ids),
  camera_sensitivity = to_site_by_round(obs_covs_raw, "Camera.Sensitivity", site_ids),
  area_behind_camera = to_site_by_round(obs_covs_raw, "Area.Behind.Camera", site_ids),
  vegetation_structure = to_site_by_round(obs_covs_raw, "Vegetation.Structure", site_ids)
)

y_missing <- is.na(y)
days_missing <- is.na(days_operational)
if (!isTRUE(all.equal(unname(y_missing), unname(days_missing)))) {
  warning(species_label, " count missingness does not exactly match Days.Operational missingness.")
}

coords_raw <- read_excel(coords_file)
names(coords_raw) <- c("site", "lat", "lon")
coords_raw$site <- as.integer(gsub("[^0-9]", "", as.character(coords_raw$site)))
coords_raw <- coords_raw[match(as.integer(site_ids), coords_raw$site), ]
if (any(is.na(coords_raw$site))) {
  stop("Camera coordinates are missing one or more sites found in the count file.")
}

coords_sf <- st_as_sf(coords_raw, coords = c("lon", "lat"), crs = 4326)
coords_proj <- st_transform(coords_sf, 2039)
coords_mat <- st_coordinates(coords_proj) / 1000
rownames(coords_mat) <- site_ids
colnames(coords_mat) <- c("x_km", "y_km")

data_list <- list(
  y = y,
  abund.covs = abund_covs,
  det.covs = det_covs,
  coords = coords_mat,
  metadata = list(
    species = species_label,
    count_file = count_file,
    projected_crs = 2039,
    coordinate_units = "km"
  )
)

saveRDS(data_list, output_file)

cat("Prepared", species_label, "single-species spAbundance data:\n")
cat("  sites:", nrow(y), "\n")
cat("  rounds:", ncol(y), "\n")
cat("  total count:", sum(y, na.rm = TRUE), "\n")
cat("  output:", output_file, "\n")
