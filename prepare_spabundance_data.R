# Prepare gazelle camera-trap data for spAbundance::msNMix.
# Run from the folder containing the CSV files.

library(readxl)

count_files <- c(
  male = "gazelle_male_abundance_60x7.csv",
  female = "gazelle_female_abundance_60x7.csv",
  juvenile = "gazelle_juvenile_abundance_60x7.csv"
)

site_cov_file <- "SiteCovs_060324.csv"
obs_cov_file <- "ObsCovs_050524.csv"
coords_file <- "Camera Coordinates.xlsx"

read_counts <- function(path) {
  x <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  stopifnot("camera_name" %in% names(x))
  y <- as.matrix(x[paste0("X", 1:7)])
  storage.mode(y) <- "integer"
  rownames(y) <- x$camera_name
  y
}

count_mats <- lapply(count_files, read_counts)

site_ids <- rownames(count_mats[[1]])
if (!all(vapply(count_mats, function(x) identical(rownames(x), site_ids), logical(1)))) {
  stop("Count files do not have identical camera ordering.")
}

y <- array(
  NA_integer_,
  dim = c(length(count_mats), length(site_ids), 7),
  dimnames = list(species = names(count_mats), site = site_ids, replicate = paste0("X", 1:7))
)

for (i in seq_along(count_mats)) {
  y[i, , ] <- count_mats[[i]]
}

site_covs_raw <- read.csv(site_cov_file, na.strings = c("NA", ""), check.names = FALSE)
site_covs_raw <- site_covs_raw[match(as.integer(site_ids), site_covs_raw$Camera.ID), ]
if (any(is.na(site_covs_raw$Camera.ID))) {
  stop("Site covariates are missing one or more camera IDs found in the count files.")
}

abund_keep <- c(
  "Predicted.Wolf.Abundance",
  "Null.Wolf.Abundance",
  "null.jackal.abundance",
  "null.boar.abundance",
  "human.abundance",
  "Nature.Reserves",
  "Minefields",
  "Built.and.Roads",
  "Firing.Zones",
  "Proximity.to.Water",
  "Wolf.pack.abundance",
  "Cattle.Abundance",
  "Cropland",
  "Animal.Structures",
  "Open.Spaces",
  "Orchards"
)

abund_covs <- site_covs_raw[abund_keep]
names(abund_covs) <- make.names(names(abund_covs))
abund_covs[] <- lapply(abund_covs, as.numeric)
abund_covs_scaled <- as.data.frame(scale(abund_covs))
names(abund_covs_scaled) <- paste0(names(abund_covs), ".z")

obs_covs_raw <- read.csv(obs_cov_file, na.strings = c("NA", ""), check.names = FALSE)
obs_covs_raw <- obs_covs_raw[order(obs_covs_raw$`Camera Number`, obs_covs_raw$`Round Number`), ]

expected <- expand.grid(`Camera Number` = as.integer(site_ids), `Round Number` = 1:7)
expected <- expected[order(expected$`Camera Number`, expected$`Round Number`), ]
if (!identical(obs_covs_raw$`Camera Number`, expected$`Camera Number`) ||
    !identical(obs_covs_raw$`Round Number`, expected$`Round Number`)) {
  stop("Observation covariates are not exactly one row per camera x round in site order.")
}

to_site_by_round <- function(column) {
  matrix(
    obs_covs_raw[[column]],
    nrow = length(site_ids),
    ncol = 7,
    byrow = TRUE,
    dimnames = list(site = site_ids, replicate = paste0("X", 1:7))
  )
}

z_matrix <- function(x) {
  (x - mean(x, na.rm = TRUE)) / stats::sd(as.vector(x), na.rm = TRUE)
}

det_covs <- list(
  days_operational = to_site_by_round("Days.Operational"),
  cattle_sightings = to_site_by_round("Cattle.Sightings")
)

det_covs$log_days_operational.z <- z_matrix(log(det_covs$days_operational))
det_covs$cattle_sightings.z <- z_matrix(det_covs$cattle_sightings)

y_missing <- apply(is.na(y), c(2, 3), any)
days_missing <- is.na(det_covs$days_operational)
if (!identical(y_missing, days_missing)) {
  warning(
    "The y missingness pattern does not match Days.Operational missingness. ",
    "Inspect camera 22 rounds 2 and 6, and cameras 22/60 round 6 before fitting."
  )
}

coords <- read_excel(coords_file)
names(coords) <- c("site", "lat", "lon")
coords$site <- as.integer(gsub("[^0-9]", "", as.character(coords$site)))
coords <- coords[match(as.integer(site_ids), coords$site), ]
if (any(is.na(coords$site))) {
  stop("Camera coordinates are missing one or more sites found in the count files.")
}

# data$coords for lfMsNMix latent factor models
coords_mat <- as.matrix(coords[c("lon", "lat")])
rownames(coords_mat) <- coords$site

data_list <- list(
  y = y,
  abund.covs = cbind(abund_covs, abund_covs_scaled),
  det.covs = det_covs,
  coords = coords_mat
)

saveRDS(data_list, "gazelle_msNMix_data_list.rds")

str(data_list)
