library(spOccupancy)

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  species = "fawn",
  habitat_cov = "Nature.Reserves",
  wolf_cov = "Null.Wolf.Abundance",
  spatial = "true",
  det_model = "effort",
  data = "gazelle_fawn_spoccupancy_data.rds",
  n_batch = "400",
  batch_length = "25",
  n_samples = "10000",
  n_burn = "",
  n_thin = "5",
  n_chains = "3",
  n_neighbors = "10",
  seed = "20260610"
)

for (arg in args) {
  pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(pieces) == 2) {
    key <- sub("^--", "", pieces[1])
    arg_map[[key]] <- pieces[2]
  }
}

as_flag <- function(x) {
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

species_label <- arg_map$species
habitat_covariates <- strsplit(arg_map$habitat_cov, ",", fixed = TRUE)[[1]]
habitat_covariates <- trimws(habitat_covariates)
habitat_covariates <- habitat_covariates[nzchar(habitat_covariates)]
wolf_covariate <- arg_map$wolf_cov
spatial_model <- as_flag(arg_map$spatial)
det_model <- arg_map$det_model
data_file <- arg_map$data
n_batch <- as.integer(arg_map$n_batch)
batch_length <- as.integer(arg_map$batch_length)
n_samples <- as.integer(arg_map$n_samples)
n_burn <- if (nzchar(arg_map$n_burn)) as.integer(arg_map$n_burn) else NULL
n_thin <- as.integer(arg_map$n_thin)
n_chains <- as.integer(arg_map$n_chains)
n_neighbors <- as.integer(arg_map$n_neighbors)
seed <- as.integer(arg_map$seed)

data_list <- readRDS(data_file)

wolf_scaled_map <- c(
  "Null.Wolf.Abundance" = "Null.Wolf.Abundance.z",
  "Wolf.pack.abundance" = "Wolf.pack.abundance.z"
)
wolf_scaled <- wolf_scaled_map[[wolf_covariate]]
if (is.null(wolf_scaled)) {
  stop("wolf_cov must be one of: Null.Wolf.Abundance, Wolf.pack.abundance")
}

habitat_scaled_map <- c(
  "Nature.Reserves" = "Nature.Reserves.z",
  "Open.Spaces" = "Open.Spaces.z",
  "Built.and.Roads" = "Built.and.Roads.z"
)
habitat_scaled <- unname(habitat_scaled_map[habitat_covariates])
if (length(habitat_scaled) == 0 || any(is.na(habitat_scaled))) {
  stop("habitat_cov must contain one or more of: Nature.Reserves, Open.Spaces, Built.and.Roads")
}

occ_formula <- stats::as.formula(
  paste("~", paste(c(habitat_scaled, "null.jackal.abundance.z", wolf_scaled), collapse = " + "))
)
det_formula_map <- list(
  effort = ~ log_days_operational.z + cattle_sightings.z,
  design = ~ path_type + camera_sensitivity + area_behind_camera + vegetation_structure
)
det_formula <- det_formula_map[[det_model]]
if (is.null(det_formula)) {
  stop("det_model must be one of: effort, design")
}

make_det_model_frame <- function(det_covs) {
  out <- lapply(det_covs, function(x) {
    if (is.matrix(x) || is.data.frame(x)) {
      as.vector(x)
    } else {
      x
    }
  })
  out <- as.data.frame(out, stringsAsFactors = FALSE)

  for (nm in names(out)) {
    if (is.character(out[[nm]])) {
      out[[nm]] <- factor(out[[nm]])
    }
  }
  out
}

X_occ <- model.matrix(occ_formula, data = data_list$occ.covs)
X_det <- model.matrix(det_formula, data = make_det_model_frame(data_list$det.covs))
p_det <- ncol(X_det)

priors <- list(
  beta.normal = list(mean = rep(0, ncol(X_occ)), var = rep(2.72, ncol(X_occ))),
  alpha.normal = list(mean = rep(0, p_det), var = rep(2.72, p_det))
)

z_init <- apply(data_list$y, 1, max, na.rm = TRUE)
z_init[!is.finite(z_init)] <- 0

set.seed(seed)

base_name <- paste0(
  gsub("[^A-Za-z0-9]+", "_", species_label),
  "_",
  if (spatial_model) "spPGOcc" else "PGOcc",
  "_",
  gsub("[^A-Za-z0-9_]+", "_", paste(habitat_covariates, collapse = "__")),
  "_",
  gsub("[^A-Za-z0-9]+", "_", wolf_covariate),
  "_det_",
  gsub("[^A-Za-z0-9]+", "_", det_model)
)
rds_file <- paste0(base_name, "_fit.rds")
summary_file <- paste0(base_name, "_summary.txt")

if (spatial_model) {
  dists <- as.matrix(stats::dist(data_list$coords))
  positive_dists <- dists[upper.tri(dists)]
  positive_dists <- positive_dists[positive_dists > 0]
  phi_lower <- 3 / max(positive_dists)
  phi_upper <- 3 / min(positive_dists)

  priors$sigma.sq.ig <- c(2, 2)
  priors$phi.unif <- c(phi_lower, phi_upper)

  inits <- list(
    alpha = 0,
    beta = 0,
    z = z_init,
    phi = mean(c(phi_lower, phi_upper)),
    sigma.sq = 1,
    w = rep(0, nrow(data_list$y))
  )

  tuning <- list(phi = 1)

  cat("Fitting spatial spPGOcc model for", species_label, "occupancy\n")
  cat("  habitat covariate(s):", paste(habitat_covariates, collapse = ", "), "\n")
  cat("  wolf covariate:", wolf_covariate, "\n")
  cat("  detection model:", det_model, "\n")
  cat("  phi prior:", sprintf("[%.4f, %.4f]", phi_lower, phi_upper), "\n")

  fit <- spPGOcc(
    occ.formula = occ_formula,
    det.formula = det_formula,
    data = data_list,
    inits = inits,
    priors = priors,
    tuning = tuning,
    cov.model = "exponential",
    NNGP = FALSE,
    n.neighbors = n_neighbors,
    n.batch = n_batch,
    batch.length = batch_length,
    n.burn = if (is.null(n_burn)) round(0.1 * n_batch * batch_length) else n_burn,
    n.thin = n_thin,
    n.chains = n_chains,
    n.omp.threads = 1,
    verbose = TRUE,
    n.report = max(1L, floor(n_batch / 10))
  )
} else {
  inits <- list(
    alpha = 0,
    beta = 0,
    z = z_init
  )

  cat("Fitting non-spatial PGOcc model for", species_label, "occupancy\n")
  cat("  habitat covariate(s):", paste(habitat_covariates, collapse = ", "), "\n")
  cat("  wolf covariate:", wolf_covariate, "\n")
  cat("  detection model:", det_model, "\n")

  fit <- PGOcc(
    occ.formula = occ_formula,
    det.formula = det_formula,
    data = data_list[c("y", "occ.covs", "det.covs")],
    inits = inits,
    priors = priors,
    n.samples = n_samples,
    n.burn = if (is.null(n_burn)) round(0.1 * n_samples) else n_burn,
    n.thin = n_thin,
    n.chains = n_chains,
    n.omp.threads = 1,
    verbose = TRUE,
    n.report = max(1L, floor(n_samples / 10))
  )
}

saveRDS(fit, rds_file)

sink(summary_file)
cat("Species:", species_label, "\n")
cat("Model:", if (spatial_model) "spPGOcc" else "PGOcc", "\n")
cat("Habitat covariate(s):", paste(habitat_covariates, collapse = ", "), "\n")
cat("Wolf covariate:", wolf_covariate, "\n\n")
cat("Detection model:", det_model, "\n\n")
print(summary(fit))
sink()

cat("Saved fit to", rds_file, "\n")
cat("Saved summary to", summary_file, "\n")
