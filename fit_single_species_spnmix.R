library(spAbundance)

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  species = "fawn",
  habitat_cov = "Open.Spaces",
  wolf_cov = "Null.Wolf.Abundance",
  include_jackal = "true",
  include_wolf = "true",
  extra_covs = "",
  det_model = "design",
  model_scope = "spatial",
  data = "gazelle_fawn_spNMix_data.rds",
  family = "NB",
  n_batch = "400",
  batch_length = "25",
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

species_label <- arg_map$species
habitat_covariate <- arg_map$habitat_cov
wolf_covariate <- arg_map$wolf_cov
include_jackal <- tolower(arg_map$include_jackal) %in% c("true", "1", "yes")
include_wolf <- tolower(arg_map$include_wolf) %in% c("true", "1", "yes")
extra_covariates_raw <- arg_map$extra_covs
det_model <- arg_map$det_model
model_scope <- arg_map$model_scope
data_file <- arg_map$data
family_name <- arg_map$family
n_batch <- as.integer(arg_map$n_batch)
batch_length <- as.integer(arg_map$batch_length)
n_burn <- if (nzchar(arg_map$n_burn)) as.integer(arg_map$n_burn) else NULL
n_thin <- as.integer(arg_map$n_thin)
n_chains <- as.integer(arg_map$n_chains)
n_neighbors <- as.integer(arg_map$n_neighbors)
seed <- as.integer(arg_map$seed)

data_list <- readRDS(data_file)

habitat_scaled_map <- c(
  "Nature.Reserves" = "Nature.Reserves.z",
  "Open.Spaces" = "Open.Spaces.z",
  "Built.and.Roads" = "Built.and.Roads.z"
)
habitat_scaled <- habitat_scaled_map[[habitat_covariate]]
if (is.null(habitat_scaled)) {
  stop("habitat_cov must be one of: Nature.Reserves, Open.Spaces, Built.and.Roads")
}

wolf_scaled_map <- c(
  "Null.Wolf.Abundance" = "Null.Wolf.Abundance.z",
  "Wolf.pack.abundance" = "Wolf.pack.abundance.z"
)
wolf_scaled <- wolf_scaled_map[[wolf_covariate]]
if (is.null(wolf_scaled)) {
  stop("wolf_cov must be one of: Null.Wolf.Abundance, Wolf.pack.abundance")
}

extra_covariates <- character(0)
if (nzchar(extra_covariates_raw)) {
  extra_covariates <- strsplit(extra_covariates_raw, ",", fixed = TRUE)[[1]]
  extra_covariates <- trimws(extra_covariates)
  extra_covariates <- extra_covariates[nzchar(extra_covariates)]
}

extra_scaled_map <- c(
  "null.boar.abundance" = "null.boar.abundance.z",
  "human.abundance" = "human.abundance.z",
  "Minefields" = "Minefields.z",
  "Firing.Zones" = "Firing.Zones.z",
  "Proximity.to.Water" = "Proximity.to.Water.z"
)
extra_scaled <- unname(extra_scaled_map[extra_covariates])
if (length(extra_covariates) > 0 && any(is.na(extra_scaled))) {
  bad <- extra_covariates[is.na(extra_scaled)]
  stop(
    "Unsupported extra_covs: ",
    paste(bad, collapse = ", "),
    ". Supported values are: ",
    paste(names(extra_scaled_map), collapse = ", ")
  )
}

abund_terms <- c(habitat_scaled)
if (include_jackal) {
  abund_terms <- c(abund_terms, "null.jackal.abundance.z")
}
abund_terms <- c(abund_terms, extra_scaled)
if (include_wolf) {
  abund_terms <- c(abund_terms, wolf_scaled)
}

abund_formula <- stats::as.formula(
  paste("~", paste(abund_terms, collapse = " + "))
)

det_formula_map <- list(
  effort = ~ log_days_operational.z + cattle_sightings.z,
  design = ~ path_type + camera_sensitivity + area_behind_camera + vegetation_structure
)
det_formula <- det_formula_map[[det_model]]
if (is.null(det_formula)) {
  stop("det_model must be one of: effort, design")
}
if (!model_scope %in% c("spatial", "nonspatial")) {
  stop("model_scope must be one of: spatial, nonspatial")
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

X_abund <- model.matrix(abund_formula, data = data_list$abund.covs)
X_det <- model.matrix(det_formula, data = make_det_model_frame(data_list$det.covs))

priors <- list(
  beta.normal = list(mean = rep(0, ncol(X_abund)), var = rep(100, ncol(X_abund))),
  alpha.normal = list(mean = rep(0, ncol(X_det)), var = rep(2.72, ncol(X_det)))
)

if (family_name == "NB") {
  priors$kappa.unif <- c(0, 10)
}

if (model_scope == "spatial") {
  priors$sigma.sq.ig <- c(2, 2)
  dists <- as.matrix(stats::dist(data_list$coords))
  positive_dists <- dists[upper.tri(dists)]
  positive_dists <- positive_dists[positive_dists > 0]
  priors$phi.unif <- c(3 / max(positive_dists), 3 / min(positive_dists))
}

N_init <- apply(data_list$y, 1, max, na.rm = TRUE)
N_init[!is.finite(N_init)] <- 0

inits <- list(
  alpha = 0,
  beta = 0,
  N = N_init
)

if (family_name == "NB") {
  inits$kappa <- 1
}

if (model_scope == "spatial") {
  inits$phi <- mean(priors$phi.unif)
  inits$sigma.sq <- 1
  inits$w <- rep(0, nrow(data_list$y))
}

tuning <- list(beta = 0.1, alpha = 0.1)
if (family_name == "NB") {
  tuning$kappa <- 0.5
}
if (model_scope == "spatial") {
  tuning$phi <- 1
  tuning$w <- 0.1
}

set.seed(seed)

base_name <- paste0(
  gsub("[^A-Za-z0-9]+", "_", species_label),
  "_",
  ifelse(model_scope == "spatial", "spNMix", "NMix"),
  "_",
  gsub("[^A-Za-z0-9]+", "_", habitat_covariate),
  "_jackal_",
  ifelse(include_jackal, "yes", "no"),
  "_wolf_",
  ifelse(include_wolf, "yes", "no"),
  "_extra_",
  ifelse(length(extra_covariates) == 0, "none", paste(gsub("[^A-Za-z0-9]+", "_", extra_covariates), collapse = "_")),
  "_",
  gsub("[^A-Za-z0-9]+", "_", wolf_covariate),
  "_det_",
  gsub("[^A-Za-z0-9]+", "_", det_model),
  "_",
  family_name
)
rds_file <- paste0(base_name, "_fit.rds")
summary_file <- paste0(base_name, "_summary.txt")

cat("Fitting single-species spNMix model for", species_label, "\n")
cat("  habitat covariate:", habitat_covariate, "\n")
cat("  wolf covariate:", wolf_covariate, "\n")
cat("  detection model:", det_model, "\n")
cat("  model scope:", model_scope, "\n")
cat("  family:", family_name, "\n")

if (model_scope == "spatial") {
  fit <- spNMix(
    abund.formula = abund_formula,
    det.formula = det_formula,
    data = data_list,
    inits = inits,
    priors = priors,
    tuning = tuning,
    cov.model = "exponential",
    NNGP = TRUE,
    n.neighbors = n_neighbors,
    n.batch = n_batch,
    batch.length = batch_length,
    n.burn = if (is.null(n_burn)) round(0.1 * n_batch * batch_length) else n_burn,
    n.thin = n_thin,
    n.chains = n_chains,
    family = family_name,
    n.omp.threads = 1,
    verbose = TRUE,
    n.report = max(1L, floor(n_batch / 10))
  )
} else {
  fit <- NMix(
    abund.formula = abund_formula,
    det.formula = det_formula,
    data = data_list,
    inits = inits,
    priors = priors,
    tuning = tuning,
    n.batch = n_batch,
    batch.length = batch_length,
    n.burn = if (is.null(n_burn)) round(0.1 * n_batch * batch_length) else n_burn,
    n.thin = n_thin,
    n.chains = n_chains,
    family = family_name,
    n.omp.threads = 1,
    verbose = TRUE,
    n.report = max(1L, floor(n_batch / 10))
  )
}

saveRDS(fit, rds_file)

sink(summary_file)
cat("Species:", species_label, "\n")
cat("Model:", ifelse(model_scope == "spatial", "spNMix", "NMix"), "\n")
cat("Habitat covariate:", habitat_covariate, "\n")
cat("Wolf covariate:", wolf_covariate, "\n")
cat("Detection model:", det_model, "\n")
cat("Model scope:", model_scope, "\n")
cat("Family:", family_name, "\n\n")
print(summary(fit))
sink()

cat("Saved fit to", rds_file, "\n")
cat("Saved summary to", summary_file, "\n")
