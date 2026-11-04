library(spAbundance)

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  data = "gazelle_multispecies_NMix_data.rds",
  family = "NB",
  n_factors = "2",
  n_batch = "400",
  batch_length = "25",
  n_burn = "",
  n_thin = "5",
  n_chains = "3",
  seed = "20260610"
)

for (arg in args) {
  pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(pieces) == 2) {
    key <- sub("^--", "", pieces[1])
    arg_map[[key]] <- pieces[2]
  }
}

data_file <- arg_map$data
family_name <- arg_map$family
n_factors <- as.integer(arg_map$n_factors)
n_batch <- as.integer(arg_map$n_batch)
batch_length <- as.integer(arg_map$batch_length)
n_burn <- if (nzchar(arg_map$n_burn)) as.integer(arg_map$n_burn) else NULL
n_thin <- as.integer(arg_map$n_thin)
n_chains <- as.integer(arg_map$n_chains)
seed <- as.integer(arg_map$seed)

data_list <- readRDS(data_file)

abund_formula <- ~ Nature.Reserves.z + Open.Spaces.z + Built.and.Roads.z +
  null.jackal.abundance.z + Null.Wolf.Abundance.z
det_formula <- ~ log_days_operational.z + cattle_sightings.z

X_abund <- model.matrix(abund_formula, data = data_list$abund.covs)

make_det_model_frame <- function(det_covs) {
  out <- lapply(det_covs, function(x) {
    if (is.matrix(x) || is.data.frame(x)) {
      as.vector(x)
    } else {
      x
    }
  })
  as.data.frame(out, stringsAsFactors = FALSE)
}

X_det <- model.matrix(det_formula, data = make_det_model_frame(data_list$det.covs))
n_species <- dim(data_list$y)[1]

priors <- list(
  beta.comm.normal = list(mean = rep(0, ncol(X_abund)), var = rep(100, ncol(X_abund))),
  alpha.comm.normal = list(mean = rep(0, ncol(X_det)), var = rep(2.72, ncol(X_det))),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig = list(a = 0.1, b = 0.1)
)
if (family_name == "NB") {
  priors$kappa.unif <- list(a = rep(0, n_species), b = rep(10, n_species))
}

N_init <- apply(data_list$y, c(1, 2), max, na.rm = TRUE)
N_init[!is.finite(N_init)] <- 0

inits <- list(
  beta.comm = 0,
  alpha.comm = 0,
  beta = 0,
  alpha = 0,
  tau.sq.beta = 0.5,
  tau.sq.alpha = 0.5,
  N = N_init
)
if (family_name == "NB") {
  inits$kappa <- rep(1, n_species)
}

tuning <- list(beta = 0.5, alpha = 0.5, lambda = 0.5, w = 0.5)
if (family_name == "NB") {
  tuning$kappa <- 0.5
}

set.seed(seed)

base_name <- paste0(
  "gazelle_lfMsNMix_",
  n_factors,
  "f_",
  family_name,
  "_shared_env_",
  n_batch,
  "b_",
  batch_length,
  "bl_",
  n_thin,
  "th_",
  n_chains,
  "ch"
)
rds_file <- paste0(base_name, "_fit.rds")
community_summary_file <- paste0(base_name, "_community_summary.txt")
species_summary_file <- paste0(base_name, "_species_summary.txt")

cat("Fitting multi-class lfMsNMix model\n")
cat("  classes:", paste(dimnames(data_list$y)[[1]], collapse = ", "), "\n")
cat("  n_factors:", n_factors, "\n")
cat("  family:", family_name, "\n")

fit <- lfMsNMix(
  abund.formula = abund_formula,
  det.formula = det_formula,
  data = data_list,
  inits = inits,
  priors = priors,
  tuning = tuning,
  n.factors = n_factors,
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

saveRDS(fit, rds_file)

sink(community_summary_file)
cat("Model: lfMsNMix\n")
cat("Family:", family_name, "\n")
cat("n_factors:", n_factors, "\n")
cat("Abundance formula: ~ Nature.Reserves.z + Open.Spaces.z + Built.and.Roads.z + null.jackal.abundance.z + Null.Wolf.Abundance.z\n")
cat("Detection formula: ~ log_days_operational.z + cattle_sightings.z\n\n")
print(summary(fit, level = "community"))
sink()

sink(species_summary_file)
cat("Model: lfMsNMix\n")
cat("Family:", family_name, "\n")
cat("n_factors:", n_factors, "\n")
cat("Abundance formula: ~ Nature.Reserves.z + Open.Spaces.z + Built.and.Roads.z + null.jackal.abundance.z + Null.Wolf.Abundance.z\n")
cat("Detection formula: ~ log_days_operational.z + cattle_sightings.z\n\n")
print(summary(fit))
sink()

cat("Saved fit to", rds_file, "\n")
cat("Saved community summary to", community_summary_file, "\n")
cat("Saved species summary to", species_summary_file, "\n")
