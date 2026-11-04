#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(spAbundance))

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  data = "gazelle_multispecies_NMix_data.rds",
  model = "lf1",
  class_order = "fawn,male,female",
  n_batch = "2000",
  batch_length = "25",
  n_burn = "",
  n_thin = "10",
  n_chains = "3",
  seed = "20260817",
  verbose = "true",
  output_tag = ""
)

for (arg in args) {
  pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(pieces) == 2) {
    arg_map[[sub("^--", "", pieces[1])]] <- pieces[2]
  }
}

model_name <- tolower(arg_map$model)
if (!model_name %in% c("ms", "lf1", "lf2")) {
  stop("--model must be one of: ms, lf1, lf2")
}

data_list <- readRDS(arg_map$data)
class_order <- trimws(strsplit(arg_map$class_order, ",", fixed = TRUE)[[1]])
available_classes <- dimnames(data_list$y)[[1]]
if (!setequal(class_order, available_classes) || length(class_order) != length(available_classes)) {
  stop("--class_order must contain each class exactly once: ",
       paste(available_classes, collapse = ","))
}
data_list$y <- data_list$y[class_order, , , drop = FALSE]

n_batch <- as.integer(arg_map$n_batch)
batch_length <- as.integer(arg_map$batch_length)
n_burn <- if (nzchar(arg_map$n_burn)) {
  as.integer(arg_map$n_burn)
} else {
  round(0.1 * n_batch * batch_length)
}
n_thin <- as.integer(arg_map$n_thin)
n_chains <- as.integer(arg_map$n_chains)
seed <- as.integer(arg_map$seed)
verbose <- tolower(arg_map$verbose) %in% c("true", "t", "1", "yes")

# The interaction is hierarchical: both predator main effects are retained.
abund_formula <- ~ Open.Spaces.z +
  Null.Wolf.Abundance.z * null.jackal.abundance.z
det_formula <- ~ log_days_operational.z + cattle_sightings.z

make_det_model_frame <- function(det_covs) {
  as.data.frame(lapply(det_covs, function(x) {
    if (is.matrix(x) || is.data.frame(x)) as.vector(x) else x
  }), stringsAsFactors = FALSE)
}

x_abund <- model.matrix(abund_formula, data = data_list$abund.covs)
x_det <- model.matrix(det_formula, data = make_det_model_frame(data_list$det.covs))
n_classes <- dim(data_list$y)[1]

priors <- list(
  beta.comm.normal = list(
    mean = rep(0, ncol(x_abund)),
    var = rep(100, ncol(x_abund))
  ),
  alpha.comm.normal = list(
    mean = rep(0, ncol(x_det)),
    var = rep(2.72, ncol(x_det))
  ),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
  kappa.unif = list(
    a = rep(0, n_classes),
    b = rep(10, n_classes)
  )
)

n_init <- apply(data_list$y, c(1, 2), max, na.rm = TRUE)
n_init[!is.finite(n_init)] <- 0
inits <- list(
  beta.comm = 0,
  alpha.comm = 0,
  beta = 0,
  alpha = 0,
  tau.sq.beta = 0.5,
  tau.sq.alpha = 0.5,
  kappa = rep(1, n_classes),
  N = n_init
)

tuning <- list(beta = 0.5, alpha = 0.5, kappa = 0.5)
if (model_name != "ms") {
  tuning$lambda <- 0.5
  tuning$w <- 0.5
}

order_tag <- paste(class_order, collapse = "-")
structure_tag <- model_name
run_tag <- if (nzchar(arg_map$output_tag)) {
  paste0("_", arg_map$output_tag)
} else {
  paste0("_", n_batch, "b_", batch_length, "bl_", n_thin, "th_", n_chains, "ch")
}
base_name <- paste0(
  "gazelle_", structure_tag,
  "_NMix_NB_open_wolf_jackal_interaction_order-", order_tag,
  run_tag
)

cat("Fitting focused multi-class N-mixture candidate\n")
cat("  structure:", model_name, "\n")
cat("  class order:", paste(class_order, collapse = ", "), "\n")
cat("  abundance:", deparse(abund_formula), "\n")
cat("  detection:", deparse(det_formula), "\n")
cat("  samples per chain:", n_batch * batch_length, "\n")
cat("  burn-in:", n_burn, "\n")
cat("  thinning:", n_thin, "\n")
cat("  chains:", n_chains, "\n")

common_args <- list(
  abund.formula = abund_formula,
  det.formula = det_formula,
  data = data_list,
  inits = inits,
  priors = priors,
  tuning = tuning,
  n.batch = n_batch,
  batch.length = batch_length,
  n.burn = n_burn,
  n.thin = n_thin,
  n.chains = n_chains,
  family = "NB",
  n.omp.threads = 1,
  verbose = verbose,
  n.report = max(1L, floor(n_batch / 10))
)

set.seed(seed)
fit <- if (model_name == "ms") {
  do.call(msNMix, common_args)
} else {
  do.call(
    lfMsNMix,
    c(common_args, list(n.factors = if (model_name == "lf1") 1L else 2L))
  )
}

fit_file <- paste0(base_name, "_fit.rds")
metadata_file <- paste0(base_name, "_metadata.rds")
community_file <- paste0(base_name, "_community_summary.txt")
class_file <- paste0(base_name, "_class_summary.txt")

saveRDS(fit, fit_file)
saveRDS(list(
  model = model_name,
  class_order = class_order,
  abundance_formula = abund_formula,
  detection_formula = det_formula,
  n_batch = n_batch,
  batch_length = batch_length,
  n_burn = n_burn,
  n_thin = n_thin,
  n_chains = n_chains,
  seed = seed,
  spAbundance_version = as.character(packageVersion("spAbundance"))
), metadata_file)

write_summary <- function(path, level) {
  sink(path)
  on.exit(sink(), add = TRUE)
  cat("Model structure:", model_name, "\n")
  cat("Class order:", paste(class_order, collapse = ", "), "\n")
  cat("Family: negative binomial\n")
  cat("Abundance formula:", paste(deparse(abund_formula), collapse = " "), "\n")
  cat("Detection formula:", paste(deparse(det_formula), collapse = " "), "\n\n")
  print(summary(fit, level = level))
}

write_summary(community_file, "community")
write_summary(class_file, "species")

cat("Saved fit:", fit_file, "\n")
cat("Saved metadata:", metadata_file, "\n")
cat("Saved summaries:", community_file, "and", class_file, "\n")
