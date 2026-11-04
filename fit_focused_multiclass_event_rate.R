#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(spAbundance))

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  data = "gazelle_multispecies_NMix_data.rds",
  model = "lf2",
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
if (!model_name %in% c("ms", "lf2")) {
  stop("--model must be one of: ms, lf2")
}

source_data <- readRDS(arg_map$data)
class_order <- trimws(strsplit(arg_map$class_order, ",", fixed = TRUE)[[1]])
available_classes <- dimnames(source_data$y)[[1]]
if (!setequal(class_order, available_classes) || length(class_order) != length(available_classes)) {
  stop("--class_order must contain each class exactly once: ",
       paste(available_classes, collapse = ","))
}

y <- source_data$y[class_order, , , drop = FALSE]
effort <- as.matrix(source_data$det.covs$days_operational)
if (!identical(dim(effort), dim(y)[2:3])) {
  stop("days_operational must be a site-by-replicate matrix")
}
missing_y <- apply(is.na(y), c(2, 3), all)
effort[missing_y] <- NA_real_
if (any(!is.na(effort) & (!is.finite(effort) | effort <= 0))) {
  stop("All observed camera efforts must be finite and positive")
}

covs <- list(
  Open.Spaces.z = source_data$abund.covs$Open.Spaces.z,
  Null.Wolf.Abundance.z = source_data$abund.covs$Null.Wolf.Abundance.z,
  null.jackal.abundance.z = source_data$abund.covs$null.jackal.abundance.z
)
data_list <- list(
  y = y,
  covs = covs,
  coords = source_data$coords,
  offset = effort
)

formula <- ~ Open.Spaces.z +
  Null.Wolf.Abundance.z * null.jackal.abundance.z
x <- model.matrix(formula, data = source_data$abund.covs)
n_classes <- dim(y)[1]

priors <- list(
  beta.comm.normal = list(mean = rep(0, ncol(x)), var = rep(100, ncol(x))),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  kappa.unif = list(a = rep(0, n_classes), b = rep(10, n_classes))
)
inits <- list(
  beta.comm = 0,
  beta = 0,
  tau.sq.beta = 0.5,
  kappa = rep(1, n_classes)
)
tuning <- list(beta = 0.5, kappa = 0.5)
if (model_name == "lf2") {
  tuning$lambda <- 0.5
  tuning$w <- 0.5
}

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

run_tag <- if (nzchar(arg_map$output_tag)) {
  paste0("_", arg_map$output_tag)
} else {
  paste0("_", n_batch, "b_", batch_length, "bl_", n_thin, "th_", n_chains, "ch")
}
base_name <- paste0(
  "gazelle_", model_name,
  "_Abund_NB_event_rate_open_wolf_jackal_interaction_order-",
  paste(class_order, collapse = "-"), run_tag
)

common_args <- list(
  formula = formula,
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
  n.report = max(1L, floor(n_batch / 10)),
  save.fitted = TRUE
)

cat("Fitting focused multivariate independent-event rate model\n")
cat("  structure:", model_name, "\n")
cat("  class order:", paste(class_order, collapse = ", "), "\n")
cat("  formula:", paste(deparse(formula), collapse = " "), "\n")
cat("  effort offset: raw days operational\n")

set.seed(seed)
fit <- if (model_name == "ms") {
  do.call(msAbund, common_args)
} else {
  do.call(lfMsAbund, c(common_args, list(n.factors = 2L)))
}

fit_file <- paste0(base_name, "_fit.rds")
summary_file <- paste0(base_name, "_summary.txt")
metadata_file <- paste0(base_name, "_metadata.rds")
saveRDS(fit, fit_file)
saveRDS(list(
  model = model_name,
  class_order = class_order,
  formula = formula,
  offset = "raw days_operational",
  n_batch = n_batch,
  batch_length = batch_length,
  n_burn = n_burn,
  n_thin = n_thin,
  n_chains = n_chains,
  seed = seed,
  spAbundance_version = as.character(packageVersion("spAbundance"))
), metadata_file)

sink(summary_file)
cat("Response: independent detection-event counts\n")
cat("Family: negative binomial\n")
cat("Model structure:", model_name, "\n")
cat("Class order:", paste(class_order, collapse = ", "), "\n")
cat("Formula:", paste(deparse(formula), collapse = " "), "\n")
cat("Offset: raw camera days operational\n\n")
print(summary(fit))
sink()

cat("Saved fit:", fit_file, "\n")
cat("Saved summary:", summary_file, "\n")
cat("Saved metadata:", metadata_file, "\n")
