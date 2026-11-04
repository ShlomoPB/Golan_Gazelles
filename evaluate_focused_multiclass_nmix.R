#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(spAbundance))

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  run_tag = "2000b_25bl_10th_3ch",
  class_order = "fawn,male,female",
  structures = "ms,lf1,lf2",
  output_prefix = "publication_focused_multiclass_nmix_comparison",
  run_waic = "false"
)
for (arg in args) {
  pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(pieces) == 2) {
    arg_map[[sub("^--", "", pieces[1])]] <- pieces[2]
  }
}

structures <- trimws(strsplit(arg_map$structures, ",", fixed = TRUE)[[1]])
if (!length(structures) || any(!structures %in% c("ms", "lf1", "lf2"))) {
  stop("--structures must contain one or more of: ms, lf1, lf2")
}
class_order <- trimws(strsplit(arg_map$class_order, ",", fixed = TRUE)[[1]])
order_tag <- paste(class_order, collapse = "-")
fit_files <- setNames(
  paste0(
    "gazelle_", structures,
    "_NMix_NB_open_wolf_jackal_interaction_order-", order_tag,
    "_", arg_map$run_tag, "_fit.rds"
  ),
  structures
)
missing_files <- fit_files[!file.exists(fit_files)]
if (length(missing_files)) {
  stop("Missing candidate fit(s): ", paste(missing_files, collapse = ", "))
}
fits <- lapply(fit_files, readRDS)
run_waic <- tolower(arg_map$run_waic) %in% c("true", "t", "1", "yes")

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}

summarize_matrix <- function(x) {
  x <- as.matrix(x)
  data.frame(
    parameter = colnames(x),
    mean = colMeans(x),
    median = apply(x, 2, median),
    lower_95 = apply(x, 2, quantile, probs = 0.025),
    upper_95 = apply(x, 2, quantile, probs = 0.975),
    posterior_probability_positive = colMeans(x > 0),
    row.names = NULL,
    check.names = FALSE
  )
}

free_loading_positions <- function(fit) {
  if (is.null(fit$lambda.samples)) return(character())
  n_classes <- length(class_order)
  n_factors <- ncol(as.matrix(fit$lambda.samples)) / n_classes
  free <- integer()
  for (factor_i in seq_len(n_factors)) {
    if (factor_i < n_classes) {
      rows <- seq.int(factor_i + 1L, n_classes)
      positions <- (factor_i - 1L) * n_classes + rows
      free <- c(free, positions)
    }
  }
  free
}

free_loading_names <- function(fit) {
  colnames(as.matrix(fit$lambda.samples))[free_loading_positions(fit)]
}

diagnostic_rows <- list()
coefficient_rows <- list()
correlation_rows <- list()

for (structure in structures) {
  fit <- fits[[structure]]
  coefficient_rows[[structure]] <- cbind(
    structure = structure,
    component = "class_abundance",
    summarize_matrix(fit$beta.samples)
  )

  sample_map <- c(
    "beta.comm" = "beta.comm.samples",
    "alpha.comm" = "alpha.comm.samples",
    "tau.sq.beta" = "tau.sq.beta.samples",
    "tau.sq.alpha" = "tau.sq.alpha.samples",
    "beta" = "beta.samples",
    "alpha" = "alpha.samples",
    "kappa" = "kappa.samples"
  )
  for (component in names(sample_map)) {
    sample_object <- fit[[sample_map[[component]]]]
    if (is.null(sample_object) || is.null(fit$rhat[[component]])) next
    parameter_names <- colnames(as.matrix(sample_object))
    diagnostic_rows[[paste(structure, component)]] <- data.frame(
      structure = structure,
      component = component,
      parameter = parameter_names,
      rhat = as.numeric(fit$rhat[[component]]),
      ess = as.numeric(fit$ESS[[component]]),
      row.names = NULL
    )
  }
  if (!is.null(fit$rhat$lambda.lower.tri)) {
    diagnostic_rows[[paste(structure, "lambda")]] <- data.frame(
      structure = structure,
      component = "lambda",
      parameter = free_loading_names(fit),
      rhat = as.numeric(fit$rhat$lambda.lower.tri),
      ess = as.numeric(fit$ESS$lambda)[free_loading_positions(fit)],
      row.names = NULL
    )
  }

  if (!is.null(fit$lambda.samples)) {
    lambda <- as.matrix(fit$lambda.samples)
    n_classes <- length(class_order)
    n_factors <- ncol(lambda) / n_classes
    pair_names <- combn(class_order, 2, simplify = FALSE)
    cor_draws <- matrix(NA_real_, nrow(lambda), length(pair_names))
    for (draw_i in seq_len(nrow(lambda))) {
      lambda_i <- matrix(lambda[draw_i, ], nrow = n_classes, ncol = n_factors)
      covariance_i <- tcrossprod(lambda_i)
      correlation_i <- cov2cor(covariance_i)
      for (pair_i in seq_along(pair_names)) {
        indices <- match(pair_names[[pair_i]], class_order)
        cor_draws[draw_i, pair_i] <- correlation_i[indices[1], indices[2]]
      }
    }
    correlation_rows[[structure]] <- data.frame(
      structure = structure,
      pair = vapply(pair_names, paste, collapse = "-", FUN.VALUE = character(1)),
      mean = colMeans(cor_draws),
      median = apply(cor_draws, 2, median),
      lower_95 = apply(cor_draws, 2, quantile, probs = 0.025),
      upper_95 = apply(cor_draws, 2, quantile, probs = 0.975),
      row.names = NULL
    )
  }
}

diagnostics <- do.call(rbind, diagnostic_rows)
coefficients <- do.call(rbind, coefficient_rows)
correlations <- do.call(rbind, correlation_rows)
write_csv(diagnostics, paste0(arg_map$output_prefix, "_diagnostics.csv"))
write_csv(coefficients, paste0(arg_map$output_prefix, "_coefficients.csv"))
write_csv(correlations, paste0(arg_map$output_prefix, "_correlations.csv"))
selected_correlations <- correlations[
  correlations$structure == "lf2",
  c("pair", "mean", "lower_95", "upper_95")
]
names(selected_correlations) <- c(
  "class_pair", "posterior_mean", "lower_95", "upper_95"
)
write_csv(
  selected_correlations,
  "publication_focused_lfmsnmix_class_correlations.csv"
)

# Use the same integration ceiling for every candidate in the WAIC comparison.
n_max <- Reduce(
  pmax,
  lapply(fits, function(fit) apply(fit$N.samples, c(2, 3), max, na.rm = TRUE))
) + 10
write_csv(
  cbind(class = class_order, as.data.frame(n_max)),
  paste0(arg_map$output_prefix, "_common_nmax.csv")
)
waic_total <- NULL
if (run_waic) {
  waic_rows <- list()
  for (structure in structures) {
    waic_i <- as.data.frame(waicAbund(
      fits[[structure]],
      N.max = n_max,
      by.species = TRUE
    ))
    waic_i$class <- class_order
    waic_i$structure <- structure
    waic_rows[[structure]] <- waic_i[, c("structure", "class", "elpd", "pD", "WAIC")]
  }
  waic_by_class <- do.call(rbind, waic_rows)
  waic_total <- aggregate(cbind(elpd, pD, WAIC) ~ structure, waic_by_class, sum)
  waic_total$delta_WAIC <- waic_total$WAIC - min(waic_total$WAIC)
  waic_total <- waic_total[order(waic_total$WAIC), ]
  write_csv(waic_by_class, paste0(arg_map$output_prefix, "_waic_by_class.csv"))
  write_csv(waic_total, paste0(arg_map$output_prefix, "_waic_total.csv"))
} else {
  message("Skipping WAIC. Use --run_waic=true for the full common-N.max calculation.")
}

ppc_settings <- data.frame(
  fit_stat = c("freeman-tukey", "freeman-tukey", "chi-squared"),
  group = c(0L, 1L, 2L),
  type = c("marginal", "marginal", "marginal"),
  stringsAsFactors = FALSE
)
ppc_rows <- list()
for (structure in structures) {
  for (setting_i in seq_len(nrow(ppc_settings))) {
    setting <- ppc_settings[setting_i, ]
    ppc <- ppcAbund(
      fits[[structure]],
      fit.stat = setting$fit_stat,
      group = setting$group,
      type = setting$type
    )
    ppc_rows[[paste(structure, setting_i)]] <- data.frame(
      structure = structure,
      class = class_order,
      fit_stat = setting$fit_stat,
      group = setting$group,
      type = setting$type,
      bayesian_p = colMeans(ppc$fit.y.rep > ppc$fit.y),
      row.names = NULL
    )
  }
}
ppc_results <- do.call(rbind, ppc_rows)
write_csv(ppc_results, paste0(arg_map$output_prefix, "_ppc.csv"))

diagnostic_summary <- aggregate(
  cbind(rhat, ess) ~ structure,
  diagnostics,
  function(x) c(min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
)
diagnostic_summary <- data.frame(
  structure = diagnostic_summary$structure,
  min_rhat = diagnostic_summary$rhat[, "min"],
  max_rhat = diagnostic_summary$rhat[, "max"],
  min_ess = diagnostic_summary$ess[, "min"],
  max_ess = diagnostic_summary$ess[, "max"],
  row.names = NULL
)
write_csv(diagnostic_summary, paste0(arg_map$output_prefix, "_diagnostic_summary.csv"))

cat("Focused candidate comparison complete.\n")
if (!is.null(waic_total)) print(waic_total)
print(diagnostic_summary)
