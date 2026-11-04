#!/usr/bin/env Rscript

# Inspect saved model objects without refitting. Detailed object structure is
# written to text so formula/settings discrepancies remain auditable.

options(stringsAsFactors = FALSE)

out <- "publication_model_object_structure_audit.txt"
con <- file(out, open = "wt")
sink(con)
on.exit({
  sink()
  close(con)
}, add = TRUE)

section <- function(title) {
  cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  cat(title, "\n")
  cat(paste(rep("=", 78), collapse = ""), "\n", sep = "")
}

show_file <- function(path) {
  cat("File:", path, "\n")
  if (!file.exists(path)) {
    cat("STATUS: MISSING\n")
    return(NULL)
  }
  info <- file.info(path)
  cat("Bytes:", info$size, "\n")
  cat("Modified:", format(info$mtime, "%Y-%m-%d %H:%M:%S %z"), "\n")
  readRDS(path)
}

section("Primary separate-class spAbundance NMix objects")
primary_file <- "publication_wolf_jackal_interaction_nmix_spabundance_fits.rds"
primary <- show_file(primary_file)
if (!is.null(primary)) {
  cat("Top-level class:", paste(class(primary), collapse = ", "), "\n")
  cat("Top-level names:", paste(names(primary), collapse = ", "), "\n")
  for (name_i in names(primary)) {
    fit <- primary[[name_i]]
    cat("\n--", name_i, "--\n")
    cat("Class:", paste(class(fit), collapse = ", "), "\n")
    cat("Names:", paste(names(fit), collapse = ", "), "\n")
    if (!is.null(fit$call)) {
      cat("Call:\n")
      print(fit$call)
    }
    for (slot_i in c("n.samples", "n.chains", "family")) {
      if (!is.null(fit[[slot_i]])) cat(slot_i, ": ", paste(fit[[slot_i]], collapse = ", "), "\n", sep = "")
    }
    sample_slots <- grep("samples$", names(fit), value = TRUE)
    for (slot_i in sample_slots) {
      cat(slot_i, " dimensions: ", paste(dim(fit[[slot_i]]), collapse = " x "), "\n", sep = "")
    }
  }
}
cat("\nAuthoritative script formula: abundance ~ Null.Wolf.Abundance.z * null.jackal.abundance.z\n")
cat("Authoritative script formula: detection ~ log_days_operational.z + cattle_sightings.z\n")
cat("Script settings: seed 20260623; 2000 batches x 25; burn 5000; thin 10; 3 chains; NB\n")
cat("Unmarked fit objects are not saved by the authoritative script; only tidy outputs are retained.\n")

section("Secondary focused two-factor lfMsNMix object")
lf_fit_file <- paste0(
  "gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-",
  "fawn-male-female_2000b_25bl_10th_3ch_fit.rds"
)
lf_metadata_file <- sub("_fit\\.rds$", "_metadata.rds", lf_fit_file)
lf_metadata <- show_file(lf_metadata_file)
if (!is.null(lf_metadata)) {
  print(lf_metadata)
}
lf_fit <- show_file(lf_fit_file)
if (!is.null(lf_fit)) {
  cat("Class:", paste(class(lf_fit), collapse = ", "), "\n")
  cat("Names:", paste(names(lf_fit), collapse = ", "), "\n")
  if (!is.null(lf_fit$call)) print(lf_fit$call)
  sample_slots <- grep("samples$", names(lf_fit), value = TRUE)
  for (slot_i in sample_slots) {
    cat(slot_i, " dimensions: ", paste(dim(lf_fit[[slot_i]]), collapse = " x "), "\n", sep = "")
  }
}

section("Primary beta-binomial model objects")
beta_file <- "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
beta <- show_file(beta_file)
if (!is.null(beta)) {
  cat("Top-level names:", paste(names(beta), collapse = ", "), "\n")
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    for (name_i in names(beta)) {
      fit <- beta[[name_i]]
      cat("\n--", name_i, "--\n")
      cat("Class:", paste(class(fit), collapse = ", "), "\n")
      if (inherits(fit, "glmmTMB")) {
        cat("Formula:\n")
        print(stats::formula(fit))
        cat("Family:", family(fit)$family, "\n")
        cat("Link:", family(fit)$link, "\n")
        cat("Observations:", stats::nobs(fit), "\n")
        cat("AIC:", stats::AIC(fit), "\n")
        cat("Optimizer convergence:", fit$fit$convergence, "\n")
        cat("Positive-definite Hessian:", fit$sdr$pdHess, "\n")
      }
    }
  } else {
    cat("glmmTMB unavailable; formulas were not extracted.\n")
  }
}

section("Authoritative derived output presence")
derived <- c(
  "publication_wolf_jackal_interaction_nmix_coefficients.csv",
  "publication_wolf_jackal_interaction_nmix_spabundance_diagnostics.csv",
  "publication_spabundance_wolf_jackal_relative_conditional_effect_ratios.csv",
  "publication_focused_lfmsnmix_class_correlations.csv",
  "publication_focused_lfmsnmix_diagnostics.csv",
  "publication_focused_lfmsnmix_ppc.csv",
  "publication_fawn_female_beta_binomial_coefficients.csv",
  "publication_fawn_female_beta_binomial_interaction_sensitivity_model_comparison.csv",
  "publication_predator_context_diel_pairwise_overlap_median_results.csv"
)
for (path in derived) {
  cat(if (file.exists(path)) "PRESENT" else "MISSING", " ", path, "\n", sep = "")
}

cat("\nAudit completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n")
