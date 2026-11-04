library(spAbundance)

# Load both models
fit_poisson <- readRDS("gazelle_msNMix_fit.rds")
fit_nb <- readRDS("gazelle_msNMix_fit_nb.rds")

# Compute WAIC
cat("Computing WAIC for Poisson model...\n")
waic_poisson <- waicAbund(fit_poisson, N.max = 150)

cat("\nComputing WAIC for Negative Binomial model...\n")
waic_nb <- waicAbund(fit_nb, N.max = 150)

# Extract key WAIC values
waic_poisson_val <- waic_poisson$waic
waic_nb_val <- waic_nb$waic
waic_diff <- waic_poisson_val - waic_nb_val

cat("\n================================================\n")
cat("WAIC COMPARISON: POISSON vs NEGATIVE BINOMIAL\n")
cat("================================================\n\n")

cat("Poisson Model:\n")
cat("  WAIC:", round(waic_poisson_val, 2), "\n")
cat("  elpd_waic:", round(waic_poisson$elpd.waic, 2), "\n")
cat("  p_waic:", round(waic_poisson$p.waic, 2), "\n")

cat("\nNegative Binomial Model:\n")
cat("  WAIC:", round(waic_nb_val, 2), "\n")
cat("  elpd_waic:", round(waic_nb$elpd.waic, 2), "\n")
cat("  p_waic:", round(waic_nb$p.waic, 2), "\n")

cat("\nDifference (Poisson - NB):\n")
cat("  WAIC difference:", round(waic_diff, 2), "\n")

if (waic_diff > 0) {
  cat("  --> NB model is BETTER (lower WAIC by", round(waic_diff, 2), ")\n")
} else {
  cat("  --> Poisson model is BETTER (lower WAIC by", round(-waic_diff, 2), ")\n")
}

cat("\nInterpretation:\n")
if (abs(waic_diff) < 4) {
  cat("  Models have similar predictive performance (difference < 4 SE)\n")
} else if (waic_diff > 4) {
  cat("  NB model is substantially better (Δ WAIC > 4 SE)\n")
} else {
  cat("  Poisson model is substantially better (Δ WAIC > 4 SE)\n")
}

cat("\n================================================\n")