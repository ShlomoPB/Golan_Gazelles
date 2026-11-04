#!/usr/bin/env Rscript

# Audit saved, superseded single-class N-mixture fits without refitting them.

options(stringsAsFactors = FALSE)

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

fmt_formula <- function(x) {
  x <- setdiff(x, "(Intercept)")
  if (!length(x)) return("~ 1")
  paste("~", paste(x, collapse = " + "))
}

named_diagnostic <- function(values, sample_matrix, prefix) {
  values <- as.numeric(values)
  term_names <- colnames(sample_matrix)
  if (is.null(term_names) || length(term_names) != length(values)) {
    term_names <- paste0(prefix, seq_along(values))
  }
  stats::setNames(values, paste0(prefix, ":", term_names))
}

files <- list.files(
  pattern = "^(fawn|female|male)_(sp)?NMix_.*_fit\\.rds$",
  full.names = FALSE
)

rows <- lapply(files, function(path) {
  fit <- readRDS(path)
  if (!inherits(fit, "NMix")) return(NULL)

  rhat <- c(
    named_diagnostic(fit$rhat$beta, fit$beta.samples, "abundance"),
    named_diagnostic(fit$rhat$alpha, fit$alpha.samples, "detection"),
    if (!is.null(fit$rhat$kappa)) {
      stats::setNames(as.numeric(fit$rhat$kappa), "overdispersion:kappa")
    }
  )
  ess <- c(
    named_diagnostic(fit$ESS$beta, fit$beta.samples, "abundance"),
    named_diagnostic(fit$ESS$alpha, fit$alpha.samples, "detection"),
    if (!is.null(fit$ESS$kappa)) {
      stats::setNames(as.numeric(fit$ESS$kappa), "overdispersion:kappa")
    }
  )

  failing <- names(rhat)[is.finite(rhat) & rhat > 1.1]
  very_low_ess <- names(ess)[is.finite(ess) & ess < 100]
  beta_names <- colnames(fit$beta.samples)
  alpha_names <- colnames(fit$alpha.samples)

  data.frame(
    class = sub("_.*$", "", path),
    fit_file = path,
    model_family = if (is.null(fit$kappa.samples)) "Poisson" else "negative binomial",
    abundance_design_terms = fmt_formula(beta_names),
    detection_design_terms = fmt_formula(alpha_names),
    n_samples_per_chain_before_burn = fit$n.samples,
    burn_in = fit$n.burn,
    thinning = fit$n.thin,
    chains = fit$n.chains,
    saved_posterior_draws = nrow(fit$beta.samples),
    rhat_min = min(rhat, na.rm = TRUE),
    rhat_max = max(rhat, na.rm = TRUE),
    ess_min = min(ess, na.rm = TRUE),
    ess_max = max(ess, na.rm = TRUE),
    n_parameters_rhat_over_1_1 = length(failing),
    parameters_rhat_over_1_1 = paste(failing, collapse = "; "),
    n_parameters_ess_under_100 = length(very_low_ess),
    parameters_ess_under_100 = paste(very_low_ess, collapse = "; "),
    posterior_predictive_check = "Not retained with this fit object",
    priors = "Not recoverable from saved fit object; inspect generating script",
    random_seed = "Not recoverable from saved fit object; inspect generating script",
    audit_status = if (length(failing) || length(very_low_ess)) {
      "REJECTED/DIAGNOSTIC: convergence or effective-sample-size failure"
    } else {
      "SUPERSEDED: no stored Rhat/ESS failure at audit threshold"
    },
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, Filter(Negate(is.null), rows))
out <- out[order(out$class, out$fit_file), ]
write_csv(out, "publication_expanded_nmix_rejection_diagnostics.csv")

cat("Audited", nrow(out), "saved superseded single-class N-mixture fits.\n")
