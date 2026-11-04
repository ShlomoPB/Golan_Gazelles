suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

fit_file <- "gazelle_lfMsNMix_2f_NB_shared_env_1200b_25bl_5th_3ch_fit.rds"
output_prefix <- "publication_multispecies_shared_env_audit"
fit <- readRDS(fit_file)

summarize_draws <- function(x) {
  tibble(
    mean = mean(x),
    median = median(x),
    conf_low = unname(quantile(x, 0.025)),
    conf_high = unname(quantile(x, 0.975))
  )
}

manual_ppc <- function(fit, stride = 50L) {
  set.seed(20260804)
  alpha_samples <- as.matrix(fit$alpha.samples)
  keep <- seq.int(1L, nrow(alpha_samples), by = stride)
  if (length(keep) > 600L) {
    keep <- keep[round(seq(1, length(keep), length.out = 600L))]
  }

  sp_names <- fit$sp.names
  n_species <- length(sp_names)
  n_coef <- length(fit$x.p.names)
  y <- fit$y
  obs_mask <- !is.na(y[1, , ])
  y_grouped <- apply(y, c(1, 2), sum, na.rm = TRUE)
  n_sites <- dim(y)[2]
  n_reps <- dim(y)[3]
  fit_y <- matrix(NA_real_, length(keep), n_species)
  fit_y_rep <- matrix(NA_real_, length(keep), n_species)

  for (ii in seq_along(keep)) {
    samp <- keep[ii]
    for (s in seq_len(n_species)) {
      cols <- seq.int(from = s, by = n_species, length.out = n_coef)
      eta <- as.numeric(fit$X.p %*% alpha_samples[samp, cols])
      p_mat <- matrix(NA_real_, n_sites, n_reps)
      p_mat[obs_mask] <- stats::plogis(eta)
      mu_vec <- fit$mu.samples[samp, s, ]
      expected_grouped <- rowSums(sweep(p_mat, 1, mu_vec, `*`), na.rm = TRUE)
      n_rep <- stats::rnbinom(n_sites, mu = mu_vec, size = fit$kappa.samples[samp, s])
      y_rep_grouped <- numeric(n_sites)

      for (j in seq_len(n_sites)) {
        reps <- which(obs_mask[j, ])
        if (length(reps)) {
          y_rep_grouped[j] <- sum(stats::rbinom(
            length(reps),
            size = n_rep[j],
            prob = p_mat[j, reps]
          ))
        }
      }

      fit_y[ii, s] <- sum((sqrt(y_grouped[s, ]) - sqrt(expected_grouped))^2)
      fit_y_rep[ii, s] <- sum((sqrt(y_rep_grouped) - sqrt(expected_grouped))^2)
    }
  }

  tibble(
    class = sp_names,
    fit_observed = colMeans(fit_y),
    fit_replicated = colMeans(fit_y_rep),
    bayesian_p = colMeans(fit_y_rep > fit_y),
    n_posterior_draws = length(keep)
  )
}

lambda <- as.matrix(fit$lambda.samples)
n_species <- length(fit$sp.names)
n_factors <- fit$q
pair_index <- combn(seq_len(n_species), 2)
correlation_draws <- matrix(NA_real_, nrow(lambda), ncol(pair_index))

for (i in seq_len(nrow(lambda))) {
  loading_i <- matrix(lambda[i, ], nrow = n_species, ncol = n_factors)
  correlation_i <- stats::cov2cor(loading_i %*% t(loading_i))
  correlation_draws[i, ] <- correlation_i[t(pair_index)]
}

pair_names <- apply(pair_index, 2, function(idx) {
  paste(fit$sp.names[idx], collapse = "-")
})

correlation_df <- bind_rows(lapply(seq_along(pair_names), function(i) {
  summarize_draws(correlation_draws[, i]) %>% mutate(pair = pair_names[i], .before = 1)
}))

free_loading <- apply(lambda, 2, stats::sd) > 1e-12
loading_names <- colnames(lambda)[free_loading]
loading_df <- bind_rows(lapply(seq_along(loading_names), function(i) {
  parameter <- loading_names[i]
  summarize_draws(lambda[, parameter]) %>%
    mutate(
      parameter = parameter,
      rhat = as.numeric(fit$rhat$lambda[i]),
      ess = as.numeric(fit$ESS$lambda[parameter]),
      .before = 1
    )
}))

beta <- as.matrix(fit$beta.samples)
beta_df <- bind_rows(lapply(seq_len(ncol(beta)), function(i) {
  summarize_draws(beta[, i]) %>%
    mutate(
      parameter = colnames(beta)[i],
      rhat = as.numeric(fit$rhat$beta[i]),
      ess = as.numeric(fit$ESS$beta[i]),
      .before = 1
    )
}))

alpha <- as.matrix(fit$alpha.samples)
alpha_df <- bind_rows(lapply(seq_len(ncol(alpha)), function(i) {
  summarize_draws(alpha[, i]) %>%
    mutate(
      parameter = colnames(alpha)[i],
      rhat = as.numeric(fit$rhat$alpha[i]),
      ess = as.numeric(fit$ESS$alpha[i]),
      .before = 1
    )
}))

ppc_df <- manual_ppc(fit)

write_csv(correlation_df, paste0(output_prefix, "_correlations.csv"))
write_csv(loading_df, paste0(output_prefix, "_loadings.csv"))
write_csv(beta_df, paste0(output_prefix, "_abundance_diagnostics.csv"))
write_csv(alpha_df, paste0(output_prefix, "_detection_diagnostics.csv"))
write_csv(ppc_df, paste0(output_prefix, "_ppc.csv"))

ecological_beta <- beta_df %>% filter(!grepl("Intercept", parameter, fixed = TRUE))
intercept_beta <- beta_df %>% filter(grepl("Intercept", parameter, fixed = TRUE))
intercept_alpha <- alpha_df %>% filter(grepl("Intercept", parameter, fixed = TRUE))

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
report <- c(
  "# Multispecies shared-environment lfMsNMix audit",
  "",
  paste0("Fit: `", fit_file, "`."),
  "",
  "## Residual latent-factor correlations",
  "",
  "Correlations are derived within each posterior draw from the standardized latent-factor covariance matrix `Lambda %*% t(Lambda)`.",
  "",
  "| Pair | Mean | 95% CrI |",
  "|---|---:|---:|"
)

for (i in seq_len(nrow(correlation_df))) {
  report <- c(report, paste0(
    "| `", correlation_df$pair[i], "` | `", fmt(correlation_df$mean[i], 2),
    "` | `[", fmt(correlation_df$conf_low[i], 2), ", ",
    fmt(correlation_df$conf_high[i], 2), "]` |"
  ))
}

report <- c(
  report,
  "",
  "## Convergence",
  "",
  paste0(
    "Free loading Rhat range: `", fmt(min(loading_df$rhat)), "-",
    fmt(max(loading_df$rhat)), "`; ecological abundance-slope Rhat range: `",
    fmt(min(ecological_beta$rhat)), "-", fmt(max(ecological_beta$rhat)), "`."
  ),
  paste0(
    "Abundance-intercept Rhat range: `", fmt(min(intercept_beta$rhat)), "-",
    fmt(max(intercept_beta$rhat)), "`; detection-intercept Rhat range: `",
    fmt(min(intercept_alpha$rhat)), "-", fmt(max(intercept_alpha$rhat)),
    "`. Absolute abundance and intercepts are not used for inference."
  ),
  "",
  "## Posterior predictive checks",
  "",
  "| Class | Bayesian p |",
  "|---|---:|"
)

for (i in seq_len(nrow(ppc_df))) {
  report <- c(report, paste0(
    "| `", ppc_df$class[i], "` | `", fmt(ppc_df$bayesian_p[i]), "` |"
  ))
}

writeLines(report, paste0(output_prefix, "_results.md"))
cat("Saved multispecies audit outputs with prefix:", output_prefix, "\n")
