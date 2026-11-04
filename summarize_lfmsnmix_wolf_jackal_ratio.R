#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(coda)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(spAbundance)
  library(tidyr)
})

parse_args <- function(args) {
  cfg <- list(mode = "full", ppc_stride = "50")
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) == 2) cfg[[pieces[1]]] <- pieces[2]
  }
  cfg
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
if (!cfg$mode %in% c("smoke", "full")) stop("--mode must be smoke or full")

fit_tag <- if (cfg$mode == "smoke") "smoke" else "2000b_25bl_5th_3ch"
output_prefix <- if (cfg$mode == "smoke") {
  "smoke_fawn_female_posterior_ratio_wolf_jackal"
} else {
  "publication_fawn_female_posterior_ratio_wolf_jackal"
}

fit_files <- c(
  "Overall" = paste0(
    "gazelle_lfMsNMix_2f_NB_wolf_jackal_interaction_overall_",
    fit_tag, "_fit.rds"
  ),
  "Early (hotter)" = paste0(
    "gazelle_lfMsNMix_2f_NB_wolf_jackal_interaction_early_hotter_",
    fit_tag, "_fit.rds"
  ),
  "Late (mild-cold)" = paste0(
    "gazelle_lfMsNMix_2f_NB_wolf_jackal_interaction_late_mild_cold_",
    fit_tag, "_fit.rds"
  )
)
missing_files <- fit_files[!file.exists(fit_files)]
if (length(missing_files)) {
  stop("Missing fitted models: ", paste(missing_files, collapse = ", "))
}

summarize_draws <- function(x) {
  x <- x[is.finite(x)]
  c(
    mean = mean(x),
    median = stats::median(x),
    lower = unname(stats::quantile(x, 0.025)),
    upper = unname(stats::quantile(x, 0.975))
  )
}

parameter_diagnostics <- function(fit, period) {
  bind_rows(
    tibble(
      parameter_group = "abundance_species",
      parameter = colnames(fit$beta.samples),
      rhat = as.vector(fit$rhat$beta),
      ess = as.vector(fit$ESS$beta)
    ),
    tibble(
      parameter_group = "detection_species",
      parameter = colnames(fit$alpha.samples),
      rhat = as.vector(fit$rhat$alpha),
      ess = as.vector(fit$ESS$alpha)
    ),
    tibble(
      parameter_group = "negative_binomial",
      parameter = colnames(fit$kappa.samples),
      rhat = as.vector(fit$rhat$kappa),
      ess = as.vector(fit$ESS$kappa)
    )
  ) %>%
    mutate(period = period, .before = 1)
}

# Site-grouped Freeman-Tukey PPC with species columns selected according to the
# term-major ordering used in alpha.samples.
manual_ppc <- function(fit, period, stride = 50L) {
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
      N_rep <- stats::rnbinom(
        n_sites,
        mu = mu_vec,
        size = fit$kappa.samples[samp, s]
      )
      y_rep_grouped <- numeric(n_sites)
      for (j in seq_len(n_sites)) {
        reps <- which(obs_mask[j, ])
        if (length(reps)) {
          y_rep_grouped[j] <- sum(stats::rbinom(
            length(reps),
            size = N_rep[j],
            prob = p_mat[j, reps]
          ))
        }
      }
      fit_y[ii, s] <- sum((sqrt(y_grouped[s, ]) - sqrt(expected_grouped))^2)
      fit_y_rep[ii, s] <- sum((sqrt(y_rep_grouped) - sqrt(expected_grouped))^2)
    }
  }

  tibble(
    period = period,
    class = sp_names,
    fit_observed = colMeans(fit_y),
    fit_replicated = colMeans(fit_y_rep),
    bayesian_p = colMeans(fit_y_rep > fit_y),
    n_posterior_draws = length(keep)
  )
}

extract_contrasts <- function(fit, period) {
  beta <- as.matrix(fit$beta.samples)
  terms <- c(
    wolf = "Null.Wolf.Abundance.z",
    jackal = "null.jackal.abundance.z",
    wolf_x_jackal = "Null.Wolf.Abundance.z:null.jackal.abundance.z"
  )
  bind_rows(lapply(names(terms), function(label) {
    term <- terms[[label]]
    fawn_col <- paste0(term, "-fawn")
    female_col <- paste0(term, "-female")
    if (!all(c(fawn_col, female_col) %in% colnames(beta))) {
      stop("Missing coefficient columns for contrast: ", label)
    }
    delta <- beta[, fawn_col] - beta[, female_col]
    sm <- summarize_draws(delta)
    tibble(
      period = period,
      contrast = label,
      mean = sm[["mean"]],
      median = sm[["median"]],
      lower = sm[["lower"]],
      upper = sm[["upper"]],
      probability_negative = mean(delta < 0),
      probability_positive = mean(delta > 0),
      probability_abs_below_0.1 = mean(abs(delta) < 0.1)
    )
  }))
}

extract_species_coefficients <- function(fit, period) {
  beta <- as.matrix(fit$beta.samples)
  keep <- grepl(
    "Null.Wolf.Abundance.z|null.jackal.abundance.z",
    colnames(beta)
  ) & grepl("-(fawn|female)$", colnames(beta))
  bind_rows(lapply(which(keep), function(i) {
    sm <- summarize_draws(beta[, i])
    tibble(
      period = period,
      parameter = colnames(beta)[i],
      mean = sm[["mean"]],
      median = sm[["median"]],
      lower = sm[["lower"]],
      upper = sm[["upper"]],
      probability_negative = mean(beta[, i] < 0),
      probability_positive = mean(beta[, i] > 0)
    )
  }))
}

ratio_chain_diagnostics <- function(fit, period, X_grid, prediction_grid) {
  n_species <- length(fit$sp.names)
  fawn_i <- match("fawn", fit$sp.names)
  female_i <- match("female", fit$sp.names)
  beta_all <- as.matrix(fit$beta.samples)
  chain_rows <- split(
    seq_len(nrow(beta_all)),
    rep(seq_len(fit$n.chains), each = fit$n.post)
  )
  chain_log_ratios <- lapply(chain_rows, function(rows) {
    beta <- beta_all[rows, , drop = FALSE]
    fawn_cols <- seq.int(from = fawn_i, by = n_species, length.out = ncol(X_grid))
    female_cols <- seq.int(from = female_i, by = n_species, length.out = ncol(X_grid))
    (beta[, fawn_cols, drop = FALSE] - beta[, female_cols, drop = FALSE]) %*% t(X_grid)
  })
  bind_rows(lapply(seq_len(nrow(prediction_grid)), function(i) {
    chains_i <- coda::mcmc.list(lapply(chain_log_ratios, function(x) coda::mcmc(x[, i])))
    tibble(
      period = period,
      wolf_z = prediction_grid$Null.Wolf.Abundance.z[i],
      jackal_z = prediction_grid$null.jackal.abundance.z[i],
      ratio_rhat = unname(coda::gelman.diag(
        chains_i,
        autoburnin = FALSE,
        multivariate = FALSE
      )$psrf[1, "Point est."]),
      ratio_ess = unname(coda::effectiveSize(chains_i)[1])
    )
  }))
}

prediction_grid <- expand.grid(
  Null.Wolf.Abundance.z = c(-1, 0, 1),
  null.jackal.abundance.z = c(-1, 0, 1),
  KEEP.OUT.ATTRS = FALSE
)
X_grid <- model.matrix(
  ~ Null.Wolf.Abundance.z * null.jackal.abundance.z,
  data = prediction_grid
)

ratio_summaries <- list()
ratio_draws <- list()
diagnostics <- list()
ppc_results <- list()
contrast_results <- list()
coefficient_results <- list()
ratio_diagnostic_results <- list()

for (period in names(fit_files)) {
  cat("Summarizing", period, "\n")
  fit <- readRDS(fit_files[[period]])
  if (!identical(colnames(X_grid), colnames(fit$X))) {
    stop("Prediction design matrix does not match fitted model for ", period)
  }
  pred <- predict(
    fit,
    X.0 = X_grid,
    type = "abundance",
    ignore.RE = TRUE,
    include.w = FALSE
  )
  mu <- pred$mu.0.samples
  fawn_i <- match("fawn", fit$sp.names)
  female_i <- match("female", fit$sp.names)
  ratio <- 100 * mu[, fawn_i, ] / mu[, female_i, ]
  colnames(ratio) <- paste0(
    "wolf_", prediction_grid$Null.Wolf.Abundance.z,
    "_jackal_", prediction_grid$null.jackal.abundance.z
  )
  ratio_draws[[period]] <- ratio

  ratio_summaries[[period]] <- bind_rows(lapply(seq_len(nrow(prediction_grid)), function(i) {
    sm <- summarize_draws(ratio[, i])
    tibble(
      period = period,
      wolf_z = prediction_grid$Null.Wolf.Abundance.z[i],
      jackal_z = prediction_grid$null.jackal.abundance.z[i],
      mean = sm[["mean"]],
      median = sm[["median"]],
      lower = sm[["lower"]],
      upper = sm[["upper"]],
      probability_below_10 = mean(ratio[, i] < 10),
      probability_above_20 = mean(ratio[, i] > 20),
      n_draws = sum(is.finite(ratio[, i]))
    )
  }))
  diagnostics[[period]] <- parameter_diagnostics(fit, period)
  ppc_results[[period]] <- manual_ppc(fit, period, as.integer(cfg$ppc_stride))
  contrast_results[[period]] <- extract_contrasts(fit, period)
  coefficient_results[[period]] <- extract_species_coefficients(fit, period)
  ratio_diagnostic_results[[period]] <- ratio_chain_diagnostics(
    fit,
    period,
    X_grid,
    prediction_grid
  )
  rm(fit, pred, mu)
  invisible(gc())
}

ratio_summary <- bind_rows(ratio_summaries) %>%
  mutate(
    period = factor(period, levels = names(fit_files)),
    wolf_level = factor(wolf_z, levels = c(-1, 0, 1)),
    jackal_level = factor(jackal_z, levels = c(-1, 0, 1)),
    label = sprintf("%.1f\n[%.1f, %.1f]", median, lower, upper)
  )
diagnostic_df <- bind_rows(diagnostics)
ppc_df <- bind_rows(ppc_results)
contrast_df <- bind_rows(contrast_results)
coefficient_df <- bind_rows(coefficient_results)
ratio_diagnostic_df <- bind_rows(ratio_diagnostic_results)
ratio_summary <- ratio_summary %>%
  left_join(ratio_diagnostic_df, by = c("period", "wolf_z", "jackal_z"))

write_csv(ratio_summary, paste0(output_prefix, "_prediction_summary.csv"))
write_csv(diagnostic_df, paste0(output_prefix, "_diagnostics.csv"))
write_csv(ppc_df, paste0(output_prefix, "_ppc.csv"))
write_csv(contrast_df, paste0(output_prefix, "_coefficient_contrasts.csv"))
write_csv(coefficient_df, paste0(output_prefix, "_species_coefficients.csv"))
write_csv(ratio_diagnostic_df, paste0(output_prefix, "_ratio_diagnostics.csv"))
saveRDS(
  list(grid = prediction_grid, ratio_draws = ratio_draws),
  paste0(output_prefix, "_posterior_draws.rds")
)

fill_upper <- max(40, max(ratio_summary$median, na.rm = TRUE))
max_ratio_rhat <- max(ratio_summary$ratio_rhat, na.rm = TRUE)
ratio_convergence_note <- if (max_ratio_rhat > 1.05) {
  paste0(
    "Exploratory sensitivity only: log-ratio chains did not converge (maximum Rhat = ",
    sprintf("%.2f", max_ratio_rhat), ")"
  )
} else {
  "Grid-cell log-ratio chains met the Rhat <= 1.05 criterion"
}
ratio_summary <- ratio_summary %>%
  mutate(text_color = if_else(median >= 0.55 * fill_upper, "white", "#263238"))

ratio_plot <- ggplot(
  ratio_summary,
  aes(x = jackal_level, y = wolf_level, fill = median)
) +
  geom_tile(color = "#fbfaf7", linewidth = 1.4) +
  geom_text(
    aes(label = label, color = text_color),
    size = 3.25,
    lineheight = 0.9,
    fontface = "bold"
  ) +
  facet_wrap(~period, nrow = 1) +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = c("#f7f0df", "#e7c283", "#d17443", "#7e231f"),
    limits = c(0, fill_upper),
    oob = scales::squish,
    name = "Posterior median\nfawns per 100 females"
  ) +
  coord_equal() +
  labs(
    x = "Jackal abundance (z)",
    y = "Wolf abundance (z)",
    title = "Exploratory latent-factor fawn:female relative-abundance ratio",
    subtitle = paste0(
      "Joint two-factor negative-binomial lfMsNMix predictions; latent site factors excluded\n",
      "Cells show posterior median [95% credible interval]\n",
      ratio_convergence_note
    ),
    caption = paste0(
      "Predictions use continuous wolf and jackal abundance at -1, 0, and +1 SD. ",
      "Ratios describe expected relative abundance, not survival or realized recruitment."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#283336"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#4c5a5e", lineheight = 1.05),
    plot.caption = element_text(color = "#5f696b", hjust = 0),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#e8e1d4", color = NA),
    legend.position = "right"
  )

ggsave(paste0(output_prefix, "_heatmap.png"), ratio_plot, width = 12.6, height = 5.3, dpi = 320)
ggsave(paste0(output_prefix, "_heatmap.pdf"), ratio_plot, width = 12.6, height = 5.3)

max_slope_rhat <- diagnostic_df %>%
  filter(
    parameter_group == "abundance_species",
    grepl("Null.Wolf.Abundance.z|null.jackal.abundance.z", parameter)
  ) %>%
  summarise(value = max(rhat, na.rm = TRUE)) %>%
  pull(value)
min_slope_ess <- diagnostic_df %>%
  filter(
    parameter_group == "abundance_species",
    grepl("Null.Wolf.Abundance.z|null.jackal.abundance.z", parameter)
  ) %>%
  summarise(value = min(ess, na.rm = TRUE)) %>%
  pull(value)
min_ratio_ess <- min(ratio_summary$ratio_ess, na.rm = TRUE)

fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
md <- c(
  "# Posterior fawn:female ratios across wolf-jackal contexts",
  "",
  "## Estimand",
  "",
  "Each cell is the posterior distribution of `100 * expected fawn abundance / expected female abundance`, calculated within the same posterior draw from a joint two-factor negative-binomial `lfMsNMix` model. Predictions exclude latent site factors (`include.w = FALSE`) and use continuous standardized predator covariates at `-1`, `0`, and `+1` SD.",
  "",
  paste0("Maximum abundance-slope Rhat: `", fmt(max_slope_rhat), "`; minimum abundance-slope ESS: `", round(min_slope_ess), "`."),
  "",
  paste0("Maximum grid-cell log-ratio Rhat: `", fmt(max_ratio_rhat), "`; minimum grid-cell log-ratio ESS: `", round(min_ratio_ess), "`."),
  "",
  "## Ratio predictions",
  "",
  "| Period | Wolf z | Jackal z | Median fawns/100 females | 95% CrI |",
  "|---|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(ratio_summary))) {
  md <- c(md, sprintf(
    "| `%s` | `%s` | `%s` | `%s` | `[%s, %s]` |",
    as.character(ratio_summary$period[i]), ratio_summary$wolf_z[i],
    ratio_summary$jackal_z[i], fmt(ratio_summary$median[i], 1),
    fmt(ratio_summary$lower[i], 1), fmt(ratio_summary$upper[i], 1)
  ))
}
md <- c(
  md,
  "",
  "## Differential predator coefficients",
  "",
  "Contrasts are fawn minus female coefficients on the log expected-abundance scale.",
  "",
  "| Period | Contrast | Median | 95% CrI | Pr(<0) |",
  "|---|---|---:|---:|---:|"
)
for (i in seq_len(nrow(contrast_df))) {
  md <- c(md, sprintf(
    "| `%s` | `%s` | `%s` | `[%s, %s]` | `%s` |",
    contrast_df$period[i], contrast_df$contrast[i],
    fmt(contrast_df$median[i]), fmt(contrast_df$lower[i]),
    fmt(contrast_df$upper[i]), fmt(contrast_df$probability_negative[i], 3)
  ))
}
md <- c(
  md,
  "",
  "## Posterior predictive checks",
  "",
  "Site-grouped Freeman-Tukey discrepancies; values near 0.5 indicate better calibration.",
  "",
  "| Period | Class | Bayesian p |",
  "|---|---|---:|"
)
for (i in seq_len(nrow(ppc_df))) {
  md <- c(md, sprintf(
    "| `%s` | `%s` | `%s` |",
    ppc_df$period[i], ppc_df$class[i], fmt(ppc_df$bayesian_p[i], 3)
  ))
}
md <- c(
  md,
  "",
  "## Interpretation constraint",
  "",
  "This is a model-predicted relative-abundance ratio. It is not a direct estimate of birth rate, fawn survival, or realized recruitment. The figure should be treated as primary only if convergence and posterior predictive diagnostics are adequate; otherwise it remains exploratory and the beta-binomial composition model should carry the ratio inference."
)
writeLines(md, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved posterior-ratio outputs with prefix:", output_prefix, "\n")
