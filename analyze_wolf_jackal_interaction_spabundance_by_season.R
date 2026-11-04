#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(spAbundance)
  library(tidyr)
})

output_prefix <- "publication_wolf_jackal_interaction_spabundance_by_season"
use_saved_fits <- "--use_saved_fits" %in% commandArgs(trailingOnly = TRUE)
fit_cache_file <- paste0(output_prefix, "_fits.rds")
seed <- 20260623
n_batch <- 2000
batch_length <- 25
n_burn <- 5000
n_thin <- 10
n_chains <- 3
contrast_draws <- 50000

species_info <- list(
  list(label = "Fawn", code = "fawn", data_rds = "gazelle_fawn_spNMix_data.rds"),
  list(label = "Female", code = "female", data_rds = "gazelle_female_spNMix_data.rds"),
  list(label = "Male", code = "male", data_rds = "gazelle_male_spNMix_data.rds")
)
period_levels <- c("Early", "Late")
focal_terms <- c(
  "Null.Wolf.Abundance.z",
  "null.jackal.abundance.z",
  "Null.Wolf.Abundance.z:null.jackal.abundance.z"
)
obs_cov_file <- "ObsCovs_050524.csv"

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_ci <- function(mean, lwr, upr, digits = 2) {
  paste0(
    fmt_num(mean, digits), " [",
    fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]"
  )
}

summarize_draws <- function(x) {
  c(
    mean = mean(x),
    median = median(x),
    lwr = unname(quantile(x, 0.025)),
    upr = unname(quantile(x, 0.975))
  )
}

make_det_model_frame <- function(det_covs) {
  out <- lapply(det_covs, function(x) {
    if (is.matrix(x) || is.data.frame(x)) as.vector(x) else x
  })
  as.data.frame(out, stringsAsFactors = FALSE)
}

mask_data_list <- function(data_list, keep_mask) {
  out <- data_list
  out$y[!keep_mask] <- NA
  for (nm in names(out$det.covs)) {
    if (is.matrix(out$det.covs[[nm]])) {
      out$det.covs[[nm]][!keep_mask] <- NA
    }
  }

  keep_sites <- rowSums(!is.na(out$y)) > 0
  out$y <- out$y[keep_sites, , drop = FALSE]
  out$abund.covs <- out$abund.covs[keep_sites, , drop = FALSE]
  out$coords <- out$coords[keep_sites, , drop = FALSE]
  for (nm in names(out$det.covs)) {
    if (is.matrix(out$det.covs[[nm]])) {
      out$det.covs[[nm]] <- out$det.covs[[nm]][keep_sites, , drop = FALSE]
    }
  }
  out
}

extract_diagnostics <- function(fit, species, period) {
  bind_rows(
    data.frame(
      species,
      period,
      component = "abundance",
      parameter = colnames(fit$beta.samples),
      rhat = as.numeric(fit$rhat$beta),
      ess = as.numeric(fit$ESS$beta)
    ),
    data.frame(
      species,
      period,
      component = "detection",
      parameter = colnames(fit$alpha.samples),
      rhat = as.numeric(fit$rhat$alpha),
      ess = as.numeric(fit$ESS$alpha)
    ),
    data.frame(
      species,
      period,
      component = "overdispersion",
      parameter = "kappa",
      rhat = as.numeric(fit$rhat$kappa),
      ess = as.numeric(fit$ESS$kappa)
    )
  )
}

summarize_relative_grid <- function(fit, grid) {
  beta <- as.matrix(fit$beta.samples)[, focal_terms, drop = FALSE]
  design <- cbind(
    grid$wolf_z,
    grid$jackal_z,
    grid$wolf_z * grid$jackal_z
  )
  relative_draws <- exp(beta %*% t(design))
  data.frame(
    grid,
    estimate = apply(relative_draws, 2, median),
    conf_low = apply(relative_draws, 2, quantile, probs = 0.025),
    conf_high = apply(relative_draws, 2, quantile, probs = 0.975)
  )
}

obs <- read.csv(
  obs_cov_file,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    camera_name = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    days_operational = as.numeric(`Days.Operational`)
  ) %>%
  arrange(camera_name, round) %>%
  mutate(
    midpoint_date = start_date + floor(as.numeric(end_date - start_date) / 2)
  )

active_obs <- obs %>%
  filter(
    !is.na(start_date),
    !is.na(end_date),
    is.finite(days_operational),
    days_operational > 0
  )
survey_start <- min(active_obs$start_date)
survey_end <- max(active_obs$end_date)
survey_split <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)

round_mid_lookup <- obs %>%
  filter(!is.na(midpoint_date)) %>%
  group_by(round) %>%
  summarise(midpoint_numeric = median(as.numeric(midpoint_date)), .groups = "drop")

obs <- obs %>%
  left_join(round_mid_lookup, by = "round") %>%
  mutate(
    midpoint_date = as.Date(
      ifelse(
        is.na(midpoint_date),
        midpoint_numeric,
        as.numeric(midpoint_date)
      ),
      origin = "1970-01-01"
    ),
    period = factor(
      ifelse(midpoint_date <= survey_split, "Early", "Late"),
      levels = period_levels
    )
  )

round_mask_summary <- obs %>%
  count(period, round, name = "n_camera_rounds")
write_csv(
  round_mask_summary,
  paste0(output_prefix, "_round_mask_summary.csv")
)

if (use_saved_fits && !file.exists(fit_cache_file)) {
  stop("Saved-fit mode requested, but fit cache was not found: ", fit_cache_file)
}
fits <- if (use_saved_fits) readRDS(fit_cache_file) else list()
coefficient_rows <- list()
diagnostic_rows <- list()
support_rows <- list()
relative_grid_rows <- list()
relative_curve_rows <- list()

heat_grid <- expand.grid(wolf_z = c(-1, 0, 1), jackal_z = c(-1, 0, 1))
line_grid <- expand.grid(
  wolf_z = c(-1, 0, 1),
  jackal_z = seq(-1.5, 1.5, length.out = 121)
)

for (sp in species_info) {
  data_full <- readRDS(sp$data_rds)
  site_ids <- as.integer(rownames(data_full$y))
  expected <- expand.grid(
    camera_name = site_ids,
    round = seq_len(ncol(data_full$y))
  ) %>%
    arrange(camera_name, round)
  obs_sp <- obs %>%
    filter(camera_name %in% site_ids) %>%
    arrange(camera_name, round)
  stopifnot(identical(obs_sp$camera_name, expected$camera_name))
  stopifnot(identical(obs_sp$round, expected$round))

  period_matrix <- matrix(
    as.character(obs_sp$period),
    nrow = length(site_ids),
    ncol = ncol(data_full$y),
    byrow = TRUE,
    dimnames = dimnames(data_full$y)
  )

  for (period_i in period_levels) {
    data_period <- mask_data_list(data_full, period_matrix == period_i)
    abund_formula <- ~ Null.Wolf.Abundance.z * null.jackal.abundance.z
    det_formula <- ~ log_days_operational.z + cattle_sightings.z
    x_abund <- model.matrix(abund_formula, data = data_period$abund.covs)
    x_det <- model.matrix(
      det_formula,
      data = make_det_model_frame(data_period$det.covs)
    )
    priors <- list(
      beta.normal = list(
        mean = rep(0, ncol(x_abund)),
        var = rep(100, ncol(x_abund))
      ),
      alpha.normal = list(
        mean = rep(0, ncol(x_det)),
        var = rep(2.72, ncol(x_det))
      ),
      kappa.unif = c(0, 10)
    )
    n_init <- apply(data_period$y, 1, max, na.rm = TRUE)
    n_init[!is.finite(n_init)] <- 0

    fit_key <- paste(sp$code, period_i, sep = "__")
    if (use_saved_fits) {
      fit <- fits[[fit_key]]
      if (is.null(fit)) {
        stop("Fit cache does not contain model: ", fit_key)
      }
    } else {
      set.seed(seed + match(sp$code, vapply(species_info, `[[`, "", "code")) * 10 + match(period_i, period_levels))
      fit <- NMix(
        abund.formula = abund_formula,
        det.formula = det_formula,
        data = data_period,
        inits = list(alpha = 0, beta = 0, N = n_init, kappa = 1),
        priors = priors,
        tuning = list(beta = 0.1, alpha = 0.1, kappa = 0.5),
        n.batch = n_batch,
        batch.length = batch_length,
        n.burn = n_burn,
        n.thin = n_thin,
        n.chains = n_chains,
        family = "NB",
        n.omp.threads = 1,
        verbose = FALSE
      )
      fits[[fit_key]] <- fit
    }
    diagnostic_rows[[length(diagnostic_rows) + 1]] <- extract_diagnostics(
      fit,
      sp$label,
      period_i
    )
    support_rows[[length(support_rows) + 1]] <- data.frame(
      species = sp$label,
      period = period_i,
      n_sites = nrow(data_period$y),
      n_nonmissing_cells = sum(!is.na(data_period$y)),
      total_count = sum(data_period$y, na.rm = TRUE)
    )

    beta <- as.matrix(fit$beta.samples)
    for (term_i in focal_terms) {
      term_index <- match(term_i, colnames(beta))
      summary_i <- summarize_draws(beta[, term_index])
      coefficient_rows[[length(coefficient_rows) + 1]] <- data.frame(
        species = sp$label,
        period = period_i,
        term = term_i,
        mean = summary_i["mean"],
        median = summary_i["median"],
        conf_low = summary_i["lwr"],
        conf_high = summary_i["upr"],
        rhat = as.numeric(fit$rhat$beta[term_index]),
        ess = as.numeric(fit$ESS$beta[term_index])
      )
    }

    relative_grid_rows[[length(relative_grid_rows) + 1]] <-
      summarize_relative_grid(fit, heat_grid) %>%
      mutate(species = sp$label, period = period_i)
    relative_curve_rows[[length(relative_curve_rows) + 1]] <-
      summarize_relative_grid(fit, line_grid) %>%
      mutate(species = sp$label, period = period_i)
  }
}

coefficient_df <- bind_rows(coefficient_rows)
diagnostic_df <- bind_rows(diagnostic_rows)
support_df <- bind_rows(support_rows)
relative_grid_df <- bind_rows(relative_grid_rows)
relative_curve_df <- bind_rows(relative_curve_rows)

set.seed(seed + 9000)
coefficient_contrast_rows <- list()
conditional_wolf_rows <- list()
conditional_wolf_contrast_rows <- list()

for (sp in species_info) {
  early_beta <- as.matrix(fits[[paste(sp$code, "Early", sep = "__")]]$beta.samples)
  late_beta <- as.matrix(fits[[paste(sp$code, "Late", sep = "__")]]$beta.samples)

  for (term_i in focal_terms) {
    early_draw <- sample(early_beta[, term_i], contrast_draws, replace = TRUE)
    late_draw <- sample(late_beta[, term_i], contrast_draws, replace = TRUE)
    difference <- late_draw - early_draw
    summary_i <- summarize_draws(difference)
    coefficient_contrast_rows[[length(coefficient_contrast_rows) + 1]] <- data.frame(
      species = sp$label,
      term = term_i,
      contrast = "Late minus Early",
      mean = summary_i["mean"],
      median = summary_i["median"],
      conf_low = summary_i["lwr"],
      conf_high = summary_i["upr"],
      posterior_probability_gt_zero = mean(difference > 0)
    )
  }

  for (period_i in period_levels) {
    beta <- if (period_i == "Early") early_beta else late_beta
    for (jackal_value in c(-1, 0, 1)) {
      log_slope <- beta[, focal_terms[1]] +
        jackal_value * beta[, focal_terms[3]]
      ratio <- exp(log_slope)
      log_summary <- summarize_draws(log_slope)
      ratio_summary <- summarize_draws(ratio)
      conditional_wolf_rows[[length(conditional_wolf_rows) + 1]] <- data.frame(
        species = sp$label,
        period = period_i,
        jackal_value_sd = jackal_value,
        log_slope_mean = log_summary["mean"],
        log_slope_conf_low = log_summary["lwr"],
        log_slope_conf_high = log_summary["upr"],
        abundance_ratio_mean = ratio_summary["mean"],
        abundance_ratio_median = ratio_summary["median"],
        abundance_ratio_conf_low = ratio_summary["lwr"],
        abundance_ratio_conf_high = ratio_summary["upr"]
      )
    }
  }

  for (jackal_value in c(-1, 0, 1)) {
    early_slope <- sample(
      early_beta[, focal_terms[1]] + jackal_value * early_beta[, focal_terms[3]],
      contrast_draws,
      replace = TRUE
    )
    late_slope <- sample(
      late_beta[, focal_terms[1]] + jackal_value * late_beta[, focal_terms[3]],
      contrast_draws,
      replace = TRUE
    )
    difference <- late_slope - early_slope
    summary_i <- summarize_draws(difference)
    conditional_wolf_contrast_rows[[length(conditional_wolf_contrast_rows) + 1]] <- data.frame(
      species = sp$label,
      jackal_value_sd = jackal_value,
      contrast = "Late minus Early wolf log-slope",
      mean = summary_i["mean"],
      median = summary_i["median"],
      conf_low = summary_i["lwr"],
      conf_high = summary_i["upr"],
      posterior_probability_gt_zero = mean(difference > 0)
    )
  }
}

coefficient_contrast_df <- bind_rows(coefficient_contrast_rows)
conditional_wolf_df <- bind_rows(conditional_wolf_rows)
conditional_wolf_contrast_df <- bind_rows(conditional_wolf_contrast_rows)

write_csv(coefficient_df, paste0(output_prefix, "_coefficients.csv"))
write_csv(diagnostic_df, paste0(output_prefix, "_mcmc_diagnostics.csv"))
write_csv(support_df, paste0(output_prefix, "_fit_support.csv"))
write_csv(relative_grid_df, paste0(output_prefix, "_relative_prediction_grid.csv"))
write_csv(relative_curve_df, paste0(output_prefix, "_relative_interaction_curves.csv"))
write_csv(
  coefficient_contrast_df,
  paste0(output_prefix, "_coefficient_season_contrasts.csv")
)
write_csv(
  conditional_wolf_df,
  paste0(output_prefix, "_conditional_wolf_effects.csv")
)
write_csv(
  conditional_wolf_contrast_df,
  paste0(output_prefix, "_conditional_wolf_season_contrasts.csv")
)
saveRDS(fits, paste0(output_prefix, "_fits.rds"))

species_levels <- vapply(species_info, `[[`, "", "label")
period_plot_labels <- c(
  Early = "Early (hotter)",
  Late = "Late (less-hot)"
)
relative_grid_df <- relative_grid_df %>%
  mutate(
    species = factor(species, levels = species_levels),
    period = factor(period, levels = period_levels)
  )
relative_curve_df <- relative_curve_df %>%
  mutate(
    species = factor(species, levels = species_levels),
    period = factor(period, levels = period_levels),
    wolf_level = factor(
      wolf_z,
      levels = c(-1, 0, 1),
      labels = c("Low wolf (-1 SD)", "Mean wolf", "High wolf (+1 SD)")
    )
  )
conditional_wolf_df <- conditional_wolf_df %>%
  mutate(
    species = factor(species, levels = species_levels),
    period = factor(period, levels = period_levels),
    jackal_level = factor(
      jackal_value_sd,
      levels = c(-1, 0, 1),
      labels = c("Low jackal (-1 SD)", "Mean jackal", "High jackal (+1 SD)")
    )
  )

heat_plot <- ggplot(
  relative_grid_df,
  aes(x = factor(jackal_z), y = factor(wolf_z), fill = log(estimate))
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", estimate)), size = 3.1) +
  facet_grid(
    period ~ species,
    labeller = labeller(period = period_plot_labels)
  ) +
  scale_fill_gradientn(
    colors = c("#F7F9FA", "#DCE6EB", "#AFC4CF", "#6F91A4", "#173F5F"),
    values = scales::rescale(
      log(c(0.15, 0.5, 1, 2, 4)),
      from = range(log(c(0.15, 4)))
    ),
    limits = log(c(0.15, 4)),
    oob = scales::squish,
    breaks = log(c(0.25, 0.5, 1, 2, 4)),
    labels = c("0.25", "0.5", "1", "2", "4"),
    name = "Relative\nabundance"
  ) +
  labs(
    x = "Jackal activity (SD)",
    y = "Wolf activity (SD)",
    title = "Class-specific relative abundance across predator contexts by survey period",
    subtitle = "Posterior medians relative to mean wolf and jackal activity within each period"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

palette <- c(
  "Low wolf (-1 SD)" = "#94AFC0",
  "Mean wolf" = "#4F7489",
  "High wolf (+1 SD)" = "#173F5F"
)
curve_plot <- ggplot(
  relative_curve_df,
  aes(x = jackal_z, y = estimate, color = wolf_level, fill = wolf_level)
) +
  geom_hline(yintercept = 1, color = "#AAB2B5", linewidth = 0.4) +
  geom_ribbon(
    aes(ymin = conf_low, ymax = conf_high),
    alpha = 0.09,
    color = NA
  ) +
  geom_line(linewidth = 0.95) +
  facet_grid(
    period ~ species,
    labeller = labeller(period = period_plot_labels)
  ) +
  scale_y_log10() +
  scale_color_manual(values = palette, name = NULL) +
  scale_fill_manual(values = palette, name = NULL) +
  labs(
    x = "Jackal activity (standardized)",
    y = "Relative expected abundance",
    title = "Wolf-jackal abundance gradients by survey period",
    subtitle = "Posterior medians and 95% credible bands; reference = mean predator activity within period",
    caption = "Relative predictions exclude the abundance intercept and do not estimate absolute population size."
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

wolf_effect_plot <- ggplot(
  conditional_wolf_df,
  aes(
    x = jackal_level,
    y = abundance_ratio_median,
    color = period,
    group = period
  )
) +
  geom_hline(yintercept = 1, color = "#879196", linetype = 2) +
  geom_errorbar(
    aes(
      ymin = abundance_ratio_conf_low,
      ymax = abundance_ratio_conf_high
    ),
    width = 0.08,
    position = position_dodge(width = 0.28)
  ) +
  geom_point(size = 2.8, position = position_dodge(width = 0.28)) +
  facet_wrap(~species, nrow = 1) +
  scale_y_log10() +
  scale_color_manual(
    values = c(Early = "#7B9BAD", Late = "#173F5F"),
    labels = period_plot_labels
  ) +
  labs(
    x = "Jackal activity",
    y = "Abundance ratio per 1-SD increase in wolf activity\n(posterior median and 95% CrI)",
    color = "Survey period",
    title = "Conditional wolf associations by class and survey period"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

ggsave(
  paste0(output_prefix, "_relative_heatmap.png"),
  heat_plot,
  width = 10.2,
  height = 7.4,
  dpi = 320
)
ggsave(
  paste0(output_prefix, "_relative_heatmap.pdf"),
  heat_plot,
  width = 10.2,
  height = 7.4
)
ggsave(
  paste0(output_prefix, "_relative_interaction_curves.png"),
  curve_plot,
  width = 10.6,
  height = 7.4,
  dpi = 320
)
ggsave(
  paste0(output_prefix, "_relative_interaction_curves.pdf"),
  curve_plot,
  width = 10.6,
  height = 7.4
)
ggsave(
  paste0(output_prefix, "_conditional_wolf_effects.png"),
  wolf_effect_plot,
  width = 9.6,
  height = 4.8,
  dpi = 320
)
ggsave(
  paste0(output_prefix, "_conditional_wolf_effects.pdf"),
  wolf_effect_plot,
  width = 9.6,
  height = 4.8
)

focal_diagnostics <- diagnostic_df %>%
  filter(component == "abundance", parameter %in% focal_terms)
intercept_diagnostics <- diagnostic_df %>%
  filter(parameter == "(Intercept)")

report <- c(
  "# Seasonal spAbundance wolf x jackal interaction models",
  "",
  "## Analysis",
  "",
  paste0(
    "The survey extended from `", survey_start, "` to `", survey_end,
    "` and was divided at its calendar midpoint (`", survey_split,
    "`), independently of gazelle detections. Camera-rounds were assigned by deployment midpoint."
  ),
  "",
  "- Family: negative binomial.",
  "- Abundance: `~ wolf_z * jackal_z`.",
  "- Detection: `~ log_days_operational.z + cattle_sightings.z`.",
  paste0(
    "- MCMC: `", n_chains, "` chains, `", n_batch * batch_length,
    "` iterations per chain, `", n_burn, "` burn-in, thinning `", n_thin, "`."
  ),
  "- Predictions are relative to mean wolf and jackal activity within each period and exclude the abundance intercept.",
  "",
  "## Fit support",
  "",
  "| Class | Period | Sites | Nonmissing counts | Total detections |",
  "|---|---|---:|---:|---:|"
)

for (i in seq_len(nrow(support_df))) {
  report <- c(
    report,
    paste0(
      "| `", support_df$species[i], "` | `", support_df$period[i], "` | `",
      support_df$n_sites[i], "` | `", support_df$n_nonmissing_cells[i], "` | `",
      support_df$total_count[i], "` |"
    )
  )
}

report <- c(
  report,
  "",
  "## Predator coefficients",
  "",
  "| Class | Period | Term | Mean [95% CrI] | Rhat | ESS |",
  "|---|---|---|---:|---:|---:|"
)
for (i in seq_len(nrow(coefficient_df))) {
  report <- c(
    report,
    paste0(
      "| `", coefficient_df$species[i], "` | `", coefficient_df$period[i],
      "` | `", coefficient_df$term[i], "` | `",
      fmt_ci(
        coefficient_df$mean[i],
        coefficient_df$conf_low[i],
        coefficient_df$conf_high[i]
      ), "` | `", fmt_num(coefficient_df$rhat[i], 3), "` | `",
      round(coefficient_df$ess[i]), "` |"
    )
  )
}

report <- c(
  report,
  "",
  "## Late-minus-Early coefficient contrasts",
  "",
  "These contrasts pair independent posterior draws from the separately fitted period models. They test differences in slopes under an approximation that does not model residual cross-period dependence at the same camera.",
  "",
  "| Class | Term | Late - Early [95% CrI] | Pr(difference > 0) |",
  "|---|---|---:|---:|"
)
for (i in seq_len(nrow(coefficient_contrast_df))) {
  report <- c(
    report,
    paste0(
      "| `", coefficient_contrast_df$species[i], "` | `",
      coefficient_contrast_df$term[i], "` | `",
      fmt_ci(
        coefficient_contrast_df$mean[i],
        coefficient_contrast_df$conf_low[i],
        coefficient_contrast_df$conf_high[i]
      ), "` | `",
      fmt_num(coefficient_contrast_df$posterior_probability_gt_zero[i], 3), "` |"
    )
  )
}

report <- c(
  report,
  "",
  "## Convergence summary",
  "",
  paste0(
    "Focal slope Rhat range: `", fmt_num(min(focal_diagnostics$rhat), 3),
    "-", fmt_num(max(focal_diagnostics$rhat), 3), "`; ESS range: `",
    round(min(focal_diagnostics$ess)), "-", round(max(focal_diagnostics$ess)), "`."
  ),
  paste0(
    "Intercept Rhat range: `", fmt_num(min(intercept_diagnostics$rhat), 3),
    "-", fmt_num(max(intercept_diagnostics$rhat), 3),
    "`. Absolute abundance should not be interpreted when intercept mixing is poor."
  ),
  "",
  "## Outputs",
  "",
  paste0("- Relative heatmap: `", output_prefix, "_relative_heatmap.png` and `.pdf`."),
  paste0("- Relative curves: `", output_prefix, "_relative_interaction_curves.png` and `.pdf`."),
  paste0("- Conditional wolf effects: `", output_prefix, "_conditional_wolf_effects.png` and `.pdf`.")
)

writeLines(report, paste0(output_prefix, "_results.md"))
cat("Saved audited seasonal spAbundance outputs with prefix:", output_prefix, "\n")
