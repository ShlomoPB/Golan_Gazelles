if (!requireNamespace("unmarked", quietly = TRUE)) {
  stop("Package 'unmarked' is required.")
}
if (!requireNamespace("spAbundance", quietly = TRUE)) {
  stop("Package 'spAbundance' is required.")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}

library(unmarked)
library(spAbundance)
library(readr)
library(dplyr)
library(ggplot2)

output_prefix <- "publication_wolf_jackal_interaction_nmix"
seed <- 20260623
n_batch <- 2000
batch_length <- 25
n_burn <- 5000
n_thin <- 10
n_chains <- 3

species_info <- list(
  list(label = "Fawn", code = "fawn", count_file = "gazelle_juvenile_abundance_60x7.csv", data_rds = "gazelle_fawn_spNMix_data.rds"),
  list(label = "Female", code = "female", count_file = "gazelle_female_abundance_60x7.csv", data_rds = "gazelle_female_spNMix_data.rds"),
  list(label = "Male", code = "male", count_file = "gazelle_male_abundance_60x7.csv", data_rds = "gazelle_male_spNMix_data.rds")
)

site_cov_file <- "SiteCovs_060324.csv"
obs_cov_file <- "ObsCovs_050524.csv"

fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_ci <- function(est, lwr, upr, digits = 2) {
  paste0(fmt_num(est, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt_num(x, 3)))
}

z_scale <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

read_count_matrix <- function(path) {
  x <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  y <- as.matrix(x[paste0("X", 1:7)])
  storage.mode(y) <- "numeric"
  rownames(y) <- x$camera_name
  colnames(y) <- paste0("X", 1:7)
  y
}

to_site_by_round <- function(df, column, site_ids) {
  matrix(
    df[[column]],
    nrow = length(site_ids),
    ncol = 7,
    byrow = TRUE,
    dimnames = list(site = site_ids, replicate = paste0("X", 1:7))
  )
}

make_unmarked_obs_covs <- function(obs_covs, site_ids) {
  list(
    Path.Type = to_site_by_round(obs_covs, "Path.Type", site_ids),
    Camera.Model = to_site_by_round(obs_covs, "Camera.Model", site_ids),
    Vegetation.Structure = to_site_by_round(obs_covs, "Vegetation.Structure", site_ids),
    Area.Behind.Camera = to_site_by_round(obs_covs, "Area.Behind.Camera", site_ids)
  )
}

make_spabund_det_frame <- function(det_covs) {
  out <- lapply(det_covs, function(x) {
    if (is.matrix(x) || is.data.frame(x)) {
      as.vector(x)
    } else {
      x
    }
  })
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  for (nm in names(out)) {
    if (is.character(out[[nm]])) {
      out[[nm]] <- factor(out[[nm]])
    }
  }
  out
}

site_covs_raw <- read.csv(site_cov_file, na.strings = c("NA", ""), check.names = FALSE)
site_covs_raw$Camera.ID <- as.integer(site_covs_raw$Camera.ID)

obs_covs_raw <- read.csv(obs_cov_file, na.strings = c("NA", ""), check.names = FALSE)

unmarked_results <- list()
spabund_results <- list()
pred_rows <- list()
coef_rows <- list()
idx_coef <- 1
diagnostic_rows <- list()

for (sp in species_info) {
  y <- read_count_matrix(sp$count_file)
  site_ids <- rownames(y)

  site_covs <- site_covs_raw[match(as.integer(site_ids), site_covs_raw$Camera.ID), ]
  site_covs <- site_covs %>%
    mutate(
      wolf_z = z_scale(Null.Wolf.Abundance),
      jackal_z = z_scale(null.jackal.abundance)
    )

  obs_covs <- obs_covs_raw[order(obs_covs_raw$`Camera Number`, obs_covs_raw$`Round Number`), ]
  expected <- expand.grid(`Camera Number` = as.integer(site_ids), `Round Number` = 1:7)
  expected <- expected[order(expected$`Camera Number`, expected$`Round Number`), ]
  stopifnot(identical(obs_covs$`Camera Number`, expected$`Camera Number`))
  stopifnot(identical(obs_covs$`Round Number`, expected$`Round Number`))

  obs_list <- make_unmarked_obs_covs(obs_covs, site_ids)

  umf <- unmarkedFramePCount(
    y = y,
    siteCovs = data.frame(
      wolf_z = site_covs$wolf_z,
      jackal_z = site_covs$jackal_z,
      row.names = site_ids
    ),
    obsCovs = obs_list
  )

  set.seed(seed)
  fit_um <- pcount(
    ~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ wolf_z * jackal_z,
    data = umf,
    mixture = "NB",
    K = max(150, max(y, na.rm = TRUE) + 50)
  )

  summ_um <- capture.output(print(summary(fit_um)))
  unmarked_results[[sp$code]] <- fit_um

  state_start <- grep("^Abundance \\(log-scale\\):", summ_um)
  det_start <- grep("^Detection \\(logit-scale\\):", summ_um)
  state_lines <- summ_um[(state_start + 2):(det_start - 2)]
  state_lines <- state_lines[nzchar(trimws(state_lines))]
  state_split <- strsplit(trimws(state_lines), "\\s+")
  state_df <- do.call(rbind, lapply(state_split, function(x) {
    data.frame(
      term = x[1],
      estimate = as.numeric(x[2]),
      std_error = as.numeric(x[3]),
      z_value = as.numeric(x[4]),
      p_value = as.numeric(x[5]),
      stringsAsFactors = FALSE
    )
  }))

  coef_rows[[idx_coef]] <- data.frame(
    species = sp$label,
    framework = "unmarked",
    term = state_df$term,
    estimate = state_df$estimate,
    std_error = state_df$std_error,
    z_value = state_df$z_value,
    p_value = state_df$p_value,
    stringsAsFactors = FALSE
  )
  idx_coef <- idx_coef + 1

  data_list <- readRDS(sp$data_rds)
  abund_formula <- ~ Null.Wolf.Abundance.z * null.jackal.abundance.z
  det_formula <- ~ log_days_operational.z + cattle_sightings.z

  X_abund <- model.matrix(abund_formula, data = data_list$abund.covs)
  X_det <- model.matrix(det_formula, data = make_spabund_det_frame(data_list$det.covs))
  priors <- list(
    beta.normal = list(mean = rep(0, ncol(X_abund)), var = rep(100, ncol(X_abund))),
    alpha.normal = list(mean = rep(0, ncol(X_det)), var = rep(2.72, ncol(X_det))),
    kappa.unif = c(0, 10)
  )
  N_init <- apply(data_list$y, 1, max, na.rm = TRUE)
  N_init[!is.finite(N_init)] <- 0
  inits <- list(alpha = 0, beta = 0, N = N_init, kappa = 1)
  tuning <- list(beta = 0.1, alpha = 0.1, kappa = 0.5)

  set.seed(seed)
  fit_sp <- NMix(
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
    verbose = FALSE
  )

  spabund_results[[sp$code]] <- fit_sp
  sum_sp <- summary(fit_sp)
  sp_text <- capture.output(print(sum_sp))
  sp_lines <- grep("Null.Wolf.Abundance.z|null.jackal.abundance.z|\\*|kappa", sp_text, value = TRUE)

  diagnostic_rows[[length(diagnostic_rows) + 1]] <- bind_rows(
    data.frame(
      species = sp$label,
      component = "abundance",
      parameter = colnames(fit_sp$beta.samples),
      rhat = as.numeric(fit_sp$rhat$beta),
      ess = as.numeric(fit_sp$ESS$beta)
    ),
    data.frame(
      species = sp$label,
      component = "detection",
      parameter = colnames(fit_sp$alpha.samples),
      rhat = as.numeric(fit_sp$rhat$alpha),
      ess = as.numeric(fit_sp$ESS$alpha)
    ),
    data.frame(
      species = sp$label,
      component = "overdispersion",
      parameter = "kappa",
      rhat = as.numeric(fit_sp$rhat$kappa),
      ess = as.numeric(fit_sp$ESS$kappa)
    )
  )

  beta_tab <- fit_sp$beta.samples
  term_map <- colnames(beta_tab)
  keep_terms <- term_map[grepl("Null.Wolf.Abundance.z|null.jackal.abundance.z", term_map)]
  for (tm in keep_terms) {
    draws <- beta_tab[, tm]
    term_index <- match(tm, term_map)
    rhat_val <- NA_real_
    ess_val <- NA_real_
    if (!is.null(fit_sp$rhat$beta) && is.finite(term_index)) {
      rhat_val <- unname(fit_sp$rhat$beta[term_index])
    }
    if (!is.null(fit_sp$ESS$beta) && is.finite(term_index)) {
      ess_val <- unname(fit_sp$ESS$beta[term_index])
    }
    coef_rows[[idx_coef]] <- data.frame(
      species = sp$label,
      framework = "spAbundance",
      term = tm,
      estimate = mean(draws),
      std_error = sd(draws),
      z_value = NA_real_,
      p_value = NA_real_,
      lwr = quantile(draws, 0.025),
      upr = quantile(draws, 0.975),
      rhat = rhat_val,
      ess = ess_val,
      stringsAsFactors = FALSE
    )
    idx_coef <- idx_coef + 1
  }

  pred_grid <- expand.grid(
    wolf_z = c(-1, 0, 1),
    jackal_z = c(-1, 0, 1),
    Path.Type = unique(as.vector(obs_list$Path.Type))[1],
    Camera.Model = unique(as.vector(obs_list$Camera.Model))[1],
    Vegetation.Structure = unique(as.vector(obs_list$Vegetation.Structure))[1],
    Area.Behind.Camera = unique(as.vector(obs_list$Area.Behind.Camera))[1]
  )
  pred_grid$framework <- "unmarked"
  pred_grid$species <- sp$label
  pred_grid$lambda <- predict(fit_um, type = "state", newdata = pred_grid)$Predicted
  pred_rows[[length(pred_rows) + 1]] <- pred_grid

  pred_grid_sp <- expand.grid(
    Null.Wolf.Abundance.z = c(-1, 0, 1),
    null.jackal.abundance.z = c(-1, 0, 1)
  )
  Xnew <- model.matrix(~ Null.Wolf.Abundance.z * null.jackal.abundance.z, data = pred_grid_sp)
  beta_draws <- fit_sp$beta.samples[, colnames(Xnew), drop = FALSE]
  eta_draws <- beta_draws %*% t(Xnew)
  mu_pred <- apply(exp(eta_draws), 2, mean)
  pred_rows[[length(pred_rows) + 1]] <- data.frame(
    wolf_z = pred_grid_sp$Null.Wolf.Abundance.z,
    jackal_z = pred_grid_sp$null.jackal.abundance.z,
    framework = "spAbundance",
    species = sp$label,
    lambda = mu_pred
  )
}

coef_df <- bind_rows(coef_rows)
pred_df <- bind_rows(pred_rows)
diagnostic_df <- bind_rows(diagnostic_rows)
write_csv(coef_df, paste0(output_prefix, "_coefficients.csv"))
write_csv(pred_df, paste0(output_prefix, "_prediction_grid.csv"))
write_csv(diagnostic_df, paste0(output_prefix, "_spabundance_diagnostics.csv"))
saveRDS(unmarked_results, paste0(output_prefix, "_unmarked_fits.rds"))
saveRDS(spabund_results, paste0(output_prefix, "_spabundance_fits.rds"))

heat_plot <- ggplot(pred_df, aes(x = factor(jackal_z), y = factor(wolf_z), fill = lambda)) +
  geom_tile(color = "white") +
  geom_text(aes(label = fmt_num(lambda, 2)), size = 3) +
  facet_grid(framework ~ species) +
  scale_fill_gradient(low = "#dceaf5", high = "#1d3557") +
  labs(
    x = "Jackal abundance (z)",
    y = "Wolf abundance (z)",
    fill = "Predicted\nabundance",
    title = "Predicted abundance across wolf x jackal interaction grid"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(paste0(output_prefix, "_heatmap.png"), heat_plot, width = 10.5, height = 5.8, dpi = 300)
ggsave(paste0(output_prefix, "_heatmap.pdf"), heat_plot, width = 10.5, height = 5.8)

pred_df_sp <- filter(pred_df, framework == "spAbundance")
heat_plot_sp <- ggplot(pred_df_sp, aes(x = factor(jackal_z), y = factor(wolf_z), fill = lambda)) +
  geom_tile(color = "white") +
  geom_text(aes(label = fmt_num(lambda, 2)), size = 3.2) +
  facet_wrap(~ species, ncol = 3) +
  scale_fill_gradient(low = "#dceaf5", high = "#1d3557") +
  labs(
    x = "Jackal abundance (z)",
    y = "Wolf abundance (z)",
    fill = "Predicted\nabundance",
    title = "spAbundance predicted abundance across wolf x jackal interaction grid"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(paste0(output_prefix, "_heatmap_spabundance_only.png"), heat_plot_sp, width = 9.4, height = 4.8, dpi = 300)
ggsave(paste0(output_prefix, "_heatmap_spabundance_only.pdf"), heat_plot_sp, width = 9.4, height = 4.8)

line_grid <- expand.grid(
  wolf_z = c(-1, 0, 1),
  jackal_z = seq(-1.5, 1.5, length.out = 121)
)
line_rows <- list()

for (sp in species_info) {
  fit_um <- unmarked_results[[sp$code]]
  um_grid <- transform(
    line_grid,
    Path.Type = unique(as.vector(make_unmarked_obs_covs(obs_covs_raw[order(obs_covs_raw$`Camera Number`, obs_covs_raw$`Round Number`), ], rownames(read_count_matrix(sp$count_file)))$Path.Type))[1],
    Camera.Model = unique(as.vector(make_unmarked_obs_covs(obs_covs_raw[order(obs_covs_raw$`Camera Number`, obs_covs_raw$`Round Number`), ], rownames(read_count_matrix(sp$count_file)))$Camera.Model))[1],
    Vegetation.Structure = unique(as.vector(make_unmarked_obs_covs(obs_covs_raw[order(obs_covs_raw$`Camera Number`, obs_covs_raw$`Round Number`), ], rownames(read_count_matrix(sp$count_file)))$Vegetation.Structure))[1],
    Area.Behind.Camera = unique(as.vector(make_unmarked_obs_covs(obs_covs_raw[order(obs_covs_raw$`Camera Number`, obs_covs_raw$`Round Number`), ], rownames(read_count_matrix(sp$count_file)))$Area.Behind.Camera))[1]
  )
  um_pred <- predict(fit_um, type = "state", newdata = um_grid)$Predicted
  line_rows[[length(line_rows) + 1]] <- data.frame(
    species = sp$label,
    framework = "unmarked",
    wolf_z = line_grid$wolf_z,
    jackal_z = line_grid$jackal_z,
    lambda = um_pred
  )

  fit_sp <- spabund_results[[sp$code]]
  sp_grid <- data.frame(
    Null.Wolf.Abundance.z = line_grid$wolf_z,
    null.jackal.abundance.z = line_grid$jackal_z
  )
  Xline <- model.matrix(~ Null.Wolf.Abundance.z * null.jackal.abundance.z, data = sp_grid)
  beta_draws <- fit_sp$beta.samples[, colnames(Xline), drop = FALSE]
  eta_draws <- beta_draws %*% t(Xline)
  mu_line <- apply(exp(eta_draws), 2, mean)
  line_rows[[length(line_rows) + 1]] <- data.frame(
    species = sp$label,
    framework = "spAbundance",
    wolf_z = line_grid$wolf_z,
    jackal_z = line_grid$jackal_z,
    lambda = mu_line
  )
}

line_df <- bind_rows(line_rows) %>%
  mutate(wolf_level = factor(wolf_z, levels = c(-1, 0, 1), labels = c("Low wolf (-1 SD)", "Mean wolf", "High wolf (+1 SD)")))

write_csv(line_df, paste0(output_prefix, "_line_predictions.csv"))

line_plot <- ggplot(line_df, aes(x = jackal_z, y = lambda, color = wolf_level)) +
  geom_line(linewidth = 1) +
  facet_grid(framework ~ species, scales = "free_y") +
  scale_color_manual(values = c("Low wolf (-1 SD)" = "#9b2226", "Mean wolf" = "#6c757d", "High wolf (+1 SD)" = "#1d3557")) +
  labs(
    x = "Jackal abundance (z)",
    y = "Predicted abundance",
    color = NULL,
    title = "Wolf x jackal interaction curves by class"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0(output_prefix, "_interaction_curves.png"), line_plot, width = 10.8, height = 5.8, dpi = 300)
ggsave(paste0(output_prefix, "_interaction_curves.pdf"), line_plot, width = 10.8, height = 5.8)

line_df_sp <- filter(line_df, framework == "spAbundance")
line_plot_sp <- ggplot(line_df_sp, aes(x = jackal_z, y = lambda, color = wolf_level)) +
  geom_line(linewidth = 1.05) +
  facet_wrap(~ species, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("Low wolf (-1 SD)" = "#9b2226", "Mean wolf" = "#6c757d", "High wolf (+1 SD)" = "#1d3557")) +
  labs(
    x = "Jackal abundance (z)",
    y = "Predicted abundance",
    color = NULL,
    title = "spAbundance interaction curves: wolf x jackal by class"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0(output_prefix, "_interaction_curves_spabundance_only.png"), line_plot_sp, width = 9.4, height = 4.8, dpi = 300)
ggsave(paste0(output_prefix, "_interaction_curves_spabundance_only.pdf"), line_plot_sp, width = 9.4, height = 4.8)

md <- c(
  "# Wolf x Jackal interaction N-mixture comparison",
  "",
  "Models were fit separately for fawns, females, and males. The `unmarked` models used the camera- and habitat-based detection model:",
  "",
  "`detection ~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera`",
  "",
  "The `spAbundance` models used the simpler effort-based detection model:",
  "",
  "`detection ~ log_days_operational.z + cattle_sightings.z`",
  "",
  "Both frameworks used the predator-only abundance model:",
  "",
  "`abundance ~ wolf * jackal`",
  "",
  "Unmarked models used `pcount(..., mixture = \"NB\")`. spAbundance models used `NMix(..., family = \"NB\")`.",
  "",
  "## Coefficient summary",
  "",
  "| Species | Framework | Term | Estimate or posterior mean | Interval / p |",
  "|---|---|---|---:|---:|"
)

for (i in seq_len(nrow(coef_df))) {
  if (coef_df$framework[i] == "unmarked") {
    md <- c(md, paste0(
      "| `", coef_df$species[i], "` | `", coef_df$framework[i], "` | `", coef_df$term[i],
      "` | `", fmt_num(coef_df$estimate[i], 3), "` | `p = ", fmt_p(coef_df$p_value[i]), "` |"
    ))
  } else {
    md <- c(md, paste0(
      "| `", coef_df$species[i], "` | `", coef_df$framework[i], "` | `", coef_df$term[i],
      "` | `", fmt_num(coef_df$estimate[i], 3), "` | `", fmt_ci(coef_df$estimate[i], coef_df$lwr[i], coef_df$upr[i], 2), "` |"
    ))
  }
}

md <- c(
  md,
  "",
  "## Interpretation",
  "",
  "Both frameworks recovered the same qualitative predator pattern: jackal effects were negative, wolf effects were positive especially for adults, and wolf x jackal interactions were positive.",
  "",
  "The adult interaction terms were supported in both `unmarked` and `spAbundance`, which makes predator context a stronger result than if it appeared in only one model family. The ecological interpretation is that the negative jackal association weakened where wolf activity was high, or equivalently that the positive wolf association was strongest where jackal activity was also high.",
  "",
  "This supports a two-predator-context hypothesis for gazelle spatial abundance. It should not be interpreted as direct evidence of wolf-mediated jackal suppression without temporally resolved predator data."
)

writeLines(md, paste0(output_prefix, "_results.md"))

cat("Saved wolf x jackal interaction N-mixture outputs with prefix:", output_prefix, "\n")
