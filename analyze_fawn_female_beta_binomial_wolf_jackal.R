#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tidyr)
})

output_prefix <- "publication_fawn_female_beta_binomial_wolf_jackal"
juvenile_file <- "gazelle_juvenile_abundance_60x7.csv"
female_file <- "gazelle_female_abundance_60x7.csv"
site_cov_file <- "SiteCovs_060324.csv"
obs_cov_file <- "ObsCovs_050524.csv"
prediction_draws <- 5000L
seed <- 20260804

z_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
}

read_count_long <- function(path, value_name) {
  read.csv(path, na.strings = c("NA", ""), check.names = FALSE) %>%
    as_tibble() %>%
    mutate(camera_id = as.integer(camera_name)) %>%
    pivot_longer(
      cols = all_of(paste0("X", 1:7)),
      names_to = "round_label",
      values_to = value_name
    ) %>%
    mutate(round = as.integer(sub("^X", "", round_label))) %>%
    select(camera_id, round, all_of(value_name))
}

build_season_lookup <- function(obs) {
  obs <- obs %>%
    mutate(
      camera_id = as.integer(`Camera Number`),
      round = as.integer(`Round Number`),
      start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
      end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
      mid_date = start_date + floor(as.numeric(end_date - start_date) / 2)
    )
  survey_start <- min(obs$start_date, na.rm = TRUE)
  survey_end <- max(obs$end_date, na.rm = TRUE)
  split_date <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)
  round_mid <- obs %>%
    filter(!is.na(mid_date)) %>%
    group_by(round) %>%
    summarise(round_mid = median(as.numeric(mid_date)), .groups = "drop")
  obs <- obs %>%
    left_join(round_mid, by = "round") %>%
    mutate(
      mid_date = as.Date(
        ifelse(is.na(mid_date), round_mid, as.numeric(mid_date)),
        origin = "1970-01-01"
      ),
      season = factor(
        if_else(mid_date <= split_date, "Early (hotter)", "Late (mild-cold)"),
        levels = c("Early (hotter)", "Late (mild-cold)")
      )
    )
  attr(obs, "split_date") <- split_date
  obs
}

fawn <- read_count_long(juvenile_file, "fawn")
female <- read_count_long(female_file, "female")
counts <- full_join(fawn, female, by = c("camera_id", "round"))
if (any(xor(is.na(counts$fawn), is.na(counts$female)))) {
  stop("Fawn and female count histories have different missing-data patterns.")
}

site_covs <- read_csv(site_cov_file, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    Polygon.Name = factor(Polygon.Name),
    wolf_z = z_scale(Null.Wolf.Abundance),
    jackal_z = z_scale(null.jackal.abundance),
    open_z = z_scale(Open.Spaces),
    built_z = z_scale(Built.and.Roads),
    reserve_z = z_scale(Nature.Reserves),
    human_z = z_scale(human.abundance)
  )
obs_raw <- read.csv(
  obs_cov_file,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
obs <- build_season_lookup(obs_raw) %>%
  transmute(
    camera_id,
    round,
    season,
    path_type = factor(Path.Type),
    camera_model = factor(Camera.Model),
    vegetation_structure = factor(Vegetation.Structure),
    area_behind_camera = factor(Area.Behind.Camera)
  )

model_data <- counts %>%
  left_join(site_covs, by = "camera_id") %>%
  left_join(obs, by = c("camera_id", "round")) %>%
  filter(
    is.finite(fawn),
    is.finite(female),
    fawn + female > 0,
    complete.cases(
      wolf_z,
      jackal_z,
      open_z,
      built_z,
      reserve_z,
      human_z,
      season,
      path_type,
      camera_model,
      vegetation_structure,
      area_behind_camera,
      Polygon.Name
    )
  ) %>%
  mutate(
    camera_id = factor(camera_id),
    Polygon.Name = droplevels(Polygon.Name),
    season = droplevels(season),
    path_type = droplevels(path_type),
    camera_model = droplevels(camera_model),
    vegetation_structure = droplevels(vegetation_structure),
    area_behind_camera = droplevels(area_behind_camera)
  )

control <- glmmTMBControl(
  optCtrl = list(iter.max = 10000, eval.max = 10000),
  profile = TRUE
)

overall_fixed <- ~ wolf_z * jackal_z + open_z + built_z + reserve_z + human_z +
  path_type + camera_model +
  vegetation_structure + area_behind_camera
seasonal_fixed <- ~ wolf_z * jackal_z * season + open_z + built_z + reserve_z +
  human_z + path_type + camera_model +
  vegetation_structure + area_behind_camera

overall_fit <- glmmTMB(
  cbind(fawn, female) ~ wolf_z * jackal_z + open_z + built_z + reserve_z +
    human_z + path_type + camera_model +
    vegetation_structure + area_behind_camera + (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)
overall_additive_fit <- update(overall_fit, . ~ . - wolf_z:jackal_z)

seasonal_fit <- glmmTMB(
  cbind(fawn, female) ~ wolf_z * jackal_z * season + open_z + built_z +
    reserve_z + human_z + path_type + camera_model +
    vegetation_structure + area_behind_camera + (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)
seasonal_no_predator_interaction_fit <- glmmTMB(
  cbind(fawn, female) ~ (wolf_z + jackal_z) * season + open_z + built_z +
    reserve_z + human_z + path_type + camera_model +
    vegetation_structure + area_behind_camera + (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)
seasonal_stable_response_fit <- glmmTMB(
  cbind(fawn, female) ~ wolf_z * jackal_z + season + open_z + built_z +
    reserve_z + human_z + path_type + camera_model +
    vegetation_structure + area_behind_camera + (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)
seasonal_additive_stable_fit <- glmmTMB(
  cbind(fawn, female) ~ wolf_z + jackal_z + season + open_z + built_z +
    reserve_z + human_z + path_type + camera_model + vegetation_structure +
    area_behind_camera + (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)

fits <- list(
  overall = overall_fit,
  overall_additive = overall_additive_fit,
  seasonal = seasonal_fit,
  seasonal_no_predator_interaction = seasonal_no_predator_interaction_fit,
  seasonal_stable_response = seasonal_stable_response_fit,
  seasonal_additive_stable = seasonal_additive_stable_fit
)
saveRDS(fits, paste0(output_prefix, "_model_fits.rds"))

fit_diagnostics <- bind_rows(lapply(names(fits), function(name) {
  fit <- fits[[name]]
  tibble(
    model = name,
    n_observations = stats::nobs(fit),
    AIC = stats::AIC(fit),
    convergence_code = fit$fit$convergence,
    positive_definite_hessian = isTRUE(fit$sdr$pdHess),
    dispersion = sigma(fit)
  )
}))

lrt_to_tibble <- function(smaller, larger, comparison) {
  tab <- as.data.frame(anova(smaller, larger))
  tibble(
    comparison = comparison,
    df_difference = tab$Df[2] - tab$Df[1],
    likelihood_ratio = tab$Chisq[2],
    p_value = tab$`Pr(>Chisq)`[2]
  )
}
lrt_results <- bind_rows(
  lrt_to_tibble(
    overall_additive_fit,
    overall_fit,
    "Overall wolf x jackal interaction"
  ),
  lrt_to_tibble(
    overall_fit,
    seasonal_stable_response_fit,
    "Common Early-Late shift"
  ),
  lrt_to_tibble(
    overall_additive_fit,
    seasonal_additive_stable_fit,
    "Common Early-Late shift in additive model"
  ),
  lrt_to_tibble(
    overall_fit,
    seasonal_fit,
    "Any seasonal structure"
  ),
  lrt_to_tibble(
    seasonal_no_predator_interaction_fit,
    seasonal_fit,
    "Seasonal wolf x jackal terms"
  ),
  lrt_to_tibble(
    seasonal_stable_response_fit,
    seasonal_fit,
    "Seasonal change in predator response"
  ),
  lrt_to_tibble(
    seasonal_additive_stable_fit,
    seasonal_no_predator_interaction_fit,
    "Seasonal change in additive predator slopes"
  )
)

extract_coefficients <- function(fit, model) {
  tab <- as.data.frame(summary(fit)$coefficients$cond)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  as_tibble(tab) %>%
    transmute(
      model = model,
      term,
      estimate = Estimate,
      std_error = `Std. Error`,
      z_value = `z value`,
      p_value = `Pr(>|z|)`,
      odds_ratio = exp(estimate),
      conf_low = exp(estimate - 1.96 * std_error),
      conf_high = exp(estimate + 1.96 * std_error)
    )
}
coefficient_df <- bind_rows(
  extract_coefficients(overall_fit, "overall"),
  extract_coefficients(seasonal_fit, "seasonal")
)

standardized_prediction <- function(
    fit,
    fixed_formula,
    grid,
    reference_data,
    n_draws,
    seed) {
  beta_hat <- glmmTMB::fixef(fit)$cond
  beta_vcov <- as.matrix(stats::vcov(fit)$cond)
  set.seed(seed)
  beta_draws <- MASS::mvrnorm(n_draws, mu = beta_hat, Sigma = beta_vcov)
  if (is.vector(beta_draws)) beta_draws <- matrix(beta_draws, nrow = 1)
  colnames(beta_draws) <- names(beta_hat)

  bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    newdata <- reference_data
    newdata$wolf_z <- grid$wolf_z[i]
    newdata$jackal_z <- grid$jackal_z[i]
    if ("season" %in% names(grid)) newdata$season <- grid$season[i]
    X <- model.matrix(fixed_formula, data = newdata)
    X <- X[, names(beta_hat), drop = FALSE]
    p_draws <- rowMeans(stats::plogis(beta_draws %*% t(X)))
    point_p <- mean(stats::plogis(as.numeric(X %*% beta_hat)))
    ratio_draws <- 100 * p_draws / (1 - p_draws)
    tibble(
      period = as.character(grid$period[i]),
      wolf_z = grid$wolf_z[i],
      jackal_z = grid$jackal_z[i],
      estimated_fawns_per_100_females = 100 * point_p / (1 - point_p),
      median = median(ratio_draws),
      lower = unname(quantile(ratio_draws, 0.025)),
      upper = unname(quantile(ratio_draws, 0.975)),
      probability_below_10 = mean(ratio_draws < 10),
      probability_above_20 = mean(ratio_draws > 20)
    )
  }))
}

# Standardize every prediction over the same observed camera-design distribution.
reference_data <- model_data %>%
  select(
    path_type,
    camera_model,
    vegetation_structure,
    area_behind_camera,
    season,
    open_z,
    built_z,
    reserve_z,
    human_z
  ) %>%
  mutate(
    open_z = 0,
    built_z = 0,
    reserve_z = 0,
    human_z = 0
  )

overall_grid <- expand.grid(
  wolf_z = c(-1, 0, 1),
  jackal_z = c(-1, 0, 1),
  KEEP.OUT.ATTRS = FALSE
) %>%
  mutate(period = "Overall")
seasonal_grid <- expand.grid(
  wolf_z = c(-1, 0, 1),
  jackal_z = c(-1, 0, 1),
  season = levels(model_data$season),
  KEEP.OUT.ATTRS = FALSE
) %>%
  mutate(period = as.character(season))

prediction_df <- bind_rows(
  standardized_prediction(
    overall_fit,
    overall_fixed,
    overall_grid,
    reference_data,
    prediction_draws,
    seed
  ),
  standardized_prediction(
    seasonal_fit,
    seasonal_fixed,
    seasonal_grid,
    reference_data,
    prediction_draws,
    seed + 1L
  )
) %>%
  mutate(
    period = factor(
      period,
      levels = c("Overall", "Early (hotter)", "Late (mild-cold)")
    ),
    wolf_level = factor(wolf_z, levels = c(-1, 0, 1)),
    jackal_level = factor(jackal_z, levels = c(-1, 0, 1)),
    label = sprintf(
      "%.1f\n[%.1f, %.1f]",
      estimated_fawns_per_100_females,
      lower,
      upper
    )
  )

write_csv(prediction_df, paste0(output_prefix, "_prediction_summary.csv"))
write_csv(coefficient_df, paste0(output_prefix, "_coefficients.csv"))
write_csv(lrt_results, paste0(output_prefix, "_likelihood_ratio_tests.csv"))
write_csv(fit_diagnostics, paste0(output_prefix, "_fit_diagnostics.csv"))
write_csv(
  model_data %>%
    select(
      camera_id, round, season, fawn, female, wolf_z, jackal_z,
      open_z, built_z, reserve_z, human_z
    ),
  paste0(output_prefix, "_analysis_data.csv")
)

fill_upper <- max(25, max(prediction_df$estimated_fawns_per_100_females))
prediction_df <- prediction_df %>%
  mutate(
    text_color = if_else(
      estimated_fawns_per_100_females >= 0.60 * fill_upper,
      "white",
      "#263238"
    )
  )

ratio_plot <- ggplot(
  prediction_df,
  aes(x = jackal_level, y = wolf_level, fill = estimated_fawns_per_100_females)
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
    name = "Fawns per\n100 females"
  ) +
  coord_equal() +
  labs(
    x = "Jackal abundance (z)",
    y = "Wolf abundance (z)",
    title = "Camera-based fawn representation across wolf-jackal context",
    subtitle = paste0(
      "Beta-binomial camera-round composition model standardized over camera-design covariates\n",
      "Cells show estimated fawns per 100 females [95% model-based interval]"
    ),
    caption = paste0(
      "Wolf and jackal abundance are continuous and evaluated at -1, 0, and +1 SD. ",
      "This is an apparent composition index, not a recruitment or survival rate."
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

fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt_p <- function(x) ifelse(x < 0.001, "<0.001", fmt(x, 3))
md <- c(
  "# Beta-binomial fawn:female composition across wolf-jackal context",
  "",
  "## Model",
  "",
  "Camera-round response: `cbind(fawn, female)`, excluding rounds with zero combined fawn and female detections. The model uses a beta-binomial logit link, continuous standardized wolf and jackal abundance, open-space, built-and-road, reserve, and human-activity controls, camera-design controls, and nested polygon/camera random intercepts.",
  "",
  paste0("Informative camera-rounds: `", nrow(model_data), "`; cameras: `", n_distinct(model_data$camera_id), "`; independent deployment-midpoint season split: `", as.character(attr(build_season_lookup(obs_raw), "split_date")), "`."),
  "",
  "Predictions are population-level fixed-effect predictions standardized over the observed distribution of path type, camera model, vegetation structure, and area behind the camera. Landscape and human-activity controls are held at their standardized means.",
  "",
  "## Model tests",
  "",
  "| Comparison | LR chi-square | df | p |",
  "|---|---:|---:|---:|"
)
for (i in seq_len(nrow(lrt_results))) {
  md <- c(md, sprintf(
    "| `%s` | `%s` | `%s` | `%s` |",
    lrt_results$comparison[i], fmt(lrt_results$likelihood_ratio[i]),
    lrt_results$df_difference[i], fmt_p(lrt_results$p_value[i])
  ))
}
md <- c(
  md,
  "",
  "## Prediction grid",
  "",
  "| Period | Wolf z | Jackal z | Fawns per 100 females | 95% interval |",
  "|---|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(prediction_df))) {
  md <- c(md, sprintf(
    "| `%s` | `%s` | `%s` | `%s` | `[%s, %s]` |",
    as.character(prediction_df$period[i]), prediction_df$wolf_z[i],
    prediction_df$jackal_z[i],
    fmt(prediction_df$estimated_fawns_per_100_females[i], 1),
    fmt(prediction_df$lower[i], 1), fmt(prediction_df$upper[i], 1)
  ))
}
md <- c(
  md,
  "",
  "## Interpretation constraint",
  "",
  "This model estimates fawn representation among classified fawn and adult-female camera counts. It accounts for repeated cameras, overdispersion, and measured camera-design covariates, but it does not fully estimate class-specific detection probability. It should be described as an apparent fawn:female composition index rather than realized recruitment, survival, or a population sex-age ratio."
)
writeLines(md, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved beta-binomial ratio outputs with prefix:", output_prefix, "\n")
