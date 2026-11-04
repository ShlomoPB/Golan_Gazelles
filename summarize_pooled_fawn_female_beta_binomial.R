#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
})

fit_file <- "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
prediction_file <- "publication_fawn_female_beta_binomial_wolf_jackal_prediction_summary.csv"
output_prefix <- "publication_fawn_female_beta_binomial_pooled"
seed <- 20260805
n_draws <- 5000L

fits <- readRDS(fit_file)
interaction_fit <- fits$overall
additive_fit <- fits$overall_additive

extract_coefficients <- function(fit, model) {
  tab <- as.data.frame(summary(fit)$coefficients$cond)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  as_tibble(tab) %>%
    transmute(
      model,
      term,
      estimate_log_odds = Estimate,
      std_error = `Std. Error`,
      odds_ratio = exp(Estimate),
      conf_low = exp(Estimate - 1.96 * `Std. Error`),
      conf_high = exp(Estimate + 1.96 * `Std. Error`),
      z_value = `z value`,
      p_value = `Pr(>|z|)`
    )
}

coefficient_df <- bind_rows(
  extract_coefficients(interaction_fit, "Pooled interaction"),
  extract_coefficients(additive_fit, "Pooled additive")
)
write_csv(coefficient_df, paste0(output_prefix, "_coefficients.csv"))

lrt_tab <- as.data.frame(anova(additive_fit, interaction_fit))
model_comparison <- tibble(
  comparison = "Pooled wolf x jackal interaction",
  additive_AIC = AIC(additive_fit),
  interaction_AIC = AIC(interaction_fit),
  delta_AIC_interaction_minus_additive = AIC(interaction_fit) - AIC(additive_fit),
  df_difference = lrt_tab$Df[2] - lrt_tab$Df[1],
  likelihood_ratio = lrt_tab$Chisq[2],
  p_value = lrt_tab$`Pr(>Chisq)`[2]
)
write_csv(model_comparison, paste0(output_prefix, "_model_comparison.csv"))

season_lrt_tab <- as.data.frame(anova(additive_fit, fits$seasonal_additive_stable))
season_comparison <- tibble(
  comparison = "Common Early-Late shift in additive model",
  pooled_AIC = AIC(additive_fit),
  season_adjusted_AIC = AIC(fits$seasonal_additive_stable),
  delta_AIC_season_adjusted_minus_pooled =
    AIC(fits$seasonal_additive_stable) - AIC(additive_fit),
  df_difference = season_lrt_tab$Df[2] - season_lrt_tab$Df[1],
  likelihood_ratio = season_lrt_tab$Chisq[2],
  p_value = season_lrt_tab$`Pr(>Chisq)`[2]
)
write_csv(season_comparison, paste0(output_prefix, "_season_comparison.csv"))

linear_contrast <- function(fit, weights, effect, other_predator_z) {
  beta <- fixef(fit)$cond
  variance <- as.matrix(vcov(fit)$cond)
  contrast <- setNames(rep(0, length(beta)), names(beta))
  contrast[names(weights)] <- weights
  estimate <- sum(contrast * beta)
  std_error <- sqrt(drop(t(contrast) %*% variance %*% contrast))
  z_value <- estimate / std_error
  tibble(
    effect,
    other_predator_z,
    estimate_log_odds = estimate,
    std_error,
    odds_ratio = exp(estimate),
    conf_low = exp(estimate - 1.96 * std_error),
    conf_high = exp(estimate + 1.96 * std_error),
    z_value,
    p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE)
  )
}

simple_effects <- bind_rows(
  lapply(c(-1, 0, 1), function(jackal) {
    linear_contrast(
      interaction_fit,
      c(wolf_z = 1, `wolf_z:jackal_z` = jackal),
      "Wolf slope",
      jackal
    )
  }),
  lapply(c(-1, 0, 1), function(wolf) {
    linear_contrast(
      interaction_fit,
      c(jackal_z = 1, `wolf_z:jackal_z` = wolf),
      "Jackal slope",
      wolf
    )
  })
)
write_csv(simple_effects, paste0(output_prefix, "_simple_effects.csv"))

# Standardize predictions over the observed camera-design distribution while
# holding landscape and human covariates at their standardized means.
model_data <- model.frame(interaction_fit)
reference_data <- model_data %>%
  select(
    path_type,
    camera_model,
    vegetation_structure,
    area_behind_camera,
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

fixed_formula <- ~ wolf_z * jackal_z + open_z + built_z + reserve_z + human_z +
  path_type + camera_model + vegetation_structure + area_behind_camera
beta_hat <- fixef(interaction_fit)$cond
beta_vcov <- as.matrix(vcov(interaction_fit)$cond)
set.seed(seed)
beta_draws <- MASS::mvrnorm(n_draws, mu = beta_hat, Sigma = beta_vcov)
if (is.vector(beta_draws)) beta_draws <- matrix(beta_draws, nrow = 1)
colnames(beta_draws) <- names(beta_hat)

prediction_grid <- expand.grid(
  jackal_z = seq(-1, 1, length.out = 101),
  wolf_z = c(-1, 0, 1),
  KEEP.OUT.ATTRS = FALSE
) %>%
  as_tibble()

curve_predictions <- bind_rows(lapply(seq_len(nrow(prediction_grid)), function(i) {
  newdata <- reference_data
  newdata$wolf_z <- prediction_grid$wolf_z[i]
  newdata$jackal_z <- prediction_grid$jackal_z[i]
  design <- model.matrix(fixed_formula, data = newdata)
  design <- design[, names(beta_hat), drop = FALSE]
  point_probability <- mean(plogis(as.numeric(design %*% beta_hat)))
  probability_draws <- rowMeans(plogis(beta_draws %*% t(design)))
  ratio_draws <- 100 * probability_draws / (1 - probability_draws)
  tibble(
    wolf_z = prediction_grid$wolf_z[i],
    jackal_z = prediction_grid$jackal_z[i],
    estimated_fawns_per_100_females = 100 * point_probability / (1 - point_probability),
    lower = unname(quantile(ratio_draws, 0.025)),
    upper = unname(quantile(ratio_draws, 0.975))
  )
})) %>%
  mutate(
    wolf_context = factor(
      wolf_z,
      levels = c(-1, 0, 1),
      labels = c("Low wolf (-1 SD)", "Mean wolf", "High wolf (+1 SD)")
    )
  )
write_csv(curve_predictions, paste0(output_prefix, "_curve_predictions.csv"))

palette <- c(
  "Low wolf (-1 SD)" = "#17807E",
  "Mean wolf" = "#D67A26",
  "High wolf (+1 SD)" = "#8B2F38"
)

curve_plot <- ggplot(
  curve_predictions,
  aes(
    x = jackal_z,
    y = estimated_fawns_per_100_females,
    color = wolf_context,
    fill = wolf_context
  )
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA) +
  geom_line(linewidth = 1.15) +
  geom_vline(xintercept = c(-1, 0, 1), color = "#859092", linewidth = 0.35, linetype = "dotted") +
  scale_color_manual(values = palette, name = NULL) +
  scale_fill_manual(values = palette, name = NULL) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(0, NA)) +
  labs(
    x = "Jackal abundance (z)",
    y = "Estimated fawns per 100 females",
    title = "Pooled fawn representation across wolf-jackal context",
    subtitle = paste0(
      "Full-survey beta-binomial model; bands are 95% model-based intervals\n",
      "Interaction p = 0.578; additive wolf p = 0.617; additive jackal p = 0.874"
    ),
    caption = paste0(
      "Predictions are standardized over camera-design covariates; landscape and human covariates are held at their means.\n",
      "This is an apparent composition index, not recruitment or survival."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#293536"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#4E5A5C", lineheight = 1.05),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "top"
  )

ggsave(paste0(output_prefix, "_interaction_curves.png"), curve_plot, width = 8.4, height = 5.7, dpi = 320)
ggsave(paste0(output_prefix, "_interaction_curves.pdf"), curve_plot, width = 8.4, height = 5.7)

# Use the supported additive structure for the primary pooled heatmap.
additive_formula <- ~ wolf_z + jackal_z + open_z + built_z + reserve_z + human_z +
  path_type + camera_model + vegetation_structure + area_behind_camera
additive_beta <- fixef(additive_fit)$cond
additive_vcov <- as.matrix(vcov(additive_fit)$cond)
set.seed(seed + 1L)
additive_draws <- MASS::mvrnorm(n_draws, mu = additive_beta, Sigma = additive_vcov)
if (is.vector(additive_draws)) additive_draws <- matrix(additive_draws, nrow = 1)
colnames(additive_draws) <- names(additive_beta)

heatmap_grid <- expand.grid(
  wolf_z = c(-1, 0, 1),
  jackal_z = c(-1, 0, 1),
  KEEP.OUT.ATTRS = FALSE
)
heatmap_data <- bind_rows(lapply(seq_len(nrow(heatmap_grid)), function(i) {
  newdata <- reference_data
  newdata$wolf_z <- heatmap_grid$wolf_z[i]
  newdata$jackal_z <- heatmap_grid$jackal_z[i]
  design <- model.matrix(additive_formula, data = newdata)
  design <- design[, names(additive_beta), drop = FALSE]
  point_probability <- mean(plogis(as.numeric(design %*% additive_beta)))
  probability_draws <- rowMeans(plogis(additive_draws %*% t(design)))
  ratio_draws <- 100 * probability_draws / (1 - probability_draws)
  tibble(
    wolf_z = heatmap_grid$wolf_z[i],
    jackal_z = heatmap_grid$jackal_z[i],
    estimated_fawns_per_100_females = 100 * point_probability / (1 - point_probability),
    lower = unname(quantile(ratio_draws, 0.025)),
    upper = unname(quantile(ratio_draws, 0.975))
  )
})) %>%
  mutate(
    wolf_level = factor(wolf_z, levels = c(-1, 0, 1)),
    jackal_level = factor(jackal_z, levels = c(-1, 0, 1)),
    label = sprintf(
      "%.1f\n[%.1f, %.1f]",
      estimated_fawns_per_100_females,
      lower,
      upper
    ),
    text_color = if_else(
      estimated_fawns_per_100_females > median(estimated_fawns_per_100_females),
      "white",
      "#263238"
    )
  )
write_csv(heatmap_data, paste0(output_prefix, "_heatmap_predictions.csv"))

heatmap_plot <- ggplot(
  heatmap_data,
  aes(x = jackal_level, y = wolf_level, fill = estimated_fawns_per_100_females)
) +
  geom_tile(color = "#FBFAF7", linewidth = 1.4) +
  geom_text(aes(label = label, color = text_color), size = 3.7, lineheight = 0.9, fontface = "bold") +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = c("#F7F0DF", "#E7C283", "#D17443", "#7E231F"),
    name = "Fawns per\n100 females"
  ) +
  coord_equal() +
  labs(
    x = "Jackal abundance (z)",
    y = "Wolf abundance (z)",
    title = "Pooled additive fawn:female composition",
    subtitle = "Preferred season-free model; cells show fawns per 100 females [95% model-based interval]",
    caption = "Wolf OR = 0.93 [0.70, 1.23]; jackal OR = 0.97 [0.65, 1.44]."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "right"
  )

ggsave(paste0(output_prefix, "_heatmap.png"), heatmap_plot, width = 7.3, height = 5.7, dpi = 320)
ggsave(paste0(output_prefix, "_heatmap.pdf"), heatmap_plot, width = 7.3, height = 5.7)

fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt_p <- function(x) ifelse(x < 0.001, "<0.001", fmt(x, 3))
predator_terms <- coefficient_df %>%
  filter(
    (model == "Pooled interaction" & term %in% c("wolf_z", "jackal_z", "wolf_z:jackal_z")) |
      (model == "Pooled additive" & term %in% c("wolf_z", "jackal_z"))
  )

md <- c(
  "# Pooled fawn:female beta-binomial analysis",
  "",
  "## Analysis",
  "",
  "Season was omitted to estimate a single full-survey predator association. The response remained at camera-round resolution rather than summing counts by camera, preserving repeated observations and their denominators. A nested polygon/camera random intercept accounted for repeated cameras, and the beta-binomial distribution accounted for extra-binomial variation.",
  "",
  "Both models adjusted for standardized open-space, built-and-road, reserve, and human covariates plus path type, camera model, vegetation structure, and area behind the camera. Continuous predator effects are per 1-SD increase.",
  "",
  "## Model comparison",
  "",
  "| Comparison | Additive AIC | Interaction AIC | Delta AIC | LR chi-square | p |",
  "|---|---:|---:|---:|---:|---:|",
  sprintf(
    "| Wolf x jackal interaction | `%s` | `%s` | `%s` | `%s` | `%s` |",
    fmt(model_comparison$additive_AIC),
    fmt(model_comparison$interaction_AIC),
    fmt(model_comparison$delta_AIC_interaction_minus_additive),
    fmt(model_comparison$likelihood_ratio),
    fmt_p(model_comparison$p_value)
  ),
  "",
  sprintf(
    "Adding the common Early-Late shift improved the additive model (pooled AIC `%s`; season-adjusted AIC `%s`; LR chi-square `%s`, p `%s`).",
    fmt(season_comparison$pooled_AIC),
    fmt(season_comparison$season_adjusted_AIC),
    fmt(season_comparison$likelihood_ratio),
    fmt_p(season_comparison$p_value)
  ),
  "",
  "## Predator coefficients",
  "",
  "| Model | Term | Odds ratio [95% CI] | p |",
  "|---|---|---:|---:|"
)
for (i in seq_len(nrow(predator_terms))) {
  md <- c(md, sprintf(
    "| %s | `%s` | `%s [%s, %s]` | `%s` |",
    predator_terms$model[i],
    predator_terms$term[i],
    fmt(predator_terms$odds_ratio[i]),
    fmt(predator_terms$conf_low[i]),
    fmt(predator_terms$conf_high[i]),
    fmt_p(predator_terms$p_value[i])
  ))
}

md <- c(
  md,
  "",
  "## Interaction-model simple slopes",
  "",
  "| Effect | Other predator z | Odds ratio [95% CI] | p |",
  "|---|---:|---:|---:|"
)
for (i in seq_len(nrow(simple_effects))) {
  md <- c(md, sprintf(
    "| %s | `%s` | `%s [%s, %s]` | `%s` |",
    simple_effects$effect[i],
    simple_effects$other_predator_z[i],
    fmt(simple_effects$odds_ratio[i]),
    fmt(simple_effects$conf_low[i]),
    fmt(simple_effects$conf_high[i]),
    fmt_p(simple_effects$p_value[i])
  ))
}

md <- c(
  md,
  "",
  "## Interpretation",
  "",
  "The interaction model was not supported over the additive model. Across the full survey, neither predator had a clear association with fawn representation, whether the other predator was low, average, or high. The weak visual tendency toward lower fawn representation when both predators were high is uncertain and should not be interpreted as an additive or synergistic predation effect.",
  "",
  "Removing season estimates a weighted average across periods. Because fawn representation increased later in the survey in the season-adjusted analysis, the pooled model answers a narrower spatial question but conceals that temporal shift. It is useful as a sensitivity analysis, not preferable to the season-adjusted model for primary demographic inference.",
  "",
  "As in the seasonal analysis, the response is an apparent fawn:female camera composition index. Predator coefficients can reflect differences in fawn detections, female detections, or residual class-specific detectability; they do not directly estimate predation, recruitment, or survival."
)
writeLines(md, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved pooled beta-binomial outputs with prefix:", output_prefix, "\n")
