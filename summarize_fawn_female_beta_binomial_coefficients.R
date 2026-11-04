#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(readr)
})

fit_file <- "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
output_prefix <- "publication_fawn_female_beta_binomial_coefficient"
fits <- readRDS(fit_file)
lrt_results <- read_csv(
  "publication_fawn_female_beta_binomial_wolf_jackal_likelihood_ratio_tests.csv",
  show_col_types = FALSE
)

fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt_p <- function(x) ifelse(x < 0.001, "<0.001", fmt(x, 3))

term_group <- function(term) {
  case_when(
    term == "(Intercept)" ~ "Intercept",
    grepl(":", term) ~ "Predator-season interactions",
    term %in% c("wolf_z", "jackal_z", "seasonLate (mild-cold)") ~ "Predators and season",
    term %in% c("open_z", "built_z", "reserve_z", "human_z") ~ "Landscape and human activity",
    grepl("^path_type", term) ~ "Path type",
    grepl("^camera_model", term) ~ "Camera model",
    grepl("^vegetation_structure", term) ~ "Vegetation structure",
    grepl("^area_behind_camera", term) ~ "Area behind camera",
    TRUE ~ "Other"
  )
}

reference_for_term <- function(term) {
  case_when(
    grepl("^path_type", term) ~ "4x4 road with off-road vehicles",
    grepl("^camera_model", term) ~ "Browning DARK OPS HD MAX",
    grepl("^vegetation_structure", term) ~ "Grassland",
    grepl("^area_behind_camera", term) ~ "Clear",
    grepl("season", term) ~ "Early (hotter)",
    TRUE ~ NA_character_
  )
}

extract_coefficients <- function(fit, model) {
  tab <- as.data.frame(summary(fit)$coefficients$cond)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  as_tibble(tab) %>%
    transmute(
      model,
      group = term_group(term),
      term,
      reference = reference_for_term(term),
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
  extract_coefficients(fits$overall, "Overall"),
  extract_coefficients(fits$seasonal, "Seasonal interaction"),
  extract_coefficients(fits$seasonal_additive_stable, "Preferred additive-season")
)
write_csv(coefficient_df, paste0(output_prefix, "_details.csv"))

linear_contrast <- function(fit, weights, label, contrast_type, season = NA_character_, context = NA_real_) {
  beta <- glmmTMB::fixef(fit)$cond
  vc <- as.matrix(stats::vcov(fit)$cond)
  L <- setNames(rep(0, length(beta)), names(beta))
  missing_terms <- setdiff(names(weights), names(L))
  if (length(missing_terms)) stop("Unknown contrast terms: ", paste(missing_terms, collapse = ", "))
  L[names(weights)] <- weights
  estimate <- sum(L * beta)
  std_error <- sqrt(drop(t(L) %*% vc %*% L))
  z_value <- estimate / std_error
  tibble(
    contrast_type,
    label,
    season,
    other_predator_z = context,
    estimate_log_odds = estimate,
    std_error,
    odds_ratio = exp(estimate),
    conf_low = exp(estimate - 1.96 * std_error),
    conf_high = exp(estimate + 1.96 * std_error),
    z_value,
    p_value = 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
  )
}

simple_effects <- list()
idx <- 1L

for (jackal in c(-1, 0, 1)) {
  simple_effects[[idx]] <- linear_contrast(
    fits$overall,
    c(wolf_z = 1, `wolf_z:jackal_z` = jackal),
    paste0("Wolf effect at jackal z = ", jackal),
    "Wolf slope",
    "Overall",
    jackal
  )
  idx <- idx + 1L
}
for (wolf in c(-1, 0, 1)) {
  simple_effects[[idx]] <- linear_contrast(
    fits$overall,
    c(jackal_z = 1, `wolf_z:jackal_z` = wolf),
    paste0("Jackal effect at wolf z = ", wolf),
    "Jackal slope",
    "Overall",
    wolf
  )
  idx <- idx + 1L
}

for (season in c("Early (hotter)", "Late (mild-cold)")) {
  late <- season == "Late (mild-cold)"
  for (jackal in c(-1, 0, 1)) {
    weights <- c(wolf_z = 1, `wolf_z:jackal_z` = jackal)
    if (late) {
      weights <- c(
        weights,
        `wolf_z:seasonLate (mild-cold)` = 1,
        `wolf_z:jackal_z:seasonLate (mild-cold)` = jackal
      )
    }
    simple_effects[[idx]] <- linear_contrast(
      fits$seasonal,
      weights,
      paste0("Wolf effect at jackal z = ", jackal),
      "Wolf slope",
      season,
      jackal
    )
    idx <- idx + 1L
  }
  for (wolf in c(-1, 0, 1)) {
    weights <- c(jackal_z = 1, `wolf_z:jackal_z` = wolf)
    if (late) {
      weights <- c(
        weights,
        `jackal_z:seasonLate (mild-cold)` = 1,
        `wolf_z:jackal_z:seasonLate (mild-cold)` = wolf
      )
    }
    simple_effects[[idx]] <- linear_contrast(
      fits$seasonal,
      weights,
      paste0("Jackal effect at wolf z = ", wolf),
      "Jackal slope",
      season,
      wolf
    )
    idx <- idx + 1L
  }
}

for (wolf in c(-1, 0, 1)) {
  for (jackal in c(-1, 0, 1)) {
    simple_effects[[idx]] <- linear_contrast(
      fits$seasonal,
      c(
        `seasonLate (mild-cold)` = 1,
        `wolf_z:seasonLate (mild-cold)` = wolf,
        `jackal_z:seasonLate (mild-cold)` = jackal,
        `wolf_z:jackal_z:seasonLate (mild-cold)` = wolf * jackal
      ),
      paste0("Late versus Early at wolf z = ", wolf, ", jackal z = ", jackal),
      "Season contrast",
      "Late versus Early",
      NA_real_
    ) %>%
      mutate(wolf_z = wolf, jackal_z = jackal)
    idx <- idx + 1L
  }
}

simple_effect_df <- bind_rows(simple_effects)
write_csv(simple_effect_df, paste0(output_prefix, "_simple_effects.csv"))

seasonal_frame <- model.frame(fits$seasonal_additive_stable)
factor_vars <- c(
  "season",
  "path_type",
  "camera_model",
  "vegetation_structure",
  "area_behind_camera"
)
factor_support <- bind_rows(lapply(factor_vars, function(variable) {
  values <- seasonal_frame[[variable]]
  camera <- seasonal_frame$camera_id
  tibble(variable, level = levels(values)) %>%
    rowwise() %>%
    mutate(
      n_camera_rounds = sum(values == level),
      n_cameras = n_distinct(camera[values == level]),
      reference = level == levels(values)[1]
    ) %>%
    ungroup()
}))
write_csv(factor_support, paste0(output_prefix, "_factor_support.csv"))

random_effect_list <- glmmTMB::VarCorr(fits$seasonal_additive_stable)$cond
random_effects <- bind_rows(lapply(names(random_effect_list), function(group) {
  component <- random_effect_list[[group]]
  std_dev <- unname(attr(component, "stddev")[1])
  tibble(group, variance = std_dev^2, std_dev)
}))
write_csv(random_effects, paste0(output_prefix, "_random_effects.csv"))

seasonal_main <- coefficient_df %>%
  filter(
    model == "Preferred additive-season",
    group %in% c(
      "Predators and season",
      "Landscape and human activity",
      "Path type",
      "Camera model",
      "Vegetation structure",
      "Area behind camera",
      "Predator-season interactions"
    )
  )

md <- c(
  "# Fawn:female beta-binomial coefficient interpretation",
  "",
  "## Response scale",
  "",
  "The response is the fawn proportion among classified fawn and adult-female counts. An odds ratio above 1 means greater fawn representation relative to females; it does not necessarily mean greater fawn abundance. Continuous covariates are standardized, so their odds ratios correspond to a 1-SD increase. Categorical effects are comparisons with the stated reference level.",
  "",
  "## Model hierarchy",
  "",
  "| Comparison | LR chi-square | df | p |",
  "|---|---:|---:|---:|"
)
for (i in seq_len(nrow(lrt_results))) {
  md <- c(md, sprintf(
    "| %s | `%s` | `%s` | `%s` |",
    lrt_results$comparison[i], fmt(lrt_results$likelihood_ratio[i]),
    lrt_results$df_difference[i], fmt_p(lrt_results$p_value[i])
  ))
}
md <- c(
  md,
  "",
  "## Preferred-model coefficients",
  "",
  "The preferred model used additive wolf and jackal slopes plus a common Early-Late shift. Models with a wolf x jackal effect or season-specific predator slopes were not supported.",
  "",
  "| Group | Term | Reference | Odds ratio [95% CI] | p |",
  "|---|---|---|---:|---:|"
)
for (i in seq_len(nrow(seasonal_main))) {
  ref <- ifelse(is.na(seasonal_main$reference[i]), "-", seasonal_main$reference[i])
  md <- c(md, sprintf(
    "| %s | `%s` | %s | `%s [%s, %s]` | `%s` |",
    seasonal_main$group[i], seasonal_main$term[i], ref,
    fmt(seasonal_main$odds_ratio[i]), fmt(seasonal_main$conf_low[i]),
    fmt(seasonal_main$conf_high[i]), fmt_p(seasonal_main$p_value[i])
  ))
}

md <- c(
  md,
  "",
  "## Conditional predator effects",
  "",
  "Each predator slope is the change associated with a 1-SD increase while holding the other predator at the stated value.",
  "",
  "| Season | Effect | Other predator z | Odds ratio [95% CI] | p |",
  "|---|---|---:|---:|---:|"
)
predator_effects <- simple_effect_df %>%
  filter(contrast_type %in% c("Wolf slope", "Jackal slope"), season != "Overall")
for (i in seq_len(nrow(predator_effects))) {
  md <- c(md, sprintf(
    "| %s | %s | `%s` | `%s [%s, %s]` | `%s` |",
    predator_effects$season[i], predator_effects$contrast_type[i],
    predator_effects$other_predator_z[i], fmt(predator_effects$odds_ratio[i]),
    fmt(predator_effects$conf_low[i]), fmt(predator_effects$conf_high[i]),
    fmt_p(predator_effects$p_value[i])
  ))
}

md <- c(
  md,
  "",
  "## Factor support",
  "",
  "| Variable | Level | Camera-rounds | Cameras | Reference |",
  "|---|---|---:|---:|---|"
)
for (i in seq_len(nrow(factor_support))) {
  md <- c(md, sprintf(
    "| %s | %s | `%s` | `%s` | %s |",
    factor_support$variable[i], factor_support$level[i],
    factor_support$n_camera_rounds[i], factor_support$n_cameras[i],
    ifelse(factor_support$reference[i], "Yes", "No")
  ))
}

md <- c(
  md,
  "",
  "## Interpretation",
  "",
  "- Open-space cover was a supported ecological covariate: a 1-SD increase was associated with approximately 55% higher odds of a classified detection being a fawn rather than an adult female.",
  "- The preferred model supported approximately 57% higher odds of fawn representation in the Late period than the Early period, after adjustment for predator, landscape, human, and camera-design variables.",
  "- No wolf or jackal effect was supported in the preferred model. In the full interaction model, their main coefficients are conditional Early-period slopes when the other predator equals its mean.",
  "- The Early wolf x jackal coefficient and its Late-period change were individually suggestive but uncertain, and the global interaction tests were not significant. Simple slopes should therefore be treated as pattern descriptions, not evidence of predator effects.",
  "- Short-shrub vegetation had a large positive coefficient relative to grassland, but only six informative camera-rounds from two cameras had this category. The estimate is sparse and should not be interpreted as a robust habitat effect.",
  "- All other camera-design intervals included an odds ratio of 1, providing no clear evidence of differential fawn-versus-female representation after adjustment.",
  "- The polygon random-intercept variance was effectively zero, while residual camera-level heterogeneity remained."
)
writeLines(md, paste0(output_prefix, "_interpretation.md"), useBytes = TRUE)

cat("Saved detailed coefficient interpretation outputs with prefix:", output_prefix, "\n")
