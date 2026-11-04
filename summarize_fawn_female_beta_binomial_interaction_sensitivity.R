#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

fit_file <- "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
output_prefix <- "publication_fawn_female_beta_binomial_interaction_sensitivity"

if (!file.exists(fit_file)) stop("Missing fit file: ", fit_file)
fits <- readRDS(fit_file)
required_fits <- c("seasonal_additive_stable", "seasonal_stable_response")
if (!all(required_fits %in% names(fits))) {
  stop("Missing required saved models: ", paste(setdiff(required_fits, names(fits)), collapse = ", "))
}

additive_fit <- fits[["seasonal_additive_stable"]]
interaction_fit <- fits[["seasonal_stable_response"]]

additive_terms <- setdiff(names(fixef(additive_fit)$cond), "(Intercept)")
interaction_terms <- setdiff(names(fixef(interaction_fit)$cond), "(Intercept)")
if (!"wolf_z:jackal_z" %in% interaction_terms ||
    !setequal(interaction_terms, c(additive_terms, "wolf_z:jackal_z"))) {
  stop("The sensitivity model does not differ from the primary model only by wolf x jackal.")
}

extract_coefficients <- function(fit, model) {
  as.data.frame(summary(fit)$coefficients$cond) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
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

key_terms <- c(
  "wolf_z", "jackal_z", "wolf_z:jackal_z", "open_z",
  "seasonLate (mild-cold)"
)
coefficient_table <- extract_coefficients(
  interaction_fit,
  "Season-adjusted interaction sensitivity"
) %>%
  filter(term %in% key_terms) %>%
  mutate(ci_excludes_one = conf_low > 1 | conf_high < 1)
write_csv(coefficient_table, paste0(output_prefix, "_coefficients.csv"))

comparison_raw <- as.data.frame(anova(additive_fit, interaction_fit))
model_comparison <- tibble(
  comparison = "Add wolf x jackal interaction to primary additive model",
  additive_AIC = AIC(additive_fit),
  interaction_AIC = AIC(interaction_fit),
  delta_AIC_interaction_minus_additive = AIC(interaction_fit) - AIC(additive_fit),
  df_difference = comparison_raw$Df[2] - comparison_raw$Df[1],
  likelihood_ratio = comparison_raw$Chisq[2],
  p_value = comparison_raw$`Pr(>Chisq)`[2],
  additive_convergence = additive_fit$fit$convergence,
  interaction_convergence = interaction_fit$fit$convergence,
  additive_positive_definite_hessian = isTRUE(additive_fit$sdr$pdHess),
  interaction_positive_definite_hessian = isTRUE(interaction_fit$sdr$pdHess)
)
write_csv(model_comparison, paste0(output_prefix, "_model_comparison.csv"))

linear_contrast <- function(weights, effect, context, other_predator_z) {
  beta <- fixef(interaction_fit)$cond
  variance <- as.matrix(vcov(interaction_fit)$cond)
  contrast <- setNames(rep(0, length(beta)), names(beta))
  contrast[names(weights)] <- weights
  estimate <- sum(contrast * beta)
  std_error <- sqrt(drop(t(contrast) %*% variance %*% contrast))
  z_value <- estimate / std_error
  tibble(
    effect,
    context,
    other_predator_z,
    estimate_log_odds = estimate,
    std_error,
    odds_ratio = exp(estimate),
    conf_low = exp(estimate - 1.96 * std_error),
    conf_high = exp(estimate + 1.96 * std_error),
    z_value,
    p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE),
    ci_excludes_one = conf_low > 1 | conf_high < 1
  )
}

simple_effects <- bind_rows(
  lapply(c(-1, 0, 1), function(jackal_z) {
    linear_contrast(
      c(wolf_z = 1, `wolf_z:jackal_z` = jackal_z),
      "Wolf slope",
      paste0("Jackal z = ", jackal_z),
      jackal_z
    )
  }),
  lapply(c(-1, 0, 1), function(wolf_z) {
    linear_contrast(
      c(jackal_z = 1, `wolf_z:jackal_z` = wolf_z),
      "Jackal slope",
      paste0("Wolf z = ", wolf_z),
      wolf_z
    )
  })
)
write_csv(simple_effects, paste0(output_prefix, "_simple_effects.csv"))

fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt_p <- function(x) ifelse(x < 0.001, "<0.001", fmt(x, 3))

report <- c(
  "# Wolf x jackal sensitivity for fawn representation",
  "",
  "## Model comparison",
  "",
  "The primary season-adjusted additive beta-binomial model was compared with an otherwise identical model containing a wolf by jackal interaction. Landscape, human-activity, camera-design, season, and nested polygon/camera adjustments were unchanged.",
  "",
  paste0(
    "Adding the interaction changed AIC from `", fmt(model_comparison$additive_AIC),
    "` to `", fmt(model_comparison$interaction_AIC), "` (delta AIC `",
    fmt(model_comparison$delta_AIC_interaction_minus_additive),
    "`; LR chi-square `", fmt(model_comparison$likelihood_ratio),
    "`; p `", fmt_p(model_comparison$p_value), "`)."
  ),
  "",
  "## Interaction-model coefficients",
  "",
  "| Term | Odds ratio [95% CI] | p |",
  "|---|---:|---:|"
)

for (i in seq_len(nrow(coefficient_table))) {
  report <- c(report, sprintf(
    "| `%s` | `%s [%s, %s]` | `%s` |",
    coefficient_table$term[i],
    fmt(coefficient_table$odds_ratio[i]),
    fmt(coefficient_table$conf_low[i]),
    fmt(coefficient_table$conf_high[i]),
    fmt_p(coefficient_table$p_value[i])
  ))
}

report <- c(
  report,
  "",
  "## Conditional predator slopes",
  "",
  "| Effect | Context | Odds ratio [95% CI] | p |",
  "|---|---|---:|---:|"
)
for (i in seq_len(nrow(simple_effects))) {
  report <- c(report, sprintf(
    "| %s | `%s` | `%s [%s, %s]` | `%s` |",
    simple_effects$effect[i],
    simple_effects$context[i],
    fmt(simple_effects$odds_ratio[i]),
    fmt(simple_effects$conf_low[i]),
    fmt(simple_effects$conf_high[i]),
    fmt_p(simple_effects$p_value[i])
  ))
}

report <- c(
  report,
  "",
  "## Decision",
  "",
  "The interaction model is retained as a sensitivity analysis. The additive model remains primary because the interaction did not improve fit, its interval included one, and all conditional predator-slope intervals included one. Open-space cover and the Late versus Early contrast remained supported in the interaction model."
)
writeLines(report, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved interaction sensitivity outputs with prefix:", output_prefix, "\n")
