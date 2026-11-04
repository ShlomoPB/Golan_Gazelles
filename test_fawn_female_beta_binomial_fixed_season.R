#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

source_fit_file <- "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
output_prefix <- "publication_fawn_female_beta_binomial_fixed_season_sensitivity"

if (!file.exists(source_fit_file)) {
  stop("Missing source fit file: ", source_fit_file)
}

source_fits <- readRDS(source_fit_file)
if (!"seasonal" %in% names(source_fits)) {
  stop("The original seasonal fit was not found.")
}

source_frame <- model.frame(source_fits[["seasonal"]])
response <- model.response(source_frame)
model_data <- tibble(
  fawn = as.numeric(response[, 1]),
  female = as.numeric(response[, 2]),
  wolf_z = source_frame[["wolf_z"]],
  jackal_z = source_frame[["jackal_z"]],
  open_z = source_frame[["open_z"]],
  season = droplevels(factor(source_frame[["season"]])),
  Polygon.Name = droplevels(factor(source_frame[["Polygon.Name"]])),
  camera_id = droplevels(factor(source_frame[["camera_id"]]))
)

if (nlevels(model_data$season) != 2L) {
  stop("Expected exactly two season levels.")
}

control <- glmmTMBControl(
  optCtrl = list(iter.max = 10000, eval.max = 10000),
  profile = TRUE
)

no_season_interaction <- glmmTMB(
  cbind(fawn, female) ~ wolf_z * jackal_z + open_z +
    (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)

fixed_season_interaction <- glmmTMB(
  cbind(fawn, female) ~ wolf_z * jackal_z + open_z + season +
    (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)

fixed_season_additive <- glmmTMB(
  cbind(fawn, female) ~ wolf_z + jackal_z + open_z + season +
    (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)

fits <- list(
  no_season_interaction = no_season_interaction,
  fixed_season_interaction = fixed_season_interaction,
  fixed_season_additive = fixed_season_additive
)
saveRDS(fits, paste0(output_prefix, "_model_fits.rds"))

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

coefficients <- bind_rows(lapply(names(fits), function(model_name) {
  extract_coefficients(fits[[model_name]], model_name)
}))
write_csv(coefficients, paste0(output_prefix, "_coefficients.csv"))

fit_diagnostics <- bind_rows(lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  tibble(
    model = model_name,
    n_observations = nobs(fit),
    AIC = AIC(fit),
    log_likelihood = as.numeric(logLik(fit)),
    convergence_code = fit$fit$convergence,
    positive_definite_hessian = isTRUE(fit$sdr$pdHess),
    dispersion = sigma(fit)
  )
}))
write_csv(fit_diagnostics, paste0(output_prefix, "_fit_diagnostics.csv"))

lrt_result <- function(smaller, larger, comparison) {
  comparison_table <- as.data.frame(anova(smaller, larger))
  tibble(
    comparison,
    smaller_model = deparse1(formula(smaller)),
    larger_model = deparse1(formula(larger)),
    smaller_AIC = AIC(smaller),
    larger_AIC = AIC(larger),
    delta_AIC_larger_minus_smaller = AIC(larger) - AIC(smaller),
    df_difference = comparison_table$Df[2] - comparison_table$Df[1],
    likelihood_ratio = comparison_table$Chisq[2],
    p_value = comparison_table$`Pr(>Chisq)`[2]
  )
}

model_tests <- bind_rows(
  lrt_result(
    no_season_interaction,
    fixed_season_interaction,
    "Add fixed Early-Late season effect"
  ),
  lrt_result(
    fixed_season_additive,
    fixed_season_interaction,
    "Wolf x jackal interaction after fixed season adjustment"
  )
)
write_csv(model_tests, paste0(output_prefix, "_model_tests.csv"))

fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt_p <- function(x) ifelse(x < 0.001, "<0.001", fmt(x, 3))
reported_terms <- c(
  "wolf_z", "jackal_z", "open_z", "seasonLate (mild-cold)",
  "wolf_z:jackal_z"
)
season_coefficients <- coefficients %>%
  filter(
    model == "fixed_season_interaction",
    term %in% reported_terms
  )

report <- c(
  "# Fixed-season sensitivity for the focused fawn:female model",
  "",
  "## Rationale",
  "",
  "Season has only two deliberately defined levels, so it was evaluated as a fixed Early-Late contrast rather than a random effect. All other fixed effects remained limited to standardized wolf activity, standardized jackal activity, their interaction, and standardized open-space cover. Nested polygon/camera random intercepts were unchanged.",
  "",
  "## Fixed-season interaction model",
  "",
  "| Term | Odds ratio [95% CI] | p |",
  "|---|---:|---:|"
)

for (i in seq_len(nrow(season_coefficients))) {
  report <- c(report, sprintf(
    "| `%s` | `%s [%s, %s]` | `%s` |",
    season_coefficients$term[i],
    fmt(season_coefficients$odds_ratio[i]),
    fmt(season_coefficients$conf_low[i]),
    fmt(season_coefficients$conf_high[i]),
    fmt_p(season_coefficients$p_value[i])
  ))
}

report <- c(
  report,
  "",
  "## Nested model tests",
  "",
  "| Comparison | LR chi-square | Delta AIC (larger - smaller) | p |",
  "|---|---:|---:|---:|"
)
for (i in seq_len(nrow(model_tests))) {
  report <- c(report, sprintf(
    "| %s | `%s` | `%s` | `%s` |",
    model_tests$comparison[i],
    fmt(model_tests$likelihood_ratio[i]),
    fmt(model_tests$delta_AIC_larger_minus_smaller[i]),
    fmt_p(model_tests$p_value[i])
  ))
}

report <- c(
  report,
  "",
  "## Interpretation",
  "",
  "This sensitivity analysis tests whether a prespecified temporal shift confounds the focused ecological associations. It should not be selected merely because individual p-values cross 0.05. Effect sizes, uncertainty, the nested season test, and the ecological rationale should be considered together."
)
writeLines(report, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved fixed-season sensitivity outputs with prefix:", output_prefix, "\n")
