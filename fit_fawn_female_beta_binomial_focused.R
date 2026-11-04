#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

source_fit_file <- "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
output_prefix <- "publication_fawn_female_beta_binomial_focused"

if (!file.exists(source_fit_file)) {
  stop("Missing source fit file: ", source_fit_file)
}

source_fits <- readRDS(source_fit_file)
if (!"overall" %in% names(source_fits) ||
    !inherits(source_fits[["overall"]], "glmmTMB")) {
  stop("The source beta-binomial fit was not found.")
}

# Recover the exact complete-case data used by the original model so the
# revised fixed effects retain the same response and sampling hierarchy.
source_frame <- model.frame(source_fits[["overall"]])
response <- model.response(source_frame)
if (!is.matrix(response) || ncol(response) != 2L) {
  stop("Expected a two-column cbind(fawn, female) response.")
}

model_data <- tibble(
  fawn = as.numeric(response[, 1]),
  female = as.numeric(response[, 2]),
  wolf_z = source_frame[["wolf_z"]],
  jackal_z = source_frame[["jackal_z"]],
  open_z = source_frame[["open_z"]],
  Polygon.Name = droplevels(factor(source_frame[["Polygon.Name"]])),
  camera_id = droplevels(factor(source_frame[["camera_id"]]))
)

if (any(!complete.cases(model_data))) {
  stop("The recovered focused-model data contain missing values.")
}
if (any(model_data$fawn + model_data$female <= 0)) {
  stop("The focused-model data contain uninformative zero-denominator rows.")
}

control <- glmmTMBControl(
  optCtrl = list(iter.max = 10000, eval.max = 10000),
  profile = TRUE
)

focused_interaction <- glmmTMB(
  cbind(fawn, female) ~ wolf_z * jackal_z + open_z +
    (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)

focused_additive <- glmmTMB(
  cbind(fawn, female) ~ wolf_z + jackal_z + open_z +
    (1 | Polygon.Name / camera_id),
  family = betabinomial(link = "logit"),
  data = model_data,
  control = control
)

fits <- list(
  focused_interaction = focused_interaction,
  focused_additive = focused_additive
)
saveRDS(fits, paste0(output_prefix, "_model_fits.rds"))

extract_coefficients <- function(fit, model) {
  coefficient_matrix <- summary(fit)$coefficients$cond
  as.data.frame(coefficient_matrix) %>%
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

coefficient_table <- bind_rows(
  extract_coefficients(focused_interaction, "Focused interaction"),
  extract_coefficients(focused_additive, "Focused additive")
)
write_csv(coefficient_table, paste0(output_prefix, "_coefficients.csv"))

fit_diagnostics <- bind_rows(lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  tibble(
    model = model_name,
    n_observations = nobs(fit),
    n_cameras = n_distinct(model_data$camera_id),
    n_polygons = n_distinct(model_data$Polygon.Name),
    AIC = AIC(fit),
    log_likelihood = as.numeric(logLik(fit)),
    convergence_code = fit$fit$convergence,
    positive_definite_hessian = isTRUE(fit$sdr$pdHess),
    dispersion = sigma(fit)
  )
}))
write_csv(fit_diagnostics, paste0(output_prefix, "_fit_diagnostics.csv"))

lrt_table <- as.data.frame(anova(focused_additive, focused_interaction))
interaction_test <- tibble(
  comparison = "Wolf x jackal interaction",
  additive_AIC = AIC(focused_additive),
  interaction_AIC = AIC(focused_interaction),
  delta_AIC_interaction_minus_additive =
    AIC(focused_interaction) - AIC(focused_additive),
  df_difference = lrt_table$Df[2] - lrt_table$Df[1],
  likelihood_ratio = lrt_table$Chisq[2],
  p_value = lrt_table$`Pr(>Chisq)`[2]
)
write_csv(interaction_test, paste0(output_prefix, "_interaction_test.csv"))

linear_contrast <- function(fit, weights, effect, context_label, context_z) {
  beta <- fixef(fit)$cond
  variance <- as.matrix(vcov(fit)$cond)
  contrast <- setNames(rep(0, length(beta)), names(beta))
  contrast[names(weights)] <- weights
  estimate <- sum(contrast * beta)
  std_error <- sqrt(drop(t(contrast) %*% variance %*% contrast))
  z_value <- estimate / std_error
  tibble(
    effect,
    context = context_label,
    other_predator_z = context_z,
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
  lapply(c(-1, 0, 1), function(jackal_z) {
    linear_contrast(
      focused_interaction,
      c(wolf_z = 1, `wolf_z:jackal_z` = jackal_z),
      "Wolf slope",
      paste0("Jackal z = ", jackal_z),
      jackal_z
    )
  }),
  lapply(c(-1, 0, 1), function(wolf_z) {
    linear_contrast(
      focused_interaction,
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

primary_coefficients <- coefficient_table %>%
  filter(model == "Focused interaction", term != "(Intercept)")

report <- c(
  "# Focused beta-binomial fawn:female composition model",
  "",
  "## Model",
  "",
  "The camera-round response was `cbind(fawn, female)`. The only fixed-effect predictors were standardized wolf activity, standardized jackal activity, their interaction, and standardized open-space cover. Nested polygon/camera random intercepts were retained to account for repeated camera-rounds.",
  "",
  paste0(
    "The analysis used `", nrow(model_data), "` informative camera-rounds from `",
    n_distinct(model_data$camera_id), "` cameras in `",
    n_distinct(model_data$Polygon.Name), "` polygons."
  ),
  "",
  "## Coefficients",
  "",
  "| Term | Odds ratio [95% CI] | p |",
  "|---|---:|---:|"
)

for (i in seq_len(nrow(primary_coefficients))) {
  report <- c(report, sprintf(
    "| `%s` | `%s [%s, %s]` | `%s` |",
    primary_coefficients$term[i],
    fmt(primary_coefficients$odds_ratio[i]),
    fmt(primary_coefficients$conf_low[i]),
    fmt(primary_coefficients$conf_high[i]),
    fmt_p(primary_coefficients$p_value[i])
  ))
}

report <- c(
  report,
  "",
  "## Interaction test",
  "",
  paste0(
    "Compared with the focused additive model, adding the wolf x jackal interaction changed AIC from `",
    fmt(interaction_test$additive_AIC), "` to `",
    fmt(interaction_test$interaction_AIC), "` (LR chi-square `",
    fmt(interaction_test$likelihood_ratio), "`, p `",
    fmt_p(interaction_test$p_value), "`)."
  ),
  "",
  "## Interpretation constraint",
  "",
  "With an interaction in the model, the wolf coefficient is the wolf association when standardized jackal activity equals zero, and the jackal coefficient is the jackal association when standardized wolf activity equals zero. The interaction coefficient estimates how each predator slope changes with a 1-SD increase in the other predator. Unsupported intervals should not be interpreted as biological effects.",
  "",
  "The response remains an apparent camera-based fawn representation index, not recruitment, survival, or a population-level demographic ratio."
)

writeLines(report, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved focused beta-binomial outputs with prefix:", output_prefix, "\n")
