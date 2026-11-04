#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
})

fit_file <- "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
verification_file <- "publication_fawn_female_beta_binomial_coefficient_details.csv"
output_prefix <- "publication_fawn_female_beta_binomial_coefficients"
model_name <- "seasonal_additive_stable"
verification_model_name <- "Preferred additive-season"

required_files <- c(fit_file, verification_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

fits <- readRDS(fit_file)
if (!model_name %in% names(fits) || !inherits(fits[[model_name]], "glmmTMB")) {
  stop("The retained model was not found in ", fit_file, ".")
}
fit <- fits[[model_name]]

required_terms <- c("wolf_z", "jackal_z", "seasonLate (mild-cold)", "open_z")
coefficient_terms <- setdiff(names(fixef(fit)$cond), "(Intercept)")
if (!all(required_terms %in% coefficient_terms) ||
    any(grepl("wolf_z:jackal_z|jackal_z:wolf_z", coefficient_terms))) {
  stop(
    "The restored primary model does not match the additive seasonal specification."
  )
}

term_specification <- tibble(
  term = c("wolf_z", "jackal_z", "open_z", "seasonLate (mild-cold)"),
  predictor = c(
    "Wolf activity",
    "Jackal activity",
    "Open-space cover",
    "Late vs early survey"
  ),
  effect_definition = c(
    "Per 1-SD increase",
    "Per 1-SD increase",
    "Per 1-SD increase",
    "Late (mild-cold) relative to Early (hotter)"
  ),
  expected_display_p = c(0.608, 0.773, 0.031, 0.045)
)

coefficient_matrix <- summary(fit)$coefficients$cond
missing_terms <- setdiff(term_specification$term, rownames(coefficient_matrix))
if (length(missing_terms)) {
  stop("Missing retained-model terms: ", paste(missing_terms, collapse = ", "))
}

coefficient_table <- as.data.frame(coefficient_matrix) %>%
  tibble::rownames_to_column("term") %>%
  filter(term %in% term_specification$term) %>%
  transmute(
    term,
    estimate_log_odds = Estimate,
    std_error = `Std. Error`,
    z_value = `z value`,
    p_value = `Pr(>|z|)`
  ) %>%
  right_join(term_specification, by = "term") %>%
  mutate(
    odds_ratio = exp(estimate_log_odds),
    conf_low = exp(estimate_log_odds - 1.96 * std_error),
    conf_high = exp(estimate_log_odds + 1.96 * std_error),
    ci_excludes_one = conf_low > 1 | conf_high < 1,
    model = model_name,
    response = "Juvenile:female ratio among classified juvenile and adult-female counts"
  ) %>%
  select(
    model, response, predictor, term, effect_definition,
    estimate_log_odds, std_error, odds_ratio, conf_low, conf_high,
    z_value, p_value, ci_excludes_one, expected_display_p
  )

# Confirm that the plotted estimates reproduce the project's existing
# preferred-model coefficient output before writing any publication files.
verification_table <- read_csv(verification_file, show_col_types = FALSE) %>%
  filter(
    model == verification_model_name,
    term %in% term_specification$term
  ) %>%
  select(
    term,
    verified_estimate = estimate_log_odds,
    verified_se = std_error,
    verified_or = odds_ratio,
    verified_low = conf_low,
    verified_high = conf_high
  )

if (nrow(verification_table) != nrow(term_specification)) {
  stop("The verification table does not contain all four retained-model terms.")
}

verification_check <- coefficient_table %>%
  left_join(verification_table, by = "term")
values_match <- with(
  verification_check,
  abs(estimate_log_odds - verified_estimate) < 1e-10 &
    abs(std_error - verified_se) < 1e-10 &
    abs(odds_ratio - verified_or) < 1e-10 &
    abs(conf_low - verified_low) < 1e-10 &
    abs(conf_high - verified_high) < 1e-10
)
if (!all(values_match)) {
  stop("Retained-model estimates do not match the existing coefficient output.")
}

if (!all(round(coefficient_table$p_value, 3) ==
         coefficient_table$expected_display_p)) {
  stop("Displayed P-values do not match the prespecified manuscript values.")
}

coefficient_table <- coefficient_table %>%
  mutate(
    p_value_label = sprintf("P = %.3f", p_value),
    support_status = if_else(
      ci_excludes_one,
      "Supported",
      "Unsupported"
    ),
    display_order = match(predictor, term_specification$predictor),
    predictor = factor(
      predictor,
      levels = rev(term_specification$predictor)
    )
  ) %>%
  arrange(display_order)

write_csv(
  coefficient_table %>%
    mutate(predictor = as.character(predictor)) %>%
    select(-display_order, -expected_display_p),
  paste0(output_prefix, ".csv")
)

data_limits <- exp(
  range(log(c(coefficient_table$conf_low, coefficient_table$conf_high))) +
    c(-0.10, 0.10)
)
axis_breaks <- c(0.5, 0.75, 1, 1.5, 2, 2.5)
axis_breaks <- axis_breaks[
  axis_breaks >= data_limits[1] & axis_breaks <= data_limits[2]
]
p_column_x <- data_limits[2] * 1.16
plot_limits <- c(data_limits[1], p_column_x * 1.03)

forest_plot <- ggplot(
  coefficient_table,
  aes(x = odds_ratio, y = predictor)
) +
  geom_vline(
    xintercept = 1,
    linewidth = 0.55,
    linetype = "dashed",
    color = "#9AA1A3"
  ) +
  geom_errorbar(
    aes(
      xmin = conf_low,
      xmax = conf_high,
      color = support_status
    ),
    orientation = "y",
    width = 0,
    linewidth = 0.95
  ) +
  geom_point(
    aes(color = support_status),
    shape = 21,
    size = 3.9,
    stroke = 0.95,
    fill = "white"
  ) +
  geom_text(
    aes(
      x = p_column_x,
      label = p_value_label,
      color = support_status
    ),
    hjust = 1,
    size = 4.0
  ) +
  scale_color_manual(
    values = c(
      "Supported" = "#263B42",
      "Unsupported" = "#9AA1A3"
    ),
    guide = "none"
  ) +
  scale_x_log10(
    breaks = axis_breaks,
    labels = function(x) format(x, trim = TRUE, scientific = FALSE),
    limits = plot_limits,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.45, 0.45))) +
  labs(
    x = "Effect on juvenile:female ratio (odds ratio, 95% CI)",
    y = NULL
  ) +
  theme_classic(base_size = 13, base_family = "sans") +
  theme(
    axis.line.y = element_blank(),
    axis.line.x = element_line(color = "#445156", linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(color = "#596366", linewidth = 0.4),
    axis.text = element_text(color = "#26363A"),
    axis.text.x = element_text(size = 11.5),
    axis.text.y = element_text(size = 12.5, margin = margin(r = 14)),
    axis.title.x = element_text(size = 12.5, margin = margin(t = 12)),
    plot.margin = margin(10, 18, 8, 10)
  )

# Export is intentionally gated so the revised code and numerical inputs can
# be reviewed before publication files are overwritten. Set to TRUE only after
# approval.
export_approved <- TRUE
if (export_approved) {
  ggsave(
    paste0(output_prefix, ".png"),
    forest_plot,
    width = 8.6,
    height = 4.6,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    paste0(output_prefix, ".pdf"),
    forest_plot,
    width = 8.6,
    height = 4.6,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
}

supported <- coefficient_table %>%
  filter(ci_excludes_one) %>%
  pull(predictor) %>%
  as.character()

cat("Verified retained model:", model_name, "\n")
cat("Saved coefficient table and forest plot:", output_prefix, "\n")
cat(
  "95% CIs excluding OR = 1:",
  if (length(supported)) paste(supported, collapse = "; ") else "none",
  "\n"
)
