#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

fit_file <- "publication_wolf_jackal_interaction_nmix_spabundance_fits.rds"
figure_grid_file <- "publication_spabundance_wolf_jackal_relative_prediction_grid.csv"
prediction_file <- "publication_spabundance_four_predator_combinations.csv"
contrast_file <- "publication_spabundance_four_predator_contrasts.csv"
weakening_file <- "publication_spabundance_jackal_weakening_contrasts.csv"
report_file <- "publication_spabundance_four_predator_combinations_results.md"

fits <- readRDS(fit_file)
expected_classes <- c("fawn", "female", "male")
if (!identical(names(fits), expected_classes)) {
  stop("Saved fit classes changed: ", paste(names(fits), collapse = ", "))
}

class_labels <- c(
  fawn = "Juvenile",
  female = "Adult female",
  male = "Adult male"
)
figure_class_labels <- c(
  fawn = "Fawn",
  female = "Female",
  male = "Male"
)

combinations <- tibble::tribble(
  ~combination_id, ~predator_combination, ~wolf_z, ~jackal_z,
  "low_low", "Low wolf / Low jackal", -1, -1,
  "low_wolf_high_jackal", "Low wolf / High jackal", -1, 1,
  "high_wolf_low_jackal", "High wolf / Low jackal", 1, -1,
  "high_high", "High wolf / High jackal", 1, 1
)

contrast_definitions <- tibble::tribble(
  ~contrast_id, ~contrast, ~numerator_id, ~denominator_id,
  "low_low_vs_high_high",
  "Low wolf / Low jackal vs High wolf / High jackal",
  "low_low", "high_high",
  "low_wolf_high_jackal_vs_high_wolf_high_jackal",
  "Low wolf / High jackal vs High wolf / High jackal",
  "low_wolf_high_jackal", "high_high",
  "low_vs_high_jackal_at_low_wolf",
  "Low jackal vs High jackal at low wolf",
  "low_low", "low_wolf_high_jackal",
  "low_vs_high_jackal_at_high_wolf",
  "Low jackal vs High jackal at high wolf",
  "high_wolf_low_jackal", "high_high",
  "low_vs_high_wolf_at_low_jackal",
  "Low wolf vs High wolf at low jackal",
  "low_low", "high_wolf_low_jackal",
  "low_vs_high_wolf_at_high_jackal",
  "Low wolf vs High wolf at high jackal",
  "low_wolf_high_jackal", "high_high"
)

coefficient_names <- c(
  "Null.Wolf.Abundance.z",
  "null.jackal.abundance.z",
  "Null.Wolf.Abundance.z:null.jackal.abundance.z"
)

prediction_rows <- list()
contrast_rows <- list()
draws_by_class <- list()

for (class_key in expected_classes) {
  beta_all <- as.matrix(fits[[class_key]]$beta.samples)
  if (nrow(beta_all) != 13500L || !all(coefficient_names %in% colnames(beta_all))) {
    stop("Unexpected posterior samples for class: ", class_key)
  }
  beta <- beta_all[, coefficient_names, drop = FALSE]
  design <- cbind(
    combinations$wolf_z,
    combinations$jackal_z,
    combinations$wolf_z * combinations$jackal_z
  )
  relative_draws <- exp(beta %*% t(design))
  colnames(relative_draws) <- combinations$combination_id
  draws_by_class[[class_key]] <- relative_draws

  prediction_rows[[class_key]] <- combinations %>%
    mutate(
      demographic_class = unname(class_labels[class_key]),
      posterior_summary = "Median",
      relative_estimate = apply(relative_draws, 2, median),
      lower_95_CrI = apply(relative_draws, 2, quantile, 0.025),
      upper_95_CrI = apply(relative_draws, 2, quantile, 0.975),
      relative_to_mean_predator_reference = TRUE,
      reference_condition = "wolf_z = 0, jackal_z = 0",
      retained_posterior_draws = nrow(relative_draws)
    ) %>%
    select(
      demographic_class, combination_id, predator_combination,
      wolf_z, jackal_z, posterior_summary, relative_estimate,
      lower_95_CrI, upper_95_CrI,
      relative_to_mean_predator_reference, reference_condition,
      retained_posterior_draws
    )

  contrast_rows[[class_key]] <- bind_rows(lapply(
    seq_len(nrow(contrast_definitions)),
    function(i) {
      definition <- contrast_definitions[i, ]
      ratio_draws <-
        relative_draws[, definition$numerator_id] /
        relative_draws[, definition$denominator_id]
      tibble(
        demographic_class = unname(class_labels[class_key]),
        contrast_id = definition$contrast_id,
        contrast = definition$contrast,
        ratio_orientation = "first condition / second condition",
        abundance_ratio = median(ratio_draws),
        lower_95_CrI = unname(quantile(ratio_draws, 0.025)),
        upper_95_CrI = unname(quantile(ratio_draws, 0.975)),
        posterior_probability_ratio_gt_1 = mean(ratio_draws > 1),
        retained_posterior_draws = length(ratio_draws)
      )
    }
  ))
}

predictions <- bind_rows(prediction_rows)
contrasts <- bind_rows(contrast_rows)
weakening <- bind_rows(lapply(expected_classes, function(class_key) {
  draws <- draws_by_class[[class_key]]
  jackal_ratio_low_wolf <-
    draws[, "low_low"] / draws[, "low_wolf_high_jackal"]
  jackal_ratio_high_wolf <-
    draws[, "high_wolf_low_jackal"] / draws[, "high_high"]
  ratio_of_ratios <- jackal_ratio_low_wolf / jackal_ratio_high_wolf
  ratio_summary <- quantile(ratio_of_ratios, c(0.025, 0.5, 0.975))
  probability_gt_1 <- mean(ratio_of_ratios > 1)
  tibble(
    demographic_class = unname(class_labels[class_key]),
    contrast = paste(
      "Jackal low/high ratio at low wolf divided by",
      "jackal low/high ratio at high wolf"
    ),
    interpretation = paste(
      "Values >1 indicate that the negative jackal association",
      "is weaker at high wolf activity"
    ),
    ratio_of_ratios = unname(ratio_summary[2]),
    lower_95_CrI = unname(ratio_summary[1]),
    upper_95_CrI = unname(ratio_summary[3]),
    posterior_probability_gt_1 = probability_gt_1
  )
}))

# Confirm exact agreement with the four matching cells in the final figure data.
figure_grid <- read_csv(figure_grid_file, show_col_types = FALSE) %>%
  filter(abs(wolf_z) == 1, abs(jackal_z) == 1) %>%
  transmute(
    figure_class = class,
    wolf_z,
    jackal_z,
    figure_estimate = estimate,
    figure_lower = conf_low,
    figure_upper = conf_high
  )
validation <- predictions %>%
  mutate(
    figure_class = unname(figure_class_labels[
      match(demographic_class, unname(class_labels))
    ])
  ) %>%
  left_join(figure_grid, by = c("figure_class", "wolf_z", "jackal_z"))
maximum_difference <- max(abs(c(
  validation$relative_estimate - validation$figure_estimate,
  validation$lower_95_CrI - validation$figure_lower,
  validation$upper_95_CrI - validation$figure_upper
)))
if (!is.finite(maximum_difference) || maximum_difference > 1e-12) {
  stop("Extracted predictions do not reproduce the final figure grid.")
}

write_csv(predictions, prediction_file)
write_csv(contrasts, contrast_file)
write_csv(weakening, weakening_file)

prediction_markdown <- c(
  "| Class | Predator combination | Median relative abundance | 95% CrI |",
  "|---|---|---:|---:|",
  vapply(seq_len(nrow(predictions)), function(i) sprintf(
    "| %s | %s | %.3f | %.3f-%.3f |",
    predictions$demographic_class[i], predictions$predator_combination[i],
    predictions$relative_estimate[i], predictions$lower_95_CrI[i],
    predictions$upper_95_CrI[i]
  ), character(1))
)

contrast_markdown <- c(
  "| Class | Contrast (first / second) | Median ratio | 95% CrI | Pr(ratio > 1) |",
  "|---|---|---:|---:|---:|",
  vapply(seq_len(nrow(contrasts)), function(i) sprintf(
    "| %s | %s | %.3f | %.3f-%.3f | %.3f |",
    contrasts$demographic_class[i], contrasts$contrast[i],
    contrasts$abundance_ratio[i], contrasts$lower_95_CrI[i],
    contrasts$upper_95_CrI[i],
    contrasts$posterior_probability_ratio_gt_1[i]
  ), character(1))
)
weakening_markdown <- c(
  "| Class | Jackal-weakening ratio of ratios | 95% CrI | Pr(>1) |",
  "|---|---:|---:|---:|",
  vapply(seq_len(nrow(weakening)), function(i) sprintf(
    "| %s | %.3f | %.3f-%.3f | %.3f |",
    weakening$demographic_class[i], weakening$ratio_of_ratios[i],
    weakening$lower_95_CrI[i], weakening$upper_95_CrI[i],
    weakening$posterior_probability_gt_1[i]
  ), character(1))
)

highest_cells <- predictions %>%
  mutate(
    demographic_class = factor(
      demographic_class, levels = unname(class_labels)
    )
  ) %>%
  group_by(demographic_class) %>%
  slice_max(relative_estimate, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(demographic_class)

high_high_comparison <- contrasts %>%
  filter(
    contrast_id ==
      "low_wolf_high_jackal_vs_high_wolf_high_jackal"
  ) %>%
  mutate(
    high_high_supported_greater = upper_95_CrI < 1
  )
high_high_supported_classes <- high_high_comparison$demographic_class[
  high_high_comparison$high_high_supported_greater
]
high_high_unsupported_classes <- high_high_comparison$demographic_class[
  !high_high_comparison$high_high_supported_greater
]

report <- c(
  "# Final spAbundance predictions at four predator combinations",
  "",
  paste(
    "Predictions use posterior medians and equal-tail 95% credible intervals",
    "from the 13,500 retained beta-coefficient draws in each final saved",
    "class-specific model. Values exclude the nonconverged abundance intercept",
    "and are relative to the mean-predator reference (wolf_z = 0, jackal_z = 0)."
  ),
  "",
  "## Predictions",
  "",
  prediction_markdown,
  "",
  "## Posterior ratios",
  "",
  contrast_markdown,
  "",
  "The second and sixth requested contrasts are mathematically identical; both",
  "are retained under the requested labels.",
  "",
  "## Direct test of jackal-association weakening",
  "",
  weakening_markdown,
  "",
  "## Factual summary",
  "",
  paste0(
    "The highest posterior-median cell was low wolf/low jackal for ",
    paste(
      highest_cells$demographic_class[
        highest_cells$combination_id == "low_low"
      ],
      collapse = " and "
    ),
    ", but not for ",
    paste(
      highest_cells$demographic_class[
        highest_cells$combination_id != "low_low"
      ],
      collapse = " and "
    ),
    ". No class was directly supported as the global maximum because at least ",
    "one required 95% pairwise ratio interval included one."
  ),
  paste(
    "High wolf/high jackal had a higher posterior median than low wolf/high",
    "jackal in all three classes. The draw-level ratio supported this direction",
    "for", paste(high_high_supported_classes, collapse = " and "),
    "but not for", paste(high_high_unsupported_classes, collapse = " and "),
    "because its 95% CrI included one."
  ),
  paste(
    "The negative jackal association was weaker at high wolf activity in all",
    "three point estimates. The draw-level ratio-of-ratios 95% CrI directly",
    "supported this weakening for adult females and adult males, but not for",
    "juveniles."
  ),
  paste(
    "Cell rankings alone do not establish pairwise differences. Direct support",
    "is determined from the draw-level ratio intervals in the contrast table."
  ),
  "",
  sprintf(
    "Validation: maximum absolute difference from the final figure grid = %.3g.",
    maximum_difference
  )
)
writeLines(report, report_file)

message("Saved: ", prediction_file)
message("Saved: ", contrast_file)
message("Saved: ", weakening_file)
message("Saved: ", report_file)
