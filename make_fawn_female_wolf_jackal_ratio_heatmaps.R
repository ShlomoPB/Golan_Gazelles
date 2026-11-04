#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

output_prefix <- "publication_fawn_female_ratio_wolf_jackal"
season_file <- "publication_wolf_jackal_interaction_spabundance_by_season_prediction_grid.csv"
overall_file <- "publication_wolf_jackal_interaction_nmix_prediction_grid.csv"

required_columns <- c("wolf_z", "jackal_z", "species", "lambda")

season_predictions <- read_csv(season_file, show_col_types = FALSE)
overall_predictions <- read_csv(overall_file, show_col_types = FALSE) %>%
  filter(framework == "spAbundance")

if (length(setdiff(c(required_columns, "season_period"), names(season_predictions)))) {
  stop("Seasonal prediction file is missing required columns.")
}
if (length(setdiff(c(required_columns, "framework"), names(overall_predictions)))) {
  stop("Overall prediction file is missing required columns.")
}

ratio_grid <- bind_rows(
  overall_predictions %>%
    transmute(
      wolf_z,
      jackal_z,
      period = "Overall",
      species,
      lambda,
      source_model = "Full seven-round spAbundance model"
    ),
  season_predictions %>%
    transmute(
      wolf_z,
      jackal_z,
      period = season_period,
      species,
      lambda,
      source_model = paste0(season_period, " spAbundance model")
    )
) %>%
  filter(species %in% c("Fawn", "Female")) %>%
  pivot_wider(names_from = species, values_from = lambda) %>%
  mutate(
    fawn_female_ratio = Fawn / Female,
    fawns_per_100_females = 100 * fawn_female_ratio,
    period = factor(
      period,
      levels = c("Overall", "Early (hotter)", "Late (mild-cold)")
    ),
    wolf_factor = factor(wolf_z, levels = c(-1, 0, 1)),
    jackal_factor = factor(jackal_z, levels = c(-1, 0, 1)),
    cell_label = sprintf("%.2f\n(%.1f/100)", fawn_female_ratio, fawns_per_100_females),
    text_color = if_else(fawns_per_100_females >= 22, "white", "#1d2528")
  ) %>%
  arrange(period, desc(wolf_z), jackal_z)

if (nrow(ratio_grid) != 27 || any(!is.finite(ratio_grid$fawn_female_ratio))) {
  stop("Expected 27 finite ratio cells: three periods x three wolf levels x three jackal levels.")
}

write_csv(
  ratio_grid %>%
    transmute(
      period = as.character(period),
      wolf_z,
      jackal_z,
      predicted_fawn_abundance = Fawn,
      predicted_female_abundance = Female,
      fawn_female_ratio,
      fawns_per_100_females,
      source_model
    ),
  paste0(output_prefix, "_prediction_grid.csv")
)

ratio_palette <- c("#f7f0df", "#e7c283", "#d17443", "#7e231f")
ratio_limits <- c(0, 40)

base_heatmap <- function(data, title, subtitle = NULL, facet = FALSE) {
  p <- ggplot(data, aes(jackal_factor, wolf_factor, fill = fawns_per_100_females)) +
    geom_tile(color = "#fbfaf7", linewidth = 1.5) +
    geom_text(
      aes(label = cell_label, color = text_color),
      size = if (facet) 3.25 else 4.4,
      lineheight = 0.92,
      fontface = "bold"
    ) +
    scale_color_identity() +
    scale_fill_gradientn(
      colors = ratio_palette,
      values = scales::rescale(c(0, 10, 20, 40)),
      limits = ratio_limits,
      oob = scales::squish,
      breaks = c(0, 10, 20, 30, 40),
      labels = c("0", "10", "20", "30", "40")
    ) +
    coord_equal() +
    labs(
      x = "Jackal abundance (z)",
      y = "Wolf abundance (z)",
      fill = "Fawns per\n100 females",
      title = title,
      subtitle = subtitle,
      caption = "Cell labels: fawn:female abundance ratio (fawns per 100 females)"
    ) +
    theme_minimal(base_size = if (facet) 10.5 else 12) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#283336"),
      plot.title = element_text(face = "bold", size = if (facet) 14 else 16),
      plot.subtitle = element_text(color = "#4c5a5e"),
      plot.caption = element_text(color = "#5f696b", hjust = 0),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "#e8e1d4", color = NA),
      legend.position = "right"
    )
  if (facet) p <- p + facet_wrap(~period, nrow = 1)
  p
}

period_file_names <- c(
  "Overall" = "overall",
  "Early (hotter)" = "early_hotter",
  "Late (mild-cold)" = "late_mild_cold"
)

for (period_i in names(period_file_names)) {
  period_data <- ratio_grid %>% filter(as.character(period) == period_i)
  period_subtitle <- if (period_i == "Overall") {
    "Full seven-round model; not an average of the two seasonal ratios"
  } else {
    paste0(period_i, " period model")
  }
  p_i <- base_heatmap(
    period_data,
    paste0("Fawn:female ratio | ", period_i),
    period_subtitle
  )
  suffix <- period_file_names[[period_i]]
  ggsave(paste0(output_prefix, "_", suffix, ".png"), p_i, width = 6.4, height = 5.6, dpi = 320)
  ggsave(paste0(output_prefix, "_", suffix, ".pdf"), p_i, width = 6.4, height = 5.6)
}

combined_plot <- base_heatmap(
  ratio_grid,
  "Predicted fawn:female abundance ratio across wolf x jackal context",
  "Overall uses the full model; seasonal panels use independently fitted Early and Late models",
  facet = TRUE
)
ggsave(paste0(output_prefix, "_combined.png"), combined_plot, width = 13.2, height = 5.0, dpi = 320)
ggsave(paste0(output_prefix, "_combined.pdf"), combined_plot, width = 13.2, height = 5.0)

fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
results_lines <- c(
  "# Predicted fawn:female ratios across wolf-jackal contexts",
  "",
  "## Definition",
  "",
  "Each cell is the ratio of the posterior mean predicted fawn abundance to the posterior mean predicted female abundance from the corresponding `spAbundance` wolf x jackal N-mixture models. Values in parentheses in the figures are the equivalent number of predicted fawns per 100 predicted females.",
  "",
  "The Overall panel was calculated from the full seven-round model. It is not an arithmetic average of the two seasonal ratios. Early and Late values come from the independently fitted seasonal models used in the original seasonal heatmap.",
  "",
  "These values describe model-based relative abundance, not observed demographic composition, recruitment, survival, or a population-level fawn:female ratio. The saved source grids contain posterior means but not joint posterior prediction draws, so uncertainty intervals cannot be reconstructed without refitting the models.",
  "",
  "## Prediction grid",
  ""
)

for (period_i in levels(ratio_grid$period)) {
  period_data <- ratio_grid %>% filter(as.character(period) == period_i)
  results_lines <- c(
    results_lines,
    paste0("### ", period_i),
    "",
    "| Wolf z | Jackal z | Predicted fawn | Predicted female | Fawn:female | Fawns per 100 females |",
    "|---:|---:|---:|---:|---:|---:|"
  )
  for (i in seq_len(nrow(period_data))) {
    results_lines <- c(results_lines, sprintf(
      "| %s | %s | %s | %s | %s | %s |",
      period_data$wolf_z[i], period_data$jackal_z[i],
      fmt(period_data$Fawn[i]), fmt(period_data$Female[i]),
      fmt(period_data$fawn_female_ratio[i]), fmt(period_data$fawns_per_100_females[i], 1)
    ))
  }
  results_lines <- c(results_lines, "")
}

overall_range <- ratio_grid %>%
  filter(period == "Overall") %>%
  summarise(minimum = min(fawn_female_ratio), maximum = max(fawn_female_ratio))
early_peak <- ratio_grid %>%
  filter(period == "Early (hotter)") %>%
  slice_max(fawn_female_ratio, n = 1, with_ties = FALSE)
late_same_cell <- ratio_grid %>%
  filter(period == "Late (mild-cold)", wolf_z == early_peak$wolf_z, jackal_z == early_peak$jackal_z)

results_lines <- c(
  results_lines,
  "## Descriptive patterns",
  "",
  paste0(
    "Across the full survey, predicted ratios ranged from `", fmt(overall_range$minimum),
    "` to `", fmt(overall_range$maximum), "` fawns per female."
  ),
  "",
  paste0(
    "The largest seasonal point estimate occurred in the Early period at wolf z = `",
    early_peak$wolf_z, "` and jackal z = `", early_peak$jackal_z, "`: `",
    fmt(early_peak$fawn_female_ratio), "` fawns per female (`",
    fmt(early_peak$fawns_per_100_females, 1), "` per 100 females). Both predicted abundances were low in that cell (fawn `",
    fmt(early_peak$Fawn), "`; female `", fmt(early_peak$Female),
    "`), so the high ratio reflects a stronger predicted female decline rather than high fawn abundance."
  ),
  "",
  paste0(
    "At that same wolf-jackal combination, the Late-period ratio was `",
    fmt(late_same_cell$fawn_female_ratio), "` (`",
    fmt(late_same_cell$fawns_per_100_females, 1),
    "` per 100 females). Because posterior ratio intervals are unavailable, this seasonal contrast is descriptive and should not be presented as statistically supported."
  ),
  "",
  "## Outputs",
  "",
  paste0("- `", output_prefix, "_combined.png`"),
  paste0("- `", output_prefix, "_overall.png`"),
  paste0("- `", output_prefix, "_early_hotter.png`"),
  paste0("- `", output_prefix, "_late_mild_cold.png`"),
  paste0("- `", output_prefix, "_prediction_grid.csv`")
)

writeLines(results_lines, paste0(output_prefix, "_results.md"), useBytes = TRUE)
