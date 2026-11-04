suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

prediction_file <- "definitive_hierarchical_diel_predictions.csv"
model_file <- "definitive_hierarchical_diel_model_summary.csv"
stability_file <- "definitive_hierarchical_diel_effect_stability.csv"
png_file <- "definitive_hierarchical_diel_curves_maintext.png"
pdf_file <- "definitive_hierarchical_diel_curves_maintext.pdf"
caption_file <- "definitive_hierarchical_diel_curves_caption.txt"

required_files <- c(prediction_file, model_file, stability_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing authoritative input file(s): ", paste(missing_files, collapse = ", "))
}
input_hashes <- tools::md5sum(required_files)

response_levels <- c("Juvenile", "Adult female", "Adult male", "Pooled species")
panel_labels <- c(
  "Juvenile" = "A. Juvenile",
  "Adult female" = "B. Adult female",
  "Adult male" = "C. Adult male",
  "Pooled species" = "D. Pooled species"
)

old_context_levels <- c(
  "Wolf -1 / Jackal -1",
  "Wolf -1 / Jackal +1",
  "Wolf +1 / Jackal -1",
  "Wolf +1 / Jackal +1"
)
context_labels <- c(
  "Wolf -1 / Jackal -1" = "Low wolf / Low jackal",
  "Wolf -1 / Jackal +1" = "Low wolf / High jackal",
  "Wolf +1 / Jackal -1" = "High wolf / Low jackal",
  "Wolf +1 / Jackal +1" = "High wolf / High jackal"
)
display_context_levels <- unname(context_labels[old_context_levels])

# Low and high predator activity correspond to predictions at −1 and +1 SD,
# respectively, along the continuous wolf and jackal activity gradients; these
# values do not represent categorical camera groups.

# Color and line type jointly encode each prediction for accessibility.
context_colors <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#6A3D9A"),
  display_context_levels
)
context_linetypes <- setNames(
  c("solid", "longdash", "dotdash", "dotted"),
  display_context_levels
)

focal_terms <- c(
  wolf = "s(Time):wolf_z",
  jackal = "s(Time):jackal_z",
  interaction = "s(Time):wolf_jackal"
)
effect_levels <- c("wolf", "jackal", "interaction")
effect_labels <- c(
  wolf = "Wolf \u00d7 time",
  jackal = "Jackal \u00d7 time",
  interaction = "Wolf \u00d7 jackal \u00d7 time"
)
model_labels <- c(
  A = "baseline",
  B = "wolf-only",
  C = "jackal-only",
  D = "additive wolf + jackal",
  E = "wolf \u00d7 jackal"
)

prediction_data <- read_csv(prediction_file, show_col_types = FALSE) %>%
  mutate(
    response_class = factor(response_class, levels = response_levels),
    context = factor(
      unname(context_labels[as.character(context)]),
      levels = display_context_levels
    )
  )

if (anyNA(prediction_data$response_class) || anyNA(prediction_data$context)) {
  stop("Prediction file contains an unexpected response class or predator context.")
}

# The selected pooled model is wolf-only. Its two jackal settings must therefore
# be exact duplicate predictions before one copy is removed for plotting.
pooled_duplicate_check <- prediction_data %>%
  filter(response_class == "Pooled species") %>%
  group_by(Time, wolf_z) %>%
  summarise(
    n_jackal = n_distinct(jackal_z),
    activity_range = diff(range(normalized_activity)),
    lower_range = diff(range(normalized_lower)),
    upper_range = diff(range(normalized_upper)),
    .groups = "drop"
  )
if (
  any(pooled_duplicate_check$n_jackal != 2L) ||
    any(as.matrix(pooled_duplicate_check[c(
      "activity_range", "lower_range", "upper_range"
    )]) > 1e-12)
) {
  stop("Pooled predictions vary with jackal activity; refusing wolf-only display.")
}

predictions <- prediction_data %>%
  filter(response_class != "Pooled species" | jackal_z == -1)

prediction_structure <- predictions %>%
  group_by(response_class) %>%
  summarise(
    contexts = n_distinct(context),
    wolf_levels = n_distinct(wolf_z),
    rows = n(),
    .groups = "drop"
  )
expected_prediction_structure <- tibble(
  response_class = factor(response_levels, levels = response_levels),
  contexts = c(4L, 4L, 4L, 2L),
  wolf_levels = c(2L, 2L, 2L, 2L),
  rows = c(960L, 960L, 960L, 480L)
)
if (!identical(prediction_structure, expected_prediction_structure)) {
  stop("The plotted prediction structure is not four A-C contexts plus two D curves.")
}

model_summary <- read_csv(model_file, show_col_types = FALSE)
stability <- read_csv(stability_file, show_col_types = FALSE)

selected_models <- model_summary %>%
  filter(stage == "Final fREML", selected) %>%
  distinct(response_class, model_id)

expected_models <- tibble(
  response_class = response_levels,
  model_id = c("E", "E", "E", "B")
)
if (!identical(
    arrange(selected_models, match(response_class, response_levels)),
    expected_models
  )) {
  stop("Selected ecological models differ from the authoritative final structure.")
}

retains_effect <- function(model_id, effect) {
  case_when(
    effect == "wolf" ~ model_id %in% c("B", "D", "E"),
    effect == "jackal" ~ model_id %in% c("C", "D", "E"),
    effect == "interaction" ~ model_id == "E",
    TRUE ~ FALSE
  )
}

format_p <- function(p) {
  ifelse(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, digits = 3, format = "f")))
}

final_tests <- model_summary %>%
  filter(stage == "Final fREML", selected, term %in% focal_terms) %>%
  transmute(
    response_class,
    effect = names(focal_terms)[match(term, focal_terms)],
    EDF,
    P
  )

pooled_wolf_test <- final_tests %>%
  filter(response_class == "Pooled species", effect == "wolf")
if (
  nrow(pooled_wolf_test) != 1L ||
    abs(pooled_wolf_test$EDF - 3.8696427845902086) > 1e-8 ||
    pooled_wolf_test$P >= 0.001
) {
  stop("The selected pooled wolf-only smooth does not match the authoritative result.")
}

plot_labels <- crossing(
  response_class = response_levels,
  effect = effect_levels
) %>%
  left_join(selected_models, by = "response_class") %>%
  left_join(final_tests, by = c("response_class", "effect")) %>%
  left_join(
    stability %>% select(response_class, effect, robustness_status),
    by = c("response_class", "effect")
  ) %>%
  mutate(
    response_class = factor(response_class, levels = response_levels),
    effect = factor(effect, levels = effect_levels),
    retained = retains_effect(model_id, as.character(effect)),
    test_label = case_when(
      !retained ~ NA_character_,
      response_class == "Pooled species" & effect == "wolf" ~ paste0(
        effect_labels[as.character(effect)], ": EDF = ",
        formatC(EDF, digits = 2, format = "f"), ", ", format_p(P), ", robust"
      ),
      robustness_status == "robust" ~ paste0(
        effect_labels[as.character(effect)], ": ", format_p(P), ", robust"
      ),
      robustness_status == "sensitivity-dependent" ~ paste0(
        effect_labels[as.character(effect)], ": ", format_p(P),
        ", sensitivity-dependent"
      ),
      TRUE ~ paste0(effect_labels[as.character(effect)], ": ", format_p(P))
    )
  ) %>%
  arrange(response_class, effect) %>%
  group_by(response_class) %>%
  summarise(
    label = paste0(
      "Selected: ", model_labels[first(model_id)], "\n",
      paste(na.omit(test_label), collapse = "\n")
    ),
    .groups = "drop"
  )

pooled_key <- tibble(
  response_class = factor("Pooled species", levels = response_levels),
  context = factor(
    c("Low wolf / Low jackal", "High wolf / Low jackal"),
    levels = display_context_levels
  ),
  label = c("Low wolf", "High wolf"),
  y = max(
    predictions$normalized_upper[predictions$response_class == "Pooled species"]
  ) * c(0.945, 0.865),
  x = 18.15,
  xend = 19.55,
  label_x = 19.85
)

maintext_plot <- ggplot(
  predictions,
  aes(
    Time, normalized_activity,
    color = context, fill = context, linetype = context
  )
) +
  geom_ribbon(
    aes(ymin = normalized_lower, ymax = normalized_upper),
    alpha = 0.11, color = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.82) +
  geom_label(
    data = plot_labels,
    inherit.aes = FALSE,
    aes(x = 0.35, y = Inf, label = label),
    hjust = 0, vjust = 1.08, size = 2.30, lineheight = 0.96,
    fill = "white", alpha = 0.91, color = "#263238", linewidth = 0.18,
    label.padding = grid::unit(0.13, "lines")
  ) +
  geom_segment(
    data = pooled_key,
    inherit.aes = FALSE,
    aes(x = x, xend = xend, y = y, yend = y, color = context, linetype = context),
    linewidth = 0.82, show.legend = FALSE
  ) +
  geom_text(
    data = pooled_key,
    inherit.aes = FALSE,
    aes(x = label_x, y = y, label = label),
    hjust = 0, vjust = 0.5, size = 2.55, color = "#263238"
  ) +
  facet_wrap(
    ~response_class, ncol = 2, scales = "free_y",
    labeller = as_labeller(panel_labels)
  ) +
  scale_color_manual(values = context_colors, drop = FALSE) +
  scale_fill_manual(values = context_colors, drop = FALSE) +
  scale_linetype_manual(values = context_linetypes, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(0, 24, 3), limits = c(0, 24), expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
  labs(
    x = "Apparent solar time",
    y = "Normalized relative activity",
    color = NULL, fill = NULL, linetype = NULL
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.margin = margin(3, 8, 3, 6),
    strip.text = element_text(face = "bold", size = 10.5, hjust = 0,
                              margin = margin(4, 4, 4, 6)),
    strip.background = element_rect(fill = "#F2F1EC", color = NA),
    panel.spacing = grid::unit(0.72, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E7E9E9", linewidth = 0.28),
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 8.8, color = "#37474F"),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 8.6),
    legend.key.width = grid::unit(1.7, "lines"),
    legend.key.height = grid::unit(0.72, "lines"),
    legend.margin = margin(0, 0, 2, 0),
    legend.box.spacing = grid::unit(0.05, "lines")
  )

ggsave(
  png_file, maintext_plot, width = 9.4, height = 6.75,
  units = "in", dpi = 500, bg = "white"
)
ggsave(
  pdf_file, maintext_plot, width = 9.4, height = 6.75,
  units = "in", device = cairo_pdf, bg = "white"
)

if (!identical(tools::md5sum(required_files), input_hashes)) {
  stop("An authoritative prediction or model-summary input changed during rendering.")
}

caption <- paste0(
  "Figure X. Predator-associated variation in diel activity of mountain gazelle ",
  "demographic classes and the pooled population. Predicted normalized diel ",
  "activity of juveniles (A), adult females (B), adult males (C), and pooled ",
  "gazelles (D) across continuous gradients of site-level wolf and jackal activity. ",
  "Low and high predator activity correspond to predictions at \u22121 and +1 SD, ",
  "respectively, along the continuous wolf and jackal activity gradients; these ",
  "values do not represent categorical camera groups. Curves are population-level ",
  "predictions from hierarchical binomial cyclic generalized additive models and ",
  "integrate to one within each plotted predator condition; shaded bands represent ",
  "simulation-based 95% intervals. Statistical boxes summarize predator-by-time ",
  "terms from the selected model and their robustness across prespecified sensitivity ",
  "analyses. The selected pooled model retained wolf only, so panel D displays only ",
  "the low- and high-wolf predictions. Apparent solar time aligns solar ",
  "noon at 12:00 but does not standardize sunrise and sunset. Panels use independent ",
  "y-axis scales to preserve class-specific curve readability."
)
writeLines(caption, caption_file, useBytes = TRUE)

cat(
  "Main-text diel figure and caption completed without model fitting.\n",
  "Verified panel D source: selected pooled wolf-only model; ",
  "wolf x time EDF = 3.87, P < 0.001.\n",
  sep = ""
)
