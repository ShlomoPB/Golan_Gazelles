#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

input_file <- "demographic_nocturnality_summary.csv"
png_file <- "demographic_diel_donut_maintext.png"
pdf_file <- "demographic_diel_donut_maintext.pdf"
caption_file <- "demographic_diel_donut_caption.txt"

if (!file.exists(input_file)) {
  stop("Missing authoritative input: ", input_file)
}

class_levels <- c("Juvenile", "Adult female", "Adult male")
expected <- tibble(
  demographic_class = class_levels,
  total_independent_detections = c(180L, 1477L, 2268L),
  daylight_detections = c(149L, 903L, 1288L),
  nighttime_detections = c(31L, 574L, 980L),
  daylight_percentage = c(82.8, 61.1, 56.8),
  nighttime_percentage = c(17.2, 38.9, 43.2)
)

authoritative <- read_csv(input_file, show_col_types = FALSE) %>%
  transmute(
    demographic_class = as.character(demographic_class),
    total_independent_detections = as.integer(total_independent_detections),
    daylight_detections = as.integer(daylight_detections),
    nighttime_detections = as.integer(nighttime_detections),
    daylight_percentage = round(100 * diurnal_proportion, 1),
    nighttime_percentage = round(100 * nocturnal_proportion, 1)
  ) %>%
  arrange(match(demographic_class, class_levels))

if (!identical(authoritative, expected)) {
  stop("Final astronomical counts or percentages differ from the figure specification.")
}

plot_data <- authoritative %>%
  select(
    demographic_class, total_independent_detections,
    daylight_detections, nighttime_detections,
    daylight_percentage
  ) %>%
  tidyr::pivot_longer(
    cols = c(daylight_detections, nighttime_detections),
    names_to = "light_period",
    values_to = "detections"
  ) %>%
  mutate(
    demographic_class = factor(demographic_class, levels = class_levels),
    light_period = factor(
      light_period,
      levels = c("daylight_detections", "nighttime_detections"),
      labels = c("Daylight", "Nighttime")
    ),
    proportion = detections / total_independent_detections
  )

period_colors <- c(
  "Daylight" = "#B9DCE8",
  "Nighttime" = "#263E4A"
)

make_donut <- function(class_name) {
  class_data <- plot_data %>%
    filter(demographic_class == class_name)
  total <- unique(class_data$total_independent_detections)
  daylight <- unique(class_data$daylight_percentage)

  ggplot(
    class_data,
    aes(x = 2, y = proportion, fill = light_period)
  ) +
    geom_col(
      width = 0.70, color = "white", linewidth = 0.75,
      show.legend = TRUE
    ) +
    coord_polar(theta = "y", start = -pi / 2, clip = "off") +
    xlim(0.55, 2.55) +
    scale_fill_manual(values = period_colors, drop = FALSE) +
    labs(
      title = class_name,
      caption = sprintf(
        "%.1f%% daylight\nn = %s",
        daylight, format(total, big.mark = ",", scientific = FALSE)
      ),
      fill = NULL
    ) +
    guides(
      fill = guide_legend(
        direction = "horizontal", byrow = TRUE,
        override.aes = list(color = NA)
      )
    ) +
    theme_void(base_family = "sans", base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold", size = 12.3, hjust = 0.5,
        color = "#1F292E", margin = margin(b = -5)
      ),
      plot.caption = element_text(
        size = 11.2, hjust = 0.5, lineheight = 1.10,
        color = "#263238", margin = margin(t = 2)
      ),
      plot.caption.position = "plot",
      plot.margin = margin(3, 8, 3, 8),
      legend.position = "none"
    )
}

donut_panels <- wrap_plots(
  lapply(class_levels, make_donut),
  nrow = 1
)

legend_data <- tibble(
  light_period = factor(
    c("Daylight", "Nighttime"),
    levels = c("Daylight", "Nighttime")
  ),
  x = c(1, 2),
  y = 1
)
legend_plot <- ggplot(
  legend_data,
  aes(x = x, y = y, fill = light_period)
) +
  geom_tile() +
  scale_fill_manual(values = period_colors, drop = FALSE) +
  guides(
    fill = guide_legend(
      title = NULL, direction = "horizontal", nrow = 1,
      override.aes = list(color = NA)
    )
  ) +
  theme_void(base_family = "sans", base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(color = "#263238", size = 9.5),
    legend.key.size = grid::unit(0.78, "lines"),
    legend.spacing.x = grid::unit(0.26, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, 0, 0)
  )
legend_grob <- cowplot::get_legend(legend_plot)

inference_plot <- ggplot() +
  annotate(
    "text", x = 0.5, y = 0.5,
    label = "All pairwise class contrasts: Holm-adjusted P < 0.001",
    family = "sans", size = 3.15, color = "#4B5A61"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

final_figure <- donut_panels / inference_plot /
  wrap_elements(full = legend_grob) +
  plot_layout(heights = c(0.84, 0.07, 0.09)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(
  png_file, final_figure,
  width = 10.5, height = 4.0, units = "in", dpi = 600,
  bg = "white"
)
ggsave(
  pdf_file, final_figure,
  width = 10.5, height = 4.0, units = "in",
  device = cairo_pdf, bg = "white"
)

caption <- paste0(
  "Figure X. Proportion of independent mountain gazelle detections occurring ",
  "during daylight and nighttime by demographic class. Nighttime was defined ",
  "astronomically as the interval between local sunset and the following ",
  "sunrise (solar elevation < 0\u00b0). Values are based on independent detections ",
  "separated by at least 30 min. Juveniles were the most strongly diurnal class, ",
  "followed by adult females and adult males. Camera-cluster bootstrap estimates ",
  "of the nocturnal proportion were 17.2% (95% CI 11.8\u201321.6%) for juveniles, ",
  "38.9% (30.8\u201345.7%) for adult females, and 43.2% (38.4\u201348.4%) for adult males; ",
  "intervals were calculated from 10,000 camera-cluster bootstrap replicates."
)
writeLines(caption, caption_file, useBytes = TRUE)

message("Saved: ", png_file)
message("Saved: ", pdf_file)
message("Saved: ", caption_file)
