#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)
})

coordinate_file <- "publication_camera_wolf_jackal_abundance_coordinate_plot_data.csv"
event_file <- "Gazella gazella_independent_events_30min_analysis_ready.csv"
output_prefix <- "publication_camera_wolf_jackal_abundance_class_detections"

required_files <- c(coordinate_file, event_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing authoritative input file(s): ", paste(missing_files, collapse = ", "))
}

class_levels <- c("Juvenile", "Adult female", "Adult male")
class_lookup <- c(
  juvenile = "Juvenile",
  female = "Adult female",
  male = "Adult male"
)
context_levels <- c(
  "Low both",
  "High jackal / Low wolf",
  "Low jackal / High wolf",
  "High both"
)
context_colors <- c(
  "Low both" = "#188977",
  "High jackal / Low wolf" = "#B12B32",
  "Low jackal / High wolf" = "#173B63",
  "High both" = "#8A5A9E"
)

coordinates <- read_csv(coordinate_file, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(camera_id),
    polygon,
    jackal_abundance_bup = as.numeric(jackal_abundance_bup),
    wolf_abundance_bup = as.numeric(wolf_abundance_bup),
    predator_context = factor(predator_context, levels = context_levels)
  ) %>%
  arrange(camera_id)

if (
  nrow(coordinates) != 60L ||
    n_distinct(coordinates$camera_id) != 60L ||
    !identical(coordinates$camera_id, 1:60) ||
    anyNA(coordinates)
) {
  stop("Coordinate input must contain complete records for cameras 1-60.")
}

events <- read_csv(event_file, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(camera_name),
    session_id = as.character(trap_station_session_id),
    gazelle_class = unname(class_lookup[as.character(gazelle_class)]),
    event_id
  )

if (anyNA(events$camera_id) || anyNA(events$session_id) || anyNA(events$gazelle_class)) {
  stop("Independent-event input contains unexpected camera or class values.")
}
if (anyDuplicated(events[c("camera_id", "session_id", "gazelle_class", "event_id")])) {
  stop("Independent event IDs are not unique within camera session and class.")
}

expected_totals <- c(
  "Juvenile" = 180L,
  "Adult female" = 1477L,
  "Adult male" = 2268L
)
observed_totals <- table(events$gazelle_class)
if (!all(observed_totals[names(expected_totals)] == expected_totals)) {
  stop("Class-event totals differ from the authoritative audit.")
}

camera_counts <- events %>%
  count(camera_id, gazelle_class, name = "independent_detections")

plot_data <- crossing(
  camera_id = 1:60,
  gazelle_class = class_levels
) %>%
  left_join(camera_counts, by = c("camera_id", "gazelle_class")) %>%
  mutate(independent_detections = replace_na(independent_detections, 0L)) %>%
  left_join(coordinates, by = "camera_id") %>%
  mutate(gazelle_class = factor(gazelle_class, levels = class_levels))

jackal_median <- median(coordinates$jackal_abundance_bup)
wolf_median <- median(coordinates$wolf_abundance_bup)
jackal_axis_max <- ceiling(max(coordinates$jackal_abundance_bup) / 25) * 25
wolf_axis_max <- ceiling(max(coordinates$wolf_abundance_bup) / 5) * 5

context_summary <- plot_data %>%
  group_by(gazelle_class, predator_context) %>%
  summarise(
    cameras = n(),
    independent_detections = sum(independent_detections),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      predator_context, "\n",
      independent_detections, " events; ", cameras, " cameras"
    ),
    x = rep(c(0.035, 0.965, 0.035, 0.965) * jackal_axis_max,
            times = length(class_levels)),
    y = rep(c(0.035, 0.035, 0.965, 0.965) * wolf_axis_max,
            times = length(class_levels)),
    hjust = rep(c(0, 1, 0, 1), times = length(class_levels)),
    vjust = rep(c(0, 0, 1, 1), times = length(class_levels))
  )

panel_labels <- plot_data %>%
  group_by(gazelle_class) %>%
  summarise(total = sum(independent_detections), .groups = "drop") %>%
  mutate(label = paste0(gazelle_class, " (n = ", format(total, big.mark = ","), ")")) %>%
  select(gazelle_class, label) %>%
  deframe()

bubble_plot <- ggplot(
  plot_data,
  aes(x = jackal_abundance_bup, y = wolf_abundance_bup)
) +
  annotate(
    "rect", xmin = 0, xmax = jackal_median, ymin = 0, ymax = wolf_median,
    fill = context_colors[["Low both"]], alpha = 0.045
  ) +
  annotate(
    "rect", xmin = jackal_median, xmax = jackal_axis_max,
    ymin = 0, ymax = wolf_median,
    fill = context_colors[["High jackal / Low wolf"]], alpha = 0.045
  ) +
  annotate(
    "rect", xmin = 0, xmax = jackal_median,
    ymin = wolf_median, ymax = wolf_axis_max,
    fill = context_colors[["Low jackal / High wolf"]], alpha = 0.045
  ) +
  annotate(
    "rect", xmin = jackal_median, xmax = jackal_axis_max,
    ymin = wolf_median, ymax = wolf_axis_max,
    fill = context_colors[["High both"]], alpha = 0.045
  ) +
  geom_vline(
    xintercept = jackal_median, color = "#5D676A",
    linewidth = 0.5, linetype = "22"
  ) +
  geom_hline(
    yintercept = wolf_median, color = "#5D676A",
    linewidth = 0.5, linetype = "22"
  ) +
  # A small baseline ring keeps zero-detection cameras visible.
  geom_point(
    aes(color = predator_context),
    size = 1.35, shape = 21, fill = "white", stroke = 0.55, alpha = 0.58,
    show.legend = FALSE
  ) +
  geom_point(
    data = filter(plot_data, independent_detections > 0),
    aes(
      size = independent_detections,
      color = predator_context,
      fill = predator_context
    ),
    shape = 21, stroke = 0.65, alpha = 0.68
  ) +
  geom_text(
    data = context_summary,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    color = "#455054", size = 2.25, lineheight = 0.94
  ) +
  facet_wrap(
    ~gazelle_class, nrow = 1,
    labeller = as_labeller(panel_labels)
  ) +
  scale_color_manual(values = context_colors, drop = FALSE, name = NULL) +
  scale_fill_manual(values = context_colors, drop = FALSE, guide = "none") +
  scale_size_area(
    max_size = 9.0,
    breaks = c(1, 10, 50, 100, 250),
    limits = c(0, max(plot_data$independent_detections)),
    name = "Independent detections"
  ) +
  scale_x_continuous(
    limits = c(0, jackal_axis_max),
    breaks = seq(0, jackal_axis_max, by = 50),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, wolf_axis_max),
    breaks = seq(0, wolf_axis_max, by = 5),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Jackal null-model abundance estimate (empirical-Bayes BUP)",
    y = "Wolf null-model abundance estimate (empirical-Bayes BUP)",
    title = "Gazelle demographic detections across wolf-jackal contexts",
    subtitle = paste0(
      "Each point is one camera; bubble area is proportional to independent ",
      "30-minute detections on a common scale"
    ),
    caption = paste0(
      "Small baseline rings denote cameras with zero detections for a class. ",
      "Dashed lines show median predator splits (jackal = ",
      sprintf("%.1f", jackal_median), "; wolf = ",
      sprintf("%.1f", wolf_median), "). Predator estimates represent relative ",
      "site-level activity, not density or direct predation risk."
    )
  ) +
  guides(
    color = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(size = 3.1, shape = 21, alpha = 0.78)
    ),
    size = guide_legend(order = 2, nrow = 1)
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E3E6E6", linewidth = 0.35),
    panel.spacing = grid::unit(0.7, "lines"),
    axis.text = element_text(color = "#344246", size = 8.5),
    axis.title = element_text(face = "bold", size = 10.5),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4D5A5E", size = 10.5),
    plot.caption = element_text(color = "#5F696B", hjust = 0, size = 8.5),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#F2F1EC", color = NA),
    legend.position = "top",
    legend.box = "vertical",
    legend.text = element_text(size = 8.6),
    legend.title = element_text(size = 9.2),
    legend.key.width = grid::unit(1.25, "lines"),
    legend.key.height = grid::unit(0.85, "lines"),
    plot.margin = margin(9, 10, 9, 10)
  )

ggsave(
  paste0(output_prefix, ".png"), bubble_plot,
  width = 15.5, height = 6.6, units = "in", dpi = 400, bg = "white"
)
ggsave(
  paste0(output_prefix, ".pdf"), bubble_plot,
  width = 15.5, height = 6.6, units = "in", device = cairo_pdf, bg = "white"
)

write_csv(
  plot_data %>%
    select(
      camera_id, polygon, gazelle_class, independent_detections,
      jackal_abundance_bup, wolf_abundance_bup, predator_context
    ),
  paste0(output_prefix, "_data.csv")
)
write_csv(
  context_summary %>%
    select(gazelle_class, predator_context, cameras, independent_detections),
  paste0(output_prefix, "_context_summary.csv")
)

cat("Class-specific camera bubble figure completed without model fitting.\n")
