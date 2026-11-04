#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  site_cov_file = "SiteCovs_060324.csv",
  output_prefix = "publication_camera_wolf_jackal_activity_coordinate_plot"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

site_cov_file <- arg_map$site_cov_file
output_prefix <- arg_map$output_prefix

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

site_data <- read_csv(site_cov_file, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    jackal_activity = as.numeric(null.jackal.abundance),
    wolf_activity = as.numeric(Null.Wolf.Abundance),
    polygon = Polygon.Name
  )

if (nrow(site_data) != 60L || n_distinct(site_data$camera_id) != 60L) {
  stop("Expected one site-covariate row for each of 60 cameras.")
}
if (any(!complete.cases(site_data[c("jackal_activity", "wolf_activity")]))) {
  stop("Wolf and jackal activity indices must be complete for all cameras.")
}

jackal_median <- median(site_data$jackal_activity)
wolf_median <- median(site_data$wolf_activity)

plot_data <- site_data %>%
  mutate(
    predator_context = case_when(
      jackal_activity <= jackal_median & wolf_activity <= wolf_median ~ "Low both",
      jackal_activity > jackal_median & wolf_activity <= wolf_median ~ "High jackal / Low wolf",
      jackal_activity <= jackal_median & wolf_activity > wolf_median ~ "Low jackal / High wolf",
      jackal_activity > jackal_median & wolf_activity > wolf_median ~ "High both"
    ),
    predator_context = factor(predator_context, levels = context_levels)
  ) %>%
  arrange(camera_id)

context_counts <- plot_data %>%
  count(predator_context, .drop = FALSE) %>%
  mutate(
    label = paste0(predator_context, "\n", "n = ", n),
    x = c(0.05, 0.95, 0.05, 0.95),
    y = c(0.05, 0.05, 0.95, 0.95),
    hjust = c(0, 1, 0, 1),
    vjust = c(0, 0, 1, 1)
  )

activity_plot <- ggplot(
  plot_data,
  aes(x = jackal_activity, y = wolf_activity, color = predator_context)
) +
  annotate(
    "rect",
    xmin = 0,
    xmax = jackal_median,
    ymin = 0,
    ymax = wolf_median,
    fill = context_colors[["Low both"]],
    alpha = 0.045
  ) +
  annotate(
    "rect",
    xmin = jackal_median,
    xmax = 1,
    ymin = 0,
    ymax = wolf_median,
    fill = context_colors[["High jackal / Low wolf"]],
    alpha = 0.045
  ) +
  annotate(
    "rect",
    xmin = 0,
    xmax = jackal_median,
    ymin = wolf_median,
    ymax = 1,
    fill = context_colors[["Low jackal / High wolf"]],
    alpha = 0.045
  ) +
  annotate(
    "rect",
    xmin = jackal_median,
    xmax = 1,
    ymin = wolf_median,
    ymax = 1,
    fill = context_colors[["High both"]],
    alpha = 0.045
  ) +
  geom_vline(
    xintercept = jackal_median,
    color = "#5D676A",
    linewidth = 0.55,
    linetype = "22"
  ) +
  geom_hline(
    yintercept = wolf_median,
    color = "#5D676A",
    linewidth = 0.55,
    linetype = "22"
  ) +
  geom_point(
    size = 3.2,
    alpha = 0.92,
    shape = 21,
    fill = "white",
    stroke = 1.15
  ) +
  geom_text_repel(
    aes(label = camera_id),
    size = 3.0,
    color = "#263238",
    seed = 20260820,
    box.padding = 0.28,
    point.padding = 0.18,
    min.segment.length = 0,
    segment.color = "#9AA2A4",
    segment.size = 0.25,
    max.overlaps = Inf
  ) +
  geom_text(
    data = context_counts,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    color = "#455054",
    size = 3.2,
    lineheight = 0.95
  ) +
  annotate(
    "label",
    x = jackal_median,
    y = 0.015,
    label = paste0("Jackal median = ", sprintf("%.3f", jackal_median)),
    hjust = 0.5,
    vjust = 0,
    size = 2.8,
    color = "#4F595C",
    fill = "white",
    linewidth = 0
  ) +
  annotate(
    "label",
    x = 0.015,
    y = wolf_median,
    label = paste0("Wolf median = ", sprintf("%.3f", wolf_median)),
    angle = 90,
    hjust = 0.5,
    vjust = 0,
    size = 2.8,
    color = "#4F595C",
    fill = "white",
    linewidth = 0
  ) +
  scale_color_manual(values = context_colors, drop = FALSE, name = NULL) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_equal(clip = "off") +
  labs(
    x = "Null-model jackal activity index",
    y = "Null-model wolf activity index",
    title = "Camera positions in wolf-jackal activity space",
    subtitle = "Numbers identify cameras; dashed lines show the median splits used for diel contexts",
    caption = "Indices represent relative site-level activity, not absolute abundance or direct predation risk."
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E3E6E6", linewidth = 0.4),
    axis.text = element_text(color = "#344246"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#4D5A5E"),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "top",
    plot.margin = margin(10, 18, 12, 14)
  )

ggsave(
  paste0(output_prefix, ".png"),
  activity_plot,
  width = 9.4,
  height = 8.8,
  dpi = 320
)
ggsave(
  paste0(output_prefix, ".pdf"),
  activity_plot,
  width = 9.4,
  height = 8.8
)

write_csv(plot_data, paste0(output_prefix, "_data.csv"))

writeLines(
  c(
    "# Camera positions in wolf-jackal activity space",
    "",
    paste0("All 60 cameras had complete predator indices. The jackal median was ", sprintf("%.6f", jackal_median), " and the wolf median was ", sprintf("%.6f", wolf_median), "."),
    "",
    capture.output(print(context_counts %>% select(predator_context, n), n = Inf)),
    "",
    "The four colors and median splits match the retained diel predator-context analysis. The axes are relative site-level activity indices and should not be interpreted as absolute animal abundance or direct predation risk."
  ),
  paste0(output_prefix, "_results.md")
)

cat("Saved wolf-jackal camera-coordinate outputs with prefix:", output_prefix, "\n")
