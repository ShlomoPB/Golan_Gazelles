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
  jackal_file = paste0(
    "C:/Users/preissbloom/OneDrive - Tel-Aviv University/Golan Wolf Study/",
    "Wolf Reports and Excels/Null Jackal abundance with SE.csv"
  ),
  wolf_file = paste0(
    "C:/Users/preissbloom/OneDrive - Tel-Aviv University/Golan Wolf Study/",
    "Wolf Reports and Excels/Null Wolf abundance with SE.csv"
  ),
  output_prefix = "publication_camera_wolf_jackal_abundance_coordinate_plot"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

read_predator_source <- function(path, standardized_column, species) {
  if (!file.exists(path)) stop("Missing ", species, " source file: ", path)

  raw <- read_csv(
    path,
    show_col_types = FALSE,
    name_repair = "unique_quiet"
  )
  required <- c("2.50%", "97.50%", "SE", "abundance", standardized_column)
  missing_columns <- setdiff(required, names(raw))
  if (length(missing_columns) > 0) {
    stop(
      species,
      " source is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  cleaned <- tibble(
    camera_id = suppressWarnings(as.integer(raw[[1]])),
    lower = as.numeric(raw[["2.50%"]]),
    upper = as.numeric(raw[["97.50%"]]),
    reported_se = as.numeric(raw[["SE"]]),
    estimate = as.numeric(raw[["abundance"]]),
    standardized = as.numeric(raw[[standardized_column]])
  ) %>%
    filter(camera_id %in% 1:60) %>%
    arrange(camera_id)

  if (
    nrow(cleaned) != 60L ||
      n_distinct(cleaned$camera_id) != 60L ||
      !identical(cleaned$camera_id, 1:60)
  ) {
    stop("Expected exactly one valid ", species, " row for each camera ID 1-60.")
  }
  if (any(!complete.cases(cleaned[c("estimate", "standardized")]))) {
    stop(species, " estimates must be complete for all 60 cameras.")
  }

  source_max <- max(cleaned$estimate)
  if (!isTRUE(all.equal(
    cleaned$standardized,
    cleaned$estimate / source_max,
    tolerance = 1e-8
  ))) {
    stop(species, " standardized values do not equal abundance / maximum abundance.")
  }

  attr(cleaned, "source_max") <- source_max
  cleaned
}

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

site_data <- read_csv(arg_map$site_cov_file, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    jackal_index_sitecovs = as.numeric(null.jackal.abundance),
    wolf_index_sitecovs = as.numeric(Null.Wolf.Abundance),
    polygon = Polygon.Name
  )

if (nrow(site_data) != 60L || n_distinct(site_data$camera_id) != 60L) {
  stop("Expected one site-covariate row for each of 60 cameras.")
}

jackal_source <- read_predator_source(
  arg_map$jackal_file,
  "jackal_abundance_stand",
  "jackal"
)
wolf_source <- read_predator_source(
  arg_map$wolf_file,
  "wolf_abundance_stand",
  "wolf"
)
jackal_max <- attr(jackal_source, "source_max")
wolf_max <- attr(wolf_source, "source_max")

plot_data <- site_data %>%
  left_join(
    jackal_source %>%
      transmute(
        camera_id,
        jackal_lower = lower,
        jackal_upper = upper,
        jackal_abundance_bup = estimate,
        jackal_index = standardized
      ),
    by = "camera_id"
  ) %>%
  left_join(
    wolf_source %>%
      transmute(
        camera_id,
        wolf_lower = lower,
        wolf_upper = upper,
        wolf_abundance_bup = estimate,
        wolf_index = standardized
      ),
    by = "camera_id"
  ) %>%
  arrange(camera_id)

if (
  max(abs(plot_data$jackal_index - plot_data$jackal_index_sitecovs)) > 1e-8 ||
    max(abs(plot_data$wolf_index - plot_data$wolf_index_sitecovs)) > 1e-8
) {
  stop("Source-file standardized estimates do not match SiteCovs_060324.csv.")
}

jackal_median <- median(plot_data$jackal_abundance_bup)
wolf_median <- median(plot_data$wolf_abundance_bup)
jackal_axis_max <- ceiling(max(plot_data$jackal_abundance_bup) / 25) * 25
wolf_axis_max <- ceiling(max(plot_data$wolf_abundance_bup) / 5) * 5

plot_data <- plot_data %>%
  mutate(
    predator_context = case_when(
      jackal_abundance_bup <= jackal_median &
        wolf_abundance_bup <= wolf_median ~ "Low both",
      jackal_abundance_bup > jackal_median &
        wolf_abundance_bup <= wolf_median ~ "High jackal / Low wolf",
      jackal_abundance_bup <= jackal_median &
        wolf_abundance_bup > wolf_median ~ "Low jackal / High wolf",
      jackal_abundance_bup > jackal_median &
        wolf_abundance_bup > wolf_median ~ "High both"
    ),
    predator_context = factor(predator_context, levels = context_levels)
  )

context_counts <- plot_data %>%
  count(predator_context, .drop = FALSE) %>%
  mutate(
    label = paste0(predator_context, "\n", "n = ", n),
    x = c(0.04, 0.96, 0.04, 0.96) * jackal_axis_max,
    y = c(0.04, 0.04, 0.96, 0.96) * wolf_axis_max,
    hjust = c(0, 1, 0, 1),
    vjust = c(0, 0, 1, 1)
  )

abundance_plot <- ggplot(
  plot_data,
  aes(
    x = jackal_abundance_bup,
    y = wolf_abundance_bup,
    color = predator_context
  )
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
    xmax = jackal_axis_max,
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
    ymax = wolf_axis_max,
    fill = context_colors[["Low jackal / High wolf"]],
    alpha = 0.045
  ) +
  annotate(
    "rect",
    xmin = jackal_median,
    xmax = jackal_axis_max,
    ymin = wolf_median,
    ymax = wolf_axis_max,
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
    y = 0.012 * wolf_axis_max,
    label = paste0("Jackal median = ", sprintf("%.1f", jackal_median)),
    hjust = 0.5,
    vjust = 0,
    size = 2.8,
    color = "#4F595C",
    fill = "white",
    linewidth = 0
  ) +
  annotate(
    "label",
    x = 0.012 * jackal_axis_max,
    y = wolf_median,
    label = paste0("Wolf median = ", sprintf("%.1f", wolf_median)),
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
    title = "Camera positions in wolf-jackal abundance-estimate space",
    subtitle = paste0(
      "Numbers identify cameras; dashed lines show the median splits used ",
      "for diel contexts"
    ),
    caption = paste0(
      "Unscaled model estimates represent relative site-level predator ",
      "activity, not density, observed counts, or direct predation risk."
    )
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    aspect.ratio = 1,
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
  paste0(arg_map$output_prefix, ".png"),
  abundance_plot,
  width = 9.4,
  height = 8.8,
  dpi = 320
)
ggsave(
  paste0(arg_map$output_prefix, ".pdf"),
  abundance_plot,
  width = 9.4,
  height = 8.8
)

write_csv(
  plot_data %>%
    select(
      camera_id,
      polygon,
      jackal_abundance_bup,
      jackal_lower,
      jackal_upper,
      wolf_abundance_bup,
      wolf_lower,
      wolf_upper,
      jackal_index,
      wolf_index,
      predator_context
    ),
  paste0(arg_map$output_prefix, "_data.csv")
)

writeLines(
  c(
    "# Camera positions in wolf-jackal abundance-estimate space",
    "",
    paste0(
      "All 60 cameras had complete estimates. The unscaled jackal median was ",
      sprintf("%.6f", jackal_median),
      " (range ",
      sprintf("%.6f", min(plot_data$jackal_abundance_bup)),
      "-",
      sprintf("%.6f", max(plot_data$jackal_abundance_bup)),
      "), and the unscaled wolf median was ",
      sprintf("%.6f", wolf_median),
      " (range ",
      sprintf("%.6f", min(plot_data$wolf_abundance_bup)),
      "-",
      sprintf("%.6f", max(plot_data$wolf_abundance_bup)),
      ")."
    ),
    "",
    paste0(
      "The source estimates reproduce the SiteCovs indices exactly after ",
      "division by the species maxima (jackal = ",
      sprintf("%.6f", jackal_max),
      "; wolf = ",
      sprintf("%.6f", wolf_max),
      ")."
    ),
    "",
    capture.output(print(context_counts %>% select(predator_context, n), n = Inf)),
    "",
    paste0(
      "The four colors and median splits match the retained diel ",
      "predator-context analysis. The axes show unscaled empirical-Bayes ",
      "BUP estimates and should not be interpreted as predator density, ",
      "observed counts, or direct predation risk."
    )
  ),
  paste0(arg_map$output_prefix, "_results.md")
)

cat(
  "Saved raw-scale wolf-jackal camera-coordinate outputs with prefix:",
  arg_map$output_prefix,
  "\n"
)
