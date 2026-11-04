#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

source_prefix <- "publication_gazelle_predator_diel_by_context_median"
output_prefix <- paste0(source_prefix, "_main_text")
gazelle_event_file <- "Gazella gazella_independent_events_30min_analysis_ready.csv"
obs_cov_file <- "ObsCovs_050524.csv"
camera_file <- "publication_predator_context_class_diel_2x2_median_camera_categories.csv"
local_timezone <- "Asia/Jerusalem"
minimum_curve_events <- 5L

density_file <- paste0(source_prefix, "_density_curves.csv")
count_file <- paste0(source_prefix, "_panel_event_counts.csv")
required_files <- c(
  density_file, count_file, gazelle_event_file, obs_cov_file, camera_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

row_levels <- c("Predators", "Fawn", "Female", "Male")
context_levels <- c(
  "Low both",
  "High jackal / Low wolf",
  "Low jackal / High wolf",
  "High both"
)
series_levels <- c("Fawn", "Female", "Male", "Wolf", "Jackal")

equation_of_time <- function(day_of_year) {
  b <- 2 * pi * (day_of_year - 81) / 364
  9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
}

parse_local_camera_time <- function(x) {
  as.POSIXct(
    sub("Z$", "", x),
    format = "%Y-%m-%dT%H:%M:%S",
    tz = local_timezone
  )
}

apparent_solar_hour <- function(timestamp, longitude) {
  day_of_year <- as.integer(format(timestamp, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(timestamp, "%H", tz = "UTC")) +
    as.numeric(format(timestamp, "%M", tz = "UTC")) / 60 +
    as.numeric(format(timestamp, "%S", tz = "UTC")) / 3600
  (utc_hour + longitude / 15 + equation_of_time(day_of_year) / 60) %% 24
}

density_raw <- read_csv(density_file, show_col_types = FALSE) %>%
  mutate(
    panel_class = as.character(panel_class),
    predator_context = factor(predator_context, levels = context_levels)
  )

panel_context_lookup <- density_raw %>%
  distinct(predator_context, panel_context) %>%
  arrange(predator_context)
panel_context_levels <- panel_context_lookup$panel_context
context_to_panel <- setNames(
  panel_context_lookup$panel_context,
  as.character(panel_context_lookup$predator_context)
)

gazelle_density <- density_raw %>%
  filter(curve_type == "Gazelle") %>%
  transmute(
    panel_row = factor(panel_class, levels = row_levels),
    predator_context,
    panel_context = factor(panel_context, levels = panel_context_levels),
    series = factor(panel_class, levels = series_levels),
    solar_hour,
    density
  )

# Predator curves were repeated down each class row in the supplementary
# overlay. Retain one copy per predator-context combination here.
predator_density <- density_raw %>%
  filter(curve_type %in% c("Wolf", "Jackal")) %>%
  distinct(predator_context, panel_context, curve_type, solar_hour, density) %>%
  transmute(
    panel_row = factor("Predators", levels = row_levels),
    predator_context,
    panel_context = factor(panel_context, levels = panel_context_levels),
    series = factor(curve_type, levels = series_levels),
    solar_hour,
    density
  )

plot_density <- bind_rows(predator_density, gazelle_density)
panel_counts <- read_csv(count_file, show_col_types = FALSE) %>%
  mutate(
    panel_class = as.character(panel_class),
    predator_context = factor(predator_context, levels = context_levels),
    panel_context = factor(panel_context, levels = panel_context_levels)
  )

predator_annotations <- panel_counts %>%
  distinct(
    predator_context, panel_context, n_cameras, wolf_events, jackal_events
  ) %>%
  mutate(
    panel_row = factor("Predators", levels = row_levels),
    label = paste0("Wolf n=", wolf_events, "\nJackal n=", jackal_events)
  )

gazelle_annotations <- panel_counts %>%
  transmute(
    panel_row = factor(panel_class, levels = row_levels),
    predator_context,
    panel_context,
    label = ifelse(
      gazelle_events < minimum_curve_events,
      paste0("Data deficient\n(n=", gazelle_events, ")"),
      paste0("n=", gazelle_events)
    )
  )

annotations <- bind_rows(predator_annotations, gazelle_annotations)

camera_context <- read_csv(camera_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(camera_name),
    predator_context = factor(predator_context, levels = context_levels)
  )

obs <- read.csv(
  obs_cov_file,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    camera_name = as.integer(`Camera Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    days_operational = as.numeric(`Days.Operational`)
  ) %>%
  filter(
    is.finite(days_operational), days_operational > 0,
    !is.na(start_date), !is.na(end_date)
  )

gazelle_events_raw <- read_csv(
  gazelle_event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  mutate(
    event_row_id = row_number(),
    camera_name = as.integer(camera_name),
    panel_row = recode(
      gazelle_class,
      juvenile = "Fawn",
      female = "Female",
      male = "Male"
    ),
    timestamp = parse_local_camera_time(event_time),
    event_date = as.Date(format(timestamp, tz = local_timezone)),
    solar_hour = apparent_solar_hour(timestamp, camelot_gps_longitude)
  )

active_event_ids <- gazelle_events_raw %>%
  select(event_row_id, camera_name, event_date) %>%
  inner_join(obs, by = "camera_name", relationship = "many-to-many") %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  group_by(event_row_id) %>%
  slice_max(start_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(event_row_id)

gazelle_events <- gazelle_events_raw %>%
  inner_join(active_event_ids, by = "event_row_id") %>%
  left_join(camera_context, by = "camera_name") %>%
  filter(
    panel_row %in% row_levels[-1],
    predator_context %in% context_levels,
    is.finite(solar_hour)
  ) %>%
  mutate(
    panel_row = factor(panel_row, levels = row_levels),
    panel_context = factor(
      unname(context_to_panel[as.character(predator_context)]),
      levels = panel_context_levels
    ),
    series = factor(as.character(panel_row), levels = series_levels)
  )

light_background <- tibble(
  xmin = c(0, 5, 7, 17, 19),
  xmax = c(5, 7, 17, 19, 24),
  light_phase = c("Night", "Dawn", "Day", "Dusk", "Night")
)
panel_background <- crossing(
  panel_row = factor(row_levels, levels = row_levels),
  panel_context = factor(panel_context_levels, levels = panel_context_levels),
  light_background
)

max_density <- max(plot_density$density, na.rm = TRUE)
annotations <- annotations %>%
  mutate(x = 23.55, y = max_density * 1.11)

series_colors <- c(
  Fawn = "#168A79",
  Female = "#B63439",
  Male = "#183B5B",
  Wolf = "#315F8C",
  Jackal = "#B45E27"
)
series_linetypes <- c(
  Fawn = "solid",
  Female = "solid",
  Male = "solid",
  Wolf = "dotted",
  Jackal = "longdash"
)
series_widths <- c(
  Fawn = 1.05,
  Female = 1.05,
  Male = 1.05,
  Wolf = 0.95,
  Jackal = 0.95
)

main_figure <- ggplot(
  plot_density,
  aes(
    x = solar_hour,
    y = density,
    color = series,
    linetype = series,
    linewidth = series
  )
) +
  geom_rect(
    data = panel_background,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
    alpha = 0.13,
    color = NA
  ) +
  geom_line() +
  geom_text(
    data = annotations,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    hjust = 1,
    vjust = 1,
    size = 2.45,
    lineheight = 0.94,
    color = "#263238"
  ) +
  facet_grid(panel_row ~ panel_context, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(0, 24, by = 6),
    limits = c(0, 24),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, max_density * 1.15),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_manual(values = series_colors, drop = FALSE, name = NULL) +
  scale_linetype_manual(values = series_linetypes, drop = FALSE, name = NULL) +
  scale_linewidth_manual(values = series_widths, drop = FALSE, name = NULL) +
  scale_fill_manual(
    values = c(
      Night = "#294C60",
      Dawn = "#E9A45D",
      Day = "#F5E6B8",
      Dusk = "#CF6A4C"
    ),
    guide = "none"
  ) +
  labs(
    x = "Apparent solar time",
    y = "Relative activity density",
    title = "Gazelle class activity relative to predator timing",
    subtitle = paste0(
      "Predator activity is shown once above each median-split camera context; gazelle classes are separated below\n",
      "All taxon-context curves are normalized independently and represent timing rather than event rate"
    ),
    caption = paste0(
      "Curves use 30-min independent events within active camera rounds. Wolf = dotted; jackal = dashed. ",
      "Background shading approximates light periods. Panels with fewer than five events are labeled data deficient and not plotted.\n",
      "Temporal alignment is descriptive and does not establish predator-induced displacement."
    )
  ) +
  theme_minimal(base_size = 9.2) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0.55, "lines"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#314043"),
    strip.text = element_text(face = "bold", size = 9.2),
    strip.background = element_rect(fill = "#E8E1D4", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#4C5A5E", lineheight = 1.05),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "top",
    legend.key.width = unit(2.2, "lines"),
    plot.margin = margin(10, 16, 10, 12)
  )

ggsave(
  paste0(output_prefix, ".png"),
  main_figure,
  width = 15.5,
  height = 12.2,
  dpi = 320
)
ggsave(
  paste0(output_prefix, ".pdf"),
  main_figure,
  width = 15.5,
  height = 12.2
)

cat("Saved main-text gazelle-predator diel figure:", output_prefix, "\n")
