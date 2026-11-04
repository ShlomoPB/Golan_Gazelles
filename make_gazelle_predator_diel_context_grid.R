#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(activity)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(overlap)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
cfg <- list(
  source_prefix = "publication_predator_context_class_diel_2x2_median",
  output_prefix = "publication_gazelle_predator_diel_by_context_median",
  full_export = file.path("..", "full-export_2026-04-15_2006.csv"),
  gazelle_event_file = "Gazella gazella_independent_events_30min.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  local_timezone = "Asia/Jerusalem",
  event_gap_minutes = "30",
  minimum_curve_events = "5",
  curve_reps = "99",
  seed = "20260810"
)

for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(cfg)) stop("Unknown argument: ", key)
    cfg[[key]] <- parts[2]
  }
}

event_gap_minutes <- as.numeric(cfg$event_gap_minutes)
minimum_curve_events <- as.integer(cfg$minimum_curve_events)
curve_reps <- as.integer(cfg$curve_reps)
seed <- as.integer(cfg$seed)
set.seed(seed)

class_levels <- c("Fawn", "Female", "Male")
context_levels <- c(
  "Low both",
  "High jackal / Low wolf",
  "Low jackal / High wolf",
  "High both"
)
predator_levels <- c("Wolf", "Jackal")
curve_levels <- c("Gazelle", "Wolf", "Jackal")

density_file <- paste0(cfg$source_prefix, "_density_curves.csv")
gazelle_count_file <- paste0(cfg$source_prefix, "_event_counts.csv")
camera_file <- paste0(cfg$source_prefix, "_camera_categories.csv")
required_files <- c(
  density_file, gazelle_count_file, camera_file, cfg$full_export,
  cfg$gazelle_event_file, cfg$obs_cov_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

equation_of_time <- function(day_of_year) {
  b <- 2 * pi * (day_of_year - 81) / 364
  9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
}

parse_local_iso_time <- function(x) {
  as.POSIXct(
    sub("Z$", "", x),
    format = "%Y-%m-%dT%H:%M:%S",
    tz = cfg$local_timezone
  )
}

parse_export_time <- function(x) {
  parsed <- as.POSIXct(
    x,
    format = "%d/%m/%Y %H:%M:%S",
    tz = cfg$local_timezone
  )
  missing <- is.na(parsed)
  parsed[missing] <- as.POSIXct(
    x[missing],
    format = "%d/%m/%Y %H:%M",
    tz = cfg$local_timezone
  )
  parsed
}

apparent_solar_hour <- function(timestamp, longitude) {
  day_of_year <- as.integer(format(timestamp, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(timestamp, "%H", tz = "UTC")) +
    as.numeric(format(timestamp, "%M", tz = "UTC")) / 60 +
    as.numeric(format(timestamp, "%S", tz = "UTC")) / 3600
  (utc_hour + longitude / 15 + equation_of_time(day_of_year) / 60) %% 24
}

to_radians <- function(hour) 2 * pi * hour / 24
to_hours <- function(radians) 24 * radians / (2 * pi)

extract_activity_density <- function(fit, taxon, context) {
  pdf <- as.data.frame(fit@pdf)
  tibble(
    predator_context = context,
    curve_type = taxon,
    solar_hour = to_hours(pdf$x),
    density = pdf$y
  )
}

camera_context <- read_csv(camera_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(camera_name),
    jackal_abundance,
    wolf_abundance,
    predator_context = factor(
      as.character(predator_context),
      levels = context_levels
    )
  )

camera_counts <- camera_context %>%
  count(predator_context, name = "n_cameras")

context_labels <- camera_counts %>%
  mutate(
    panel_context = paste0(
      as.character(predator_context), "\n", n_cameras, " cameras"
    )
  )
context_label_lookup <- setNames(
  context_labels$panel_context,
  as.character(context_labels$predator_context)
)
panel_context_levels <- unname(context_label_lookup[context_levels])

obs <- read.csv(
  cfg$obs_cov_file,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    camera_name = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    days_operational = as.numeric(`Days.Operational`)
  ) %>%
  filter(
    is.finite(days_operational), days_operational > 0,
    !is.na(start_date), !is.na(end_date)
  )

message("Reading wolf and jackal records from: ", normalizePath(cfg$full_export))
raw_predator <- fread(
  cfg$full_export,
  select = c(
    "Camera Name", "Date/Time", "Species", "Trap Station Session ID",
    "Camelot GPS Latitude", "Camelot GPS Longitude"
  ),
  na.strings = c("", "NA"),
  showProgress = FALSE
)
setnames(
  raw_predator,
  c(
    "Camera Name", "Date/Time", "Species", "Trap Station Session ID",
    "Camelot GPS Latitude", "Camelot GPS Longitude"
  ),
  c("camera_name", "event_time_raw", "species", "session_id", "latitude", "longitude")
)

raw_predator <- raw_predator[
  species %chin% c("Canis lupus", "Canis aureus")
]
raw_predator[, `:=`(
  camera_name = as.integer(camera_name),
  session_id = as.character(session_id),
  timestamp = parse_export_time(event_time_raw),
  taxon = fifelse(species == "Canis lupus", "Wolf", "Jackal"),
  latitude = as.numeric(latitude),
  longitude = as.numeric(longitude)
)]
raw_predator <- raw_predator[
  !is.na(camera_name) & !is.na(session_id) & !is.na(timestamp)
]
raw_predator[, event_date := as.Date(timestamp, tz = cfg$local_timezone)]
raw_predator[, raw_row_id := .I]

obs_dt <- as.data.table(obs)
active_predator_rows <- merge(
  raw_predator,
  obs_dt,
  by = "camera_name",
  allow.cartesian = TRUE
)[
  event_date >= start_date & event_date <= end_date
]
setorder(active_predator_rows, raw_row_id, -start_date)
active_predator_rows <- active_predator_rows[, .SD[1], by = raw_row_id]

# Exact duplicate media annotations are collapsed before applying the
# 30-minute independence interval within camera, session, and taxon.
predator_deduplicated <- active_predator_rows[
  , .(
    latitude = median(latitude, na.rm = TRUE),
    longitude = median(longitude, na.rm = TRUE),
    raw_records = .N
  ),
  by = .(camera_name, session_id, taxon, timestamp)
]
setorder(predator_deduplicated, camera_name, session_id, taxon, timestamp)
predator_deduplicated[, gap_minutes := as.numeric(difftime(
  timestamp, shift(timestamp), units = "mins"
)), by = .(camera_name, session_id, taxon)]
predator_deduplicated[, independent_event_id := cumsum(
  is.na(gap_minutes) | gap_minutes > event_gap_minutes
), by = .(camera_name, session_id, taxon)]

predator_events_dt <- predator_deduplicated[
  , .(
    timestamp = min(timestamp),
    latitude = median(latitude, na.rm = TRUE),
    longitude = median(longitude, na.rm = TRUE),
    exact_timestamps = .N,
    raw_records = sum(raw_records)
  ),
  by = .(camera_name, session_id, taxon, independent_event_id)
]

predator_events <- as_tibble(predator_events_dt) %>%
  mutate(
    solar_hour = apparent_solar_hour(timestamp, longitude),
    time_rad = to_radians(solar_hour)
  ) %>%
  left_join(camera_context, by = "camera_name") %>%
  filter(
    taxon %in% predator_levels,
    predator_context %in% context_levels,
    is.finite(solar_hour)
  ) %>%
  mutate(
    taxon = factor(taxon, levels = predator_levels),
    predator_context = factor(predator_context, levels = context_levels)
  )

write_csv(
  predator_events %>%
    select(
      camera_name, session_id, taxon, timestamp, solar_hour,
      predator_context, exact_timestamps, raw_records
    ),
  paste0(cfg$output_prefix, "_predator_independent_events.csv")
)

predator_audit <- tibble(
  stage = c(
    "Raw full-export target records",
    "Target records in active formal rounds",
    paste0(event_gap_minutes, "-min independent predator events")
  ),
  wolf = c(
    raw_predator[taxon == "Wolf", .N],
    active_predator_rows[taxon == "Wolf", .N],
    sum(predator_events$taxon == "Wolf")
  ),
  jackal = c(
    raw_predator[taxon == "Jackal", .N],
    active_predator_rows[taxon == "Jackal", .N],
    sum(predator_events$taxon == "Jackal")
  )
)
write_csv(predator_audit, paste0(cfg$output_prefix, "_predator_event_audit.csv"))

predator_counts <- predator_events %>%
  group_by(predator_context, taxon) %>%
  summarise(
    n_events = n(),
    n_cameras_with_events = n_distinct(camera_name),
    .groups = "drop"
  ) %>%
  complete(
    predator_context = factor(context_levels, levels = context_levels),
    taxon = factor(predator_levels, levels = predator_levels),
    fill = list(n_events = 0L, n_cameras_with_events = 0L)
  )

predator_density_rows <- list()
predator_density_index <- 1L
for (context_i in context_levels) {
  for (taxon_i in predator_levels) {
    times_i <- predator_events$time_rad[
      predator_events$predator_context == context_i &
        predator_events$taxon == taxon_i
    ]
    if (length(times_i) < minimum_curve_events) next
    fit_i <- fitact(
      times_i,
      reps = curve_reps,
      sample = "data",
      show = FALSE
    )
    predator_density_rows[[predator_density_index]] <-
      extract_activity_density(fit_i, taxon_i, context_i)
    predator_density_index <- predator_density_index + 1L
  }
}
predator_density <- bind_rows(predator_density_rows) %>%
  mutate(
    predator_context = factor(predator_context, levels = context_levels),
    curve_type = factor(curve_type, levels = curve_levels)
  )

gazelle_density <- read_csv(density_file, show_col_types = FALSE) %>%
  transmute(
    panel_class = factor(class, levels = class_levels),
    predator_context = factor(
      as.character(predator_context),
      levels = context_levels
    ),
    curve_type = factor("Gazelle", levels = curve_levels),
    solar_hour,
    density
  )

predator_density_panels <- crossing(
  panel_class = factor(class_levels, levels = class_levels),
  predator_density
)

plot_density <- bind_rows(gazelle_density, predator_density_panels) %>%
  mutate(
    panel_context = factor(
      unname(context_label_lookup[as.character(predator_context)]),
      levels = panel_context_levels
    ),
    curve_type = factor(curve_type, levels = curve_levels)
  )

gazelle_counts <- read_csv(gazelle_count_file, show_col_types = FALSE) %>%
  transmute(
    predator_context = factor(
      as.character(predator_context),
      levels = context_levels
    ),
    panel_class = factor(class, levels = class_levels),
    gazelle_events = n_events,
    gazelle_cameras = n_cameras
  )

panel_counts <- gazelle_counts %>%
  left_join(
    predator_counts %>%
      select(predator_context, taxon, n_events) %>%
      pivot_wider(
        names_from = taxon,
        values_from = n_events,
        names_glue = "{tolower(taxon)}_events"
      ),
    by = "predator_context"
  ) %>%
  left_join(camera_counts, by = "predator_context") %>%
  mutate(
    panel_context = factor(
      unname(context_label_lookup[as.character(predator_context)]),
      levels = panel_context_levels
    ),
    annotation = paste0(
      ifelse(
        gazelle_events < minimum_curve_events,
        paste0("Gazelle: data deficient (n=", gazelle_events, ")"),
        paste0("Gazelle n=", gazelle_events)
      ),
      "\nWolf n=", wolf_events,
      ifelse(wolf_events < minimum_curve_events, " (ticks)", ""),
      "\nJackal n=", jackal_events,
      ifelse(jackal_events < minimum_curve_events, " (ticks)", "")
    )
  )

gazelle_events_raw <- read_csv(
  cfg$gazelle_event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  mutate(
    event_row_id = row_number(),
    camera_name = as.integer(camera_name),
    panel_class = recode(
      gazelle_class,
      juvenile = "Fawn",
      female = "Female",
      male = "Male"
    ),
    timestamp = parse_local_iso_time(event_time),
    event_date = as.Date(format(timestamp, tz = cfg$local_timezone)),
    solar_hour = apparent_solar_hour(timestamp, camelot_gps_longitude),
    time_rad = to_radians(solar_hour)
  )

active_gazelle_ids <- gazelle_events_raw %>%
  select(event_row_id, camera_name, event_date) %>%
  inner_join(
    obs %>% select(camera_name, start_date, end_date),
    by = "camera_name",
    relationship = "many-to-many"
  ) %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  group_by(event_row_id) %>%
  slice_max(start_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(event_row_id)

gazelle_events <- gazelle_events_raw %>%
  inner_join(active_gazelle_ids, by = "event_row_id") %>%
  left_join(camera_context, by = "camera_name") %>%
  filter(
    panel_class %in% class_levels,
    predator_context %in% context_levels,
    is.finite(solar_hour)
  ) %>%
  mutate(
    panel_class = factor(panel_class, levels = class_levels),
    predator_context = factor(predator_context, levels = context_levels),
    panel_context = factor(
      unname(context_label_lookup[as.character(predator_context)]),
      levels = panel_context_levels
    )
  )

sparse_gazelle_groups <- panel_counts %>%
  filter(gazelle_events < minimum_curve_events) %>%
  select(panel_class, predator_context)
sparse_gazelle_events <- gazelle_events %>%
  inner_join(
    sparse_gazelle_groups,
    by = c("panel_class", "predator_context")
  ) %>%
  mutate(curve_type = factor("Gazelle", levels = curve_levels))

sparse_predator_groups <- predator_counts %>%
  filter(n_events < minimum_curve_events) %>%
  select(predator_context, taxon)
sparse_predator_events <- predator_events %>%
  inner_join(
    sparse_predator_groups,
    by = c("predator_context", "taxon")
  ) %>%
  crossing(panel_class = factor(class_levels, levels = class_levels)) %>%
  mutate(
    panel_context = factor(
      unname(context_label_lookup[as.character(predator_context)]),
      levels = panel_context_levels
    ),
    curve_type = factor(as.character(taxon), levels = curve_levels)
  )

light_background <- tibble(
  xmin = c(0, 5, 7, 17, 19),
  xmax = c(5, 7, 17, 19, 24),
  light_phase = c("Night", "Dawn", "Day", "Dusk", "Night")
)
panel_background <- crossing(
  panel_class = factor(class_levels, levels = class_levels),
  panel_context = factor(panel_context_levels, levels = panel_context_levels),
  light_background
)

max_density <- max(plot_density$density, na.rm = TRUE)
panel_counts <- panel_counts %>%
  mutate(x = 23.65, y = max_density * 1.13)

curve_colors <- c(
  Gazelle = "#147D70",
  Wolf = "#315F8C",
  Jackal = "#B45E27"
)
curve_linetypes <- c(
  Gazelle = "solid",
  Wolf = "dotted",
  Jackal = "longdash"
)
curve_widths <- c(Gazelle = 1.05, Wolf = 0.9, Jackal = 0.9)

activity_grid <- ggplot(
  plot_density,
  aes(
    x = solar_hour,
    y = density,
    color = curve_type,
    linetype = curve_type,
    linewidth = curve_type
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
  geom_rug(
    data = sparse_predator_events,
    inherit.aes = FALSE,
    aes(x = solar_hour, color = curve_type),
    sides = "b",
    linewidth = 0.7,
    alpha = 0.9
  ) +
  geom_text(
    data = panel_counts,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = annotation),
    hjust = 1,
    vjust = 1,
    size = 2.45,
    lineheight = 0.94,
    color = "#263238"
  ) +
  facet_grid(panel_class ~ panel_context, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(0, 24, by = 6),
    limits = c(0, 24),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, max_density * 1.16),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_manual(values = curve_colors, drop = FALSE, name = NULL) +
  scale_linetype_manual(values = curve_linetypes, drop = FALSE, name = NULL) +
  scale_linewidth_manual(values = curve_widths, drop = FALSE, name = NULL) +
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
    title = "Gazelle and predator diel activity across wolf-jackal contexts",
    subtitle = paste0(
      "One gazelle class-context curve per panel; wolf and jackal curves use the same cameras and active survey rounds\n",
      "Gazelle = solid; wolf = dotted; jackal = dashed. Each taxon-context density is normalized independently"
    ),
    caption = paste0(
      "All curves use 30-min independent events and apparent solar time. Background shading approximates night, dawn, day, and dusk.\n",
      "Curve alignment indicates temporal coincidence, not predator-induced displacement; habitat and camera-level spatial selection remain potential confounders."
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
  paste0(cfg$output_prefix, ".png"),
  activity_grid,
  width = 15.5,
  height = 10.5,
  dpi = 320
)
ggsave(
  paste0(cfg$output_prefix, ".pdf"),
  activity_grid,
  width = 15.5,
  height = 10.5
)

write_csv(plot_density, paste0(cfg$output_prefix, "_density_curves.csv"))
write_csv(panel_counts, paste0(cfg$output_prefix, "_panel_event_counts.csv"))

overlap_rows <- list()
overlap_index <- 1L
for (class_i in class_levels) {
  for (context_i in context_levels) {
    gazelle_times <- gazelle_events$time_rad[
      gazelle_events$panel_class == class_i &
        gazelle_events$predator_context == context_i
    ]
    for (predator_i in predator_levels) {
      predator_times <- predator_events$time_rad[
        predator_events$taxon == predator_i &
          predator_events$predator_context == context_i
      ]
      estimator <- if (min(length(gazelle_times), length(predator_times)) < 75) {
        "Dhat1"
      } else {
        "Dhat4"
      }
      estimable <- min(length(gazelle_times), length(predator_times)) >=
        minimum_curve_events
      overlap_value <- if (estimable) {
        tryCatch(
          unname(overlapEst(gazelle_times, predator_times, type = estimator)),
          error = function(e) NA_real_
        )
      } else {
        NA_real_
      }
      overlap_rows[[overlap_index]] <- tibble(
        class = class_i,
        predator_context = context_i,
        predator = predator_i,
        gazelle_events = length(gazelle_times),
        predator_events = length(predator_times),
        estimator,
        temporal_overlap = overlap_value,
        status = ifelse(
          estimable,
          "Descriptive estimate",
          paste0("Fewer than ", minimum_curve_events, " events")
        )
      )
      overlap_index <- overlap_index + 1L
    }
  }
}
overlap_results <- bind_rows(overlap_rows) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(predator_context, levels = context_levels),
    predator = factor(predator, levels = predator_levels)
  ) %>%
  arrange(class, predator_context, predator)
write_csv(
  overlap_results,
  paste0(cfg$output_prefix, "_gazelle_predator_overlap.csv")
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

report <- c(
  "# Gazelle-predator diel activity grid",
  "",
  "## Analysis",
  "",
  paste0(
    "The figure splits the previous class-oriented diel figure into 12 panels: three gazelle classes by four median-split wolf-jackal camera contexts. Gazelle curves are the existing 30-min independent-event curves. Wolf (`Canis lupus`) and jackal (`Canis aureus`) records were extracted from `",
    cfg$full_export,
    "`, restricted to the same documented active camera rounds, collapsed to 30-min independent events within camera session and taxon, and converted from Israel local civil time to apparent solar time."
  ),
  "",
  "Each curve is normalized independently to integrate to one. Curve height therefore represents relative timing within a taxon-context combination, not abundance or event rate.",
  "",
  "## Predator event audit",
  "",
  "| Processing stage | Wolf | Jackal |",
  "|---|---:|---:|"
)
for (i in seq_len(nrow(predator_audit))) {
  report <- c(report, sprintf(
    "| %s | %s | %s |",
    predator_audit$stage[i], predator_audit$wolf[i], predator_audit$jackal[i]
  ))
}

report <- c(
  report,
  "",
  "## Context event counts",
  "",
  "| Gazelle class | Predator context | Cameras | Gazelle events | Wolf events | Jackal events |",
  "|---|---|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(panel_counts))) {
  report <- c(report, sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    panel_counts$panel_class[i], panel_counts$predator_context[i],
    panel_counts$n_cameras[i], panel_counts$gazelle_events[i],
    panel_counts$wolf_events[i], panel_counts$jackal_events[i]
  ))
}

report <- c(
  report,
  "",
  "## Descriptive gazelle-predator temporal overlap",
  "",
  "These point estimates quantify curve similarity but are not formal tests of predator-linked activity shifts. `Dhat1` was used when the smaller sample was below 75 events and `Dhat4` otherwise.",
  "",
  "| Class | Context | Predator | Events | Estimator | Overlap | Status |",
  "|---|---|---|---:|---|---:|---|"
)
for (i in seq_len(nrow(overlap_results))) {
  report <- c(report, sprintf(
    "| %s | %s | %s | `%s / %s` | `%s` | `%s` | %s |",
    overlap_results$class[i], overlap_results$predator_context[i],
    overlap_results$predator[i], overlap_results$gazelle_events[i],
    overlap_results$predator_events[i], overlap_results$estimator[i],
    fmt(overlap_results$temporal_overlap[i]), overlap_results$status[i]
  ))
}

fawn_overlap <- overlap_results %>%
  filter(class == "Fawn", is.finite(temporal_overlap))
female_max <- overlap_results %>%
  filter(class == "Female", is.finite(temporal_overlap)) %>%
  slice_max(temporal_overlap, n = 1, with_ties = FALSE)
male_max <- overlap_results %>%
  filter(class == "Male", is.finite(temporal_overlap)) %>%
  slice_max(temporal_overlap, n = 1, with_ties = FALSE)

report <- c(
  report,
  "",
  "## Descriptive patterns",
  "",
  paste0(
    "- Wolf and jackal activity was predominantly nocturnal in every context, whereas the gazelle curves were more concentrated in daylight. Fawn-predator overlap was consistently low where estimable (`",
    fmt(min(fawn_overlap$temporal_overlap)), "-",
    fmt(max(fawn_overlap$temporal_overlap)), "`)."
  ),
  paste0(
    "- Female activity had its greatest predator overlap with `",
    female_max$predator, "` in the `", female_max$predator_context,
    "` context (`Dhat = ", fmt(female_max$temporal_overlap), "`)."
  ),
  paste0(
    "- Male activity had its greatest predator overlap with `",
    male_max$predator, "` in the `", male_max$predator_context,
    "` context (`Dhat = ", fmt(male_max$temporal_overlap),
    "`). The relatively nocturnal male curve in that context visually aligns with both predator curves more than the male curves in the other contexts."
  ),
  "- The class-specific shifts do not show a single consistent rule in which all gazelle classes move toward or away from predator activity. The strongest visual evidence is therefore context- and class-specific, not a general predator-tracking or predator-avoidance pattern."
)

report <- c(
  report,
  "",
  "## Interpretation constraints",
  "",
  "- Similar or opposing peaks are descriptive. Shared responses to daylight, temperature, habitat, prey behavior, and camera placement can produce apparent temporal alignment or separation.",
  "- Predator curves are pooled across the cameras in each context and repeated across the three gazelle-class panels. They do not represent predator activity conditional on a particular gazelle event.",
  "- The high-jackal/low-wolf fawn panel contains only three fawn events and is labeled data deficient rather than displaying a density curve.",
  "- A formal test of linked temporal shifts would require a camera-aware model with predator activity as a time-varying predictor or a prespecified randomization analysis; visual curve displacement alone is not causal evidence.",
  "",
  "## Outputs",
  "",
  paste0("- Figure: `", cfg$output_prefix, ".png` and `.pdf`."),
  paste0("- Density curves: `", cfg$output_prefix, "_density_curves.csv`."),
  paste0("- Panel counts: `", cfg$output_prefix, "_panel_event_counts.csv`."),
  paste0("- Predator event audit: `", cfg$output_prefix, "_predator_event_audit.csv`."),
  paste0("- Independent predator events: `", cfg$output_prefix, "_predator_independent_events.csv`."),
  paste0("- Descriptive overlap: `", cfg$output_prefix, "_gazelle_predator_overlap.csv`.")
)
writeLines(report, paste0(cfg$output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved gazelle-predator diel grid with prefix:", cfg$output_prefix, "\n")
