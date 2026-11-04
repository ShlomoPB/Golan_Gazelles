#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(activity)
  library(dplyr)
  library(ggplot2)
  library(overlap)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  event_file = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  site_cov_file = "SiteCovs_060324.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  output_prefix = "publication_predator_context_class_diel_2x2",
  local_timezone = "Asia/Jerusalem",
  lower_quantile = "0.33",
  upper_quantile = "0.67",
  minimum_curve_events = "5",
  minimum_overlap_events = "20",
  curve_reps = "99",
  overlap_bootstrap_reps = "500",
  seed = "20260809"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

event_file <- arg_map$event_file
site_cov_file <- arg_map$site_cov_file
obs_cov_file <- arg_map$obs_cov_file
output_prefix <- arg_map$output_prefix
local_timezone <- arg_map$local_timezone
lower_quantile <- as.numeric(arg_map$lower_quantile)
upper_quantile <- as.numeric(arg_map$upper_quantile)
minimum_curve_events <- as.integer(arg_map$minimum_curve_events)
minimum_overlap_events <- as.integer(arg_map$minimum_overlap_events)
curve_reps <- as.integer(arg_map$curve_reps)
overlap_bootstrap_reps <- as.integer(arg_map$overlap_bootstrap_reps)
seed <- as.integer(arg_map$seed)

if (
  !is.finite(lower_quantile) || !is.finite(upper_quantile) ||
    lower_quantile < 0 || upper_quantile > 1 || lower_quantile > upper_quantile
) {
  stop("Quantiles must satisfy 0 <= lower_quantile <= upper_quantile <= 1.")
}

classification_label <- if (
  isTRUE(all.equal(lower_quantile, 0.5)) &&
    isTRUE(all.equal(upper_quantile, 0.5))
) {
  "median split"
} else {
  paste0(
    round(100 * lower_quantile), "th/",
    round(100 * upper_quantile), "th-percentile extremes"
  )
}
assignment_description <- if (classification_label == "median split") {
  "All cameras were assigned to one of the four predator contexts."
} else {
  "Cameras between the low and high cutoffs were classified as intermediate and excluded from this four-corner figure."
}

set.seed(seed)

class_levels <- c("Fawn", "Female", "Male")
context_levels <- c(
  "Low both",
  "High jackal / Low wolf",
  "Low jackal / High wolf",
  "High both"
)
class_labels <- c(
  juvenile = "Fawn",
  female = "Female",
  male = "Male"
)
class_colors <- c(
  Fawn = "#078C83",
  Female = "#C52E36",
  Male = "#173B63"
)

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

to_radians <- function(hour) 2 * pi * hour / 24
to_hours <- function(radians) 24 * radians / (2 * pi)

extract_activity_density <- function(fit, context, class_label) {
  density <- as.data.frame(fit@pdf)
  tibble(
    context,
    class = class_label,
    solar_hour = to_hours(density$x),
    density = density$y
  )
}

pick_overlap_estimator <- function(n_a, n_b) {
  if (min(n_a, n_b) < 75) "Dhat1" else "Dhat4"
}

clustered_overlap <- function(data, context, class_a, class_b, camera_ids, reps) {
  context_data <- filter(data, predator_context == context)
  time_a <- context_data$time_rad[context_data$class == class_a]
  time_b <- context_data$time_rad[context_data$class == class_b]
  n_a <- length(time_a)
  n_b <- length(time_b)
  cameras_a <- n_distinct(context_data$camera_name[context_data$class == class_a])
  cameras_b <- n_distinct(context_data$camera_name[context_data$class == class_b])
  estimator <- pick_overlap_estimator(n_a, n_b)

  result <- tibble(
    predator_context = context,
    class_a,
    class_b,
    n_events_a = n_a,
    n_events_b = n_b,
    n_cameras_a = cameras_a,
    n_cameras_b = cameras_b,
    estimator,
    overlap = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    bootstrap_successes = 0L,
    estimable = FALSE,
    reason = NA_character_
  )

  if (min(n_a, n_b) < minimum_overlap_events) {
    result$reason <- paste0("Fewer than ", minimum_overlap_events, " events in at least one class")
    return(result)
  }

  result$overlap <- unname(overlapEst(time_a, time_b, type = estimator))
  boot_values <- replicate(reps, {
    sampled_cameras <- sample(camera_ids, length(camera_ids), replace = TRUE)
    boot_a <- unlist(lapply(sampled_cameras, function(camera_i) {
      context_data$time_rad[
        context_data$camera_name == camera_i & context_data$class == class_a
      ]
    }), use.names = FALSE)
    boot_b <- unlist(lapply(sampled_cameras, function(camera_i) {
      context_data$time_rad[
        context_data$camera_name == camera_i & context_data$class == class_b
      ]
    }), use.names = FALSE)
    if (length(boot_a) < 5 || length(boot_b) < 5) {
      return(NA_real_)
    }
    tryCatch(
      unname(overlapEst(boot_a, boot_b, type = estimator)),
      error = function(e) NA_real_
    )
  })
  boot_values <- boot_values[is.finite(boot_values)]
  result$bootstrap_successes <- length(boot_values)
  if (length(boot_values) < ceiling(0.8 * reps)) {
    result$reason <- "Too few successful camera-cluster bootstrap replicates"
    return(result)
  }

  result$conf_low <- unname(quantile(boot_values, 0.025))
  result$conf_high <- unname(quantile(boot_values, 0.975))
  result$estimable <- TRUE
  result$reason <- NA_character_
  result
}

obs <- read.csv(
  obs_cov_file,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    camera_name = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y")
  ) %>%
  filter(!is.na(camera_name), !is.na(start_date), !is.na(end_date))

site_covariates <- read_csv(site_cov_file, show_col_types = FALSE)
jackal_low <- quantile(
  site_covariates$null.jackal.abundance,
  lower_quantile,
  na.rm = TRUE
)
jackal_high <- quantile(
  site_covariates$null.jackal.abundance,
  upper_quantile,
  na.rm = TRUE
)
wolf_low <- quantile(
  site_covariates$Null.Wolf.Abundance,
  lower_quantile,
  na.rm = TRUE
)
wolf_high <- quantile(
  site_covariates$Null.Wolf.Abundance,
  upper_quantile,
  na.rm = TRUE
)

camera_context <- site_covariates %>%
  transmute(
    camera_name = as.integer(Camera.ID),
    jackal_abundance = null.jackal.abundance,
    wolf_abundance = Null.Wolf.Abundance,
    predator_context = case_when(
      jackal_abundance <= jackal_low & wolf_abundance <= wolf_low ~ "Low both",
      jackal_abundance >= jackal_high & wolf_abundance <= wolf_low ~ "High jackal / Low wolf",
      jackal_abundance <= jackal_low & wolf_abundance >= wolf_high ~ "Low jackal / High wolf",
      jackal_abundance >= jackal_high & wolf_abundance >= wolf_high ~ "High both",
      TRUE ~ "Intermediate / mixed"
    )
  )

raw_events <- read_csv(
  event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  mutate(
    event_row_id = row_number(),
    camera_name = as.integer(camera_name),
    class = unname(class_labels[as.character(gazelle_class)]),
    timestamp = parse_local_camera_time(event_time),
    event_date = as.Date(format(timestamp, tz = local_timezone)),
    solar_hour = apparent_solar_hour(timestamp, camelot_gps_longitude),
    time_rad = to_radians(solar_hour)
  )

# Keep only detections within a documented active camera round. Boundary dates
# can belong to adjacent rounds, so assign the event to the later deployment.
active_matches <- raw_events %>%
  select(event_row_id, camera_name, event_date) %>%
  inner_join(obs, by = "camera_name", relationship = "many-to-many") %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  group_by(event_row_id) %>%
  slice_max(start_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(event_row_id, round)

events <- raw_events %>%
  inner_join(active_matches, by = "event_row_id") %>%
  left_join(camera_context, by = "camera_name") %>%
  filter(predator_context %in% context_levels, class %in% class_levels) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(predator_context, levels = context_levels)
  )

camera_counts <- camera_context %>%
  filter(predator_context %in% context_levels) %>%
  count(predator_context, name = "n_cameras") %>%
  mutate(predator_context = factor(predator_context, levels = context_levels)) %>%
  arrange(predator_context)

event_counts <- events %>%
  count(predator_context, class, name = "n_events") %>%
  complete(
    predator_context = factor(context_levels, levels = context_levels),
    class = factor(class_levels, levels = class_levels),
    fill = list(n_events = 0L)
  ) %>%
  left_join(camera_counts, by = "predator_context") %>%
  arrange(predator_context, class)

density_rows <- list()
density_index <- 1L
for (context_i in context_levels) {
  for (class_i in class_levels) {
    class_times <- events$time_rad[
      events$predator_context == context_i & events$class == class_i
    ]
    if (length(class_times) < minimum_curve_events) {
      next
    }
    activity_fit <- fitact(
      class_times,
      reps = curve_reps,
      sample = "data",
      show = FALSE
    )
    density_rows[[density_index]] <- extract_activity_density(
      activity_fit,
      context_i,
      class_i
    )
    density_index <- density_index + 1L
  }
}
density_curves <- bind_rows(density_rows) %>%
  mutate(
    predator_context = factor(context, levels = context_levels),
    class = factor(class, levels = class_levels)
  )

sparse_fawn_events <- events %>%
  filter(class == "Fawn") %>%
  add_count(predator_context, class, name = "class_context_n") %>%
  filter(class_context_n < minimum_curve_events)

facet_lookup <- camera_counts %>%
  transmute(
    predator_context,
    facet_label = paste0(as.character(predator_context), " (", n_cameras, " cameras)")
  )
facet_levels <- facet_lookup$facet_label

density_curves <- density_curves %>%
  left_join(facet_lookup, by = "predator_context") %>%
  mutate(facet_label = factor(facet_label, levels = facet_levels))
sparse_fawn_events <- sparse_fawn_events %>%
  left_join(facet_lookup, by = "predator_context") %>%
  mutate(facet_label = factor(facet_label, levels = facet_levels))

panel_annotations <- event_counts %>%
  group_by(predator_context, n_cameras) %>%
  summarise(
    label = paste0(
      "Fawn n=", n_events[class == "Fawn"],
      ifelse(n_events[class == "Fawn"] < minimum_curve_events, " (ticks only)", ""),
      "\nFemale n=", n_events[class == "Female"],
      "\nMale n=", n_events[class == "Male"]
    ),
    .groups = "drop"
  ) %>%
  left_join(facet_lookup, by = "predator_context") %>%
  mutate(
    facet_label = factor(facet_label, levels = facet_levels),
    x = 23.6,
    y = Inf
  )

light_background <- tibble(
  xmin = c(0, 5, 7, 17, 19),
  xmax = c(5, 7, 17, 19, 24),
  light_phase = c("Night", "Dawn", "Day", "Dusk", "Night")
)
panel_background <- crossing(
  facet_label = factor(facet_levels, levels = facet_levels),
  light_background
)

activity_plot <- ggplot(
  density_curves,
  aes(x = solar_hour, y = density, color = class)
) +
  geom_rect(
    data = panel_background,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_phase),
    alpha = 0.16,
    color = NA
  ) +
  geom_line(linewidth = 1.05) +
  geom_rug(
    data = sparse_fawn_events,
    inherit.aes = FALSE,
    aes(x = solar_hour, color = class),
    sides = "b",
    linewidth = 0.9,
    alpha = 0.9
  ) +
  geom_text(
    data = panel_annotations,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    hjust = 1,
    vjust = 1.15,
    size = 3.0,
    lineheight = 0.98,
    color = "#263238"
  ) +
  facet_wrap(~facet_label, ncol = 2) +
  scale_x_continuous(
    breaks = seq(0, 24, by = 3),
    limits = c(0, 24),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_manual(values = class_colors, drop = FALSE, name = NULL) +
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
    title = "Gazelle diel activity across wolf-jackal contexts",
    subtitle = paste0(
      "Independent 30-min events pooled across the survey; low/high cameras defined by a ", classification_label, "\n",
      "Curves are normalized within class and context; class-context samples below five events are shown as ticks only"
    ),
    caption = paste0(
      "Camera clocks were interpreted as Israel local civil time and converted to apparent solar time. ",
      "Background shading indicates approximate night, dawn, day, and dusk periods."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = grid::unit(1.0, "lines"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#314043"),
    strip.text = element_text(face = "bold", size = 10.5),
    strip.background = element_rect(fill = "#E8E1D4", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#4C5A5E", lineheight = 1.05),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "top"
  )

ggsave(
  paste0(output_prefix, ".png"),
  activity_plot,
  width = 11.5,
  height = 8.2,
  dpi = 320
)
ggsave(
  paste0(output_prefix, ".pdf"),
  activity_plot,
  width = 11.5,
  height = 8.2
)

pair_matrix <- t(combn(class_levels, 2))
overlap_rows <- list()
overlap_index <- 1L
for (context_i in context_levels) {
  camera_ids <- camera_context$camera_name[camera_context$predator_context == context_i]
  for (pair_i in seq_len(nrow(pair_matrix))) {
    overlap_rows[[overlap_index]] <- clustered_overlap(
      events,
      context_i,
      pair_matrix[pair_i, 1],
      pair_matrix[pair_i, 2],
      camera_ids,
      overlap_bootstrap_reps
    )
    overlap_index <- overlap_index + 1L
  }
}
overlap_results <- bind_rows(overlap_rows) %>%
  mutate(
    predator_context = factor(predator_context, levels = context_levels),
    pair = paste(class_a, "vs", class_b)
  ) %>%
  arrange(predator_context, pair)

write_csv(camera_context, paste0(output_prefix, "_camera_categories.csv"))
write_csv(event_counts, paste0(output_prefix, "_event_counts.csv"))
write_csv(density_curves, paste0(output_prefix, "_density_curves.csv"))
write_csv(overlap_results, paste0(output_prefix, "_pairwise_overlap.csv"))

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

md <- c(
  "# Four-context gazelle diel activity figure",
  "",
  "## Design",
  "",
  paste0(
    "The figure uses independent 30-min events within documented active camera rounds. Camera timestamps were interpreted as Israel local civil time and converted to apparent solar time. Low and high predator contexts were defined using a ", classification_label, ": jackal cutoffs `",
    fmt(jackal_low), "` and `", fmt(jackal_high), "`; wolf cutoffs `",
    fmt(wolf_low), "` and `", fmt(wolf_high), "`. ", assignment_description
  ),
  "",
  "Kernel curves are normalized within each class and predator context. They describe temporal distribution, not abundance or event rate.",
  "",
  "## Panel support",
  "",
  "| Predator context | Cameras | Fawn events | Female events | Male events |",
  "|---|---:|---:|---:|---:|"
)
for (context_i in context_levels) {
  context_counts <- filter(event_counts, predator_context == context_i)
  md <- c(md, sprintf(
    "| %s | `%s` | `%s` | `%s` | `%s` |",
    context_i,
    unique(context_counts$n_cameras),
    context_counts$n_events[context_counts$class == "Fawn"],
    context_counts$n_events[context_counts$class == "Female"],
    context_counts$n_events[context_counts$class == "Male"]
  ))
}

md <- c(
  md,
  "",
  paste0(
    "Class-context groups with fewer than ", minimum_curve_events,
    " independent events were displayed as event-time ticks rather than assigned a kernel-density curve. Pairwise overlap was not estimated when either class had fewer than ",
    minimum_overlap_events, " events."
  ),
  "",
  "## Pairwise temporal overlap",
  "",
  "Intervals use a camera-cluster bootstrap. Estimates requiring fewer than 20 events in either class were not calculated.",
  "",
  "| Predator context | Pair | Events | Cameras with detections | Overlap [95% CI] | Status |",
  "|---|---|---:|---:|---:|---|"
)
for (i in seq_len(nrow(overlap_results))) {
  overlap_text <- if (overlap_results$estimable[i]) {
    paste0(
      fmt(overlap_results$overlap[i]), " [",
      fmt(overlap_results$conf_low[i]), ", ",
      fmt(overlap_results$conf_high[i]), "]"
    )
  } else {
    "NA"
  }
  status <- ifelse(
    overlap_results$estimable[i],
    "Estimated",
    overlap_results$reason[i]
  )
  md <- c(md, sprintf(
    "| %s | %s | `%s / %s` | `%s / %s` | `%s` | %s |",
    as.character(overlap_results$predator_context[i]),
    overlap_results$pair[i],
    overlap_results$n_events_a[i],
    overlap_results$n_events_b[i],
    overlap_results$n_cameras_a[i],
    overlap_results$n_cameras_b[i],
    overlap_text,
    status
  ))
}

writeLines(md, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved 2x2 predator-context diel outputs with prefix:", output_prefix, "\n")
