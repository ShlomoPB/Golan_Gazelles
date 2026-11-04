#!/usr/bin/env Rscript

# Overall gazelle-class activity distributions using the retained local-clock
# interpretation and date-resolved active camera-rounds.

suppressPackageStartupMessages({
  library(activity)
  library(dplyr)
  library(ggplot2)
  library(overlap)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  event_file = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  output_prefix = "publication_overall_gazelle_diel_activity",
  local_timezone = "Asia/Jerusalem",
  curve_reps = "199",
  bootstrap_reps = "500",
  seed = "20260820"
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
obs_cov_file <- arg_map$obs_cov_file
output_prefix <- arg_map$output_prefix
local_timezone <- arg_map$local_timezone
curve_reps <- as.integer(arg_map$curve_reps)
bootstrap_reps <- as.integer(arg_map$bootstrap_reps)
seed <- as.integer(arg_map$seed)

set.seed(seed)

class_levels <- c("Fawn", "Female", "Male")
class_labels <- c(juvenile = "Fawn", female = "Female", male = "Male")
class_colors <- c(Fawn = "#138F7A", Female = "#B72A32", Male = "#173F66")

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

assign_rounds <- function(events, active_obs) {
  event_index <- events %>%
    transmute(event_row_id = row_number(), camera_name, event_date)
  matches <- event_index %>%
    inner_join(
      active_obs %>% select(camera_name, round, start_date, end_date),
      by = "camera_name",
      relationship = "many-to-many"
    ) %>%
    filter(event_date >= start_date, event_date <= end_date) %>%
    group_by(event_row_id) %>%
    # A shared boundary date is assigned to the later-starting round.
    slice_max(start_date, n = 1, with_ties = FALSE) %>%
    ungroup()
  matches$round[match(seq_len(nrow(events)), matches$event_row_id)]
}

extract_density <- function(fit, class_label) {
  density <- as.data.frame(fit@pdf)
  tibble(
    class = class_label,
    solar_hour = to_hours(density$x),
    density = density$y
  )
}

pick_overlap_estimator <- function(n_a, n_b) {
  if (min(n_a, n_b) < 75) "Dhat1" else "Dhat4"
}

clustered_overlap <- function(events, class_a, class_b, cameras, reps) {
  time_a <- events$time_rad[events$class == class_a]
  time_b <- events$time_rad[events$class == class_b]
  estimator <- pick_overlap_estimator(length(time_a), length(time_b))
  point <- unname(overlapEst(time_a, time_b, type = estimator))

  boot <- replicate(reps, {
    sampled_cameras <- sample(cameras, length(cameras), replace = TRUE)
    boot_a <- unlist(lapply(sampled_cameras, function(camera_i) {
      events$time_rad[events$camera_name == camera_i & events$class == class_a]
    }), use.names = FALSE)
    boot_b <- unlist(lapply(sampled_cameras, function(camera_i) {
      events$time_rad[events$camera_name == camera_i & events$class == class_b]
    }), use.names = FALSE)
    if (length(boot_a) < 5 || length(boot_b) < 5) return(NA_real_)
    tryCatch(
      unname(overlapEst(boot_a, boot_b, type = estimator)),
      error = function(e) NA_real_
    )
  })
  boot <- boot[is.finite(boot)]
  if (length(boot) < ceiling(0.8 * reps)) {
    stop("Too few successful bootstrap replicates for ", class_a, " and ", class_b)
  }

  tibble(
    class_1 = class_a,
    class_2 = class_b,
    n_events_1 = length(time_a),
    n_events_2 = length(time_b),
    cameras_with_events_1 = n_distinct(
      events$camera_name[events$class == class_a]
    ),
    cameras_with_events_2 = n_distinct(
      events$camera_name[events$class == class_b]
    ),
    estimator,
    overlap = point,
    conf_low = unname(quantile(boot, 0.025)),
    conf_high = unname(quantile(boot, 0.975)),
    bootstrap_successes = length(boot)
  )
}

obs <- read_csv(obs_cov_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    days_operational = as.numeric(Days.Operational)
  )

active_obs <- obs %>%
  filter(
    is.finite(days_operational),
    days_operational > 0,
    !is.na(start_date),
    !is.na(end_date)
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

raw_events$round <- assign_rounds(raw_events, active_obs)
events <- raw_events %>%
  filter(!is.na(round), class %in% class_levels) %>%
  mutate(class = factor(class, levels = class_levels))

event_audit <- tibble(
  metric = c(
    "prepared_event_rows",
    "retained_active_round_event_rows",
    "excluded_event_rows",
    "active_camera_rounds",
    "date_resolved_camera_days",
    "timestamp_interpretation"
  ),
  value = c(
    nrow(raw_events),
    nrow(events),
    nrow(raw_events) - nrow(events),
    nrow(active_obs),
    sum(active_obs$days_operational),
    paste0("Camera clock parsed as ", local_timezone)
  )
)

event_counts <- events %>%
  group_by(class) %>%
  summarise(
    independent_events = n(),
    summed_classified_animals = sum(event_count, na.rm = TRUE),
    cameras_with_events = n_distinct(camera_name),
    .groups = "drop"
  )

density_results <- bind_rows(lapply(class_levels, function(class_i) {
  class_times <- events$time_rad[events$class == class_i]
  fit <- fitact(class_times, reps = curve_reps, sample = "data", show = FALSE)
  extract_density(fit, class_i)
}))

class_pairs <- combn(class_levels, 2, simplify = FALSE)
camera_ids <- sort(unique(active_obs$camera_name))
overlap_results <- bind_rows(lapply(class_pairs, function(pair_i) {
  clustered_overlap(
    events,
    class_a = pair_i[1],
    class_b = pair_i[2],
    cameras = camera_ids,
    reps = bootstrap_reps
  )
}))

shade <- tibble(
  xmin = c(0, 5, 7, 17, 19),
  xmax = c(5, 7, 17, 19, 24),
  light_period = c("Night", "Twilight", "Daylight", "Twilight", "Night")
)

sample_label <- paste(
  paste0(event_counts$class, ": n=", event_counts$independent_events),
  collapse = "   "
)

plot <- ggplot(density_results, aes(solar_hour, density, color = class)) +
  geom_rect(
    data = shade,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = light_period),
    alpha = 0.16,
    color = NA
  ) +
  geom_line(linewidth = 1.15) +
  scale_color_manual(values = class_colors, name = NULL) +
  scale_fill_manual(
    values = c(Night = "#607482", Twilight = "#E6A66C", Daylight = "#F4ECD8"),
    guide = "none"
  ) +
  scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24)) +
  labs(
    x = "Apparent solar time",
    y = "Relative activity density",
    title = "Diel activity of gazelle demographic classes",
    subtitle = sample_label,
    caption = paste0(
      "Curves are normalized distributions, not event rates.\n",
      "Camera clocks were interpreted as Israel local civil time and converted to apparent solar time; ",
      "shading indicates approximate light periods."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#46555B"),
    plot.caption = element_text(hjust = 0, color = "#58676C", lineheight = 1.05),
    legend.position = "top"
  )

write_csv(event_audit, paste0(output_prefix, "_event_audit.csv"))
write_csv(event_counts, paste0(output_prefix, "_event_counts.csv"))
write_csv(density_results, paste0(output_prefix, "_density_curves.csv"))
write_csv(overlap_results, paste0(output_prefix, "_class_overlap.csv"))

ggsave(
  paste0(output_prefix, ".png"),
  plot,
  width = 9.2,
  height = 5.5,
  dpi = 320
)
ggsave(
  paste0(output_prefix, ".pdf"),
  plot,
  width = 9.2,
  height = 5.5
)

cat("Saved corrected overall gazelle diel outputs with prefix:", output_prefix, "\n")
