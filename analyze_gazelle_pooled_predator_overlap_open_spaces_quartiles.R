#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(overlap)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
cfg <- list(
  site_cov_file = "SiteCovs_060324.csv",
  gazelle_event_file = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  predator_event_file = "publication_gazelle_predator_diel_by_context_median_predator_independent_events.csv",
  output_prefix = "publication_gazelle_pooled_predator_overlap_by_open_spaces_quartile",
  site_covariate = "Open.Spaces",
  stratum_label = "open-space cover",
  stratum_slug = "open-space",
  local_timezone = "Asia/Jerusalem",
  bootstrap_reps = "500",
  permutation_reps = "999",
  workers = "3",
  minimum_events = "20",
  small_sample_threshold = "75",
  guild_gap_minutes = "30",
  seed = "20260819"
)

for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(cfg)) stop("Unknown argument: ", key)
    cfg[[key]] <- parts[2]
  }
}

bootstrap_reps <- as.integer(cfg$bootstrap_reps)
permutation_reps <- as.integer(cfg$permutation_reps)
workers <- as.integer(cfg$workers)
minimum_events <- as.integer(cfg$minimum_events)
small_sample_threshold <- as.integer(cfg$small_sample_threshold)
guild_gap_minutes <- as.numeric(cfg$guild_gap_minutes)
seed <- as.integer(cfg$seed)
set.seed(seed)

class_levels <- c("Male", "Female", "Fawn")
open_space_levels <- c(
  paste0("Q1: low ", cfg$stratum_label),
  paste0("Q4: high ", cfg$stratum_label)
)

required_files <- c(
  cfg$site_cov_file, cfg$gazelle_event_file, cfg$predator_event_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

equation_of_time <- function(day_of_year) {
  b <- 2 * pi * (day_of_year - 81) / 364
  9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
}

parse_local_camera_time <- function(x) {
  as.POSIXct(
    sub("Z$", "", x),
    format = "%Y-%m-%dT%H:%M:%S",
    tz = cfg$local_timezone
  )
}

parse_utc_time <- function(x) {
  as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

apparent_solar_hour <- function(timestamp, longitude) {
  day_of_year <- as.integer(format(timestamp, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(timestamp, "%H", tz = "UTC")) +
    as.numeric(format(timestamp, "%M", tz = "UTC")) / 60 +
    as.numeric(format(timestamp, "%S", tz = "UTC")) / 3600
  (utc_hour + longitude / 15 + equation_of_time(day_of_year) / 60) %% 24
}

to_radians <- function(hour) 2 * pi * hour / 24
to_hours <- function(radians) pmin(24, pmax(0, 24 * radians / (2 * pi)))

safe_overlap <- function(gazelle_times, predator_times, estimator) {
  if (min(length(gazelle_times), length(predator_times)) < 2) {
    return(NA_real_)
  }
  tryCatch(
    unname(overlapEst(gazelle_times, predator_times, type = estimator)),
    error = function(e) NA_real_
  )
}

resample_camera_times <- function(data, sampled_cameras) {
  unlist(lapply(sampled_cameras, function(camera_i) {
    data$time_rad[data$camera_name == camera_i]
  }), use.names = FALSE)
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt(x, 3)))
}

site <- read.csv(
  cfg$site_cov_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    camera_name = as.integer(Camera.ID),
    open_space_cover = as.numeric(.data[[cfg$site_covariate]])
  )

q25 <- unname(quantile(site$open_space_cover, 0.25, type = 7, na.rm = TRUE))
q75 <- unname(quantile(site$open_space_cover, 0.75, type = 7, na.rm = TRUE))

camera_groups <- site %>%
  mutate(
    open_space_group = case_when(
      open_space_cover <= q25 ~ open_space_levels[1],
      open_space_cover >= q75 ~ open_space_levels[2],
      TRUE ~ NA_character_
    ),
    open_space_group = factor(open_space_group, levels = open_space_levels)
  )
quartile_cameras <- camera_groups %>% filter(!is.na(open_space_group))

write_csv(
  camera_groups %>%
    rename(
      site_covariate_value = open_space_cover,
      stratum_group = open_space_group
    ),
  paste0(cfg$output_prefix, "_camera_groups.csv")
)

gazelle_events <- read_csv(
  cfg$gazelle_event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  transmute(
    camera_name = as.integer(camera_name),
    class = recode(
      gazelle_class,
      male = "Male",
      female = "Female",
      juvenile = "Fawn"
    ),
    timestamp = parse_local_camera_time(event_time),
    longitude = as.numeric(camelot_gps_longitude)
  ) %>%
  mutate(
    solar_hour = apparent_solar_hour(timestamp, longitude),
    time_rad = to_radians(solar_hour)
  ) %>%
  inner_join(
    quartile_cameras %>% select(camera_name, open_space_group),
    by = "camera_name"
  ) %>%
  filter(class %in% class_levels, is.finite(time_rad)) %>%
  mutate(
    class = factor(class, levels = class_levels),
    open_space_group = factor(open_space_group, levels = open_space_levels)
  )

predator_species_events <- read_csv(
  cfg$predator_event_file,
  col_types = cols(timestamp = col_character()),
  show_col_types = FALSE
) %>%
  transmute(
    camera_name = as.integer(camera_name),
    session_id = as.character(session_id),
    taxon = as.character(taxon),
    timestamp = parse_utc_time(timestamp),
    solar_hour = as.numeric(solar_hour)
  ) %>%
  filter(taxon %in% c("Wolf", "Jackal"), !is.na(timestamp), is.finite(solar_hour))

# Reapply the independence interval after pooling taxa to define guild events.
predator_guild_events <- predator_species_events %>%
  arrange(camera_name, session_id, timestamp) %>%
  group_by(camera_name, session_id) %>%
  mutate(
    gap_minutes = as.numeric(difftime(timestamp, lag(timestamp), units = "mins")),
    guild_event_id = cumsum(is.na(gap_minutes) | gap_minutes > guild_gap_minutes)
  ) %>%
  group_by(camera_name, session_id, guild_event_id) %>%
  summarise(
    timestamp = first(timestamp),
    solar_hour = first(solar_hour),
    contains_wolf = any(taxon == "Wolf"),
    contains_jackal = any(taxon == "Jackal"),
    source_events = n(),
    .groups = "drop"
  )

guild_audit <- predator_guild_events %>%
  summarise(
    guild_events = n(),
    wolf_only = sum(contains_wolf & !contains_jackal),
    jackal_only = sum(contains_jackal & !contains_wolf),
    both_taxa = sum(contains_wolf & contains_jackal)
  ) %>%
  mutate(
    source_wolf_events = sum(predator_species_events$taxon == "Wolf"),
    source_jackal_events = sum(predator_species_events$taxon == "Jackal"),
    source_total_events = nrow(predator_species_events),
    jackal_source_fraction = source_jackal_events / source_total_events,
    events_removed_by_guild_decluster = source_total_events - guild_events
  )

predator_events <- predator_guild_events %>%
  transmute(
    camera_name,
    session_id,
    guild_event_id,
    timestamp,
    solar_hour,
    time_rad = to_radians(solar_hour),
    contains_wolf,
    contains_jackal,
    source_events
  ) %>%
  inner_join(
    quartile_cameras %>% select(camera_name, open_space_group),
    by = "camera_name"
  ) %>%
  mutate(open_space_group = factor(open_space_group, levels = open_space_levels))

write_csv(guild_audit, paste0(cfg$output_prefix, "_predator_guild_audit.csv"))
write_csv(
  predator_events %>% rename(stratum_group = open_space_group),
  paste0(cfg$output_prefix, "_predator_guild_events.csv")
)

gazelle_counts <- gazelle_events %>%
  group_by(class, open_space_group) %>%
  summarise(
    gazelle_events = n(),
    gazelle_cameras = n_distinct(camera_name),
    .groups = "drop"
  ) %>%
  complete(
    class = factor(class_levels, levels = class_levels),
    open_space_group = factor(open_space_levels, levels = open_space_levels),
    fill = list(gazelle_events = 0L, gazelle_cameras = 0L)
  )

predator_counts <- predator_events %>%
  group_by(open_space_group) %>%
  summarise(
    predator_events = n(),
    predator_cameras = n_distinct(camera_name),
    .groups = "drop"
  )

event_counts <- crossing(
  class = factor(class_levels, levels = class_levels),
  open_space_group = factor(open_space_levels, levels = open_space_levels)
) %>%
  left_join(gazelle_counts, by = c("class", "open_space_group")) %>%
  left_join(predator_counts, by = "open_space_group")

write_csv(
  event_counts %>% rename(stratum_group = open_space_group),
  paste0(cfg$output_prefix, "_event_counts.csv")
)

tasks <- event_counts %>%
  group_by(class) %>%
  summarise(
    minimum_sample = min(c(gazelle_events, predator_events)),
    .groups = "drop"
  ) %>%
  mutate(
    class = as.character(class),
    estimator = ifelse(
      minimum_sample < small_sample_threshold,
      "Dhat1",
      "Dhat4"
    ),
    eligible = minimum_sample >= minimum_events
  ) %>%
  arrange(match(class, class_levels))

camera_ids <- lapply(open_space_levels, function(group_i) {
  quartile_cameras$camera_name[quartile_cameras$open_space_group == group_i]
})
names(camera_ids) <- open_space_levels
extreme_cameras <- quartile_cameras$camera_name
n_low <- length(camera_ids[[open_space_levels[1]]])

permutation_groups <- replicate(
  permutation_reps,
  sample(extreme_cameras, length(extreme_cameras), replace = FALSE)
)
bootstrap_cameras <- lapply(open_space_levels, function(group_i) {
  replicate(
    bootstrap_reps,
    sample(
      camera_ids[[group_i]],
      length(camera_ids[[group_i]]),
      replace = TRUE
    )
  )
})
names(bootstrap_cameras) <- open_space_levels

analyze_task <- function(task_i) {
  class_i <- tasks$class[task_i]
  estimator_i <- tasks$estimator[task_i]
  gazelle_i <- gazelle_events %>% filter(class == class_i)

  observed <- sapply(open_space_levels, function(group_i) {
    safe_overlap(
      gazelle_i$time_rad[gazelle_i$open_space_group == group_i],
      predator_events$time_rad[predator_events$open_space_group == group_i],
      estimator_i
    )
  })
  observed_delta <- observed[open_space_levels[2]] - observed[open_space_levels[1]]

  bootstrap_matrix <- matrix(
    NA_real_, nrow = bootstrap_reps, ncol = length(open_space_levels),
    dimnames = list(NULL, open_space_levels)
  )
  if (tasks$eligible[task_i]) {
    for (rep_i in seq_len(bootstrap_reps)) {
      for (group_i in open_space_levels) {
        sampled_cameras <- bootstrap_cameras[[group_i]][, rep_i]
        bootstrap_matrix[rep_i, group_i] <- safe_overlap(
          resample_camera_times(
            gazelle_i[gazelle_i$open_space_group == group_i, ],
            sampled_cameras
          ),
          resample_camera_times(
            predator_events[predator_events$open_space_group == group_i, ],
            sampled_cameras
          ),
          estimator_i
        )
      }
    }
  }

  bootstrap_delta <-
    bootstrap_matrix[, open_space_levels[2]] -
    bootstrap_matrix[, open_space_levels[1]]
  successful_bootstrap <- is.finite(bootstrap_delta)
  enough_bootstrap <- sum(successful_bootstrap) >= ceiling(0.8 * bootstrap_reps)

  if (enough_bootstrap) {
    low_ci <- unname(quantile(
      bootstrap_matrix[successful_bootstrap, open_space_levels[1]],
      c(0.025, 0.975)
    ))
    high_ci <- unname(quantile(
      bootstrap_matrix[successful_bootstrap, open_space_levels[2]],
      c(0.025, 0.975)
    ))
    delta_ci <- unname(quantile(
      bootstrap_delta[successful_bootstrap],
      c(0.025, 0.975)
    ))
  } else {
    low_ci <- high_ci <- delta_ci <- c(NA_real_, NA_real_)
  }

  permutation_delta <- rep(NA_real_, permutation_reps)
  if (tasks$eligible[task_i]) {
    for (rep_i in seq_len(permutation_reps)) {
      perm_low <- permutation_groups[seq_len(n_low), rep_i]
      perm_high <- permutation_groups[(n_low + 1):length(extreme_cameras), rep_i]
      low_perm <- safe_overlap(
        gazelle_i$time_rad[gazelle_i$camera_name %in% perm_low],
        predator_events$time_rad[predator_events$camera_name %in% perm_low],
        estimator_i
      )
      high_perm <- safe_overlap(
        gazelle_i$time_rad[gazelle_i$camera_name %in% perm_high],
        predator_events$time_rad[predator_events$camera_name %in% perm_high],
        estimator_i
      )
      permutation_delta[rep_i] <- high_perm - low_perm
    }
  }

  successful_permutations <- permutation_delta[is.finite(permutation_delta)]
  permutation_p <- if (
    is.finite(observed_delta) &&
      length(successful_permutations) >= ceiling(0.8 * permutation_reps)
  ) {
    (sum(abs(successful_permutations) >= abs(observed_delta)) + 1) /
      (length(successful_permutations) + 1)
  } else {
    NA_real_
  }

  list(
    result = tibble(
      class = class_i,
      estimator = estimator_i,
      eligible = tasks$eligible[task_i],
      overlap_q1_low = observed[open_space_levels[1]],
      overlap_q4_high = observed[open_space_levels[2]],
      delta_high_minus_low = observed_delta,
      overlap_q1_conf_low = low_ci[1],
      overlap_q1_conf_high = low_ci[2],
      overlap_q4_conf_low = high_ci[1],
      overlap_q4_conf_high = high_ci[2],
      delta_conf_low = delta_ci[1],
      delta_conf_high = delta_ci[2],
      bootstrap_successes = sum(successful_bootstrap),
      permutation_successes = length(successful_permutations),
      permutation_p = permutation_p
    ),
    bootstrap = tibble(
      class = class_i,
      replicate = seq_len(bootstrap_reps),
      overlap_q1_low = bootstrap_matrix[, open_space_levels[1]],
      overlap_q4_high = bootstrap_matrix[, open_space_levels[2]],
      delta_high_minus_low = bootstrap_delta
    ),
    permutation = tibble(
      class = class_i,
      replicate = seq_len(permutation_reps),
      delta_high_minus_low = permutation_delta
    )
  )
}

worker_count <- max(1L, min(workers, nrow(tasks), parallel::detectCores() - 1L))
message("Running ", nrow(tasks), " pooled-predator contrasts on ", worker_count, " workers.")
cluster <- parallel::makeCluster(worker_count)
invisible(parallel::clusterEvalQ(cluster, {
  suppressPackageStartupMessages({
    library(dplyr)
    library(overlap)
    library(tibble)
  })
}))
parallel::clusterExport(
  cluster,
  c(
    "tasks", "gazelle_events", "predator_events", "open_space_levels",
    "bootstrap_reps", "bootstrap_cameras", "permutation_reps",
    "permutation_groups", "extreme_cameras", "n_low", "safe_overlap",
    "resample_camera_times"
  ),
  envir = environment()
)
analysis_rows <- tryCatch(
  parallel::parLapply(cluster, seq_len(nrow(tasks)), analyze_task),
  finally = parallel::stopCluster(cluster)
)

results <- bind_rows(lapply(analysis_rows, `[[`, "result")) %>%
  mutate(
    p_fdr = p.adjust(permutation_p, method = "BH"),
    supported = eligible & is.finite(p_fdr) & p_fdr < 0.05 &
      ((delta_conf_low > 0 & delta_conf_high > 0) |
        (delta_conf_low < 0 & delta_conf_high < 0)),
    class = factor(class, levels = class_levels)
  ) %>%
  arrange(class)

write_csv(results, paste0(cfg$output_prefix, "_results.csv"))
write_csv(
  bind_rows(lapply(analysis_rows, `[[`, "bootstrap")),
  paste0(cfg$output_prefix, "_bootstrap_draws.csv")
)
write_csv(
  bind_rows(lapply(analysis_rows, `[[`, "permutation")),
  paste0(cfg$output_prefix, "_permutation_draws.csv")
)

density_grid <- seq(0, 2 * pi, length.out = 512)
line_rows <- list()
overlap_rows <- list()
annotation_rows <- list()
row_i <- 1L

for (task_i in seq_len(nrow(tasks))) {
  class_i <- tasks$class[task_i]
  estimator_i <- tasks$estimator[task_i]
  adjust_i <- if (estimator_i == "Dhat1") 0.8 else 1
  result_i <- results %>% filter(class == class_i)

  for (group_i in open_space_levels) {
    gazelle_times <- gazelle_events$time_rad[
      gazelle_events$class == class_i &
        gazelle_events$open_space_group == group_i
    ]
    predator_times <- predator_events$time_rad[
      predator_events$open_space_group == group_i
    ]
    gazelle_density <- densityFit(
      gazelle_times,
      density_grid,
      getBandWidth(gazelle_times) / adjust_i
    )
    predator_density <- densityFit(
      predator_times,
      density_grid,
      getBandWidth(predator_times) / adjust_i
    )
    overlap_density <- pmin(gazelle_density, predator_density)

    line_rows[[row_i]] <- bind_rows(
      tibble(
        class = class_i,
        open_space_group = group_i,
        solar_hour = to_hours(density_grid),
        curve = "Gazelle class",
        density = gazelle_density
      ),
      tibble(
        class = class_i,
        open_space_group = group_i,
        solar_hour = to_hours(density_grid),
        curve = "Pooled predators",
        density = predator_density
      )
    )
    overlap_rows[[row_i]] <- tibble(
      class = class_i,
      open_space_group = group_i,
      solar_hour = to_hours(density_grid),
      overlap_density = overlap_density
    )

    count_i <- event_counts %>%
      filter(class == class_i, open_space_group == group_i)
    overlap_i <- if (group_i == open_space_levels[1]) {
      result_i$overlap_q1_low
    } else {
      result_i$overlap_q4_high
    }
    panel_max <- max(c(gazelle_density, predator_density), na.rm = TRUE)
    annotation_rows[[row_i]] <- tibble(
      class = class_i,
      open_space_group = group_i,
      x = 23.5,
      y = panel_max * 0.98,
      label = paste0(
        estimator_i, " = ", fmt(overlap_i, 2),
        "\nGazelle n = ", count_i$gazelle_events,
        "; predator n = ", count_i$predator_events,
        ifelse(!result_i$eligible, "\nQ1-Q4 test data deficient", "")
      )
    )
    row_i <- row_i + 1L
  }
}

line_data <- bind_rows(line_rows) %>%
  mutate(
    class = factor(class, levels = class_levels),
    open_space_group = factor(open_space_group, levels = open_space_levels),
    curve = factor(curve, levels = c("Gazelle class", "Pooled predators"))
  )
overlap_data <- bind_rows(overlap_rows) %>%
  mutate(
    class = factor(class, levels = class_levels),
    open_space_group = factor(open_space_group, levels = open_space_levels)
  )
annotation_data <- bind_rows(annotation_rows) %>%
  mutate(
    class = factor(class, levels = class_levels),
    open_space_group = factor(open_space_group, levels = open_space_levels)
  )

write_csv(
  line_data %>% rename(stratum_group = open_space_group),
  paste0(cfg$output_prefix, "_density_curves.csv")
)

curve_colors <- c(
  "Gazelle class" = "#087F6D",
  "Pooled predators" = "#A33E2E"
)
curve_linetypes <- c(
  "Gazelle class" = "solid",
  "Pooled predators" = "longdash"
)

activity_plot <- ggplot() +
  annotate(
    "rect", xmin = 0, xmax = 5.4, ymin = -Inf, ymax = Inf,
    fill = "#DDE4E6", alpha = 0.55
  ) +
  annotate(
    "rect", xmin = 5.4, xmax = 7.0, ymin = -Inf, ymax = Inf,
    fill = "#F3E6D7", alpha = 0.40
  ) +
  annotate(
    "rect", xmin = 17.0, xmax = 18.8, ymin = -Inf, ymax = Inf,
    fill = "#F0DCD7", alpha = 0.40
  ) +
  annotate(
    "rect", xmin = 18.8, xmax = 24, ymin = -Inf, ymax = Inf,
    fill = "#DDE4E6", alpha = 0.55
  ) +
  geom_area(
    data = line_data,
    aes(x = solar_hour, y = density, fill = curve, group = curve),
    stat = "identity",
    position = "identity",
    alpha = 0.12,
    show.legend = FALSE
  ) +
  geom_ribbon(
    data = overlap_data,
    aes(x = solar_hour, ymin = 0, ymax = overlap_density),
    fill = "#D9A441",
    alpha = 0.30
  ) +
  geom_line(
    data = line_data,
    aes(
      x = solar_hour,
      y = density,
      color = curve,
      linetype = curve,
      group = curve
    ),
    linewidth = 1.05
  ) +
  geom_text(
    data = annotation_data,
    aes(x = x, y = y, label = label),
    hjust = 1,
    vjust = 1,
    size = 3.0,
    color = "#26383A"
  ) +
  # A common y-scale is essential because overlap is an area comparison.
  facet_grid(open_space_group ~ class) +
  scale_color_manual(values = curve_colors, name = NULL) +
  scale_fill_manual(values = curve_colors, name = NULL) +
  scale_linetype_manual(values = curve_linetypes, name = NULL) +
  scale_x_continuous(
    breaks = seq(0, 24, 3),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  labs(
    x = "Apparent solar time",
    y = "Activity density",
    title = paste0(
      "Gazelle overlap with pooled predator activity across ",
      cfg$stratum_slug,
      " extremes"
    ),
    subtitle = paste0(
      "Top row: Q1 (low ", cfg$stratum_label,
      "); bottom row: Q4 (high ", cfg$stratum_label, ")\n",
      "Shaded gold area is the shared portion of the two activity distributions"
    ),
    caption = paste0(
      "Predator activity pools wolf and jackal detections and reapplies a 30-min guild-level independence interval. ",
      "Jackals comprise ", fmt(100 * guild_audit$jackal_source_fraction, 1),
      "% of source predator events. Curves are normalized within each panel and describe timing, not event rates."
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#D5DCDA", linewidth = 0.3),
    panel.background = element_rect(fill = "#FCFBF7", color = NA),
    strip.background = element_rect(fill = "#E8E0D2", color = NA),
    strip.text = element_text(face = "bold", size = 10.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#344346"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4D5B5E", lineheight = 1.05),
    plot.caption = element_text(color = "#5E686A", hjust = 0, size = 8.2),
    legend.position = "top",
    legend.justification = "center",
    panel.spacing = grid::unit(0.65, "lines"),
    plot.margin = margin(10, 14, 10, 10)
  )

ggsave(
  paste0(cfg$output_prefix, ".png"),
  activity_plot,
  width = 13.0,
  height = 8.3,
  dpi = 400,
  bg = "white"
)
ggsave(
  paste0(cfg$output_prefix, ".pdf"),
  activity_plot,
  width = 13.0,
  height = 8.3,
  bg = "white"
)

fawn_result <- results %>% filter(class == "Fawn")
fawn_minimum <- event_counts %>%
  filter(class == "Fawn") %>%
  summarise(minimum = min(gazelle_events)) %>%
  pull(minimum)
sample_note <- if (fawn_result$eligible) {
  "All three class comparisons met the prespecified minimum event threshold."
} else {
  paste0(
    "The fawn comparison is descriptive only because its smallest stratum contained ",
    fawn_minimum,
    " fawn events."
  )
}

report <- c(
  paste0(
    "# Gazelle overlap with pooled predator activity across ",
    cfg$stratum_slug,
    " quartiles"
  ),
  "",
  "## Design",
  "",
  paste0(
    "The top row represents the ", length(camera_ids[[open_space_levels[1]]]),
    " cameras in Q1 (`", cfg$site_covariate, " <= ",
    fmt(q25), "`) and the bottom row represents the ",
    length(camera_ids[[open_space_levels[2]]]),
    " cameras in Q4 (`", cfg$site_covariate, " >= ",
    fmt(q75), "`). The middle ", nrow(site) - length(extreme_cameras),
    " cameras were excluded."
  ),
  "",
  paste0(
    "Wolf and jackal records were pooled and then declustered within camera session using a ",
    guild_gap_minutes,
    "-min interval across taxa. This produced `",
    guild_audit$guild_events,
    "` guild events from `",
    guild_audit$source_total_events,
    "` species-specific independent events. Jackals comprised `",
    fmt(100 * guild_audit$jackal_source_fraction, 1),
    "%` of the source events, so the pooled curve is frequency weighted and primarily reflects jackal timing."
  ),
  "",
  paste0(
    "Overlap was estimated with `Dhat1` for fawns and `Dhat4` for adults. Q4-Q1 uncertainty used ",
    bootstrap_reps,
    " camera-cluster bootstrap replicates, and null differences used ",
    permutation_reps,
    " camera-label permutations. FDR correction was applied across eligible class comparisons."
  ),
  "",
  "## Results",
  "",
  "| Class | Estimator | Q1 overlap [95% CI] | Q4 overlap [95% CI] | Q4-Q1 [95% CI] | Permutation P | FDR P |",
  "|---|---|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(results))) {
  report <- c(report, sprintf(
    "| %s | `%s` | `%s [%s, %s]` | `%s [%s, %s]` | `%s [%s, %s]` | `%s` | `%s` |",
    results$class[i], results$estimator[i],
    fmt(results$overlap_q1_low[i]),
    fmt(results$overlap_q1_conf_low[i]),
    fmt(results$overlap_q1_conf_high[i]),
    fmt(results$overlap_q4_high[i]),
    fmt(results$overlap_q4_conf_low[i]),
    fmt(results$overlap_q4_conf_high[i]),
    fmt(results$delta_high_minus_low[i]),
    fmt(results$delta_conf_low[i]),
    fmt(results$delta_conf_high[i]),
    fmt_p(results$permutation_p[i]),
    fmt_p(results$p_fdr[i])
  ))
}

supported_n <- sum(results$supported, na.rm = TRUE)
report <- c(
  report,
  "",
  paste0(
    "Supported Q4-Q1 differences after FDR correction: `",
    supported_n,
    "` of `",
    sum(results$eligible),
    "` eligible comparisons."
  ),
  "",
  sample_note,
  "",
  paste0(
    "Overlap describes temporal coincidence at sampled cameras. It does not isolate ",
    cfg$stratum_label,
    " effects from differences in predator composition, abundance, camera placement, or other site characteristics."
  ),
  "",
  "## Outputs",
  "",
  paste0("- Figure: `", cfg$output_prefix, ".png` and `.pdf`."),
  paste0("- Results: `", cfg$output_prefix, "_results.csv`."),
  paste0("- Density curves: `", cfg$output_prefix, "_density_curves.csv`."),
  paste0("- Event audit: `", cfg$output_prefix, "_event_counts.csv` and `_predator_guild_audit.csv`.")
)

writeLines(report, paste0(cfg$output_prefix, "_results.md"))
cat("Saved pooled-predator ", cfg$stratum_slug, " overlap analysis: ", cfg$output_prefix, "\n", sep = "")
