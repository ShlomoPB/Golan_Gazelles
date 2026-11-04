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
  output_prefix = "publication_gazelle_predator_overlap_by_open_spaces_quartile",
  local_timezone = "Asia/Jerusalem",
  bootstrap_reps = "500",
  permutation_reps = "999",
  workers = "6",
  minimum_events = "20",
  small_sample_threshold = "75",
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
seed <- as.integer(cfg$seed)
set.seed(seed)

class_levels <- c("Fawn", "Female", "Male")
predator_levels <- c("Wolf", "Jackal")
open_space_levels <- c("Q1: low open-space cover", "Q4: high open-space cover")

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

apparent_solar_hour <- function(timestamp, longitude) {
  day_of_year <- as.integer(format(timestamp, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(timestamp, "%H", tz = "UTC")) +
    as.numeric(format(timestamp, "%M", tz = "UTC")) / 60 +
    as.numeric(format(timestamp, "%S", tz = "UTC")) / 3600
  (utc_hour + longitude / 15 + equation_of_time(day_of_year) / 60) %% 24
}

to_radians <- function(hour) 2 * pi * hour / 24
to_hours <- function(radians) 24 * radians / (2 * pi)

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
    open_space_cover = as.numeric(Open.Spaces)
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
camera_summary <- quartile_cameras %>%
  group_by(open_space_group) %>%
  summarise(
    n_cameras = n(),
    min_open_space_cover = min(open_space_cover),
    max_open_space_cover = max(open_space_cover),
    .groups = "drop"
  ) %>%
  mutate(q25_threshold = q25, q75_threshold = q75)

write_csv(
  camera_groups,
  paste0(cfg$output_prefix, "_camera_groups.csv")
)
write_csv(
  camera_summary,
  paste0(cfg$output_prefix, "_camera_summary.csv")
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
      juvenile = "Fawn",
      female = "Female",
      male = "Male"
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

predator_events <- read_csv(
  cfg$predator_event_file,
  col_types = cols(timestamp = col_character()),
  show_col_types = FALSE
) %>%
  transmute(
    camera_name = as.integer(camera_name),
    predator = taxon,
    solar_hour = as.numeric(solar_hour),
    time_rad = to_radians(solar_hour)
  ) %>%
  inner_join(
    quartile_cameras %>% select(camera_name, open_space_group),
    by = "camera_name"
  ) %>%
  filter(predator %in% predator_levels, is.finite(time_rad)) %>%
  mutate(
    predator = factor(predator, levels = predator_levels),
    open_space_group = factor(open_space_group, levels = open_space_levels)
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
  group_by(predator, open_space_group) %>%
  summarise(
    predator_events = n(),
    predator_cameras = n_distinct(camera_name),
    .groups = "drop"
  ) %>%
  complete(
    predator = factor(predator_levels, levels = predator_levels),
    open_space_group = factor(open_space_levels, levels = open_space_levels),
    fill = list(predator_events = 0L, predator_cameras = 0L)
  )

event_counts <- crossing(
  class = factor(class_levels, levels = class_levels),
  predator = factor(predator_levels, levels = predator_levels),
  open_space_group = factor(open_space_levels, levels = open_space_levels)
) %>%
  left_join(gazelle_counts, by = c("class", "open_space_group")) %>%
  left_join(predator_counts, by = c("predator", "open_space_group"))

write_csv(
  event_counts,
  paste0(cfg$output_prefix, "_event_counts.csv")
)

tasks <- event_counts %>%
  group_by(class, predator) %>%
  summarise(
    minimum_sample = min(c(gazelle_events, predator_events)),
    .groups = "drop"
  ) %>%
  mutate(
    class = as.character(class),
    predator = as.character(predator),
    estimator = ifelse(
      minimum_sample < small_sample_threshold,
      "Dhat1",
      "Dhat4"
    ),
    eligible = minimum_sample >= minimum_events
  )

camera_ids <- lapply(open_space_levels, function(group_i) {
  quartile_cameras$camera_name[quartile_cameras$open_space_group == group_i]
})
names(camera_ids) <- open_space_levels
extreme_cameras <- quartile_cameras$camera_name
n_low <- length(camera_ids[[open_space_levels[1]]])

# Use the same camera-label permutations for all six tests.
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
  predator_i <- tasks$predator[task_i]
  estimator_i <- tasks$estimator[task_i]

  gazelle_i <- gazelle_events %>% filter(class == class_i)
  predator_i_data <- predator_events %>% filter(predator == predator_i)

  observed <- sapply(open_space_levels, function(group_i) {
    safe_overlap(
      gazelle_i$time_rad[gazelle_i$open_space_group == group_i],
      predator_i_data$time_rad[predator_i_data$open_space_group == group_i],
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
            predator_i_data[predator_i_data$open_space_group == group_i, ],
            sampled_cameras
          ),
          estimator_i
        )
      }
    }
  }

  bootstrap_delta <-
    bootstrap_matrix[, open_space_levels[2]] - bootstrap_matrix[, open_space_levels[1]]
  successful_bootstrap <- is.finite(bootstrap_delta)
  enough_bootstrap <- sum(successful_bootstrap) >= ceiling(0.8 * bootstrap_reps)

  if (enough_bootstrap) {
    overlap_low_ci <- unname(quantile(
      bootstrap_matrix[successful_bootstrap, open_space_levels[1]],
      c(0.025, 0.975)
    ))
    overlap_high_ci <- unname(quantile(
      bootstrap_matrix[successful_bootstrap, open_space_levels[2]],
      c(0.025, 0.975)
    ))
    delta_ci <- unname(quantile(
      bootstrap_delta[successful_bootstrap],
      c(0.025, 0.975)
    ))
  } else {
    overlap_low_ci <- overlap_high_ci <- delta_ci <- c(NA_real_, NA_real_)
  }

  permutation_delta <- rep(NA_real_, permutation_reps)
  if (tasks$eligible[task_i]) {
    for (rep_i in seq_len(permutation_reps)) {
      perm_low <- permutation_groups[seq_len(n_low), rep_i]
      perm_high <- permutation_groups[(n_low + 1):length(extreme_cameras), rep_i]
      overlap_low_perm <- safe_overlap(
        gazelle_i$time_rad[gazelle_i$camera_name %in% perm_low],
        predator_i_data$time_rad[predator_i_data$camera_name %in% perm_low],
        estimator_i
      )
      overlap_high_perm <- safe_overlap(
        gazelle_i$time_rad[gazelle_i$camera_name %in% perm_high],
        predator_i_data$time_rad[predator_i_data$camera_name %in% perm_high],
        estimator_i
      )
      permutation_delta[rep_i] <- overlap_high_perm - overlap_low_perm
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

  result_row <- tibble(
    class = class_i,
    predator = predator_i,
    estimator = estimator_i,
    eligible = tasks$eligible[task_i],
    overlap_q1_low_open_spaces = observed[open_space_levels[1]],
    overlap_q4_high_open_spaces = observed[open_space_levels[2]],
    delta_high_minus_low = observed_delta,
    overlap_q1_conf_low = overlap_low_ci[1],
    overlap_q1_conf_high = overlap_low_ci[2],
    overlap_q4_conf_low = overlap_high_ci[1],
    overlap_q4_conf_high = overlap_high_ci[2],
    delta_conf_low = delta_ci[1],
    delta_conf_high = delta_ci[2],
    bootstrap_successes = sum(successful_bootstrap),
    permutation_successes = length(successful_permutations),
    permutation_p = permutation_p
  )

  bootstrap_row <- tibble(
    class = class_i,
    predator = predator_i,
    replicate = seq_len(bootstrap_reps),
    overlap_q1_low_open_spaces = bootstrap_matrix[, open_space_levels[1]],
    overlap_q4_high_open_spaces = bootstrap_matrix[, open_space_levels[2]],
    delta_high_minus_low = bootstrap_delta
  )

  permutation_row <- tibble(
    class = class_i,
    predator = predator_i,
    replicate = seq_len(permutation_reps),
    delta_high_minus_low = permutation_delta
  )

  list(
    result = result_row,
    bootstrap = bootstrap_row,
    permutation = permutation_row
  )
}

worker_count <- max(1L, min(workers, nrow(tasks), parallel::detectCores() - 1L))
message("Running ", nrow(tasks), " contrasts on ", worker_count, " workers.")
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
    "bootstrap_reps", "bootstrap_cameras", "camera_ids", "permutation_reps",
    "permutation_groups", "extreme_cameras", "n_low", "safe_overlap",
    "resample_camera_times"
  ),
  envir = environment()
)
analysis_rows <- tryCatch(
  parallel::parLapply(cluster, seq_len(nrow(tasks)), analyze_task),
  finally = parallel::stopCluster(cluster)
)

result_rows <- lapply(analysis_rows, `[[`, "result")
bootstrap_rows <- lapply(analysis_rows, `[[`, "bootstrap")
permutation_rows <- lapply(analysis_rows, `[[`, "permutation")

for (task_i in seq_along(analysis_rows)) {
  result_i <- analysis_rows[[task_i]]$result
  message(
    result_i$class, " x ", result_i$predator, ": ", result_i$estimator,
    "; low=", fmt(result_i$overlap_q1_low_open_spaces),
    "; high=", fmt(result_i$overlap_q4_high_open_spaces),
    "; high-low=", fmt(result_i$delta_high_minus_low),
    "; permutation p=", fmt_p(result_i$permutation_p)
  )
}

results <- bind_rows(result_rows) %>%
  mutate(
    p_fdr = p.adjust(permutation_p, method = "BH"),
    supported = eligible & is.finite(p_fdr) & p_fdr < 0.05 &
      ((delta_conf_low > 0 & delta_conf_high > 0) |
        (delta_conf_low < 0 & delta_conf_high < 0)),
    direction = case_when(
      !eligible ~ "Data deficient",
      delta_high_minus_low > 0 ~ "Greater overlap at high-open-space cameras",
      delta_high_minus_low < 0 ~ "Lower overlap at high-open-space cameras",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(predator, class)

write_csv(results, paste0(cfg$output_prefix, "_results.csv"))
write_csv(
  bind_rows(bootstrap_rows),
  paste0(cfg$output_prefix, "_bootstrap_draws.csv")
)
write_csv(
  bind_rows(permutation_rows),
  paste0(cfg$output_prefix, "_permutation_draws.csv")
)

# The requested two curves are pointwise shared-density profiles for Q1 and Q4.
# The annotated overlap coefficients remain the inferential summaries.
density_grid <- seq(0, 2 * pi, length.out = 512)
profile_rows <- list()
annotation_rows <- list()
profile_i <- 1L
annotation_i <- 1L

for (task_i in seq_len(nrow(tasks))) {
  class_i <- tasks$class[task_i]
  predator_i <- tasks$predator[task_i]
  estimator_i <- tasks$estimator[task_i]
  adjust_i <- if (estimator_i == "Dhat1") 0.8 else 1
  result_i <- results %>%
    filter(class == class_i, predator == predator_i)

  panel_profiles <- list()
  for (group_i in open_space_levels) {
    gazelle_times <- gazelle_events$time_rad[
      gazelle_events$class == class_i &
        gazelle_events$open_space_group == group_i
    ]
    predator_times <- predator_events$time_rad[
      predator_events$predator == predator_i &
        predator_events$open_space_group == group_i
    ]

    gazelle_bw <- getBandWidth(gazelle_times)
    predator_bw <- getBandWidth(predator_times)
    gazelle_density <- densityFit(
      gazelle_times, density_grid, gazelle_bw / adjust_i
    )
    predator_density <- densityFit(
      predator_times, density_grid, predator_bw / adjust_i
    )
    shared_density <- pmin(gazelle_density, predator_density)

    panel_profiles[[group_i]] <- shared_density
    profile_rows[[profile_i]] <- tibble(
      class = class_i,
      predator = predator_i,
      open_space_group = group_i,
      solar_hour = to_hours(density_grid),
      shared_density = shared_density
    )
    profile_i <- profile_i + 1L
  }

  panel_max <- max(unlist(panel_profiles), na.rm = TRUE)
  counts_i <- event_counts %>%
    filter(class == class_i, predator == predator_i)
  for (group_index in seq_along(open_space_levels)) {
    group_i <- open_space_levels[group_index]
    overlap_i <- if (group_index == 1) {
      result_i$overlap_q1_low_open_spaces
    } else {
      result_i$overlap_q4_high_open_spaces
    }
    count_i <- counts_i %>% filter(open_space_group == group_i)
    annotation_rows[[annotation_i]] <- tibble(
      class = class_i,
      predator = predator_i,
      open_space_group = group_i,
      x = 0.45,
      y = panel_max * (1.00 - 0.14 * (group_index - 1)),
      label = paste0(
        ifelse(group_index == 1, "Q1", "Q4"),
        ": ", estimator_i, " = ", fmt(overlap_i, 2),
        " (n=", count_i$gazelle_events, "/", count_i$predator_events, ")"
      )
    )
    annotation_i <- annotation_i + 1L
  }

  comparison_label <- if (!result_i$eligible) {
    low_count <- counts_i %>%
      filter(open_space_group == open_space_levels[1]) %>%
      pull(gazelle_events)
    paste0("Data deficient for comparison\n(Q1 gazelle n = ", low_count, ")")
  } else {
    paste0(
      "High-low = ", fmt(result_i$delta_high_minus_low, 2),
      "\n95% CI ", fmt(result_i$delta_conf_low, 2), " to ",
      fmt(result_i$delta_conf_high, 2),
      "\nFDR p = ", fmt_p(result_i$p_fdr)
    )
  }

  annotation_rows[[annotation_i]] <- tibble(
    class = class_i,
    predator = predator_i,
    open_space_group = NA_character_,
    x = 23.5,
    y = panel_max,
    label = comparison_label
  )
  annotation_i <- annotation_i + 1L
}

profile_data <- bind_rows(profile_rows) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator = factor(predator, levels = predator_levels),
    open_space_group = factor(open_space_group, levels = open_space_levels)
  )
annotation_data <- bind_rows(annotation_rows) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator = factor(predator, levels = predator_levels),
    open_space_group = factor(open_space_group, levels = open_space_levels)
  )

write_csv(
  profile_data,
  paste0(cfg$output_prefix, "_shared_density_curves.csv")
)

open_space_colors <- c(
  "Q1: low open-space cover" = "#526D82",
  "Q4: high open-space cover" = "#C07A2D"
)
open_space_linetypes <- c(
  "Q1: low open-space cover" = "longdash",
  "Q4: high open-space cover" = "solid"
)

overlap_plot <- ggplot(profile_data) +
  annotate(
    "rect", xmin = 0, xmax = 5.4, ymin = -Inf, ymax = Inf,
    fill = "#DCE3E5", alpha = 0.72
  ) +
  annotate(
    "rect", xmin = 5.4, xmax = 7.0, ymin = -Inf, ymax = Inf,
    fill = "#F5E8D8", alpha = 0.58
  ) +
  annotate(
    "rect", xmin = 17.0, xmax = 18.8, ymin = -Inf, ymax = Inf,
    fill = "#F1DAD5", alpha = 0.58
  ) +
  annotate(
    "rect", xmin = 18.8, xmax = 24, ymin = -Inf, ymax = Inf,
    fill = "#DCE3E5", alpha = 0.72
  ) +
  geom_line(
    aes(
      x = solar_hour, y = shared_density,
      color = open_space_group, linetype = open_space_group
    ),
    linewidth = 1.05
  ) +
  geom_text(
    data = annotation_data %>% filter(!is.na(open_space_group)),
    aes(x = x, y = y, label = label, color = open_space_group),
    hjust = 0, vjust = 1, size = 3.0, show.legend = FALSE
  ) +
  geom_text(
    data = annotation_data %>% filter(is.na(open_space_group)),
    aes(x = x, y = y, label = label),
    hjust = 1, vjust = 1, size = 2.8, color = "#26383A"
  ) +
  facet_grid(predator ~ class, scales = "free_y") +
  scale_color_manual(values = open_space_colors, name = NULL) +
  scale_linetype_manual(values = open_space_linetypes, name = NULL) +
  scale_x_continuous(
    breaks = seq(0, 24, 3),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  labs(
    x = "Apparent solar time",
    y = "Shared gazelle-predator kernel density",
    title = "Gazelle-predator diel overlap across open-space extremes",
    subtitle = paste0(
      "Curves show the pointwise shared density at Q1 and Q4 cameras; ",
      "middle 50% excluded\n",
      "Annotations report overlap coefficients and camera-level high-minus-low tests"
    ),
    caption = paste0(
      "Q1 <= ", fmt(q25), "; Q4 >= ", fmt(q75),
      " for site-level Open.Spaces. Event counts are gazelle/predator. ",
      "Shading indicates approximate night and twilight. Curves describe timing, not event rates."
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#D7DEDC", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "#FBFAF6", color = NA),
    strip.background = element_rect(fill = "#E8E1D4", color = NA),
    strip.text = element_text(face = "bold", size = 10.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#344346"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#4D5B5E", lineheight = 1.05),
    plot.caption = element_text(color = "#5E686A", hjust = 0, size = 8.5),
    legend.position = "top",
    legend.justification = "center",
    panel.spacing = grid::unit(0.7, "lines"),
    plot.margin = margin(10, 14, 10, 10)
  )

ggsave(
  paste0(cfg$output_prefix, ".png"),
  overlap_plot,
  width = 13.2,
  height = 8.4,
  dpi = 400,
  bg = "white"
)
ggsave(
  paste0(cfg$output_prefix, ".pdf"),
  overlap_plot,
  width = 13.2,
  height = 8.4,
  bg = "white"
)

report <- c(
  "# Gazelle-predator diel overlap across open-space quartiles",
  "",
  "## Design and inference",
  "",
  paste0(
    "Cameras were assigned to the bottom (`Q1 <= ", fmt(q25),
    "`) or top (`Q4 >= ", fmt(q75),
    "`) quartile of the site-level `Open.Spaces` cover covariate. The middle 50% of cameras was excluded. Q1 contained `",
    camera_summary$n_cameras[camera_summary$open_space_group == open_space_levels[1]],
    "` cameras and Q4 contained `",
    camera_summary$n_cameras[camera_summary$open_space_group == open_space_levels[2]],
    "` cameras."
  ),
  "",
  paste0(
    "For each gazelle class-predator pair, overlap was estimated separately in Q1 and Q4. A common estimator was used for both strata: `Dhat1` when the smallest of the four contributing event samples contained fewer than ",
    small_sample_threshold, " events and `Dhat4` otherwise. Formal inference required at least ",
    minimum_events, " events in every contributing sample."
  ),
  "",
  paste0(
    "Uncertainty in each overlap and in `Q4 - Q1` was estimated with ",
    bootstrap_reps,
    " camera-cluster bootstrap replicates, retaining all events from each resampled camera. The null distribution for the quartile difference was obtained from ",
    permutation_reps,
    " permutations of the open-space-group labels among the 30 extreme cameras while preserving the 15/15 group sizes. Permutation p-values were adjusted across the eligible tests with the Benjamini-Hochberg false-discovery-rate procedure."
  ),
  "",
  "## Event counts",
  "",
  "| Class | Predator | Open-space group | Gazelle events | Predator events | Gazelle cameras | Predator cameras |",
  "|---|---|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(event_counts))) {
  report <- c(report, sprintf(
    "| %s | %s | %s | %d | %d | %d | %d |",
    event_counts$class[i], event_counts$predator[i],
    event_counts$open_space_group[i], event_counts$gazelle_events[i],
    event_counts$predator_events[i], event_counts$gazelle_cameras[i],
    event_counts$predator_cameras[i]
  ))
}

report <- c(
  report,
  "",
  "## Results",
  "",
  "| Class | Predator | Estimator | Q1 overlap [95% CI] | Q4 overlap [95% CI] | Q4-Q1 [95% CI] | Permutation P | FDR P | Direction |",
  "|---|---|---|---:|---:|---:|---:|---:|---|"
)

for (i in seq_len(nrow(results))) {
  report <- c(report, sprintf(
    "| %s | %s | `%s` | `%s [%s, %s]` | `%s [%s, %s]` | `%s [%s, %s]` | `%s` | `%s` | %s |",
    results$class[i], results$predator[i], results$estimator[i],
    fmt(results$overlap_q1_low_open_spaces[i]),
    fmt(results$overlap_q1_conf_low[i]),
    fmt(results$overlap_q1_conf_high[i]),
    fmt(results$overlap_q4_high_open_spaces[i]),
    fmt(results$overlap_q4_conf_low[i]),
    fmt(results$overlap_q4_conf_high[i]),
    fmt(results$delta_high_minus_low[i]),
    fmt(results$delta_conf_low[i]),
    fmt(results$delta_conf_high[i]),
    fmt_p(results$permutation_p[i]),
    fmt_p(results$p_fdr[i]),
    results$direction[i]
  ))
}

supported_n <- sum(results$supported, na.rm = TRUE)
report <- c(
  report,
  "",
  "## Publication assessment",
  "",
  paste0(
    "Supported high-versus-low open-space differences after FDR correction: `",
    supported_n, "` of `", sum(results$eligible), "` eligible tests."
  ),
  "",
  if (supported_n > 0) {
    "The analysis provides evidence that temporal overlap differs between open-space extremes for at least one class-predator pair. It remains an observational, post hoc quartile analysis and is best treated as a sensitivity or Supporting Information result unless it was an a priori hypothesis."
  } else {
    "The analysis does not provide multiplicity-adjusted evidence that gazelle-predator temporal overlap differs between open-space extremes. The figure may be useful descriptively, but adding it as a main inferential result would overstate the evidence."
  },
  "",
  "Fawn contrasts were retained as descriptive curves only because Q1 contained five fawn events from two cameras; they were excluded from bootstrap and permutation inference.",
  "",
  "The quartile comparison discards the middle 50% of cameras and does not isolate open-space cover from predator activity, camera placement, or other site differences. Overlap coefficients quantify similarity in detection timing and do not demonstrate behavioral avoidance or attraction.",
  "",
  "Overlap coefficients follow Ridout and Linkie (2009) and were calculated with the `overlap` R package.",
  "",
  "## Outputs",
  "",
  paste0("- Main 2 x 3 figure: `", cfg$output_prefix, ".png` and `.pdf`."),
  paste0("- Results: `", cfg$output_prefix, "_results.csv`."),
  paste0("- Event and camera audit: `", cfg$output_prefix, "_event_counts.csv`, `_camera_summary.csv`, and `_camera_groups.csv`."),
  paste0("- Reproducibility draws: `", cfg$output_prefix, "_bootstrap_draws.csv` and `_permutation_draws.csv`.")
)

writeLines(report, paste0(cfg$output_prefix, "_results.md"))
cat("Saved open-space-quartile gazelle-predator overlap analysis:", cfg$output_prefix, "\n")
