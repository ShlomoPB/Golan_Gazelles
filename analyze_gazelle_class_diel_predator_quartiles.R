#!/usr/bin/env Rscript

# Compare class-specific gazelle diel activity between the lowest and highest
# quartiles of site-level wolf and jackal activity.

suppressPackageStartupMessages({
  library(activity)
  library(dplyr)
  library(ggplot2)
  library(overlap)
  library(readr)
  library(tidyr)
})

set.seed(20260819)

event_file <- "Gazella gazella_independent_events_30min_analysis_ready.csv"
site_file <- "SiteCovs_060324.csv"
output_dir <- file.path("results", "H3_predator_quartiles")
local_timezone <- "Asia/Jerusalem"
bootstrap_reps <- 500L
randomization_reps <- 999L
minimum_test_events <- 20L
small_sample_threshold <- 75L
activity_fit_reps <- 99L
parallel_workers <- 3L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(event_file, site_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required inputs: ", paste(missing_files, collapse = ", "))
}

class_levels <- c("Fawn", "Female", "Male")
predator_levels <- c("Wolf", "Jackal")
quartile_levels <- c("Q1: low activity", "Q4: high activity")

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

safe_overlap <- function(times_q1, times_q4, estimator) {
  if (min(length(times_q1), length(times_q4)) < 2) return(NA_real_)
  tryCatch(
    unname(overlap::overlapEst(times_q1, times_q4, type = estimator)),
    error = function(e) NA_real_
  )
}

resample_camera_values <- function(data, sampled_cameras, variable) {
  unlist(lapply(sampled_cameras, function(camera_i) {
    data[[variable]][data$camera_name == camera_i]
  }), use.names = FALSE)
}

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt(x, 3)))
}

site <- read_csv(site_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(Camera.ID),
    Wolf = as.numeric(Null.Wolf.Abundance),
    Jackal = as.numeric(null.jackal.abundance)
  )

if (nrow(site) != 60L || n_distinct(site$camera_name) != 60L) {
  stop("Site covariate file does not contain 60 unique cameras.")
}

site_long <- site %>%
  pivot_longer(
    cols = all_of(predator_levels),
    names_to = "predator",
    values_to = "activity_index"
  ) %>%
  group_by(predator) %>%
  mutate(
    q25_threshold = quantile(activity_index, 0.25, type = 7, na.rm = TRUE),
    q75_threshold = quantile(activity_index, 0.75, type = 7, na.rm = TRUE),
    quartile = case_when(
      activity_index <= q25_threshold ~ quartile_levels[1],
      activity_index >= q75_threshold ~ quartile_levels[2],
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    predator = factor(predator, levels = predator_levels),
    quartile = factor(quartile, levels = quartile_levels)
  )

camera_groups <- site_long %>% filter(!is.na(quartile))
camera_summary <- camera_groups %>%
  group_by(predator, quartile) %>%
  summarise(
    cameras = n(),
    activity_min = min(activity_index),
    activity_max = max(activity_index),
    q25_threshold = first(q25_threshold),
    q75_threshold = first(q75_threshold),
    .groups = "drop"
  )

if (any(camera_summary$cameras != 15L)) {
  stop("Quartile thresholds did not produce exactly 15 Q1 and 15 Q4 cameras per predator.")
}

write_csv(site_long, file.path(output_dir, "H3_predator_quartile_camera_groups.csv"))
write_csv(camera_summary, file.path(output_dir, "H3_predator_quartile_camera_summary.csv"))

events <- read_csv(
  event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  transmute(
    camera_name = as.integer(camera_name),
    class = recode(
      as.character(gazelle_class),
      juvenile = "Fawn",
      female = "Female",
      male = "Male"
    ),
    timestamp = parse_local_camera_time(event_time),
    longitude = as.numeric(camelot_gps_longitude)
  ) %>%
  mutate(
    solar_hour = apparent_solar_hour(timestamp, longitude),
    time_rad = to_radians(solar_hour),
    nocturnality_hours = abs(solar_hour - 12),
    class = factor(class, levels = class_levels)
  ) %>%
  filter(class %in% class_levels, is.finite(time_rad), is.finite(nocturnality_hours))

expected_counts <- c(Fawn = 180L, Female = 1477L, Male = 2268L)
observed_counts <- table(events$class)
if (
  nrow(events) != sum(expected_counts) ||
    any(as.integer(observed_counts[names(expected_counts)]) != expected_counts)
) {
  stop("Final gazelle event totals do not match the validated analysis-ready data.")
}

analysis_events <- events %>%
  inner_join(
    camera_groups %>% select(camera_name, predator, quartile),
    by = "camera_name",
    relationship = "many-to-many"
  ) %>%
  mutate(
    predator = factor(predator, levels = predator_levels),
    quartile = factor(quartile, levels = quartile_levels)
  )

event_counts <- analysis_events %>%
  group_by(predator, class, quartile) %>%
  summarise(
    events = n(),
    cameras_with_events = n_distinct(camera_name),
    mean_nocturnality_hours = mean(nocturnality_hours),
    median_nocturnality_hours = median(nocturnality_hours),
    .groups = "drop"
  ) %>%
  complete(
    predator = factor(predator_levels, levels = predator_levels),
    class = factor(class_levels, levels = class_levels),
    quartile = factor(quartile_levels, levels = quartile_levels),
    fill = list(events = 0L, cameras_with_events = 0L)
  ) %>%
  left_join(
    camera_summary %>% select(predator, quartile, cameras_assigned = cameras),
    by = c("predator", "quartile")
  )

write_csv(event_counts, file.path(output_dir, "H3_predator_quartile_event_counts.csv"))

tasks <- crossing(
  predator = predator_levels,
  class = class_levels
) %>%
  left_join(
    event_counts %>%
      group_by(predator, class) %>%
      summarise(minimum_events = min(events), .groups = "drop"),
    by = c("predator", "class")
  ) %>%
  mutate(
    estimator = ifelse(minimum_events < small_sample_threshold, "Dhat1", "Dhat4"),
    formal_test_eligible = minimum_events >= minimum_test_events
  ) %>%
  arrange(match(predator, predator_levels), match(class, class_levels))

activity_fits <- vector("list", nrow(tasks))
fit_index <- 1L
for (task_i in seq_len(nrow(tasks))) {
  predator_i <- tasks$predator[task_i]
  class_i <- tasks$class[task_i]
  fit_groups <- vector("list", length(quartile_levels))
  names(fit_groups) <- quartile_levels
  for (quartile_i in quartile_levels) {
    times_i <- analysis_events$time_rad[
      analysis_events$predator == predator_i &
        analysis_events$class == class_i &
        analysis_events$quartile == quartile_i
    ]
    fit_groups[[quartile_i]] <- if (length(times_i) >= 5L) {
      activity::fitact(
        times_i,
        reps = activity_fit_reps,
        sample = "data",
        show = FALSE
      )
    } else {
      NULL
    }
  }
  activity_fits[[fit_index]] <- fit_groups
  fit_index <- fit_index + 1L
}

analyze_task <- function(task_i) {
  predator_i <- tasks$predator[task_i]
  class_i <- tasks$class[task_i]
  estimator_i <- tasks$estimator[task_i]
  eligible_i <- tasks$formal_test_eligible[task_i]
  event_i <- analysis_events %>%
    filter(predator == predator_i, class == class_i)

  group_data <- lapply(quartile_levels, function(quartile_i) {
    event_i %>% filter(quartile == quartile_i)
  })
  names(group_data) <- quartile_levels
  camera_ids <- lapply(quartile_levels, function(quartile_i) {
    camera_groups$camera_name[
      camera_groups$predator == predator_i & camera_groups$quartile == quartile_i
    ]
  })
  names(camera_ids) <- quartile_levels

  times_q1 <- group_data[[quartile_levels[1]]]$time_rad
  times_q4 <- group_data[[quartile_levels[2]]]$time_rad
  observed_overlap <- safe_overlap(times_q1, times_q4, estimator_i)
  mean_q1 <- mean(group_data[[quartile_levels[1]]]$nocturnality_hours)
  mean_q4 <- mean(group_data[[quartile_levels[2]]]$nocturnality_hours)
  nocturnality_delta <- mean_q4 - mean_q1

  bootstrap_overlap <- rep(NA_real_, bootstrap_reps)
  bootstrap_mean_q1 <- rep(NA_real_, bootstrap_reps)
  bootstrap_mean_q4 <- rep(NA_real_, bootstrap_reps)
  bootstrap_nocturnality_delta <- rep(NA_real_, bootstrap_reps)

  for (rep_i in seq_len(bootstrap_reps)) {
    sampled_q1 <- sample(camera_ids[[quartile_levels[1]]], 15L, replace = TRUE)
    sampled_q4 <- sample(camera_ids[[quartile_levels[2]]], 15L, replace = TRUE)
    times_q1_i <- resample_camera_values(
      group_data[[quartile_levels[1]]], sampled_q1, "time_rad"
    )
    times_q4_i <- resample_camera_values(
      group_data[[quartile_levels[2]]], sampled_q4, "time_rad"
    )
    noct_q1_i <- resample_camera_values(
      group_data[[quartile_levels[1]]], sampled_q1, "nocturnality_hours"
    )
    noct_q4_i <- resample_camera_values(
      group_data[[quartile_levels[2]]], sampled_q4, "nocturnality_hours"
    )
    bootstrap_overlap[rep_i] <- safe_overlap(times_q1_i, times_q4_i, estimator_i)
    if (length(noct_q1_i)) bootstrap_mean_q1[rep_i] <- mean(noct_q1_i)
    if (length(noct_q4_i)) bootstrap_mean_q4[rep_i] <- mean(noct_q4_i)
    bootstrap_nocturnality_delta[rep_i] <-
      bootstrap_mean_q4[rep_i] - bootstrap_mean_q1[rep_i]
  }

  overlap_success <- is.finite(bootstrap_overlap)
  nocturnality_success <- is.finite(bootstrap_nocturnality_delta)
  enough_overlap <- sum(overlap_success) >= ceiling(0.8 * bootstrap_reps)
  enough_nocturnality <- sum(nocturnality_success) >= ceiling(0.8 * bootstrap_reps)
  overlap_ci <- if (enough_overlap) {
    unname(quantile(bootstrap_overlap[overlap_success], c(0.025, 0.975)))
  } else {
    c(NA_real_, NA_real_)
  }
  nocturnality_ci <- if (enough_nocturnality) {
    unname(quantile(
      bootstrap_nocturnality_delta[nocturnality_success],
      c(0.025, 0.975)
    ))
  } else {
    c(NA_real_, NA_real_)
  }

  compare_result <- if (
    eligible_i &&
      !is.null(activity_fits[[task_i]][[quartile_levels[1]]]) &&
      !is.null(activity_fits[[task_i]][[quartile_levels[2]]])
  ) {
    tryCatch(
      activity::compareCkern(
        activity_fits[[task_i]][[quartile_levels[1]]],
        activity_fits[[task_i]][[quartile_levels[2]]],
        reps = randomization_reps
      ),
      error = function(e) c(obs = NA, null = NA, seNull = NA, pNull = NA)
    )
  } else {
    c(obs = NA, null = NA, seNull = NA, pNull = NA)
  }
  compare_p <- if (is.finite(compare_result["pNull"])) {
    (round(compare_result["pNull"] * randomization_reps) + 1) /
      (randomization_reps + 1)
  } else {
    NA_real_
  }

  # Permute camera labels, rather than event rows, for the nocturnality test.
  extreme_cameras <- unlist(camera_ids, use.names = FALSE)
  nocturnality_permutation_delta <- rep(NA_real_, randomization_reps)
  if (eligible_i) {
    for (rep_i in seq_len(randomization_reps)) {
      permuted_cameras <- sample(extreme_cameras, length(extreme_cameras))
      permuted_q1 <- permuted_cameras[seq_len(15L)]
      permuted_q4 <- permuted_cameras[16:30]
      noct_q1_i <- event_i$nocturnality_hours[event_i$camera_name %in% permuted_q1]
      noct_q4_i <- event_i$nocturnality_hours[event_i$camera_name %in% permuted_q4]
      if (length(noct_q1_i) && length(noct_q4_i)) {
        nocturnality_permutation_delta[rep_i] <- mean(noct_q4_i) - mean(noct_q1_i)
      }
    }
  }
  valid_nocturnality_permutations <-
    nocturnality_permutation_delta[is.finite(nocturnality_permutation_delta)]
  nocturnality_p_two_sided <- if (
    eligible_i &&
      length(valid_nocturnality_permutations) >= ceiling(0.8 * randomization_reps)
  ) {
    (sum(abs(valid_nocturnality_permutations) >= abs(nocturnality_delta)) + 1) /
      (length(valid_nocturnality_permutations) + 1)
  } else {
    NA_real_
  }
  nocturnality_p_decrease <- if (is.finite(nocturnality_p_two_sided)) {
    (sum(valid_nocturnality_permutations <= nocturnality_delta) + 1) /
      (length(valid_nocturnality_permutations) + 1)
  } else {
    NA_real_
  }

  result <- tibble(
    predator = predator_i,
    class = class_i,
    n_q1 = nrow(group_data[[quartile_levels[1]]]),
    n_q4 = nrow(group_data[[quartile_levels[2]]]),
    cameras_with_events_q1 = n_distinct(group_data[[quartile_levels[1]]]$camera_name),
    cameras_with_events_q4 = n_distinct(group_data[[quartile_levels[2]]]$camera_name),
    estimator = estimator_i,
    formal_test_eligible = eligible_i,
    overlap = observed_overlap,
    overlap_conf_low = overlap_ci[1],
    overlap_conf_high = overlap_ci[2],
    camera_bootstrap_overlap_successes = sum(overlap_success),
    compareCkern_overlap = unname(compare_result["obs"]),
    compareCkern_null_overlap = unname(compare_result["null"]),
    compareCkern_p = compare_p,
    mean_nocturnality_q1 = mean_q1,
    mean_nocturnality_q4 = mean_q4,
    nocturnality_delta_q4_minus_q1 = nocturnality_delta,
    nocturnality_delta_conf_low = nocturnality_ci[1],
    nocturnality_delta_conf_high = nocturnality_ci[2],
    camera_bootstrap_nocturnality_successes = sum(nocturnality_success),
    nocturnality_permutation_successes = length(valid_nocturnality_permutations),
    nocturnality_p_two_sided = nocturnality_p_two_sided,
    nocturnality_p_decrease = nocturnality_p_decrease
  )

  bootstrap <- tibble(
    predator = predator_i,
    class = class_i,
    replicate = seq_len(bootstrap_reps),
    overlap = bootstrap_overlap,
    mean_nocturnality_q1 = bootstrap_mean_q1,
    mean_nocturnality_q4 = bootstrap_mean_q4,
    nocturnality_delta_q4_minus_q1 = bootstrap_nocturnality_delta
  )

  permutation <- tibble(
    predator = predator_i,
    class = class_i,
    replicate = seq_len(randomization_reps),
    nocturnality_delta_q4_minus_q1 = nocturnality_permutation_delta
  )

  list(result = result, bootstrap = bootstrap, permutation = permutation)
}

# Seed each worker reproducibly because bootstrap and compareCkern are stochastic.
worker_count <- max(1L, min(parallel_workers, nrow(tasks)))
if (worker_count > 1L) {
  cluster <- parallel::makeCluster(worker_count)
  parallel::clusterEvalQ(cluster, {
    suppressPackageStartupMessages({
      library(activity)
      library(dplyr)
      library(overlap)
      library(tibble)
    })
  })
  parallel::clusterExport(
    cluster,
    c(
      "tasks", "analysis_events", "camera_groups", "activity_fits",
      "quartile_levels", "bootstrap_reps", "randomization_reps",
      "safe_overlap", "resample_camera_values"
    ),
    envir = environment()
  )
  parallel::clusterSetRNGStream(cluster, iseed = 20260819L)
  analysis_rows <- tryCatch(
    parallel::parLapply(cluster, seq_len(nrow(tasks)), analyze_task),
    finally = parallel::stopCluster(cluster)
  )
} else {
  analysis_rows <- lapply(seq_len(nrow(tasks)), analyze_task)
}

results <- bind_rows(lapply(analysis_rows, `[[`, "result")) %>%
  mutate(
    distributional_difference = case_when(
      !formal_test_eligible ~ "Data deficient",
      compareCkern_p < 0.05 ~ "P < 0.05",
      TRUE ~ "Not supported"
    ),
    nocturnality_test = case_when(
      !formal_test_eligible ~ "Data deficient",
      nocturnality_p_two_sided < 0.05 ~ "P < 0.05",
      TRUE ~ "Not supported"
    ),
    nocturnality_direction = case_when(
      nocturnality_delta_q4_minus_q1 < 0 ~ "More diurnal at Q4",
      nocturnality_delta_q4_minus_q1 > 0 ~ "More nocturnal at Q4",
      TRUE ~ "No mean change"
    ),
    nocturnality_interval_excludes_zero =
      nocturnality_delta_conf_low > 0 | nocturnality_delta_conf_high < 0
  ) %>%
  arrange(match(predator, predator_levels), match(class, class_levels))

bootstrap_draws <- bind_rows(lapply(analysis_rows, `[[`, "bootstrap"))
permutation_draws <- bind_rows(lapply(analysis_rows, `[[`, "permutation"))
write_csv(results, file.path(output_dir, "H3_predator_quartile_results.csv"))
write_csv(bootstrap_draws, file.path(output_dir, "H3_predator_quartile_bootstrap_draws.csv"))
write_csv(
  permutation_draws,
  file.path(output_dir, "H3_predator_quartile_nocturnality_permutation_draws.csv")
)

density_rows <- list()
density_i <- 1L
for (task_i in seq_len(nrow(tasks))) {
  for (quartile_i in quartile_levels) {
    fit_i <- activity_fits[[task_i]][[quartile_i]]
    if (is.null(fit_i)) next
    pdf_i <- as.data.frame(fit_i@pdf)
    density_rows[[density_i]] <- tibble(
      predator = tasks$predator[task_i],
      class = tasks$class[task_i],
      quartile = quartile_i,
      solar_hour = to_hours(pdf_i$x),
      density = pdf_i$y
    )
    density_i <- density_i + 1L
  }
}

density_curves <- bind_rows(density_rows) %>%
  mutate(
    predator = factor(predator, levels = predator_levels),
    class = factor(class, levels = class_levels),
    quartile = factor(quartile, levels = quartile_levels)
  )

shared_density <- density_curves %>%
  select(predator, class, quartile, solar_hour, density) %>%
  pivot_wider(names_from = quartile, values_from = density) %>%
  mutate(
    shared_density = pmin(
      .data[[quartile_levels[1]]],
      .data[[quartile_levels[2]]]
    )
  )

write_csv(density_curves, file.path(output_dir, "H3_predator_quartile_density_curves.csv"))

annotations <- results %>%
  mutate(
    predator = factor(predator, levels = predator_levels),
    class = factor(class, levels = class_levels),
    kernel_label = ifelse(
      formal_test_eligible,
      paste0("compareCkern P = ", fmt_p(compareCkern_p)),
      "Kernel test: data deficient (n < 20)"
    ),
    nocturnality_test_label = ifelse(
      formal_test_eligible,
      paste0("Nocturnality P = ", fmt_p(nocturnality_p_two_sided)),
      "Nocturnality test: data deficient"
    ),
    label = paste0(
      estimator, " = ", fmt(overlap),
      " [", fmt(overlap_conf_low), ", ", fmt(overlap_conf_high), "]\n",
      "Delta nocturnality = ",
      ifelse(nocturnality_delta_q4_minus_q1 >= 0, "+", ""),
      fmt(nocturnality_delta_q4_minus_q1), " h",
      " [", fmt(nocturnality_delta_conf_low), ", ",
      fmt(nocturnality_delta_conf_high), "]\n",
      kernel_label, " | ", nocturnality_test_label,
      " | n = ", n_q1, "/", n_q4
    )
  )

quartile_colors <- c(
  "Q1: low activity" = "#315D6D",
  "Q4: high activity" = "#B76B32"
)
quartile_linetypes <- c(
  "Q1: low activity" = "longdash",
  "Q4: high activity" = "solid"
)

quartile_plot <- ggplot() +
  annotate(
    "rect", xmin = 0, xmax = 5.4, ymin = -Inf, ymax = Inf,
    fill = "#DCE3E5", alpha = 0.65
  ) +
  annotate(
    "rect", xmin = 5.4, xmax = 7.0, ymin = -Inf, ymax = Inf,
    fill = "#F5E8D8", alpha = 0.50
  ) +
  annotate(
    "rect", xmin = 17.0, xmax = 18.8, ymin = -Inf, ymax = Inf,
    fill = "#F1DAD5", alpha = 0.50
  ) +
  annotate(
    "rect", xmin = 18.8, xmax = 24, ymin = -Inf, ymax = Inf,
    fill = "#DCE3E5", alpha = 0.65
  ) +
  geom_area(
    data = shared_density,
    aes(x = solar_hour, y = shared_density),
    fill = "#8C9698",
    alpha = 0.35
  ) +
  geom_line(
    data = density_curves,
    aes(x = solar_hour, y = density, color = quartile, linetype = quartile),
    linewidth = 1.0
  ) +
  geom_label(
    data = annotations,
    aes(x = 23.6, y = Inf, label = label),
    hjust = 1,
    vjust = 1.1,
    size = 3.0,
    lineheight = 0.95,
    linewidth = 0,
    fill = scales::alpha("white", 0.78)
  ) +
  facet_grid(
    predator ~ class,
    scales = "free_y",
    labeller = labeller(
      predator = c(Wolf = "Wolf activity", Jackal = "Jackal activity")
    )
  ) +
  scale_color_manual(values = quartile_colors, name = NULL) +
  scale_linetype_manual(values = quartile_linetypes, name = NULL) +
  scale_x_continuous(
    breaks = seq(0, 24, 3),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  labs(
    x = "Apparent solar time",
    y = "Activity density",
    title = "Gazelle diel activity at low and high predator-activity sites",
    subtitle = paste0(
      "Q1 and Q4 contain 15 cameras per predator; gray fill is shared density. ",
      "Intervals use camera-cluster bootstrap."
    ),
    caption = paste0(
      "Delta nocturnality is Q4 - Q1 mean hours from solar noon: negative values indicate greater diurnality at high-activity sites. ",
      "P values are unadjusted; nocturnality P values use two-sided camera-label permutation tests."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(hjust = 0, size = 8.5),
    legend.position = "top",
    legend.justification = "center"
  )

ggsave(
  file.path(output_dir, "publication_gazelle_class_diel_predator_quartiles.png"),
  quartile_plot,
  width = 13,
  height = 8.3,
  dpi = 400,
  bg = "white"
)
ggsave(
  file.path(output_dir, "publication_gazelle_class_diel_predator_quartiles.pdf"),
  quartile_plot,
  width = 13,
  height = 8.3,
  bg = "white"
)

result_table_rows <- vapply(seq_len(nrow(results)), function(i) {
  result_i <- results[i, ]
  kernel_result <- if (result_i$formal_test_eligible) {
    paste0(fmt_p(result_i$compareCkern_p), " (", result_i$distributional_difference, ")")
  } else {
    "Data deficient"
  }
  nocturnality_result <- if (result_i$formal_test_eligible) {
    paste0(fmt_p(result_i$nocturnality_p_two_sided), " (", result_i$nocturnality_test, ")")
  } else {
    "Data deficient"
  }
  paste0(
    "| ", result_i$predator,
    " | ", result_i$class,
    " | ", result_i$n_q1, " / ", result_i$n_q4,
    " | ", result_i$estimator, " = ", fmt(result_i$overlap),
    " [", fmt(result_i$overlap_conf_low), ", ", fmt(result_i$overlap_conf_high), "]",
    " | ", kernel_result,
    " | ", fmt(result_i$mean_nocturnality_q1), " / ",
    fmt(result_i$mean_nocturnality_q4),
    " | ", ifelse(result_i$nocturnality_delta_q4_minus_q1 >= 0, "+", ""),
    fmt(result_i$nocturnality_delta_q4_minus_q1),
    " [", fmt(result_i$nocturnality_delta_conf_low), ", ",
    fmt(result_i$nocturnality_delta_conf_high), "]",
    " | ", nocturnality_result, " |"
  )
}, character(1))

supported_distributions <- results %>%
  filter(distributional_difference == "P < 0.05") %>%
  transmute(label = paste0(tolower(class), " at ", tolower(predator), "-quartile sites")) %>%
  pull(label)
supported_nocturnality <- results %>%
  filter(nocturnality_test == "P < 0.05") %>%
  transmute(label = paste0(
    tolower(class), " at ", tolower(predator), "-quartile sites (Delta = ",
    fmt(nocturnality_delta_q4_minus_q1), " h)"
  )) %>%
  pull(label)

summary_lines <- c(
  "# Gazelle diel activity across predator-activity quartiles",
  "",
  paste0(
    "The analysis used ", nrow(events),
    " validated class-specific 30-minute gazelle events. For each predator, Q1 and Q4 each contained 15 of the 60 cameras."
  ),
  "",
  "Nocturnality was calculated as absolute apparent-solar hours from solar noon (0 = noon; 12 = midnight). Negative Q4-Q1 changes indicate greater diurnality at high predator-activity sites.",
  "",
  "## Quartile thresholds",
  "",
  capture.output(print(camera_summary, n = Inf)),
  "",
  "## Results",
  "",
  "| Predator | Class | Events Q1 / Q4 | Q1-Q4 overlap [95% CI] | compareCkern P | Mean nocturnality Q1 / Q4 (h) | Delta nocturnality [95% CI] (h) | Nocturnality P |",
  "|---|---|---:|---:|---:|---:|---:|---:|",
  result_table_rows,
  "",
  paste0(
    "Using unadjusted compareCkern P values, whole-distribution differences had P < 0.05 for ",
    ifelse(
      length(supported_distributions),
      paste(supported_distributions, collapse = ", "),
      "no class-predator comparison"
    ), "."
  ),
  paste0(
    "Using unadjusted two-sided camera-label permutation tests, nocturnality changes had P < 0.05 for ",
    ifelse(
      length(supported_nocturnality),
      paste(supported_nocturnality, collapse = ", "),
      "no comparison"
    ), "."
  ),
  "The fawn contrast across jackal quartiles was data deficient because Q4 contained 14 events from four cameras.",
  "",
  "## Interpretation constraints",
  "",
  "Quartiles are descriptive site categories. The compareCkern test resamples event times rather than camera clusters, so its unadjusted P values should be interpreted alongside camera-bootstrap intervals and sample sizes. Nocturnality P values use camera-label permutations. Predator indices are relative camera-level activity measures, not direct predation-risk estimates."
)
writeLines(summary_lines, file.path(output_dir, "H3_predator_quartile_results.md"))

session_lines <- c(
  R.version.string,
  paste("activity", as.character(packageVersion("activity"))),
  paste("overlap", as.character(packageVersion("overlap"))),
  paste("ggplot2", as.character(packageVersion("ggplot2")))
)
writeLines(session_lines, file.path(output_dir, "H3_predator_quartile_session_info.txt"))

cat("\nPredator quartile camera summary\n")
print(camera_summary, n = Inf)
cat("\nClass-specific Q1-Q4 diel results\n")
print(results, n = Inf)
cat("\nOutputs saved to: ", normalizePath(output_dir), "\n", sep = "")
