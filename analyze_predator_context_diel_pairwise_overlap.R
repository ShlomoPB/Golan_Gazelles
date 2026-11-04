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
  source_prefix = "publication_predator_context_class_diel_2x2_median",
  output_prefix = "publication_predator_context_diel_pairwise_overlap_median",
  event_file = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  local_timezone = "Asia/Jerusalem",
  small_sample_threshold = "75",
  formal_minimum_events = "20",
  compare_reps = "999",
  camera_bootstrap_reps = "500",
  parallel_workers = "3",
  seed = "20260811"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

source_prefix <- arg_map$source_prefix
output_prefix <- arg_map$output_prefix
event_file <- arg_map$event_file
obs_cov_file <- arg_map$obs_cov_file
local_timezone <- arg_map$local_timezone
small_sample_threshold <- as.integer(arg_map$small_sample_threshold)
formal_minimum_events <- as.integer(arg_map$formal_minimum_events)
compare_reps <- as.integer(arg_map$compare_reps)
camera_bootstrap_reps <- as.integer(arg_map$camera_bootstrap_reps)
parallel_workers <- as.integer(arg_map$parallel_workers)
seed <- as.integer(arg_map$seed)

set.seed(seed)

class_levels <- c("Fawn", "Female", "Male")
context_levels <- c(
  "Low both",
  "High jackal / Low wolf",
  "Low jackal / High wolf",
  "High both"
)
context_short <- c(
  "Low both" = "Low both",
  "High jackal / Low wolf" = "High jackal / low wolf",
  "Low jackal / High wolf" = "Low jackal / high wolf",
  "High both" = "High both"
)
class_labels <- c(
  juvenile = "Fawn",
  female = "Female",
  male = "Male"
)

camera_file <- paste0(source_prefix, "_camera_categories.csv")
required_files <- c(camera_file, event_file, obs_cov_file)
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

camera_context <- read_csv(camera_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(camera_name),
    predator_context = factor(as.character(predator_context), levels = context_levels)
  )

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

active_matches <- raw_events %>%
  select(event_row_id, camera_name, event_date) %>%
  inner_join(obs, by = "camera_name", relationship = "many-to-many") %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  group_by(event_row_id) %>%
  slice_max(start_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(event_row_id)

events <- raw_events %>%
  inner_join(active_matches, by = "event_row_id") %>%
  left_join(camera_context, by = "camera_name") %>%
  filter(class %in% class_levels, predator_context %in% context_levels) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(predator_context, levels = context_levels)
  )

camera_ids_by_context <- lapply(context_levels, function(context_i) {
  camera_context$camera_name[camera_context$predator_context == context_i]
})
names(camera_ids_by_context) <- context_levels

event_times <- list()
activity_fits <- list()
for (class_i in class_levels) {
  event_times[[class_i]] <- list()
  activity_fits[[class_i]] <- list()
  for (context_i in context_levels) {
    times_i <- events$time_rad[
      events$class == class_i & events$predator_context == context_i
    ]
    event_times[[class_i]][[context_i]] <- times_i
    activity_fits[[class_i]][[context_i]] <- tryCatch(
      fitact(times_i, sample = "none", show = FALSE),
      error = function(e) NULL
    )
  }
}

camera_cluster_interval <- function(class_i, context_a, context_b, estimator) {
  data_a <- filter(events, class == class_i, predator_context == context_a)
  data_b <- filter(events, class == class_i, predator_context == context_b)
  cameras_a <- camera_ids_by_context[[context_a]]
  cameras_b <- camera_ids_by_context[[context_b]]

  bootstrap_values <- replicate(camera_bootstrap_reps, {
    sampled_a <- sample(cameras_a, length(cameras_a), replace = TRUE)
    sampled_b <- sample(cameras_b, length(cameras_b), replace = TRUE)
    times_a <- unlist(lapply(sampled_a, function(camera_i) {
      data_a$time_rad[data_a$camera_name == camera_i]
    }), use.names = FALSE)
    times_b <- unlist(lapply(sampled_b, function(camera_i) {
      data_b$time_rad[data_b$camera_name == camera_i]
    }), use.names = FALSE)
    if (length(times_a) < 2 || length(times_b) < 2) return(NA_real_)
    tryCatch(
      unname(overlapEst(times_a, times_b, type = estimator)),
      error = function(e) NA_real_
    )
  })
  bootstrap_values <- bootstrap_values[is.finite(bootstrap_values)]
  if (length(bootstrap_values) < ceiling(0.8 * camera_bootstrap_reps)) {
    return(list(
      conf_low = NA_real_,
      conf_high = NA_real_,
      successes = length(bootstrap_values)
    ))
  }
  list(
    conf_low = unname(quantile(bootstrap_values, 0.025)),
    conf_high = unname(quantile(bootstrap_values, 0.975)),
    successes = length(bootstrap_values)
  )
}

context_pairs <- combn(context_levels, 2, simplify = FALSE)
comparison_rows <- list()
compare_tasks <- list()
comparison_index <- 1L

for (class_i in class_levels) {
  for (pair_i in context_pairs) {
    context_a <- pair_i[1]
    context_b <- pair_i[2]
    times_a <- event_times[[class_i]][[context_a]]
    times_b <- event_times[[class_i]][[context_b]]
    n_a <- length(times_a)
    n_b <- length(times_b)
    estimator <- if (min(n_a, n_b) < small_sample_threshold) "Dhat1" else "Dhat4"
    dhat <- tryCatch(
      unname(overlapEst(times_a, times_b, type = estimator)),
      error = function(e) NA_real_
    )
    interval <- camera_cluster_interval(class_i, context_a, context_b, estimator)

    fit_a <- activity_fits[[class_i]][[context_a]]
    fit_b <- activity_fits[[class_i]][[context_b]]
    eligible <- min(n_a, n_b) >= formal_minimum_events
    compare_tasks[[comparison_index]] <- list(
      fit_a = if (eligible) fit_a else NULL,
      fit_b = if (eligible) fit_b else NULL,
      reps = compare_reps
    )

    comparison_rows[[comparison_index]] <- tibble(
      class = class_i,
      context_a,
      context_b,
      pair = paste(
        unname(context_short[context_a]),
        "vs",
        unname(context_short[context_b])
      ),
      n_a,
      n_b,
      cameras_with_events_a = n_distinct(
        events$camera_name[events$class == class_i & events$predator_context == context_a]
      ),
      cameras_with_events_b = n_distinct(
        events$camera_name[events$class == class_i & events$predator_context == context_b]
      ),
      estimator,
      overlap = dhat,
      conf_low = interval$conf_low,
      conf_high = interval$conf_high,
      camera_bootstrap_successes = interval$successes,
      formal_test_eligible = eligible
    )
    comparison_index <- comparison_index + 1L
  }
}

run_compare_task <- function(task) {
  if (is.null(task$fit_a) || is.null(task$fit_b)) {
    return(c(obs = NA_real_, null = NA_real_, seNull = NA_real_, pNull = NA_real_))
  }
  tryCatch(
    activity::compareCkern(task$fit_a, task$fit_b, reps = task$reps),
    error = function(e) c(obs = NA_real_, null = NA_real_, seNull = NA_real_, pNull = NA_real_)
  )
}

worker_count <- max(1L, min(parallel_workers, length(compare_tasks)))
if (worker_count > 1L) {
  cluster <- parallel::makeCluster(worker_count)
  parallel::clusterSetRNGStream(cluster, iseed = seed + 100000L)
  compare_outputs <- tryCatch(
    parallel::parLapply(cluster, compare_tasks, run_compare_task),
    finally = parallel::stopCluster(cluster)
  )
} else {
  compare_outputs <- lapply(compare_tasks, run_compare_task)
}

compare_matrix <- do.call(rbind, lapply(compare_outputs, function(result) {
  result[c("obs", "null", "seNull", "pNull")]
}))
raw_p <- compare_matrix[, "pNull"]
monte_carlo_p <- ifelse(
  is.finite(raw_p),
  (round(raw_p * compare_reps) + 1) / (compare_reps + 1),
  NA_real_
)

comparison_results <- bind_rows(comparison_rows) %>%
  mutate(
    compareCkern_overlap = compare_matrix[, "obs"],
    compareCkern_null_overlap = compare_matrix[, "null"],
    compareCkern_p_raw = raw_p,
    compareCkern_p = monte_carlo_p
  ) %>%
  group_by(class) %>%
  mutate(
    overlap_rank_all = rank(overlap, ties.method = "min", na.last = "keep")
  ) %>%
  ungroup()

comparison_results <- comparison_results %>%
  group_by(class) %>%
  mutate(
    overlap_rank_eligible = ifelse(
      formal_test_eligible,
      rank(ifelse(formal_test_eligible, overlap, NA_real_), ties.method = "min", na.last = "keep"),
      NA_real_
    ),
    significance = case_when(
      !formal_test_eligible ~ "Descriptive only",
      compareCkern_p < 0.001 ~ "***",
      compareCkern_p < 0.01 ~ "**",
      compareCkern_p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup() %>%
  arrange(class, overlap)

write_csv(comparison_results, paste0(output_prefix, "_results.csv"))

pair_levels <- unique(unlist(lapply(context_pairs, function(pair_i) {
  paste(
    unname(context_short[pair_i[1]]),
    "vs",
    unname(context_short[pair_i[2]])
  )
})))

plot_data <- comparison_results %>%
  mutate(
    class = factor(class, levels = class_levels),
    pair = factor(pair, levels = rev(pair_levels)),
    inference = case_when(
      !formal_test_eligible ~ "Descriptive only (<20 events)",
      compareCkern_p < 0.05 ~ "Raw P < 0.05",
      TRUE ~ "Not supported"
    ),
    estimate_label = paste0(
      sprintf("%.2f", overlap),
      " ",
      significance
    )
  )

inference_colors <- c(
  "Raw P < 0.05" = "#A52A2A",
  "Not supported" = "#245B78",
  "Descriptive only (<20 events)" = "#8A8F91"
)

overlap_plot <- ggplot(
  plot_data,
  aes(x = overlap, y = pair, color = inference)
) +
  geom_segment(
    aes(x = conf_low, xend = conf_high, yend = pair),
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = estimate_label),
    hjust = -0.12,
    size = 3.0,
    color = "#263238"
  ) +
  facet_wrap(~class, ncol = 1) +
  scale_color_manual(values = inference_colors, name = NULL) +
  scale_x_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Temporal overlap (Dhat)",
    y = NULL,
    title = "Within-class temporal overlap among wolf-jackal contexts",
    subtitle = paste0(
      "Intervals use camera-cluster bootstrap; tests use unadjusted compareCkern P values\n",
      "Dhat1 used when the smaller sample was <", small_sample_threshold, " events; otherwise Dhat4"
    ),
    caption = paste0(
      "* P<0.05; ** P<0.01; *** P<0.001. P values are unadjusted. Tests with fewer than ",
      formal_minimum_events,
      " events in either context are descriptive only. Lower overlap indicates a larger temporal shift."
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(color = "#304044"),
    axis.title.x = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#E8E1D4", color = NA),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4C5A5E", lineheight = 1.05),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "top"
  )

ggsave(
  paste0(output_prefix, ".png"),
  overlap_plot,
  width = 10.6,
  height = 10.2,
  dpi = 320
)
ggsave(
  paste0(output_prefix, ".pdf"),
  overlap_plot,
  width = 10.6,
  height = 10.2
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt(x, 3)))
}

md <- c(
  "# Pairwise predator-context temporal overlap within gazelle classes",
  "",
  "## Methods",
  "",
  paste0(
    "All six context pairs were evaluated within each class using the all-camera median classification. `Dhat1` was used when the smaller sample contained fewer than ",
    small_sample_threshold,
    " events; otherwise `Dhat4` was used. Intervals were obtained from ",
    camera_bootstrap_reps,
    " camera-cluster bootstrap replicates. Distributional differences were tested with `activity::compareCkern` using ",
    compare_reps,
    " null replicates and a +1 Monte Carlo correction. We report unadjusted P-values."
  ),
  "",
  paste0(
    "Because `compareCkern` resamples events rather than cameras, its p-values do not fully account for camera-level clustering. Comparisons with fewer than ",
    formal_minimum_events,
    " events in either context were calculated descriptively but excluded from formal inference."
  ),
  "",
  "## Results",
  "",
  "| Class | Context pair | Events | Estimator | Dhat [camera-bootstrap 95% CI] | compareCkern P |",
  "|---|---|---:|---|---:|---:|"
)
for (i in seq_len(nrow(comparison_results))) {
  interval_text <- if (is.finite(comparison_results$conf_low[i])) {
    paste0(
      fmt(comparison_results$overlap[i]), " [",
      fmt(comparison_results$conf_low[i]), ", ",
      fmt(comparison_results$conf_high[i]), "]"
    )
  } else {
    paste0(fmt(comparison_results$overlap[i]), " [NA]")
  }
  md <- c(md, sprintf(
    "| %s | %s | `%s / %s` | `%s` | `%s` | `%s` |",
    comparison_results$class[i],
    comparison_results$pair[i],
    comparison_results$n_a[i],
    comparison_results$n_b[i],
    comparison_results$estimator[i],
    interval_text,
    fmt_p(comparison_results$compareCkern_p[i])
  ))
}

md <- c(md, "", "## Lowest eligible overlap by class", "")
for (class_i in class_levels) {
  lowest <- comparison_results %>%
    filter(class == class_i, formal_test_eligible, is.finite(overlap)) %>%
    slice_min(overlap, n = 1, with_ties = FALSE)
  if (nrow(lowest)) {
    md <- c(md, paste0(
      "- **", class_i, ":** ", lowest$pair,
      ", `Dhat = ", fmt(lowest$overlap),
      "` (95% camera-bootstrap CI `", fmt(lowest$conf_low), "-", fmt(lowest$conf_high),
      "`); unadjusted `compareCkern P = ",
      fmt_p(lowest$compareCkern_p), "`."
    ))
  }
}

writeLines(md, paste0(output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved pairwise diel overlap outputs with prefix:", output_prefix, "\n")
