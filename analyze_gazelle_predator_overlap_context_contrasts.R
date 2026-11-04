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
  source_prefix = "publication_gazelle_predator_diel_by_context_median",
  output_prefix = "publication_gazelle_predator_overlap_context_contrasts",
  gazelle_event_file = "Gazella gazella_independent_events_30min.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  camera_file = "publication_predator_context_class_diel_2x2_median_camera_categories.csv",
  local_timezone = "Asia/Jerusalem",
  bootstrap_reps = "999",
  minimum_events = "5",
  small_sample_threshold = "75",
  seed = "20260811"
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
minimum_events <- as.integer(cfg$minimum_events)
small_sample_threshold <- as.integer(cfg$small_sample_threshold)
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

predator_event_file <- paste0(
  cfg$source_prefix,
  "_predator_independent_events.csv"
)
required_files <- c(
  predator_event_file, cfg$gazelle_event_file, cfg$obs_cov_file,
  cfg$camera_file
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

camera_context <- read_csv(cfg$camera_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(camera_name),
    predator_context = factor(predator_context, levels = context_levels)
  )

obs <- read.csv(
  cfg$obs_cov_file,
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

gazelle_raw <- read_csv(
  cfg$gazelle_event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  mutate(
    event_row_id = row_number(),
    camera_name = as.integer(camera_name),
    class = recode(
      gazelle_class,
      juvenile = "Fawn",
      female = "Female",
      male = "Male"
    ),
    timestamp = parse_local_camera_time(event_time),
    event_date = as.Date(format(timestamp, tz = cfg$local_timezone)),
    solar_hour = apparent_solar_hour(timestamp, camelot_gps_longitude),
    time_rad = to_radians(solar_hour)
  )

active_gazelle_ids <- gazelle_raw %>%
  select(event_row_id, camera_name, event_date) %>%
  inner_join(obs, by = "camera_name", relationship = "many-to-many") %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  group_by(event_row_id) %>%
  slice_max(start_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(event_row_id)

gazelle_events <- gazelle_raw %>%
  inner_join(active_gazelle_ids, by = "event_row_id") %>%
  left_join(camera_context, by = "camera_name") %>%
  filter(
    class %in% class_levels,
    predator_context %in% context_levels,
    is.finite(time_rad)
  ) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(predator_context, levels = context_levels)
  )

predator_events <- read_csv(
  predator_event_file,
  col_types = cols(timestamp = col_character()),
  show_col_types = FALSE
) %>%
  transmute(
    camera_name = as.integer(camera_name),
    predator = factor(taxon, levels = predator_levels),
    predator_context = factor(predator_context, levels = context_levels),
    solar_hour,
    time_rad = to_radians(solar_hour)
  ) %>%
  filter(
    predator %in% predator_levels,
    predator_context %in% context_levels,
    is.finite(time_rad)
  )

event_counts <- crossing(
  class = factor(class_levels, levels = class_levels),
  predator = factor(predator_levels, levels = predator_levels),
  predator_context = factor(context_levels, levels = context_levels)
) %>%
  rowwise() %>%
  mutate(
    gazelle_events = sum(
      gazelle_events$class == class &
        gazelle_events$predator_context == predator_context
    ),
    predator_events = sum(
      predator_events$predator == predator &
        predator_events$predator_context == predator_context
    )
  ) %>%
  ungroup()

context_pairs <- combn(context_levels, 2, simplify = FALSE)
contrast_tasks <- list()
task_index <- 1L
for (class_i in class_levels) {
  for (predator_i in predator_levels) {
    for (pair_i in context_pairs) {
      context_a <- pair_i[1]
      context_b <- pair_i[2]
      counts_a <- event_counts %>%
        filter(
          class == class_i,
          predator == predator_i,
          predator_context == context_a
        )
      counts_b <- event_counts %>%
        filter(
          class == class_i,
          predator == predator_i,
          predator_context == context_b
        )
      minimum_pair_sample <- min(
        counts_a$gazelle_events,
        counts_a$predator_events,
        counts_b$gazelle_events,
        counts_b$predator_events
      )
      estimator <- if (minimum_pair_sample < small_sample_threshold) {
        "Dhat1"
      } else {
        "Dhat4"
      }
      contrast_tasks[[task_index]] <- tibble(
        class = class_i,
        predator = predator_i,
        context_a,
        context_b,
        estimator,
        eligible = minimum_pair_sample >= minimum_events
      )
      task_index <- task_index + 1L
    }
  }
}
contrast_tasks <- bind_rows(contrast_tasks)

required_overlap_keys <- bind_rows(
  contrast_tasks %>%
    filter(eligible) %>%
    transmute(class, predator, predator_context = context_a, estimator),
  contrast_tasks %>%
    filter(eligible) %>%
    transmute(class, predator, predator_context = context_b, estimator)
) %>%
  distinct() %>%
  arrange(class, predator, predator_context, estimator) %>%
  mutate(
    key = paste(class, predator, predator_context, estimator, sep = "||")
  )

camera_ids_by_context <- lapply(context_levels, function(context_i) {
  camera_context$camera_name[camera_context$predator_context == context_i]
})
names(camera_ids_by_context) <- context_levels

bootstrap_overlap_key <- function(class_i, predator_i, context_i, estimator) {
  cameras <- camera_ids_by_context[[context_i]]
  gazelle_i <- gazelle_events %>%
    filter(class == class_i, predator_context == context_i)
  predator_i_data <- predator_events %>%
    filter(predator == predator_i, predator_context == context_i)

  observed_gazelle <- gazelle_i$time_rad
  observed_predator <- predator_i_data$time_rad
  observed_overlap <- tryCatch(
    unname(overlapEst(
      observed_gazelle,
      observed_predator,
      type = estimator
    )),
    error = function(e) NA_real_
  )

  bootstrap_values <- replicate(bootstrap_reps, {
    sampled_cameras <- sample(cameras, length(cameras), replace = TRUE)
    gazelle_times <- unlist(lapply(sampled_cameras, function(camera_i) {
      gazelle_i$time_rad[gazelle_i$camera_name == camera_i]
    }), use.names = FALSE)
    predator_times <- unlist(lapply(sampled_cameras, function(camera_i) {
      predator_i_data$time_rad[predator_i_data$camera_name == camera_i]
    }), use.names = FALSE)
    if (min(length(gazelle_times), length(predator_times)) < 2) {
      return(NA_real_)
    }
    tryCatch(
      unname(overlapEst(gazelle_times, predator_times, type = estimator)),
      error = function(e) NA_real_
    )
  })
  list(observed = observed_overlap, bootstrap = bootstrap_values)
}

message(
  "Bootstrapping ", nrow(required_overlap_keys),
  " class-predator-context-estimator combinations with ",
  bootstrap_reps, " camera replicates each."
)
overlap_store <- list()
bootstrap_draw_rows <- list()
for (i in seq_len(nrow(required_overlap_keys))) {
  key_i <- required_overlap_keys[i, ]
  result_i <- bootstrap_overlap_key(
    key_i$class,
    key_i$predator,
    key_i$predator_context,
    key_i$estimator
  )
  overlap_store[[key_i$key]] <- result_i
  bootstrap_draw_rows[[i]] <- tibble(
    class = key_i$class,
    predator = key_i$predator,
    predator_context = key_i$predator_context,
    estimator = key_i$estimator,
    replicate = seq_len(bootstrap_reps),
    overlap = result_i$bootstrap
  )
}
bootstrap_draws <- bind_rows(bootstrap_draw_rows)
write_csv(
  bootstrap_draws,
  paste0(cfg$output_prefix, "_bootstrap_draws.csv")
)

key_name <- function(class_i, predator_i, context_i, estimator_i) {
  paste(class_i, predator_i, context_i, estimator_i, sep = "||")
}

contrast_rows <- list()
for (i in seq_len(nrow(contrast_tasks))) {
  task <- contrast_tasks[i, ]
  if (!task$eligible) {
    contrast_rows[[i]] <- task %>%
      mutate(
        overlap_a = NA_real_,
        overlap_b = NA_real_,
        delta_overlap = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        bootstrap_se = NA_real_,
        bootstrap_successes = 0L,
        p_bootstrap = NA_real_,
        reason = paste0("Fewer than ", minimum_events, " events in one context")
      )
    next
  }

  result_a <- overlap_store[[key_name(
    task$class, task$predator, task$context_a, task$estimator
  )]]
  result_b <- overlap_store[[key_name(
    task$class, task$predator, task$context_b, task$estimator
  )]]
  delta_bootstrap <- result_a$bootstrap - result_b$bootstrap
  delta_bootstrap <- delta_bootstrap[is.finite(delta_bootstrap)]
  sufficient_bootstrap <- length(delta_bootstrap) >= ceiling(0.8 * bootstrap_reps)
  observed_delta <- result_a$observed - result_b$observed

  if (sufficient_bootstrap) {
    interval <- unname(quantile(delta_bootstrap, c(0.025, 0.975)))
    lower_tail <- (sum(delta_bootstrap <= 0) + 1) /
      (length(delta_bootstrap) + 1)
    upper_tail <- (sum(delta_bootstrap >= 0) + 1) /
      (length(delta_bootstrap) + 1)
    p_bootstrap <- min(1, 2 * min(lower_tail, upper_tail))
  } else {
    interval <- c(NA_real_, NA_real_)
    p_bootstrap <- NA_real_
  }

  contrast_rows[[i]] <- task %>%
    mutate(
      overlap_a = result_a$observed,
      overlap_b = result_b$observed,
      delta_overlap = observed_delta,
      conf_low = interval[1],
      conf_high = interval[2],
      bootstrap_se = ifelse(
        sufficient_bootstrap,
        sd(delta_bootstrap),
        NA_real_
      ),
      bootstrap_successes = length(delta_bootstrap),
      p_bootstrap,
      reason = ifelse(
        sufficient_bootstrap,
        NA_character_,
        "Too few successful camera-bootstrap contrasts"
      )
    )
}

contrasts <- bind_rows(contrast_rows) %>%
  group_by(class, predator) %>%
  mutate(
    p_fdr_within_class_predator = p.adjust(p_bootstrap, method = "BH")
  ) %>%
  ungroup()

contrasts$p_fdr_global <- NA_real_
eligible_index <- which(is.finite(contrasts$p_bootstrap))
contrasts$p_fdr_global[eligible_index] <- p.adjust(
  contrasts$p_bootstrap[eligible_index],
  method = "BH"
)
contrasts <- contrasts %>%
  mutate(
    supported = is.finite(p_fdr_within_class_predator) &
      p_fdr_within_class_predator < 0.05 &
      ((conf_low > 0 & conf_high > 0) | (conf_low < 0 & conf_high < 0)),
    context_contrast = paste(context_a, "minus", context_b)
  ) %>%
  arrange(class, predator, context_a, context_b)

write_csv(contrasts, paste0(cfg$output_prefix, "_results.csv"))

pair_levels <- unique(unlist(lapply(context_pairs, function(pair_i) {
  paste(pair_i[1], "minus", pair_i[2])
})))
plot_data <- contrasts %>%
  filter(is.finite(delta_overlap)) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator = factor(predator, levels = predator_levels),
    context_contrast = factor(context_contrast, levels = rev(pair_levels)),
    inference = ifelse(
      supported,
      "Supported after within-family FDR",
      "Not supported"
    )
  )

contrast_colors <- c(
  "Supported after within-family FDR" = "#A54231",
  "Not supported" = "#285F73"
)

contrast_plot <- ggplot(
  plot_data,
  aes(x = delta_overlap, y = context_contrast, color = inference)
) +
  geom_vline(xintercept = 0, color = "#6F797C", linewidth = 0.55) +
  geom_segment(
    aes(x = conf_low, xend = conf_high, yend = context_contrast),
    linewidth = 0.85
  ) +
  geom_point(size = 2.7) +
  facet_grid(class ~ predator, scales = "free_y", space = "free_y") +
  scale_color_manual(values = contrast_colors, name = NULL) +
  labs(
    x = "Difference in temporal overlap (first context minus second)",
    y = NULL,
    title = "Context contrasts in gazelle-predator temporal overlap",
    subtitle = paste0(
      "Intervals use ", bootstrap_reps,
      " camera-cluster bootstrap replicates; colors use within-class/predator FDR\n",
      "Positive values indicate greater overlap in the first named context"
    ),
    caption = paste0(
      "A common estimator was used within each contrast: Dhat1 when any contributing sample was below ",
      small_sample_threshold, " events, otherwise Dhat4. ",
      "Fawn contrasts involving high-jackal/low-wolf cameras were excluded (n=3)."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "#E8E1D4", color = NA),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(color = "#314043"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4C5A5E"),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "top"
  )

ggsave(
  paste0(cfg$output_prefix, ".png"),
  contrast_plot,
  width = 13.5,
  height = 10.5,
  dpi = 320
)
ggsave(
  paste0(cfg$output_prefix, ".pdf"),
  contrast_plot,
  width = 13.5,
  height = 10.5
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt(x, 3)))
}

report <- c(
  "# Camera-bootstrap context contrasts in gazelle-predator overlap",
  "",
  "## Methods",
  "",
  paste0(
    "For each gazelle class and predator, temporal-overlap coefficients were contrasted among all six pairs of median-split wolf-jackal contexts. Cameras were resampled with replacement within each context, and all gazelle and predator events from a sampled camera were retained together. This preserves within-camera dependence between the two taxa. Each contrast used ",
    bootstrap_reps, " replicates."
  ),
  "",
  paste0(
    "A common overlap estimator was used on both sides of each contrast: `Dhat1` when any of the four contributing samples contained fewer than ",
    small_sample_threshold,
    " events and `Dhat4` otherwise. Contrasts were excluded when either taxon had fewer than ",
    minimum_events, " events in either context."
  ),
  "",
  "Approximate two-sided p-values were calculated from the proportion of camera-bootstrap differences on either side of zero, with a +1 correction. P-values were adjusted across the context comparisons within each class-predator family and, secondarily, across all eligible comparisons. A contrast was labeled supported only when its percentile interval excluded zero and its within-family FDR p-value was below 0.05.",
  "",
  "## Results",
  "",
  "| Class | Predator | Context contrast | Estimator | Overlap | Difference [95% camera-bootstrap CI] | Within-family FDR p | Global FDR p | Status |",
  "|---|---|---|---|---:|---:|---:|---:|---|"
)

for (i in seq_len(nrow(contrasts))) {
  overlap_text <- if (is.finite(contrasts$overlap_a[i])) {
    paste0(fmt(contrasts$overlap_a[i]), " / ", fmt(contrasts$overlap_b[i]))
  } else {
    "NA"
  }
  difference_text <- if (is.finite(contrasts$delta_overlap[i])) {
    paste0(
      fmt(contrasts$delta_overlap[i]), " [",
      fmt(contrasts$conf_low[i]), ", ", fmt(contrasts$conf_high[i]), "]"
    )
  } else {
    "NA"
  }
  status <- if (!contrasts$eligible[i]) {
    contrasts$reason[i]
  } else if (contrasts$supported[i]) {
    "Supported"
  } else {
    "Not supported"
  }
  report <- c(report, sprintf(
    "| %s | %s | %s | `%s` | `%s` | `%s` | `%s` | `%s` | %s |",
    contrasts$class[i], contrasts$predator[i],
    contrasts$context_contrast[i], contrasts$estimator[i],
    overlap_text, difference_text,
    fmt_p(contrasts$p_fdr_within_class_predator[i]),
    fmt_p(contrasts$p_fdr_global[i]), status
  ))
}

supported_results <- contrasts %>% filter(supported)
report <- c(
  report,
  "",
  "## Multiplicity-aware summary",
  "",
  paste0(
    "Supported within-family context contrasts: `",
    nrow(supported_results), "` of `",
    sum(is.finite(contrasts$p_bootstrap)), "` eligible comparisons."
  )
)
if (nrow(supported_results)) {
  for (i in seq_len(nrow(supported_results))) {
    report <- c(report, paste0(
      "- **", supported_results$class[i], " x ",
      supported_results$predator[i], ":** ",
      supported_results$context_contrast[i], ", difference `",
      fmt(supported_results$delta_overlap[i]), "` (95% CI `",
      fmt(supported_results$conf_low[i]), "-",
      fmt(supported_results$conf_high[i]), "`; within-family FDR `p = ",
      fmt_p(supported_results$p_fdr_within_class_predator[i]), "`)."
    ))
  }
} else {
  report <- c(report, "- No context difference met both support criteria.")
}

report <- c(
  report,
  "",
  "## Interpretation constraints",
  "",
  "- These are camera-level contrasts in temporal overlap, not tests of individual behavioral avoidance or attraction.",
  "- Contexts were defined by camera-level predator abundance and were not experimentally assigned. Habitat and other site differences may therefore confound overlap contrasts.",
  "- Bootstrap p-values are approximate. The percentile confidence intervals and effect magnitudes should be emphasized over threshold-based significance.",
  "- The fawn high-jackal/low-wolf sample contained only three events and was excluded from formal contrasts.",
  "",
  "## Outputs",
  "",
  paste0("- Contrast table: `", cfg$output_prefix, "_results.csv`."),
  paste0("- Figure: `", cfg$output_prefix, ".png` and `.pdf`."),
  paste0("- Reproducible bootstrap draws: `", cfg$output_prefix, "_bootstrap_draws.csv`.")
)
writeLines(report, paste0(cfg$output_prefix, "_results.md"), useBytes = TRUE)

cat("Saved gazelle-predator overlap context contrasts:", cfg$output_prefix, "\n")
