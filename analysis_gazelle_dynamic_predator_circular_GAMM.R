#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(mvtnorm)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
cfg <- list(
  full_export = file.path("..", "full-export_2026-04-15_2006.csv"),
  gazelle_events = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  obs_covs = "ObsCovs_050524.csv",
  site_covs = "SiteCovs_060324.csv",
  output_prefix = "dynamic_predator_GAMM",
  timezone = "Asia/Jerusalem",
  event_gap_minutes = "30",
  simulations = "1200",
  prediction_step = "0.1",
  run_site_shape = "true",
  run_3h_sensitivity = "true",
  seed = "20260820"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2L) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(cfg)) stop("Unknown argument: ", key)
    cfg[[key]] <- parts[2]
  }
}

event_gap_minutes <- as.numeric(cfg$event_gap_minutes)
simulations <- as.integer(cfg$simulations)
prediction_step <- as.numeric(cfg$prediction_step)
run_site_shape <- tolower(cfg$run_site_shape) %in% c("true", "1", "yes")
run_3h_sensitivity <- tolower(cfg$run_3h_sensitivity) %in% c("true", "1", "yes")
seed <- as.integer(cfg$seed)
set.seed(seed)

out <- function(suffix) paste0(cfg$output_prefix, "_", suffix)
required_files <- c(cfg$full_export, cfg$gazelle_events, cfg$obs_covs, cfg$site_covs)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

class_keys <- c("juvenile", "female", "male")
class_labels <- c(
  juvenile = "Juvenile",
  female = "Adult female",
  male = "Adult male"
)
response_columns <- c(
  juvenile = "juvenile_events",
  female = "female_events",
  male = "male_events"
)
basis_k <- c(juvenile = 12L, female = 16L, male = 16L)
expected_gazelle_events <- c(juvenile = 180L, female = 1477L, male = 2268L)
context_levels <- c(
  "Low dynamic wolf / Low dynamic jackal",
  "Low dynamic wolf / High dynamic jackal",
  "High dynamic wolf / Low dynamic jackal",
  "High dynamic wolf / High dynamic jackal"
)
context_labels <- c(
  "Low wolf / Low jackal",
  "Low wolf / High jackal",
  "High wolf / Low jackal",
  "High wolf / High jackal"
)
context_colors <- c(
  "Low dynamic wolf / Low dynamic jackal" = "#188977",
  "Low dynamic wolf / High dynamic jackal" = "#B12B32",
  "High dynamic wolf / Low dynamic jackal" = "#173B63",
  "High dynamic wolf / High dynamic jackal" = "#8A5A9E"
)
context_linetypes <- c(
  "Low dynamic wolf / Low dynamic jackal" = "solid",
  "Low dynamic wolf / High dynamic jackal" = "22",
  "High dynamic wolf / Low dynamic jackal" = "42",
  "High dynamic wolf / High dynamic jackal" = "dotdash"
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f", big.mark = ","))
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", formatC(x, 3, format = "f")))
}

parse_export_time <- function(x) {
  parsed <- as.POSIXct(x, format = "%d/%m/%Y %H:%M:%S", tz = cfg$timezone)
  missing <- is.na(parsed)
  parsed[missing] <- as.POSIXct(
    x[missing], format = "%d/%m/%Y %H:%M", tz = cfg$timezone
  )
  parsed
}

parse_gazelle_time <- function(x) {
  as.POSIXct(
    sub("Z$", "", x),
    format = "%Y-%m-%dT%H:%M:%S",
    tz = cfg$timezone
  )
}

equation_of_time <- function(day_of_year) {
  b <- 2 * pi * (day_of_year - 81) / 364
  9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
}

apparent_solar_hour <- function(timestamp, longitude) {
  day_of_year <- as.integer(format(timestamp, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(timestamp, "%H", tz = "UTC")) +
    as.numeric(format(timestamp, "%M", tz = "UTC")) / 60 +
    as.numeric(format(timestamp, "%S", tz = "UTC")) / 3600
  (utc_hour + longitude / 15 + equation_of_time(day_of_year) / 60) %% 24
}

effort <- read.csv(
  cfg$obs_covs,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    camera_id = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    effort_days = as.numeric(Days.Operational)
  ) %>%
  filter(
    !is.na(camera_id), !is.na(round),
    !is.na(start_date), !is.na(end_date),
    is.finite(effort_days), effort_days > 0
  ) %>%
  arrange(camera_id, round)

if (nrow(effort) != 344L || n_distinct(effort$camera_id, effort$round) != 344L) {
  stop("Expected 344 unique positive-effort camera-rounds; found ", nrow(effort))
}

site <- read_csv(cfg$site_covs, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    polygon = as.character(Polygon.Name)
  ) %>%
  arrange(camera_id)
if (nrow(site) != 60L || n_distinct(site$camera_id) != 60L) {
  stop("Expected 60 unique cameras in the site-covariate file.")
}

# Predator event extraction ---------------------------------------------------

message("Reading predator records from ", normalizePath(cfg$full_export))
predator_raw <- fread(
  cfg$full_export,
  select = c(
    "Camera Name", "Date/Time", "Species", "Trap Station Session ID",
    "Camelot GPS Latitude", "Camelot GPS Longitude"
  ),
  na.strings = c("", "NA"),
  showProgress = FALSE
)
setnames(
  predator_raw,
  c(
    "Camera Name", "Date/Time", "Species", "Trap Station Session ID",
    "Camelot GPS Latitude", "Camelot GPS Longitude"
  ),
  c("camera_id", "event_time_raw", "species", "session_id", "latitude", "longitude")
)
predator_raw <- predator_raw[species %chin% c("Canis lupus", "Canis aureus")]
predator_raw[, `:=`(
  raw_row_id = .I,
  camera_id = as.integer(camera_id),
  session_id = as.character(session_id),
  timestamp_uncorrected = parse_export_time(event_time_raw),
  taxon = fifelse(species == "Canis lupus", "Wolf", "Jackal"),
  latitude = as.numeric(latitude),
  longitude = as.numeric(longitude)
)]
predator_raw[, clock_offset_days := fcase(
  camera_id == 22L, 148,
  camera_id == 24L, 238,
  default = 0
)]
predator_raw[, timestamp := timestamp_uncorrected + clock_offset_days * 86400]
predator_raw[, event_date := as.Date(timestamp, tz = cfg$timezone)]
predator_raw <- predator_raw[
  !is.na(camera_id) & !is.na(session_id) & !is.na(timestamp)
]

predator_clock_audit <- predator_raw[
  camera_id %in% c(22L, 24L),
  .(
    raw_records = .N,
    offset_days = unique(clock_offset_days),
    uncorrected_first = min(timestamp_uncorrected),
    uncorrected_last = max(timestamp_uncorrected),
    corrected_first = min(timestamp),
    corrected_last = max(timestamp)
  ),
  by = .(camera_id, taxon)
]
write_csv(as_tibble(predator_clock_audit), out("predator_clock_corrections.csv"))

predator_round_matches <- as_tibble(predator_raw) %>%
  select(
    raw_row_id, camera_id, session_id, taxon, timestamp,
    event_date, latitude, longitude
  ) %>%
  inner_join(effort, by = "camera_id", relationship = "many-to-many") %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  arrange(raw_row_id, desc(start_date), desc(round)) %>%
  group_by(raw_row_id) %>%
  slice_head(n = 1) %>%
  ungroup()

predator_deduplicated <- predator_round_matches %>%
  group_by(camera_id, round, session_id, taxon, timestamp) %>%
  summarise(
    latitude = median(latitude, na.rm = TRUE),
    longitude = median(longitude, na.rm = TRUE),
    raw_records = n(),
    .groups = "drop"
  ) %>%
  arrange(camera_id, round, session_id, taxon, timestamp) %>%
  group_by(camera_id, round, session_id, taxon) %>%
  mutate(
    gap_minutes = as.numeric(difftime(timestamp, lag(timestamp), units = "mins")),
    independent_event_id = cumsum(is.na(gap_minutes) | gap_minutes > event_gap_minutes)
  ) %>%
  ungroup()

predator_events <- predator_deduplicated %>%
  group_by(camera_id, round, session_id, taxon, independent_event_id) %>%
  summarise(
    timestamp = min(timestamp),
    latitude = median(latitude, na.rm = TRUE),
    longitude = median(longitude, na.rm = TRUE),
    exact_timestamps = n(),
    raw_records = sum(raw_records),
    .groups = "drop"
  ) %>%
  mutate(
    solar_hour = apparent_solar_hour(timestamp, longitude),
    solar_bin = floor(solar_hour)
  ) %>%
  filter(is.finite(solar_hour), solar_bin %in% 0:23)

predator_event_audit <- bind_rows(
  as_tibble(predator_raw)[, c("taxon")] %>%
    count(taxon, name = "n") %>%
    mutate(stage = "Raw full-export annotations"),
  predator_round_matches %>%
    count(taxon, name = "n") %>%
    mutate(stage = "Clock-corrected annotations in active rounds"),
  predator_events %>%
    count(taxon, name = "n") %>%
    mutate(stage = paste0(event_gap_minutes, "-min independent events"))
) %>%
  select(stage, taxon, n) %>%
  arrange(stage, taxon)
write_csv(predator_event_audit, out("predator_event_audit.csv"))
write_csv(
  predator_events %>%
    select(
      camera_id, round, session_id, taxon, timestamp, solar_hour,
      solar_bin, exact_timestamps, raw_records
    ),
  out("predator_independent_events.csv")
)

# Matched hourly table --------------------------------------------------------

gazelle_raw <- read_csv(
  cfg$gazelle_events,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  mutate(
    event_row_id = row_number(),
    camera_id = as.integer(camera_name),
    timestamp = parse_gazelle_time(event_time),
    event_date = as.Date(format(timestamp, tz = cfg$timezone)),
    solar_hour = apparent_solar_hour(timestamp, camelot_gps_longitude),
    solar_bin = floor(solar_hour)
  )

gazelle_round_matches <- gazelle_raw %>%
  select(
    event_row_id, camera_id, gazelle_class, timestamp,
    event_date, solar_hour, solar_bin
  ) %>%
  inner_join(effort, by = "camera_id", relationship = "many-to-many") %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  arrange(event_row_id, desc(start_date), desc(round)) %>%
  group_by(event_row_id) %>%
  slice_head(n = 1) %>%
  ungroup()

gazelle_totals <- gazelle_round_matches %>% count(gazelle_class, name = "n")
observed_gazelle <- setNames(gazelle_totals$n, gazelle_totals$gazelle_class)
if (
  nrow(gazelle_round_matches) != sum(expected_gazelle_events) ||
    !all(observed_gazelle[names(expected_gazelle_events)] == expected_gazelle_events)
) {
  stop("Frozen gazelle event totals were not preserved.")
}

gazelle_counts <- gazelle_round_matches %>%
  count(camera_id, round, solar_bin, gazelle_class, name = "n_events") %>%
  pivot_wider(
    names_from = gazelle_class,
    values_from = n_events,
    values_fill = 0,
    names_glue = "{gazelle_class}_events"
  )

predator_counts <- predator_events %>%
  count(camera_id, round, solar_bin, taxon, name = "n_events") %>%
  pivot_wider(
    names_from = taxon,
    values_from = n_events,
    values_fill = 0,
    names_glue = "{tolower(taxon)}_events"
  )

hourly <- effort %>%
  select(camera_id, round, start_date, end_date, effort_days) %>%
  crossing(solar_bin = 0:23) %>%
  mutate(solar_hour = solar_bin + 0.5) %>%
  left_join(site, by = "camera_id") %>%
  left_join(gazelle_counts, by = c("camera_id", "round", "solar_bin")) %>%
  left_join(predator_counts, by = c("camera_id", "round", "solar_bin")) %>%
  mutate(
    across(
      c(juvenile_events, female_events, male_events, wolf_events, jackal_events),
      ~ replace_na(.x, 0L)
    )
  ) %>%
  arrange(camera_id, round, solar_bin) %>%
  group_by(camera_id, round) %>%
  mutate(
    wolf_events_3h = lag(wolf_events, default = last(wolf_events)) +
      wolf_events + lead(wolf_events, default = first(wolf_events)),
    jackal_events_3h = lag(jackal_events, default = last(jackal_events)) +
      jackal_events + lead(jackal_events, default = first(jackal_events))
  ) %>%
  ungroup()

if (nrow(hourly) != 344L * 24L) stop("Expected 8,256 matched hourly rows.")

standardize_unique <- function(values, ids) {
  unique_values <- tibble(id = ids, value = values) %>% distinct(id, .keep_all = TRUE)
  center <- mean(unique_values$value)
  scale <- sd(unique_values$value)
  list(z = (values - center) / scale, center = center, scale = scale)
}

add_exposure_decomposition <- function(data, count_column, prefix, window_hours) {
  count_values <- data[[count_column]]
  log_rate <- log1p(count_values / (data$effort_days * window_hours))
  result <- data %>%
    mutate(.log_rate = log_rate) %>%
    group_by(camera_id, round) %>%
    mutate(
      .round_mean = mean(.log_rate),
      .hour_deviation = .log_rate - .round_mean
    ) %>%
    ungroup() %>%
    group_by(camera_id) %>%
    mutate(
      .camera_mean = mean(.log_rate),
      .round_deviation = .round_mean - .camera_mean
    ) %>%
    ungroup()

  between <- standardize_unique(result$.camera_mean, result$camera_id)
  round_component <- standardize_unique(
    result$.round_deviation,
    paste(result$camera_id, result$round)
  )
  hour_scale <- sd(result$.hour_deviation)
  if (!is.finite(hour_scale) || hour_scale <= 0) stop("Invalid exposure scale: ", prefix)

  result[[paste0(prefix, "_log_rate")]] <- result$.log_rate
  result[[paste0(prefix, "_between_z")]] <- between$z
  result[[paste0(prefix, "_round_z")]] <- result$.round_deviation / round_component$scale
  result[[paste0(prefix, "_hour_z")]] <- result$.hour_deviation / hour_scale
  result <- result %>%
    select(-.log_rate, -.round_mean, -.hour_deviation, -.camera_mean, -.round_deviation)

  attributes <- tibble(
    predictor = prefix,
    component = c("between_camera", "between_round_within_camera", "hour_within_camera_round"),
    center = c(between$center, 0, 0),
    scale = c(between$scale, round_component$scale, hour_scale),
    window_hours = window_hours
  )
  list(data = result, attributes = attributes)
}

decompositions <- list(
  list(column = "wolf_events", prefix = "wolf_1h", window = 1),
  list(column = "jackal_events", prefix = "jackal_1h", window = 1),
  list(column = "wolf_events_3h", prefix = "wolf_3h", window = 3),
  list(column = "jackal_events_3h", prefix = "jackal_3h", window = 3)
)
scale_rows <- list()
for (i in seq_along(decompositions)) {
  definition <- decompositions[[i]]
  transformed <- add_exposure_decomposition(
    hourly,
    definition$column,
    definition$prefix,
    definition$window
  )
  hourly <- transformed$data
  scale_rows[[i]] <- transformed$attributes
}
exposure_scales <- bind_rows(scale_rows)

hourly <- hourly %>%
  mutate(
    camera_round = paste(camera_id, round, sep = "_"),
    camera_factor = factor(camera_id),
    camera_round_factor = factor(camera_round),
    polygon_factor = factor(polygon),
    round_factor = factor(round)
  )

write_csv(exposure_scales, out("exposure_scaling.csv"))
write_csv(
  hourly %>%
    select(
      camera_id, polygon, round, start_date, end_date, effort_days,
      solar_bin, solar_hour, juvenile_events, female_events, male_events,
      wolf_events, jackal_events, wolf_events_3h, jackal_events_3h,
      contains("_log_rate"), contains("_between_z"),
      contains("_round_z"), contains("_hour_z"), camera_round
    ),
  out("hourly_data.csv")
)

support_summary <- bind_rows(lapply(c("1h", "3h"), function(window) {
  wolf_hour <- hourly[[paste0("wolf_", window, "_hour_z")]]
  jackal_hour <- hourly[[paste0("jackal_", window, "_hour_z")]]
  tibble(
    window = window,
    statistic = c(
      "wolf_hour_z_min", "wolf_hour_z_q25", "wolf_hour_z_median",
      "wolf_hour_z_q75", "wolf_hour_z_max",
      "jackal_hour_z_min", "jackal_hour_z_q25", "jackal_hour_z_median",
      "jackal_hour_z_q75", "jackal_hour_z_max",
      "hourly_wolf_jackal_correlation", "cells_with_both_predators"
    ),
    value = c(
      min(wolf_hour), quantile(wolf_hour, 0.25), median(wolf_hour),
      quantile(wolf_hour, 0.75), max(wolf_hour),
      min(jackal_hour), quantile(jackal_hour, 0.25), median(jackal_hour),
      quantile(jackal_hour, 0.75), max(jackal_hour),
      cor(wolf_hour, jackal_hour),
      sum(hourly$wolf_events > 0 & hourly$jackal_events > 0)
    )
  )
}))
write_csv(support_summary, out("predictor_support.csv"))

# Dynamic circular GAMMs ------------------------------------------------------

dynamic_formula <- function(k, site_shape = FALSE) {
  site_term <- if (site_shape) {
    "s(solar_hour, by = camera_factor, bs = 'cc', k = 4, m = 1, id = 20) + "
  } else {
    ""
  }
  as.formula(paste0(
    "n_events ~ wolf_between_z * jackal_between_z + ",
    "wolf_round_z * jackal_round_z + ",
    "wolf_hour_z * jackal_hour_z + round_factor + ",
    "s(solar_hour, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_hour_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = jackal_hour_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_jackal_hour_z, bs = 'cc', k = ", k, ") + ",
    site_term,
    "s(camera_round_factor, bs = 're') + ",
    "s(camera_factor, bs = 're') + ",
    "s(polygon_factor, bs = 're') + offset(log(effort_days))"
  ))
}

prepare_model_data <- function(class_key, window = "1h") {
  wolf_prefix <- paste0("wolf_", window)
  jackal_prefix <- paste0("jackal_", window)
  hourly %>%
    transmute(
      n_events = .data[[response_columns[[class_key]]]],
      solar_hour,
      effort_days,
      wolf_between_z = .data[[paste0(wolf_prefix, "_between_z")]],
      jackal_between_z = .data[[paste0(jackal_prefix, "_between_z")]],
      wolf_round_z = .data[[paste0(wolf_prefix, "_round_z")]],
      jackal_round_z = .data[[paste0(jackal_prefix, "_round_z")]],
      wolf_hour_z = .data[[paste0(wolf_prefix, "_hour_z")]],
      jackal_hour_z = .data[[paste0(jackal_prefix, "_hour_z")]],
      wolf_jackal_hour_z = wolf_hour_z * jackal_hour_z,
      round_factor,
      camera_round_factor,
      camera_factor,
      polygon_factor,
      camera_id,
      round,
      solar_bin
    )
}

fit_model <- function(data, k, site_shape = FALSE) {
  bam(
    dynamic_formula(k, site_shape),
    data = data,
    family = nb(link = "log"),
    method = "fREML",
    discrete = TRUE,
    select = TRUE,
    knots = list(solar_hour = c(0, 24)),
    control = gam.control(trace = FALSE)
  )
}

convergence_label <- function(fit) {
  if (!is.null(fit$outer.info$conv)) return(as.character(fit$outer.info$conv))
  if (!is.null(fit$mgcv.conv)) {
    if (is.logical(fit$mgcv.conv) && length(fit$mgcv.conv) == 1L) {
      return(ifelse(fit$mgcv.conv, "full convergence", "not converged"))
    }
    if (is.list(fit$mgcv.conv) && !is.null(fit$mgcv.conv$fully.converged)) {
      return(ifelse(fit$mgcv.conv$fully.converged, "full convergence", "not converged"))
    }
  }
  "No convergence warning returned"
}

extract_summary <- function(fit, class_key, window, model_structure) {
  fit_summary <- summary(fit)
  p_table <- as.data.frame(fit_summary$p.table)
  s_table <- as.data.frame(fit_summary$s.table)
  smooth_terms <- rownames(fit_summary$s.table)
  bind_rows(
    tibble(
      class = class_labels[[class_key]],
      exposure_window = window,
      model_structure,
      term = rownames(fit_summary$p.table),
      term_type = "parametric",
      estimate = p_table[[1]],
      SE = p_table[[2]],
      EDF = NA_real_,
      reference_df = NA_real_,
      test_statistic = p_table[[3]],
      P = p_table[[4]]
    ),
    tibble(
      class = class_labels[[class_key]],
      exposure_window = window,
      model_structure,
      term = smooth_terms,
      term_type = case_when(
        grepl("camera|polygon", smooth_terms) ~ "random smooth",
        smooth_terms == "s(solar_hour)" ~ "baseline cyclic smooth",
        TRUE ~ "dynamic-predator cyclic smooth"
      ),
      estimate = NA_real_,
      SE = NA_real_,
      EDF = s_table[[1]],
      reference_df = s_table[[2]],
      test_statistic = s_table[[3]],
      P = s_table[[4]]
    )
  ) %>%
    mutate(supported = P < 0.05)
}

extract_diagnostics <- function(fit, data, class_key, window, structure, k) {
  fit_summary <- summary(fit)
  theta <- fit$family$getTheta(trans = TRUE)
  pearson <- residuals(fit, type = "pearson")
  fitted_counts <- fitted(fit)
  predicted_zero <- (theta / (theta + fitted_counts))^theta
  lag_summary <- data %>%
    mutate(pearson = pearson) %>%
    arrange(camera_round_factor, solar_bin) %>%
    group_by(camera_round_factor) %>%
    summarise(
      lag1 = suppressWarnings(cor(pearson[-n()], pearson[-1], use = "complete.obs")),
      .groups = "drop"
    )
  # Use every model row so the basis check is reproducible rather than based on
  # mgcv's default random 5,000-row subsample.
  set.seed(
    seed + 100L * match(class_key, class_keys) + ifelse(window == "1h", 1L, 3L)
  )
  k_table <- as.data.frame(k.check(fit, subsample = nrow(data), n.rep = 500))
  k_terms <- rownames(k_table)
  ecological_k <- k_table[
    grepl("solar_hour", k_terms) & !grepl("camera_factor", k_terms),
    , drop = FALSE
  ]
  concurvity_matrix <- concurvity(fit, full = TRUE)
  ecological_columns <- grepl("solar_hour", colnames(concurvity_matrix)) &
    !grepl("camera_factor", colnames(concurvity_matrix))
  tibble(
    class = class_labels[[class_key]],
    exposure_window = window,
    model_structure = structure,
    rows = nrow(data),
    events = sum(data$n_events),
    positive_cells = sum(data$n_events > 0),
    k = k,
    theta = theta,
    deviance_explained = fit_summary$dev.expl,
    convergence = convergence_label(fit),
    AIC = AIC(fit),
    pearson_dispersion = sum(pearson^2, na.rm = TRUE) / fit$df.residual,
    observed_zero_fraction = mean(data$n_events == 0),
    predicted_zero_fraction = mean(predicted_zero),
    mean_lag1_pearson = mean(lag_summary$lag1, na.rm = TRUE),
    minimum_ecological_k_index = min(ecological_k[[3]], na.rm = TRUE),
    minimum_ecological_k_P = min(ecological_k[[4]], na.rm = TRUE),
    maximum_ecological_concurvity = max(
      concurvity_matrix["estimate", ecological_columns], na.rm = TRUE
    )
  )
}

primary_models <- list()
candidate_models <- list()
model_data <- list()
site_shape_selection <- tibble()

for (class_key in class_keys) {
  k <- basis_k[[class_key]]
  data_i <- prepare_model_data(class_key, "1h")
  model_data[[class_key]] <- data_i
  message("Fitting 1-h dynamic GAMM for ", class_labels[[class_key]], "...")
  random_intercept_fit <- fit_model(data_i, k, site_shape = FALSE)
  candidate_models[[class_key]] <- random_intercept_fit
  selected_fit <- random_intercept_fit
  selected_structure <- "random intercepts"

  if (run_site_shape && class_key != "juvenile") {
    message("Fitting camera-specific curve sensitivity for ", class_labels[[class_key]], "...")
    site_shape_fit <- fit_model(data_i, k, site_shape = TRUE)
    delta_aic <- AIC(site_shape_fit) - AIC(random_intercept_fit)
    select_site_shape <- convergence_label(site_shape_fit) == "full convergence" &&
      is.finite(delta_aic) && delta_aic < -10
    site_shape_selection <- bind_rows(
      site_shape_selection,
      tibble(
        class = class_labels[[class_key]],
        random_intercept_AIC = AIC(random_intercept_fit),
        site_shape_AIC = AIC(site_shape_fit),
        delta_AIC = delta_aic,
        site_shape_convergence = convergence_label(site_shape_fit),
        selected = select_site_shape
      )
    )
    candidate_models[[paste0(class_key, "_site_shape")]] <- site_shape_fit
    if (select_site_shape) {
      selected_fit <- site_shape_fit
      selected_structure <- "camera-specific cyclic deviations"
    }
  }
  attr(selected_fit, "selected_structure") <- selected_structure
  primary_models[[class_key]] <- selected_fit
}

three_hour_models <- list()
if (run_3h_sensitivity) {
  for (class_key in class_keys) {
    k <- basis_k[[class_key]]
    data_i <- prepare_model_data(class_key, "3h")
    use_site_shape <- identical(
      attr(primary_models[[class_key]], "selected_structure"),
      "camera-specific cyclic deviations"
    )
    message("Fitting 3-h exposure sensitivity for ", class_labels[[class_key]], "...")
    three_hour_models[[class_key]] <- fit_model(data_i, k, use_site_shape)
    attr(three_hour_models[[class_key]], "selected_structure") <- ifelse(
      use_site_shape,
      "camera-specific cyclic deviations",
      "random intercepts"
    )
  }
}

model_summary <- bind_rows(lapply(class_keys, function(class_key) {
  extract_summary(
    primary_models[[class_key]],
    class_key,
    "1h",
    attr(primary_models[[class_key]], "selected_structure")
  )
}))
if (run_3h_sensitivity) {
  model_summary <- bind_rows(
    model_summary,
    bind_rows(lapply(class_keys, function(class_key) {
      extract_summary(
        three_hour_models[[class_key]],
        class_key,
        "3h",
        attr(three_hour_models[[class_key]], "selected_structure")
      )
    }))
  )
}

diagnostics <- bind_rows(lapply(class_keys, function(class_key) {
  extract_diagnostics(
    primary_models[[class_key]],
    model_data[[class_key]],
    class_key,
    "1h",
    attr(primary_models[[class_key]], "selected_structure"),
    basis_k[[class_key]]
  )
}))
if (run_3h_sensitivity) {
  diagnostics <- bind_rows(
    diagnostics,
    bind_rows(lapply(class_keys, function(class_key) {
      extract_diagnostics(
        three_hour_models[[class_key]],
        prepare_model_data(class_key, "3h"),
        class_key,
        "3h",
        attr(three_hour_models[[class_key]], "selected_structure"),
        basis_k[[class_key]]
      )
    }))
  )
}

write_csv(model_summary, out("model_summary.csv"))
write_csv(diagnostics, out("diagnostics.csv"))
write_csv(site_shape_selection, out("site_shape_selection.csv"))
saveRDS(primary_models, out("models.rds"))
saveRDS(candidate_models, out("candidate_models.rds"))
if (run_3h_sensitivity) saveRDS(three_hour_models, out("3h_sensitivity_models.rds"))

dynamic_coefficients <- model_summary %>%
  filter(
    term %in% c(
      "wolf_hour_z", "jackal_hour_z", "wolf_hour_z:jackal_hour_z"
    )
  ) %>%
  mutate(
    effect = case_when(
      term == "wolf_hour_z" ~ "Wolf activity",
      term == "jackal_hour_z" ~ "Jackal activity",
      TRUE ~ "Wolf x jackal interaction"
    ),
    beta_lower = estimate - qnorm(0.975) * SE,
    beta_upper = estimate + qnorm(0.975) * SE,
    incidence_rate_ratio = exp(estimate),
    IRR_lower = exp(beta_lower),
    IRR_upper = exp(beta_upper)
  ) %>%
  select(
    class, exposure_window, model_structure, effect, term,
    estimate, SE, beta_lower, beta_upper,
    incidence_rate_ratio, IRR_lower, IRR_upper, P, supported
  )
write_csv(dynamic_coefficients, out("dynamic_effect_coefficients.csv"))

# Population-level predictions ----------------------------------------------

simulate_curves <- function(fit, data, class_key) {
  grid <- expand_grid(
    solar_hour = seq(prediction_step / 2, 24 - prediction_step / 2, by = prediction_step),
    wolf_hour_z = c(-1, 1),
    jackal_hour_z = c(-1, 1)
  ) %>%
    mutate(
      wolf_jackal_hour_z = wolf_hour_z * jackal_hour_z,
      wolf_between_z = 0,
      jackal_between_z = 0,
      wolf_round_z = 0,
      jackal_round_z = 0,
      effort_days = 1,
      round_factor = factor(4, levels = levels(data$round_factor)),
      camera_round_factor = factor(
        levels(data$camera_round_factor)[1], levels = levels(data$camera_round_factor)
      ),
      camera_factor = factor(
        levels(data$camera_factor)[1], levels = levels(data$camera_factor)
      ),
      polygon_factor = factor(
        levels(data$polygon_factor)[1], levels = levels(data$polygon_factor)
      ),
      context = case_when(
        wolf_hour_z < 0 & jackal_hour_z < 0 ~ context_levels[1],
        wolf_hour_z < 0 & jackal_hour_z > 0 ~ context_levels[2],
        wolf_hour_z > 0 & jackal_hour_z < 0 ~ context_levels[3],
        TRUE ~ context_levels[4]
      ),
      context = factor(context, levels = context_levels)
    ) %>%
    arrange(context, solar_hour)

  excluded <- vapply(fit$smooth, function(x) x$label, character(1))
  excluded <- excluded[grepl("camera|polygon", excluded)]
  X <- predict(fit, newdata = grid, type = "lpmatrix", exclude = excluded)
  covariance <- tryCatch(vcov(fit, unconditional = TRUE), error = function(e) vcov(fit))
  draws <- rmvnorm(simulations, mean = coef(fit), sigma = covariance, method = "svd")
  point_raw <- exp(drop(X %*% coef(fit)))
  draw_raw <- exp(pmin(X %*% t(draws), 700))
  point_normalized <- point_raw
  draw_normalized <- draw_raw

  for (context_i in context_levels) {
    index <- which(grid$context == context_i)
    point_area <- sum(point_raw[index]) * prediction_step
    draw_area <- colSums(draw_raw[index, , drop = FALSE]) * prediction_step
    point_normalized[index] <- point_raw[index] / point_area
    draw_normalized[index, ] <- sweep(draw_raw[index, , drop = FALSE], 2, draw_area, "/")
  }

  grid %>%
    transmute(
      class = class_labels[[class_key]],
      solar_hour,
      wolf_dynamic_level = ifelse(wolf_hour_z < 0, "-1 SD", "+1 SD"),
      jackal_dynamic_level = ifelse(jackal_hour_z < 0, "-1 SD", "+1 SD"),
      wolf_hour_z,
      jackal_hour_z,
      context,
      predicted_raw = point_raw,
      raw_lower = apply(draw_raw, 1, quantile, 0.025),
      raw_upper = apply(draw_raw, 1, quantile, 0.975),
      predicted_normalized = point_normalized,
      lower = apply(draw_normalized, 1, quantile, 0.025),
      upper = apply(draw_normalized, 1, quantile, 0.975)
    )
}

predictions <- bind_rows(lapply(class_keys, function(class_key) {
  message("Simulating dynamic exposure curves for ", class_labels[[class_key]], "...")
  simulate_curves(primary_models[[class_key]], model_data[[class_key]], class_key)
})) %>%
  mutate(
    class = factor(class, levels = unname(class_labels[class_keys])),
    context = factor(context, levels = context_levels)
  )
write_csv(predictions, out("predictions.csv"))

prediction_plot <- ggplot(
  predictions,
  aes(
    x = solar_hour,
    y = predicted_normalized,
    color = context,
    fill = context,
    linetype = context
  )
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.08, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ class, ncol = 1, scales = "free_y") +
  scale_color_manual(values = context_colors, labels = context_labels) +
  scale_fill_manual(values = context_colors, labels = context_labels) +
  scale_linetype_manual(values = context_linetypes, labels = context_labels) +
  scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24), expand = c(0, 0)) +
  labs(
    title = "Gazelle diel activity under dynamic wolf-jackal exposure",
    subtitle = "Predictions use within-camera-round hourly predator anomalies at -1 and +1 SD",
    x = "Apparent solar time",
    y = "Normalized predicted activity density",
    color = NULL,
    fill = NULL,
    linetype = NULL,
    caption = paste0(
      "Predator predictors are clock-corrected 30-min independent event rates, not site assignments. ",
      "Shading shows simulation-based 95% intervals."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(color = "#4F5B60"),
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.key.width = grid::unit(1.5, "cm"),
    plot.caption = element_text(hjust = 0, color = "#5C666A")
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

ggsave(out("predicted_curves.png"), prediction_plot, width = 11, height = 10, dpi = 320)
ggsave(out("predicted_curves.pdf"), prediction_plot, width = 11, height = 10)

exposure_plot_data <- hourly %>%
  select(
    camera_id, round, solar_hour, wolf_events, jackal_events,
    wolf_1h_hour_z, jackal_1h_hour_z
  ) %>%
  pivot_longer(
    cols = c(wolf_1h_hour_z, jackal_1h_hour_z),
    names_to = "predator",
    values_to = "dynamic_z"
  ) %>%
  mutate(predator = ifelse(grepl("wolf", predator), "Wolf", "Jackal")) %>%
  group_by(predator, solar_hour) %>%
  summarise(
    mean_dynamic_z = mean(dynamic_z),
    q25 = quantile(dynamic_z, 0.25),
    q75 = quantile(dynamic_z, 0.75),
    .groups = "drop"
  )

exposure_plot <- ggplot(
  exposure_plot_data,
  aes(x = solar_hour, y = mean_dynamic_z, color = predator, fill = predator)
) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, color = "#626B70", linewidth = 0.4) +
  scale_color_manual(values = c(Wolf = "#173B63", Jackal = "#B12B32")) +
  scale_fill_manual(values = c(Wolf = "#173B63", Jackal = "#B12B32")) +
  scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24), expand = c(0, 0)) +
  labs(
    title = "Observed dynamic predator exposure by solar hour",
    subtitle = "Mean and interquartile range of within-camera-round standardized event rates",
    x = "Apparent solar time",
    y = "Dynamic predator activity (SD)",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )
ggsave(out("observed_predator_exposure.png"), exposure_plot, width = 9, height = 5.5, dpi = 320)
ggsave(out("observed_predator_exposure.pdf"), exposure_plot, width = 9, height = 5.5)

coefficient_plot_data <- dynamic_coefficients %>%
  filter(effect != "Wolf x jackal interaction") %>%
  mutate(
    class = factor(class, levels = unname(class_labels[class_keys])),
    effect = factor(effect, levels = c("Jackal activity", "Wolf activity")),
    exposure_window = factor(
      exposure_window,
      levels = c("1h", "3h"),
      labels = c("1-hour bins", "3-hour sensitivity")
    )
  )

coefficient_plot <- ggplot(
  coefficient_plot_data,
  aes(
    x = incidence_rate_ratio,
    y = effect,
    xmin = IRR_lower,
    xmax = IRR_upper,
    color = exposure_window
  )
) +
  geom_vline(xintercept = 1, color = "#5F6669", linewidth = 0.55) +
  geom_errorbar(
    orientation = "y",
    width = 0.12,
    linewidth = 0.65,
    position = position_dodge(width = 0.42)
  ) +
  geom_point(size = 2.5, position = position_dodge(width = 0.42)) +
  facet_wrap(~ class, ncol = 1) +
  scale_x_log10() +
  scale_color_manual(values = c("1-hour bins" = "#173B63", "3-hour sensitivity" = "#B36A2E")) +
  labs(
    title = "Gazelle event-rate associations with dynamic predator activity",
    subtitle = "Within-camera-round effects; interaction estimates are reported separately",
    x = "Incidence-rate ratio per 1-SD increase (95% CI)",
    y = NULL,
    color = NULL,
    caption = "Positive associations indicate temporal concordance, not attraction or abundance responses."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(color = "#4F5B60"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    plot.caption = element_text(hjust = 0, color = "#5C666A")
  )
ggsave(out("dynamic_effect_coefficients.png"), coefficient_plot, width = 8, height = 8, dpi = 320)
ggsave(out("dynamic_effect_coefficients.pdf"), coefficient_plot, width = 8, height = 8)

# Scientific report ----------------------------------------------------------

focal_smooths <- model_summary %>%
  filter(
    exposure_window %in% c("1h", "3h"),
    term %in% c(
      "s(solar_hour):wolf_hour_z",
      "s(solar_hour):jackal_hour_z",
      "s(solar_hour):wolf_jackal_hour_z"
    )
  ) %>%
  mutate(
    effect = case_when(
      grepl("wolf_jackal", term) ~ "Dynamic wolf x jackal by time",
      grepl("wolf_hour", term) ~ "Dynamic wolf by time",
      TRUE ~ "Dynamic jackal by time"
    )
  )

diagnostic_lines <- c(
  "# Dynamic predator circular GAMM diagnostics",
  "",
  "## Predictor construction",
  "",
  paste0(
    "Predator timestamps from cameras 22 and 24 were corrected by +148 and +238 days, respectively, before active-round assignment and apparent-solar-time conversion. Records were collapsed to ",
    event_gap_minutes,
    "-min independent events within camera, Camelot session, and taxon."
  ),
  paste0(
    "One-hour predator event rates were decomposed into between-camera means, between-round deviations within cameras, and hourly deviations within camera-rounds. The cyclic predator smooths used the final hourly component."
  ),
  "A three-hour circular moving-window rate was fitted as a sensitivity to sparse one-hour predator counts.",
  "",
  "## Model diagnostics",
  "",
  "| Class | Window | Structure | Theta | Deviance explained | Dispersion | Lag-1 residual r | Min k-index | Min k P | Max concurvity |",
  "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
  paste0(
    "| ", diagnostics$class, " | ", diagnostics$exposure_window, " | ",
    diagnostics$model_structure, " | ", fmt(diagnostics$theta), " | ",
    fmt(100 * diagnostics$deviance_explained, 1), "% | ",
    fmt(diagnostics$pearson_dispersion), " | ", fmt(diagnostics$mean_lag1_pearson),
    " | ", fmt(diagnostics$minimum_ecological_k_index), " | ",
    fmt_p(diagnostics$minimum_ecological_k_P), " | ",
    fmt(diagnostics$maximum_ecological_concurvity), " |"
  ),
  "",
  "## Dynamic predator smooths",
  "",
  "| Class | Window | Effect | EDF | Approximate P |",
  "|---|---|---|---:|---:|",
  paste0(
    "| ", focal_smooths$class, " | ", focal_smooths$exposure_window, " | ",
    focal_smooths$effect, " | ", fmt(focal_smooths$EDF, 2), " | ",
    fmt_p(focal_smooths$P), " |"
  ),
  "",
  "Smooth-term P values are approximate. Interpretation requires agreement across temporal windows, acceptable concurvity and basis checks, and biologically coherent prediction curves."
)
writeLines(diagnostic_lines, out("diagnostics.md"))

dynamic_table_lines <- c(
  "| Class | Window | Effect | Beta | 95% CI | IRR | 95% CI | P |",
  "|---|---|---|---:|---:|---:|---:|---:|",
  paste0(
    "| ", dynamic_coefficients$class, " | ",
    dynamic_coefficients$exposure_window, " | ",
    dynamic_coefficients$effect, " | ",
    fmt(dynamic_coefficients$estimate), " | ",
    fmt(dynamic_coefficients$beta_lower), " to ",
    fmt(dynamic_coefficients$beta_upper), " | ",
    fmt(dynamic_coefficients$incidence_rate_ratio), " | ",
    fmt(dynamic_coefficients$IRR_lower), " to ",
    fmt(dynamic_coefficients$IRR_upper), " | ",
    fmt_p(dynamic_coefficients$P), " |"
  )
)

male_1h <- dynamic_coefficients %>%
  filter(class == "Adult male", exposure_window == "1h")
male_3h <- dynamic_coefficients %>%
  filter(class == "Adult male", exposure_window == "3h")
juvenile_interaction_1h <- dynamic_coefficients %>%
  filter(
    class == "Juvenile", exposure_window == "1h",
    effect == "Wolf x jackal interaction"
  )
juvenile_interaction_3h <- dynamic_coefficients %>%
  filter(
    class == "Juvenile", exposure_window == "3h",
    effect == "Wolf x jackal interaction"
  )

results_lines <- c(
  "# Dynamic predator circular GAMM results",
  "",
  "## Analysis scope",
  "",
  paste0(
    "The analysis used ", format(nrow(hourly), big.mark = ","),
    " camera-round-hour cells, ", format(sum(expected_gazelle_events), big.mark = ","),
    " independent gazelle class-events, ",
    predator_event_audit$n[
      predator_event_audit$stage == paste0(event_gap_minutes, "-min independent events") &
        predator_event_audit$taxon == "Wolf"
    ],
    " independent wolf events, and ",
    format(
      predator_event_audit$n[
        predator_event_audit$stage == paste0(event_gap_minutes, "-min independent events") &
          predator_event_audit$taxon == "Jackal"
      ],
      big.mark = ","
    ),
    " independent jackal events."
  ),
  "Unlike the spatial-gradient GAMM, these models used predator rates that changed among camera-round-hour cells. Camera assignments to low or high site-level predator contexts were not used.",
  "",
  "## Dynamic event-rate associations",
  "",
  dynamic_table_lines,
  "",
  paste0(
    "Adult-male coefficients for hourly wolf and jackal activity were positive after accounting for their camera means, round-level deviations, broad diel shape, survey round, exposure, and spatial clustering. For the 1-hour predictors, the wolf IRR was ",
    fmt(male_1h$incidence_rate_ratio[male_1h$effect == "Wolf activity"]),
    " (95% CI ", fmt(male_1h$IRR_lower[male_1h$effect == "Wolf activity"]),
    "-", fmt(male_1h$IRR_upper[male_1h$effect == "Wolf activity"]),
    "; P ", fmt_p(male_1h$P[male_1h$effect == "Wolf activity"]),
    ") and the jackal IRR was ",
    fmt(male_1h$incidence_rate_ratio[male_1h$effect == "Jackal activity"]),
    " (", fmt(male_1h$IRR_lower[male_1h$effect == "Jackal activity"]),
    "-", fmt(male_1h$IRR_upper[male_1h$effect == "Jackal activity"]),
    "; P ", fmt_p(male_1h$P[male_1h$effect == "Jackal activity"]), ")."
  ),
  paste0(
    "Adult-male coefficients remained positive in the 3-hour sensitivity (wolf IRR ",
    fmt(male_3h$incidence_rate_ratio[male_3h$effect == "Wolf activity"]),
    ", P ", fmt_p(male_3h$P[male_3h$effect == "Wolf activity"]),
    "; jackal IRR ",
    fmt(male_3h$incidence_rate_ratio[male_3h$effect == "Jackal activity"]),
    ", P ", fmt_p(male_3h$P[male_3h$effect == "Jackal activity"]),
    "). However, the wolf interval excluded one only for the 1-hour definition, whereas the jackal interval excluded one only for the 3-hour definition; neither predator effect was independently supported in both windows. Adult-female hourly effects were unsupported in both windows."
  ),
  paste0(
    "Juvenile hourly wolf, jackal, and wolf x jackal effects were unsupported in both exposure definitions. The 1-hour interaction estimate was beta = ",
    fmt(juvenile_interaction_1h$estimate), " (P ",
    fmt_p(juvenile_interaction_1h$P),
    "), and the 3-hour estimate was beta = ",
    fmt(juvenile_interaction_3h$estimate), " (P ",
    fmt_p(juvenile_interaction_3h$P),
    "). Only seven 1-hour cells contained a juvenile event and both predator types, so the juvenile interaction remains weakly informed."
  ),
  "",
  "## Diel-shape modification",
  ""
)
for (class_key in class_keys) {
  class_label <- class_labels[[class_key]]
  smooth_i <- focal_smooths %>%
    filter(class == class_label, exposure_window == "1h")
  results_lines <- c(
    results_lines,
    paste0("### ", class_label),
    "",
    paste0(
      "Dynamic wolf-by-time: EDF = ",
      fmt(smooth_i$EDF[smooth_i$effect == "Dynamic wolf by time"], 2),
      ", P ", fmt_p(smooth_i$P[smooth_i$effect == "Dynamic wolf by time"]),
      "; dynamic jackal-by-time: EDF = ",
      fmt(smooth_i$EDF[smooth_i$effect == "Dynamic jackal by time"], 2),
      ", P ", fmt_p(smooth_i$P[smooth_i$effect == "Dynamic jackal by time"]),
      "; dynamic wolf x jackal by time: EDF = ",
      fmt(smooth_i$EDF[smooth_i$effect == "Dynamic wolf x jackal by time"], 2),
      ", P ", fmt_p(smooth_i$P[smooth_i$effect == "Dynamic wolf x jackal by time"]),
      "."
    ),
    ""
  )
}
results_lines <- c(
  results_lines,
  "No predator-by-time smooth was supported consistently across the 1-hour and 3-hour analyses. The models therefore provided no robust evidence that short-term wolf or jackal activity altered the shape or timing of gazelle diel distributions. Adult males showed a positive but window-dependent event-rate signal, not a shift to a different part of the day.",
  "",
  "## Interpretive limits",
  "",
  "The hourly predictors are average event rates for a solar-hour bin across each camera round, not exact predator presence on the same date as each gazelle event. Associations therefore measure fine-scale temporal concordance or segregation within camera-rounds, not immediate behavioral responses to individual encounters.",
  paste0(
    "All fits converged and the final baseline smooths used less than their maximum available degrees of freedom after increasing the basis sizes to ",
    basis_k[["juvenile"]], " for juveniles and ", basis_k[["female"]],
    " for adults. However, full-data k-index tests remained low (all P < 0.01) and residual lag-1 correlation remained especially high for juveniles. The dynamic analysis should therefore be treated as exploratory rather than replacing the primary camera-context analysis."
  ),
  "Observed predator rates are estimated with error. The three-hour sensitivity evaluates whether conclusions survive modest temporal smoothing, but a fully joint multivariate model would be required to propagate predator-process uncertainty exactly.",
  "Positive associations indicate temporal concordance and negative associations indicate segregation after controlling broad diel shape, survey round, exposure, and spatial clustering. Neither establishes causal attraction or avoidance."
)
writeLines(results_lines, out("results.md"))
writeLines(capture.output(sessionInfo()), out("session_info.txt"))

cat("Dynamic predator circular GAMM analysis completed.\n")
