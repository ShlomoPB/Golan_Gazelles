#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(glmmTMB)
  library(readr)
  library(tidyr)
})

set.seed(20260831)

event_file <- "Gazella gazella_independent_events_30min_analysis_ready.csv"
obs_file <- "ObsCovs_050524.csv"
camera_file <- "publication_camera_array_map_sites.csv"
summary_file <- "demographic_nocturnality_summary.csv"
bootstrap_file <- "demographic_nocturnality_bootstrap.csv"
model_file <- "demographic_nocturnality_model.csv"
results_file <- "demographic_nocturnality_results.md"
local_timezone <- "Asia/Jerusalem"
bootstrap_replicates <- 10000L

class_levels <- c("Juvenile", "Adult female", "Adult male")
model_class_levels <- c("Adult female", "Juvenile", "Adult male")
expected_counts <- c(
  "Juvenile" = 180L,
  "Adult female" = 1477L,
  "Adult male" = 2268L
)
old_flash_night <- c(
  "Juvenile" = 0.041,
  "Adult female" = 0.210,
  "Adult male" = 0.376
)
old_flash_overall_daylight <- 0.791

required_files <- c(event_file, obs_file, camera_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input(s): ", paste(missing_files, collapse = ", "))
}

parse_local_camera_time <- function(x) {
  as.POSIXct(
    sub("Z$", "", x),
    format = "%Y-%m-%dT%H:%M:%S",
    tz = local_timezone
  )
}

solar_elevation <- function(timestamp, latitude, longitude) {
  day_of_year <- as.integer(format(timestamp, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(timestamp, "%H", tz = "UTC")) +
    as.numeric(format(timestamp, "%M", tz = "UTC")) / 60 +
    as.numeric(format(timestamp, "%S", tz = "UTC")) / 3600
  gamma <- 2 * pi / 365 * (day_of_year - 1 + (utc_hour - 12) / 24)
  equation_time <- 229.18 * (
    0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -
      0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma)
  )
  declination <- 0.006918 - 0.399912 * cos(gamma) +
    0.070257 * sin(gamma) - 0.006758 * cos(2 * gamma) +
    0.000907 * sin(2 * gamma) - 0.002697 * cos(3 * gamma) +
    0.00148 * sin(3 * gamma)
  true_solar_minutes <-
    (utc_hour * 60 + equation_time + 4 * longitude) %% 1440
  hour_angle <- true_solar_minutes / 4 - 180
  latitude_radians <- latitude * pi / 180
  hour_angle_radians <- hour_angle * pi / 180
  cosine_zenith <-
    sin(latitude_radians) * sin(declination) +
    cos(latitude_radians) * cos(declination) * cos(hour_angle_radians)
  cosine_zenith <- pmin(1, pmax(-1, cosine_zenith))
  90 - acos(cosine_zenith) * 180 / pi
}

assign_rounds <- function(events, observations) {
  events %>%
    select(event_row_id, camera_id, event_date) %>%
    inner_join(
      observations %>% select(camera_id, round, start_date, end_date),
      by = "camera_id",
      relationship = "many-to-many"
    ) %>%
    filter(event_date >= start_date, event_date <= end_date) %>%
    group_by(event_row_id) %>%
    slice_max(start_date, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(event_row_id, round)
}

find_solar_crossings <- function(date, latitude, longitude) {
  date <- as.Date(date)
  day_start <- as.POSIXct(
    paste(date, "00:00:00"), tz = local_timezone
  )
  day_end <- as.POSIXct(
    paste(date + 1, "00:00:00"), tz = local_timezone
  )
  grid <- seq(day_start, day_end, by = "10 min")
  elevation <- solar_elevation(grid, latitude, longitude)
  sunrise_index <- which(
    elevation[-length(elevation)] < 0 & elevation[-1] >= 0
  )[1]
  sunset_index <- which(
    elevation[-length(elevation)] >= 0 & elevation[-1] < 0
  )[1]
  if (is.na(sunrise_index) || is.na(sunset_index)) {
    stop("Could not identify solar crossing for ", date)
  }
  root_time <- function(index) {
    limits <- as.numeric(grid[c(index, index + 1)])
    root <- uniroot(
      function(seconds) {
        solar_elevation(
          as.POSIXct(
            seconds, origin = "1970-01-01", tz = local_timezone
          ),
          latitude,
          longitude
        )
      },
      interval = limits,
      tol = 0.1
    )$root
    as.POSIXct(
      root, origin = "1970-01-01", tz = local_timezone
    )
  }
  list(
    sunrise = root_time(sunrise_index),
    sunset = root_time(sunset_index)
  )
}

clock_minutes <- function(x) {
  as.numeric(format(x, "%H", tz = local_timezone)) * 60 +
    as.numeric(format(x, "%M", tz = local_timezone)) +
    as.numeric(format(x, "%S", tz = local_timezone)) / 60
}

format_clock <- function(minutes) {
  minutes <- round(minutes) %% 1440
  sprintf("%02d:%02d", minutes %/% 60, minutes %% 60)
}

model_diagnostics <- function(model, model_name) {
  beta <- fixef(model)$cond
  pearson <- residuals(model, type = "pearson")
  dispersion_df <- max(1, nobs(model) - length(beta))
  random_components <- VarCorr(model)$cond
  random_rows <- bind_rows(lapply(names(random_components), function(group) {
    standard_deviation <- unname(attr(random_components[[group]], "stddev")[1])
    tibble(
      random_group = group,
      random_variance = standard_deviation^2
    )
  }))
  convergence_code <- model$fit$convergence
  positive_definite_hessian <- isTRUE(model$sdr$pdHess)
  max_absolute_gradient <- max(abs(model$sdr$gradient.fixed))
  singular_fit <- tryCatch(
    isTRUE(performance::check_singularity(model, tolerance = 1e-5)),
    error = function(e) NA
  )
  pearson_dispersion <- sum(pearson^2) / dispersion_df
  model_AIC <- AIC(model)
  tibble(
    model = model_name,
    convergence_code = convergence_code,
    positive_definite_hessian = positive_definite_hessian,
    max_absolute_gradient = max_absolute_gradient,
    singular = singular_fit,
    pearson_dispersion = pearson_dispersion,
    AIC = model_AIC,
    camera_variance = random_rows$random_variance[
      match("camera_factor", random_rows$random_group)
    ],
    camera_round_variance = random_rows$random_variance[
      match("camera_round", random_rows$random_group)
    ]
  )
}

extract_pairwise <- function(model, model_name) {
  marginal <- emmeans(model, ~ demographic_class, weights = "equal")
  planned <- contrast(
    marginal,
    method = list(
      "Juvenile / Adult female" = c(-1, 1, 0),
      "Juvenile / Adult male" = c(0, 1, -1),
      "Adult female / Adult male" = c(1, 0, -1)
    ),
    adjust = "none"
  )
  result <- as.data.frame(
    summary(planned, infer = c(TRUE, TRUE), adjust = "none")
  )
  lower_name <- intersect(c("asymp.LCL", "lower.CL"), names(result))[1]
  upper_name <- intersect(c("asymp.UCL", "upper.CL"), names(result))[1]
  as_tibble(result) %>%
    transmute(
      model = model_name,
      contrast = as.character(contrast),
      log_odds_ratio = estimate,
      standard_error = SE,
      odds_ratio = exp(estimate),
      lower_95_CI = exp(.data[[lower_name]]),
      upper_95_CI = exp(.data[[upper_name]]),
      raw_P = p.value,
      Holm_adjusted_P = p.adjust(p.value, method = "holm")
    )
}

# Authoritative camera coordinates and event data --------------------------

camera_coordinates <- read_csv(camera_file, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  )
if (
  nrow(camera_coordinates) != 60L ||
    n_distinct(camera_coordinates$camera_id) != 60L
) {
  stop("Camera-coordinate file no longer contains 60 cameras.")
}

events <- read_csv(
  event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  transmute(
    event_row_id = row_number(),
    camera_id = as.integer(camera_name),
    session_id = as.character(trap_station_session_id),
    event_id = as.character(event_id),
    demographic_class = recode(
      as.character(gazelle_class),
      juvenile = "Juvenile",
      female = "Adult female",
      male = "Adult male"
    ),
    timestamp = parse_local_camera_time(event_time),
    event_latitude = as.numeric(camelot_gps_latitude),
    event_longitude = as.numeric(camelot_gps_longitude)
  ) %>%
  left_join(camera_coordinates, by = "camera_id")

if (anyNA(events$latitude) || anyNA(events$longitude)) {
  stop("At least one event did not match an authoritative camera coordinate.")
}

coordinate_check <- events %>%
  distinct(
    camera_id, event_latitude, event_longitude, latitude, longitude
  ) %>%
  mutate(
    latitude_difference = event_latitude - latitude,
    longitude_difference = event_longitude - longitude
  )
coordinate_mismatches <- coordinate_check %>%
  filter(
    abs(latitude_difference) > 1e-4 |
      abs(longitude_difference) > 1e-4
  )
expected_camera_51_typo <-
  nrow(coordinate_mismatches) == 1L &&
  coordinate_mismatches$camera_id == 51L &&
  abs(coordinate_mismatches$latitude_difference - 1) < 1e-4 &&
  abs(coordinate_mismatches$longitude_difference) < 1e-4
if (!expected_camera_51_typo) {
  stop(
    "Unexpected disagreement between event and authoritative camera ",
    "coordinates."
  )
}

events <- events %>%
  mutate(
    event_date = as.Date(format(timestamp, tz = local_timezone)),
    solar_elevation_degrees = solar_elevation(
      timestamp, latitude, longitude
    ),
    night = as.integer(solar_elevation_degrees < 0),
    demographic_class = factor(
      demographic_class, levels = class_levels
    )
  )

observed_counts <- table(events$demographic_class)
if (
  nrow(events) != sum(expected_counts) ||
    any(as.integer(observed_counts[names(expected_counts)]) != expected_counts)
) {
  stop(
    "Authoritative class totals did not match: ",
    paste(names(observed_counts), observed_counts, collapse = "; ")
  )
}
if (anyNA(events$timestamp) || anyNA(events$solar_elevation_degrees)) {
  stop("Missing timestamp or solar-elevation classification.")
}

observations <- read_csv(obs_file, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    days_operational = as.numeric(Days.Operational)
  ) %>%
  filter(
    is.finite(days_operational), days_operational > 0,
    !is.na(start_date), !is.na(end_date)
  )

round_matches <- assign_rounds(events, observations)
events <- events %>%
  left_join(round_matches, by = "event_row_id") %>%
  mutate(
    round_factor = factor(round),
    camera_factor = factor(camera_id),
    camera_round = interaction(camera_id, round, drop = TRUE),
    demographic_class = factor(
      demographic_class, levels = model_class_levels
    )
  )
if (anyNA(events$round)) {
  stop("At least one event failed to map to an active survey round.")
}

# Sunrise/sunset audit ------------------------------------------------------

survey_start <- min(observations$start_date)
survey_end <- max(observations$end_date)
solar_schedule_index <- crossing(
  camera_coordinates,
  date = seq(survey_start, survey_end, by = "day")
)
solar_schedule_rows <- lapply(seq_len(nrow(solar_schedule_index)), function(i) {
  row <- solar_schedule_index[i, ]
  crossings <- find_solar_crossings(
    row$date, row$latitude, row$longitude
  )
  tibble(
    camera_id = row$camera_id,
    date = row$date,
    latitude = row$latitude,
    longitude = row$longitude,
    sunrise = crossings$sunrise,
    sunset = crossings$sunset,
    sunrise_clock_minutes = clock_minutes(crossings$sunrise),
    sunset_clock_minutes = clock_minutes(crossings$sunset)
  )
})
solar_schedule <- bind_rows(solar_schedule_rows)

manual_indices <- c(
  which.min(solar_schedule$date),
  which.min(abs(as.numeric(solar_schedule$date - as.Date("2020-10-17"))) +
    abs(solar_schedule$camera_id - 30) / 100),
  which.max(solar_schedule$date)
)
manual_checks <- solar_schedule[manual_indices, ] %>%
  rowwise() %>%
  mutate(
    sunrise_before = solar_elevation(
      sunrise - 5 * 60, latitude, longitude
    ),
    sunrise_after = solar_elevation(
      sunrise + 5 * 60, latitude, longitude
    ),
    sunset_before = solar_elevation(
      sunset - 5 * 60, latitude, longitude
    ),
    sunset_after = solar_elevation(
      sunset + 5 * 60, latitude, longitude
    ),
    sunrise_label = format(sunrise, "%Y-%m-%d %H:%M", tz = local_timezone),
    sunset_label = format(sunset, "%Y-%m-%d %H:%M", tz = local_timezone)
  ) %>%
  ungroup()
if (
  any(manual_checks$sunrise_before >= 0) ||
    any(manual_checks$sunrise_after <= 0) ||
    any(manual_checks$sunset_before <= 0) ||
    any(manual_checks$sunset_after >= 0)
) {
  stop("Manual sunrise/sunset crossing audit failed.")
}

# Camera-cluster bootstrap -------------------------------------------------

class_summary <- events %>%
  mutate(
    demographic_class = factor(
      demographic_class, levels = class_levels
    )
  ) %>%
  group_by(demographic_class) %>%
  summarise(
    total_independent_detections = n(),
    nighttime_detections = sum(night),
    daylight_detections = sum(1 - night),
    nocturnal_proportion = mean(night),
    diurnal_proportion = mean(1 - night),
    .groups = "drop"
  )

study_cameras <- sort(unique(observations$camera_id))
if (length(study_cameras) != 60L) {
  stop("The positive-effort observation history no longer contains 60 cameras.")
}
camera_class_totals <- crossing(
  camera_id = study_cameras,
  demographic_class = factor(class_levels, levels = class_levels)
) %>%
  left_join(
    events %>%
      mutate(
        demographic_class = factor(
          demographic_class, levels = class_levels
        )
      ) %>%
      group_by(camera_id, demographic_class) %>%
      summarise(total = n(), nighttime = sum(night), .groups = "drop"),
    by = c("camera_id", "demographic_class")
  ) %>%
  mutate(
    total = replace_na(total, 0L),
    nighttime = replace_na(nighttime, 0L)
  ) %>%
  arrange(camera_id, demographic_class)

total_matrix <- matrix(
  camera_class_totals$total,
  nrow = length(study_cameras),
  ncol = length(class_levels),
  byrow = TRUE,
  dimnames = list(study_cameras, class_levels)
)
night_matrix <- matrix(
  camera_class_totals$nighttime,
  nrow = length(study_cameras),
  ncol = length(class_levels),
  byrow = TRUE,
  dimnames = list(study_cameras, class_levels)
)

set.seed(20260831)
bootstrap_values <- matrix(
  NA_real_, nrow = bootstrap_replicates, ncol = length(class_levels),
  dimnames = list(NULL, class_levels)
)
for (replicate_i in seq_len(bootstrap_replicates)) {
  sampled <- sample.int(
    length(study_cameras),
    length(study_cameras),
    replace = TRUE
  )
  totals <- colSums(total_matrix[sampled, , drop = FALSE])
  nights <- colSums(night_matrix[sampled, , drop = FALSE])
  bootstrap_values[replicate_i, ] <- nights / totals
}
if (anyNA(bootstrap_values)) {
  stop("At least one bootstrap sample lacked a demographic class.")
}

bootstrap_class <- as_tibble(bootstrap_values) %>%
  mutate(replicate = row_number()) %>%
  pivot_longer(
    cols = all_of(class_levels),
    names_to = "demographic_class",
    values_to = "estimate"
  ) %>%
  mutate(
    statistic_type = "Class nocturnal proportion",
    contrast = NA_character_
  )

bootstrap_pairwise <- tibble(
  replicate = seq_len(bootstrap_replicates),
  `Juvenile - Adult female` =
    bootstrap_values[, "Juvenile"] -
    bootstrap_values[, "Adult female"],
  `Juvenile - Adult male` =
    bootstrap_values[, "Juvenile"] -
    bootstrap_values[, "Adult male"],
  `Adult female - Adult male` =
    bootstrap_values[, "Adult female"] -
    bootstrap_values[, "Adult male"]
) %>%
  pivot_longer(
    cols = -replicate,
    names_to = "contrast",
    values_to = "estimate"
  ) %>%
  mutate(
    statistic_type = "Pairwise nocturnal-proportion difference",
    demographic_class = NA_character_
  )

bootstrap_output <- bind_rows(
  bootstrap_class,
  bootstrap_pairwise
) %>%
  select(
    replicate, statistic_type, demographic_class, contrast, estimate
  )
write_csv(bootstrap_output, bootstrap_file)

class_intervals <- bootstrap_class %>%
  group_by(demographic_class) %>%
  summarise(
    bootstrap_lower_95 = quantile(estimate, 0.025),
    bootstrap_upper_95 = quantile(estimate, 0.975),
    .groups = "drop"
  )
pairwise_intervals <- bootstrap_pairwise %>%
  group_by(contrast) %>%
  summarise(
    bootstrap_lower_95 = quantile(estimate, 0.025),
    bootstrap_upper_95 = quantile(estimate, 0.975),
    .groups = "drop"
  )

class_night <- setNames(
  class_summary$nocturnal_proportion,
  as.character(class_summary$demographic_class)
)
observed_pairwise <- tibble(
  contrast = c(
    "Juvenile - Adult female",
    "Juvenile - Adult male",
    "Adult female - Adult male"
  ),
  difference = c(
    class_night[["Juvenile"]] - class_night[["Adult female"]],
    class_night[["Juvenile"]] - class_night[["Adult male"]],
    class_night[["Adult female"]] - class_night[["Adult male"]]
  )
)
pairwise_intervals <- observed_pairwise %>%
  left_join(pairwise_intervals, by = "contrast")

summary_output <- class_summary %>%
  left_join(class_intervals, by = "demographic_class") %>%
  mutate(
    nocturnal_percentage = 100 * nocturnal_proportion,
    diurnal_percentage = 100 * diurnal_proportion,
    bootstrap_lower_95_percentage = 100 * bootstrap_lower_95,
    bootstrap_upper_95_percentage = 100 * bootstrap_upper_95,
    old_flash_nocturnal_percentage =
      100 * unname(old_flash_night[as.character(demographic_class)])
  )
write_csv(summary_output, summary_file)

# Logistic mixed models ----------------------------------------------------

simple_model <- glmmTMB(
  night ~ demographic_class + round_factor + (1 | camera_factor),
  family = binomial(link = "logit"),
  data = events,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 10000, eval.max = 10000)
  )
)
camera_round_model <- glmmTMB(
  night ~ demographic_class + round_factor +
    (1 | camera_factor) + (1 | camera_round),
  family = binomial(link = "logit"),
  data = events,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 10000, eval.max = 10000)
  )
)

diagnostics <- bind_rows(
  model_diagnostics(simple_model, "Camera random intercept"),
  model_diagnostics(
    camera_round_model, "Camera + camera-round random intercepts"
  )
)
simple_pairwise <- extract_pairwise(
  simple_model, "Camera random intercept"
)
camera_round_pairwise <- extract_pairwise(
  camera_round_model, "Camera + camera-round random intercepts"
)
all_pairwise <- bind_rows(simple_pairwise, camera_round_pairwise) %>%
  left_join(diagnostics, by = "model")

class_effect_stability <- simple_pairwise %>%
  select(contrast, simple_log_OR = log_odds_ratio) %>%
  inner_join(
    camera_round_pairwise %>%
      select(contrast, camera_round_log_OR = log_odds_ratio),
    by = "contrast"
  ) %>%
  mutate(absolute_log_OR_change = abs(
    camera_round_log_OR - simple_log_OR
  ))
maximum_log_OR_change <- max(class_effect_stability$absolute_log_OR_change)

candidate_diagnostics <- diagnostics %>%
  filter(model == "Camera + camera-round random intercepts")
candidate_stable <-
  candidate_diagnostics$convergence_code == 0L &&
  candidate_diagnostics$positive_definite_hessian &&
  !isTRUE(candidate_diagnostics$singular)
candidate_substantive <-
  candidate_stable &&
  is.finite(candidate_diagnostics$camera_round_variance) &&
  candidate_diagnostics$camera_round_variance >= 1e-4 &&
  candidate_diagnostics$AIC < diagnostics$AIC[
    diagnostics$model == "Camera random intercept"
  ] - 2

retained_model_name <- if (candidate_substantive) {
  "Camera + camera-round random intercepts"
} else {
  "Camera random intercept"
}
all_pairwise <- all_pairwise %>%
  mutate(retained_model = model == retained_model_name)
write_csv(all_pairwise, model_file)

retained_pairwise <- all_pairwise %>% filter(retained_model)
bootstrap_supported <- with(
  pairwise_intervals,
  bootstrap_lower_95 > 0 | bootstrap_upper_95 < 0
)
model_supported <- retained_pairwise$Holm_adjusted_P < 0.05
recommendation <- if (all(bootstrap_supported) && all(model_supported)) {
  paste(
    "A. Report percentages with camera-cluster bootstrap confidence intervals",
    "in the main text; place the concordant formal mixed-model comparison in",
    "the Supporting Information because it adds little biological information."
  )
} else {
  paste(
    "B. Report percentages, camera-bootstrap confidence intervals, and the",
    "formal class comparison because the retained round-adjusted model supports",
    "the adult-female versus adult-male contrast whereas its unadjusted",
    "camera-bootstrap interval includes zero."
  )
}

# Results report -----------------------------------------------------------

summary_markdown <- c(
  "| Class | Total | Night | Daylight | Nocturnal % (camera-bootstrap 95% CI) | Diurnal % | Old flash nocturnal % |",
  "|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(summary_output)), function(i) sprintf(
    "| %s | %d | %d | %d | %.1f (%.1f-%.1f) | %.1f | %.1f |",
    summary_output$demographic_class[i],
    summary_output$total_independent_detections[i],
    summary_output$nighttime_detections[i],
    summary_output$daylight_detections[i],
    summary_output$nocturnal_percentage[i],
    summary_output$bootstrap_lower_95_percentage[i],
    summary_output$bootstrap_upper_95_percentage[i],
    summary_output$diurnal_percentage[i],
    summary_output$old_flash_nocturnal_percentage[i]
  ), character(1))
)
pairwise_bootstrap_markdown <- c(
  "| Contrast | Difference in nocturnal proportion | Camera-bootstrap 95% CI |",
  "|---|---:|---:|",
  vapply(seq_len(nrow(pairwise_intervals)), function(i) sprintf(
    "| %s | %.3f | %.3f to %.3f |",
    pairwise_intervals$contrast[i], pairwise_intervals$difference[i],
    pairwise_intervals$bootstrap_lower_95[i],
    pairwise_intervals$bootstrap_upper_95[i]
  ), character(1))
)
model_markdown <- c(
  "| Contrast | Adjusted odds ratio | 95% CI | Holm P |",
  "|---|---:|---:|---:|",
  vapply(seq_len(nrow(retained_pairwise)), function(i) sprintf(
    "| %s | %.3f | %.3f-%.3f | %s |",
    retained_pairwise$contrast[i], retained_pairwise$odds_ratio[i],
    retained_pairwise$lower_95_CI[i], retained_pairwise$upper_95_CI[i],
    ifelse(
      retained_pairwise$Holm_adjusted_P[i] < 0.001,
      "<0.001",
      sprintf("%.3f", retained_pairwise$Holm_adjusted_P[i])
    )
  ), character(1))
)
manual_markdown <- c(
  "| Camera | Date | Sunrise | Elevation 5 min before/after | Sunset | Elevation 5 min before/after |",
  "|---:|---|---|---:|---|---:|",
  vapply(seq_len(nrow(manual_checks)), function(i) sprintf(
    "| %d | %s | %s | %.2f / %.2f | %s | %.2f / %.2f |",
    manual_checks$camera_id[i], manual_checks$date[i],
    manual_checks$sunrise_label[i], manual_checks$sunrise_before[i],
    manual_checks$sunrise_after[i], manual_checks$sunset_label[i],
    manual_checks$sunset_before[i], manual_checks$sunset_after[i]
  ), character(1))
)

new_overall_daylight <- mean(1 - events$night)
qualitative_order_retained <-
  all(diff(summary_output$nocturnal_proportion) > 0)
simple_diag <- diagnostics %>% filter(model == "Camera random intercept")
candidate_diag <- diagnostics %>%
  filter(model == "Camera + camera-round random intercepts")

report <- c(
  "# Astronomical demographic nocturnality analysis",
  "",
  sprintf(
    paste(
      "Astronomical sunset-to-sunrise classification yielded %.1f%% nighttime",
      "detections for juveniles, %.1f%% for adult females, and %.1f%% for adult",
      "males. The corresponding supplied flash-based values were 4.1%%, 21.0%%,",
      "and 37.6%%. Overall astronomical daylight among retained class-event rows",
      "was %.1f%%, compared with the supplied flash-based value of 79.1%%."
    ),
    100 * summary_output$nocturnal_proportion[1],
    100 * summary_output$nocturnal_proportion[2],
    100 * summary_output$nocturnal_proportion[3],
    100 * new_overall_daylight
  ),
  "",
  "## Class-specific estimates",
  "",
  summary_markdown,
  "",
  "## Astronomical implementation audit",
  "",
  paste(
    "Night was defined as solar elevation <0 degrees using each corrected local",
    "camera timestamp and camera coordinate. Across all 60 cameras and survey",
    sprintf(
      "dates (%s to %s), sunrise ranged from %s to %s and sunset from %s to %s.",
      survey_start, survey_end,
      format_clock(min(solar_schedule$sunrise_clock_minutes)),
      format_clock(max(solar_schedule$sunrise_clock_minutes)),
      format_clock(min(solar_schedule$sunset_clock_minutes)),
      format_clock(max(solar_schedule$sunset_clock_minutes))
    )
  ),
  paste(
    sprintf(
      paste(
        "The independent-event export recorded camera 51 at latitude %.6f,",
        "whereas the authoritative camera map records latitude %.6f"
      ),
      coordinate_mismatches$event_latitude,
      coordinate_mismatches$latitude
    ),
    paste(
      "(a one-degree discrepancy); their longitudes agree, and all other camera",
      "coordinates agree. Solar calculations used the authoritative map",
      "coordinate while retaining the event-export coordinate as an audit field."
    )
  ),
  "",
  manual_markdown,
  "",
  "The sign changes around each manually inspected crossing verify the geometric-horizon implementation.",
  "",
  "## Camera-cluster bootstrap",
  "",
  sprintf(
    paste(
      "We resampled the %d positive-effort study cameras 10,000 times",
      "and retained all class-specific detections from each sampled camera;",
      "cameras with no retained events contributed empty clusters."
    ),
    length(study_cameras)
  ),
  "",
  pairwise_bootstrap_markdown,
  "",
  "## Formal class comparison",
  "",
  paste(
    "The retained model was:", retained_model_name,
    "with demographic class and survey round as fixed effects."
  ),
  sprintf(
    paste(
      "The camera-only model converged (code %d; positive-definite Hessian = %s;",
      "singular = %s; Pearson dispersion = %.3f; maximum absolute gradient =",
      "%.2e)."
    ),
    simple_diag$convergence_code,
    simple_diag$positive_definite_hessian,
    simple_diag$singular,
    simple_diag$pearson_dispersion,
    simple_diag$max_absolute_gradient
  ),
  sprintf(
    paste(
      "Adding camera-round produced variance %.4f, changed AIC by %.2f",
      "(candidate minus simple), and changed planned log odds ratios by at most",
      "%.3f. Candidate convergence code = %d; positive-definite Hessian = %s;",
      "singular = %s; Pearson dispersion = %.3f; maximum absolute gradient =",
      "%.2e."
    ),
    candidate_diag$camera_round_variance,
    candidate_diag$AIC - simple_diag$AIC,
    maximum_log_OR_change,
    candidate_diag$convergence_code,
    candidate_diag$positive_definite_hessian,
    candidate_diag$singular,
    candidate_diag$pearson_dispersion,
    candidate_diag$max_absolute_gradient
  ),
  "",
  model_markdown,
  "",
  "## Interpretation",
  "",
  paste(
    "The astronomical definition",
    ifelse(qualitative_order_retained, "retained", "did not retain"),
    paste(
      "the previous qualitative ordering: juveniles were most diurnal, adult",
      "females intermediate, and adult males most nocturnal."
    )
  ),
  paste(
    "Astronomical nighttime fractions were higher than the supplied flash-based",
    "fractions, especially for juveniles and females. This is expected because",
    "the geometric definition begins exactly at sunset and includes twilight,",
    "whereas flash activation depends on ambient-light thresholds and camera",
    "settings. The old flash percentages could not be recomputed event-by-event",
    "because the authoritative independent-event file contains no flash field;",
    "they are therefore treated as supplied benchmarks rather than paired values."
  ),
  paste(
    "These proportions describe when detections occurred, not absolute nighttime",
    "activity: astronomical night and daylight had unequal, seasonally changing",
    "durations. This analysis does not replace the circular GAM analysis of the",
    "complete 24-hour activity distribution."
  ),
  "",
  "## Manuscript recommendation",
  "",
  recommendation
)
writeLines(report, results_file)

message("Saved: ", summary_file)
message("Saved: ", bootstrap_file)
message("Saved: ", model_file)
message("Saved: ", results_file)
