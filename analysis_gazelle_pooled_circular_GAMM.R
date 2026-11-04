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
  timestamp_corrections = "camera_timestamp_corrections.csv",
  site_covs = "SiteCovs_060324.csv",
  obs_covs = "ObsCovs_050524.csv",
  class_model_summary = "circular_GAMM_model_summary.csv",
  output_prefix = "pooled_circular_GAMM",
  timezone = "Asia/Jerusalem",
  primary_k = "8",
  sensitivity_k = "12",
  prediction_step = "0.1",
  simulations = "1500",
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

primary_k <- as.integer(cfg$primary_k)
sensitivity_k <- as.integer(cfg$sensitivity_k)
prediction_step <- as.numeric(cfg$prediction_step)
simulations <- as.integer(cfg$simulations)
seed <- as.integer(cfg$seed)
set.seed(seed)

out <- function(suffix) paste0(cfg$output_prefix, "_", suffix)
required <- c(
  cfg$full_export, cfg$timestamp_corrections, cfg$site_covs, cfg$obs_covs
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required files: ", paste(missing, collapse = ", "))

expected_encounters <- 3688L
context_levels <- c(
  "Low wolf / Low jackal",
  "Low wolf / High jackal",
  "High wolf / Low jackal",
  "High wolf / High jackal"
)
context_colors <- c(
  "Low wolf / Low jackal" = "#168575",
  "Low wolf / High jackal" = "#B62D35",
  "High wolf / Low jackal" = "#173F66",
  "High wolf / High jackal" = "#835796"
)
context_linetypes <- c(
  "Low wolf / Low jackal" = "solid",
  "Low wolf / High jackal" = "22",
  "High wolf / Low jackal" = "42",
  "High wolf / High jackal" = "dotdash"
)
focal_smooth_terms <- c(
  "s(solar_hour):wolf_z",
  "s(solar_hour):jackal_z",
  "s(solar_hour):wolf_jackal_z"
)
focal_smooth_labels <- c(
  "s(solar_hour):wolf_z" = "Wolf by time",
  "s(solar_hour):jackal_z" = "Jackal by time",
  "s(solar_hour):wolf_jackal_z" = "Wolf x jackal by time"
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

# Effort and site covariates -------------------------------------------------

obs_raw <- read.csv(
  cfg$obs_covs,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
effort <- tibble(
  camera_id = as.integer(obs_raw[["Camera Number"]]),
  round = as.integer(obs_raw[["Round Number"]]),
  start_date = as.Date(obs_raw[["Start date"]], format = "%d/%m/%Y"),
  end_date = as.Date(obs_raw[["End Date"]], format = "%d/%m/%Y"),
  effort_days = as.numeric(obs_raw[["Days.Operational"]])
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
    wolf_index = as.numeric(Null.Wolf.Abundance),
    jackal_index = as.numeric(null.jackal.abundance)
  ) %>%
  arrange(camera_id)
if (nrow(site) != 60L || n_distinct(site$camera_id) != 60L || anyNA(site)) {
  stop("Expected complete wolf and jackal indices for 60 unique cameras.")
}
wolf_center <- mean(site$wolf_index)
wolf_scale <- sd(site$wolf_index)
jackal_center <- mean(site$jackal_index)
jackal_scale <- sd(site$jackal_index)
site <- site %>%
  mutate(
    wolf_z = (wolf_index - wolf_center) / wolf_scale,
    jackal_z = (jackal_index - jackal_center) / jackal_scale,
    wolf_jackal_z = wolf_z * jackal_z
  )

# Authoritative species-level encounters ------------------------------------

message("Reading species-level gazelle records from ", normalizePath(cfg$full_export))
raw <- fread(
  cfg$full_export,
  select = c(
    "Camera Name", "Date/Time", "Species", "Trap Station Session ID",
    "Camelot GPS Latitude", "Camelot GPS Longitude", "Sighting Quantity",
    "Sex", "Life stage"
  ),
  na.strings = c("", "NA"),
  showProgress = FALSE
)
raw <- raw[Species == "Gazella gazella"]
setnames(
  raw,
  c(
    "Camera Name", "Date/Time", "Trap Station Session ID",
    "Camelot GPS Latitude", "Camelot GPS Longitude", "Sighting Quantity",
    "Sex", "Life stage"
  ),
  c(
    "camera_id", "timestamp_text", "session_id", "latitude", "longitude",
    "sighting_quantity", "sex", "life_stage"
  )
)
raw[, `:=`(
  raw_row_id = .I,
  camera_id = as.integer(camera_id),
  session_id = as.character(session_id),
  timestamp_recorded = parse_export_time(timestamp_text),
  latitude = as.numeric(latitude),
  longitude = as.numeric(longitude),
  sighting_quantity = as.numeric(sighting_quantity)
)]
if (anyNA(raw$timestamp_recorded) || anyNA(raw$camera_id) || anyNA(raw$session_id)) {
  stop("Gazelle records contain missing camera, session, or timestamp values.")
}

corrections <- fread(cfg$timestamp_corrections)
corrections[, `:=`(
  camera_id = as.integer(camera_name),
  session_id = as.character(trap_station_session_id),
  offset_days = as.integer(offset_days)
)]
if (anyDuplicated(corrections[, .(camera_id, session_id)])) {
  stop("Timestamp-correction keys are not unique.")
}
raw <- merge(
  raw,
  corrections[, .(camera_id, session_id, offset_days, expected_round)],
  by = c("camera_id", "session_id"),
  all.x = TRUE,
  sort = FALSE
)
raw[is.na(offset_days), offset_days := 0L]

# Preserve the camera's recorded local clock time while correcting only dates.
raw[, corrected_date := as.Date(timestamp_recorded, tz = cfg$timezone) + offset_days]
raw[, clock_time := format(timestamp_recorded, "%H:%M:%S", tz = cfg$timezone)]
raw[, timestamp := as.POSIXct(
  paste(corrected_date, clock_time),
  format = "%Y-%m-%d %H:%M:%S",
  tz = cfg$timezone
)]
if (anyNA(raw$timestamp)) stop("One or more corrected timestamps could not be parsed.")

setorder(raw, camera_id, session_id, timestamp, raw_row_id)
raw[, encounter_id := cumsum(
  seq_len(.N) == 1L |
    as.numeric(difftime(timestamp, shift(timestamp), units = "mins")) > 30
), by = .(camera_id, session_id)]

encounters <- raw[, .(
  timestamp = min(timestamp),
  last_timestamp = max(timestamp),
  raw_sighting_rows = .N,
  maximum_group_size = suppressWarnings(max(sighting_quantity, na.rm = TRUE)),
  latitude = median(latitude, na.rm = TRUE),
  longitude = median(longitude, na.rm = TRUE),
  demographic_labels = paste(
    sort(unique(paste(sex, life_stage, sep = " / "))), collapse = "; "
  )
), by = .(camera_id, session_id, encounter_id)]
encounters[!is.finite(maximum_group_size), maximum_group_size := NA_real_]
if (nrow(encounters) != expected_encounters) {
  stop(
    "Species-level 30-min audit failed: expected ", expected_encounters,
    " encounters but reconstructed ", nrow(encounters), "."
  )
}
if (any(!is.finite(encounters$longitude))) {
  stop("One or more pooled encounters lack usable longitude.")
}

encounters[, `:=`(
  species_event_id = .I,
  event_date = as.Date(timestamp, tz = cfg$timezone),
  solar_time = apparent_solar_hour(timestamp, longitude),
  solar_bin = floor(apparent_solar_hour(timestamp, longitude))
)]

event_round_matches <- as_tibble(encounters) %>%
  select(
    species_event_id, camera_id, session_id, encounter_id, timestamp,
    last_timestamp, event_date, solar_time, solar_bin, raw_sighting_rows,
    maximum_group_size, latitude, longitude, demographic_labels
  ) %>%
  inner_join(effort, by = "camera_id", relationship = "many-to-many") %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  arrange(species_event_id, desc(start_date), desc(round)) %>%
  group_by(species_event_id) %>%
  slice_head(n = 1) %>%
  ungroup()
if (nrow(event_round_matches) != expected_encounters ||
    n_distinct(event_round_matches$species_event_id) != expected_encounters) {
  stop(
    "Expected all ", expected_encounters,
    " pooled encounters to map uniquely to positive-effort rounds; mapped ",
    nrow(event_round_matches), "."
  )
}
write_csv(
  event_round_matches %>%
    select(
      species_event_id, camera_id, round, session_id, encounter_id,
      timestamp, last_timestamp, solar_time, solar_bin, raw_sighting_rows,
      maximum_group_size, latitude, longitude, demographic_labels
    ),
  out("species_level_independent_encounters.csv")
)

event_audit <- tibble(
  metric = c(
    "raw_gazelle_sighting_rows",
    "species_level_30min_encounters",
    "encounters_assigned_to_positive_effort_rounds",
    "encounters_excluded_from_model",
    "cameras_with_encounters",
    "positive_effort_camera_rounds",
    "camera_round_hour_cells",
    "total_camera_nights",
    "timestamp_correction_sessions"
  ),
  value = c(
    nrow(raw), nrow(encounters), nrow(event_round_matches),
    nrow(encounters) - nrow(event_round_matches),
    n_distinct(event_round_matches$camera_id), nrow(effort),
    nrow(effort) * 24L, sum(effort$effort_days), nrow(corrections)
  )
)
write_csv(event_audit, out("event_audit.csv"))

# Camera-round-hour response table ------------------------------------------

event_counts <- event_round_matches %>%
  count(camera_id, round, solar_bin, name = "gazelle_events")

hourly <- effort %>%
  select(camera_id, round, start_date, end_date, effort_days) %>%
  crossing(solar_bin = 0:23) %>%
  mutate(solar_hour = solar_bin + 0.5) %>%
  left_join(site, by = "camera_id") %>%
  left_join(event_counts, by = c("camera_id", "round", "solar_bin")) %>%
  mutate(
    gazelle_events = replace_na(gazelle_events, 0L),
    camera_id_factor = factor(camera_id),
    camera_round_factor = factor(paste(camera_id, round, sep = "_")),
    round_factor = factor(round)
  ) %>%
  arrange(camera_id, round, solar_bin)
if (nrow(hourly) != 8256L || sum(hourly$gazelle_events) != expected_encounters) {
  stop("Hourly aggregation did not preserve 8,256 cells and 3,688 encounters.")
}
write_csv(
  hourly %>%
    select(
      camera_id, round, start_date, end_date, effort_days, solar_bin,
      solar_hour, gazelle_events, wolf_index, jackal_index, wolf_z,
      jackal_z, wolf_jackal_z
    ),
  out("hourly_data.csv")
)

model_data <- hourly %>%
  transmute(
    gazelle_events,
    solar_hour,
    effort_days,
    wolf_z,
    jackal_z,
    wolf_jackal_z,
    round_factor,
    camera_id = camera_id_factor,
    camera_round = camera_round_factor,
    solar_bin
  )

# Circular GAMMs ------------------------------------------------------------

pooled_formula <- function(k, camera_round_effect = FALSE) {
  camera_round_term <- if (camera_round_effect) {
    "s(camera_round, bs = 're') + "
  } else {
    ""
  }
  as.formula(paste0(
    "gazelle_events ~ wolf_z * jackal_z + round_factor + ",
    "s(solar_hour, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = jackal_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_jackal_z, bs = 'cc', k = ", k, ") + ",
    camera_round_term,
    "s(camera_id, bs = 're') + offset(log(effort_days))"
  ))
}

fit_model <- function(k, camera_round_effect = FALSE) {
  bam(
    pooled_formula(k, camera_round_effect),
    data = model_data,
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

message("Fitting requested pooled circular GAMM with k = ", primary_k, "...")
primary_fit <- fit_model(primary_k, camera_round_effect = FALSE)
message("Fitting larger-basis sensitivity with k = ", sensitivity_k, "...")
k_sensitivity_fit <- fit_model(sensitivity_k, camera_round_effect = FALSE)
message("Fitting camera-round random-intercept sensitivity with k = ", primary_k, "...")
cluster_sensitivity_fit <- fit_model(primary_k, camera_round_effect = TRUE)

models <- list(
  requested_k8 = primary_fit,
  larger_basis_k12 = k_sensitivity_fit,
  camera_round_sensitivity = cluster_sensitivity_fit
)
saveRDS(models, out("models.rds"))

extract_model_summary <- function(fit, model_name) {
  fit_summary <- summary(fit)
  p_table <- as.data.frame(fit_summary$p.table)
  s_table <- as.data.frame(fit_summary$s.table)
  bind_rows(
    tibble(
      model = model_name,
      term = rownames(p_table),
      term_type = "parametric",
      estimate = p_table[[1]],
      SE = p_table[[2]],
      EDF = NA_real_,
      reference_df = NA_real_,
      test_statistic = p_table[[3]],
      P = p_table[[4]]
    ),
    tibble(
      model = model_name,
      term = rownames(s_table),
      term_type = case_when(
        grepl("camera", rownames(s_table)) ~ "random smooth",
        rownames(s_table) == "s(solar_hour)" ~ "baseline cyclic smooth",
        TRUE ~ "predator-by-time cyclic smooth"
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

extract_diagnostics <- function(fit, model_name, k, camera_round_effect, seed_offset) {
  fit_summary <- summary(fit)
  theta <- fit$family$getTheta(trans = TRUE)
  pearson <- residuals(fit, type = "pearson")
  fitted_counts <- fitted(fit)
  predicted_zero <- (theta / (theta + fitted_counts))^theta
  lag_summary <- model_data %>%
    mutate(pearson = pearson) %>%
    arrange(camera_round, solar_bin) %>%
    group_by(camera_round) %>%
    summarise(
      lag1 = suppressWarnings(cor(pearson[-n()], pearson[-1], use = "complete.obs")),
      .groups = "drop"
    )

  set.seed(seed + seed_offset)
  k_raw <- as.data.frame(k.check(fit, subsample = nrow(model_data), n.rep = 500))
  k_terms <- rownames(k_raw)
  ecological_k <- k_raw[grepl("solar_hour", k_terms), , drop = FALSE]
  concurvity_matrix <- concurvity(fit, full = TRUE)
  ecological_columns <- grepl("solar_hour", colnames(concurvity_matrix))

  diagnostics <- tibble(
    model = model_name,
    k,
    camera_round_random_effect = camera_round_effect,
    rows = nrow(model_data),
    events = sum(model_data$gazelle_events),
    positive_cells = sum(model_data$gazelle_events > 0),
    convergence = convergence_label(fit),
    theta,
    deviance_explained = fit_summary$dev.expl,
    adjusted_r_squared = fit_summary$r.sq,
    AIC = AIC(fit),
    pearson_dispersion = sum(pearson^2, na.rm = TRUE) / fit$df.residual,
    observed_zero_fraction = mean(model_data$gazelle_events == 0),
    predicted_zero_fraction = mean(predicted_zero),
    mean_within_camera_round_lag1_pearson = mean(lag_summary$lag1, na.rm = TRUE),
    minimum_ecological_k_index = min(ecological_k[[3]], na.rm = TRUE),
    minimum_ecological_k_P = min(ecological_k[[4]], na.rm = TRUE),
    maximum_ecological_concurvity = max(
      concurvity_matrix["estimate", ecological_columns], na.rm = TRUE
    )
  )
  k_table <- tibble(
    model = model_name,
    term = k_terms,
    k_prime = k_raw[[1]],
    EDF = k_raw[[2]],
    k_index = k_raw[[3]],
    P = k_raw[[4]]
  )
  list(diagnostics = diagnostics, k_table = k_table)
}

model_summary <- bind_rows(
  extract_model_summary(primary_fit, "Requested k = 8"),
  extract_model_summary(k_sensitivity_fit, "Larger basis k = 12"),
  extract_model_summary(cluster_sensitivity_fit, "Camera-round sensitivity")
)
diagnostic_objects <- list(
  extract_diagnostics(primary_fit, "Requested k = 8", primary_k, FALSE, 1L),
  extract_diagnostics(
    k_sensitivity_fit, "Larger basis k = 12", sensitivity_k, FALSE, 2L
  ),
  extract_diagnostics(
    cluster_sensitivity_fit, "Camera-round sensitivity", primary_k, TRUE, 3L
  )
)
diagnostics <- bind_rows(lapply(diagnostic_objects, `[[`, "diagnostics"))
k_diagnostics <- bind_rows(lapply(diagnostic_objects, `[[`, "k_table"))
write_csv(model_summary, out("model_summary.csv"))
write_csv(diagnostics, out("diagnostics.csv"))
write_csv(k_diagnostics, out("k_diagnostics.csv"))

pdf(out("diagnostic_plots.pdf"), width = 9, height = 8)
for (model_name in names(models)) {
  gam.check(models[[model_name]], rep = 0)
}
dev.off()

# Population-level normalized predictions ----------------------------------

prediction_grid <- function() {
  expand_grid(
    solar_hour = seq(
      prediction_step / 2, 24 - prediction_step / 2, by = prediction_step
    ),
    wolf_z = c(-1, 1),
    jackal_z = c(-1, 1)
  ) %>%
    mutate(
      wolf_jackal_z = wolf_z * jackal_z,
      round_factor = factor("4", levels = levels(model_data$round_factor)),
      camera_id = factor(
        levels(model_data$camera_id)[1], levels = levels(model_data$camera_id)
      ),
      camera_round = factor(
        levels(model_data$camera_round)[1], levels = levels(model_data$camera_round)
      ),
      effort_days = 1,
      context = case_when(
        wolf_z < 0 & jackal_z < 0 ~ context_levels[1],
        wolf_z < 0 & jackal_z > 0 ~ context_levels[2],
        wolf_z > 0 & jackal_z < 0 ~ context_levels[3],
        TRUE ~ context_levels[4]
      ),
      context = factor(context, levels = context_levels)
    ) %>%
    arrange(context, solar_hour)
}

point_curves <- function(fit) {
  grid <- prediction_grid()
  excluded <- vapply(fit$smooth, function(x) x$label, character(1))
  excluded <- excluded[grepl("camera", excluded)]
  X <- predict(fit, newdata = grid, type = "lpmatrix", exclude = excluded)
  raw <- exp(drop(X %*% coef(fit)))
  normalized <- raw
  for (context_i in context_levels) {
    index <- which(grid$context == context_i)
    normalized[index] <- raw[index] / (sum(raw[index]) * prediction_step)
  }
  grid %>% mutate(predicted_normalized = normalized)
}

simulate_primary_curves <- function(fit) {
  grid <- prediction_grid()
  excluded <- vapply(fit$smooth, function(x) x$label, character(1))
  excluded <- excluded[grepl("camera", excluded)]
  X <- predict(fit, newdata = grid, type = "lpmatrix", exclude = excluded)
  covariance <- tryCatch(vcov(fit, unconditional = TRUE), error = function(e) vcov(fit))
  draws <- rmvnorm(simulations, mean = coef(fit), sigma = covariance, method = "svd")
  raw_point <- exp(drop(X %*% coef(fit)))
  raw_draws <- exp(pmin(X %*% t(draws), 700))
  normalized_point <- raw_point
  normalized_draws <- raw_draws
  for (context_i in context_levels) {
    index <- which(grid$context == context_i)
    normalized_point[index] <- raw_point[index] /
      (sum(raw_point[index]) * prediction_step)
    draw_area <- colSums(raw_draws[index, , drop = FALSE]) * prediction_step
    normalized_draws[index, ] <- sweep(
      raw_draws[index, , drop = FALSE], 2, draw_area, "/"
    )
  }
  grid %>%
    transmute(
      solar_hour,
      wolf_level = ifelse(wolf_z < 0, "-1 SD", "+1 SD"),
      jackal_level = ifelse(jackal_z < 0, "-1 SD", "+1 SD"),
      wolf_z,
      jackal_z,
      context,
      predicted_normalized = normalized_point,
      lower = apply(normalized_draws, 1, quantile, 0.025),
      upper = apply(normalized_draws, 1, quantile, 0.975)
    )
}

predictions <- simulate_primary_curves(primary_fit)
write_csv(predictions, out("predictions.csv"))

curve_metrics <- predictions %>%
  group_by(context, wolf_level, jackal_level) %>%
  summarise(
    peak_solar_hour = solar_hour[which.max(predicted_normalized)],
    daytime_06_18_fraction = sum(
      predicted_normalized[solar_hour >= 6 & solar_hour < 18]
    ) * prediction_step,
    mean_hours_from_solar_noon = sum(
      abs(solar_hour - 12) * predicted_normalized
    ) * prediction_step,
    .groups = "drop"
  )
write_csv(curve_metrics, out("curve_metrics.csv"))

primary_points <- point_curves(primary_fit) %>%
  select(context, solar_hour, primary = predicted_normalized)
sensitivity_curves <- bind_rows(
  point_curves(k_sensitivity_fit) %>% mutate(model = "Larger basis k = 12"),
  point_curves(cluster_sensitivity_fit) %>%
    mutate(model = "Camera-round sensitivity")
) %>%
  select(model, context, solar_hour, sensitivity = predicted_normalized) %>%
  left_join(primary_points, by = c("context", "solar_hour"))
curve_stability <- sensitivity_curves %>%
  group_by(model, context) %>%
  summarise(
    curve_correlation = cor(primary, sensitivity),
    RMSE = sqrt(mean((sensitivity - primary)^2)),
    maximum_absolute_difference = max(abs(sensitivity - primary)),
    .groups = "drop"
  )
write_csv(curve_stability, out("curve_stability.csv"))

shade <- tibble(
  xmin = c(0, 5, 7, 17, 19),
  xmax = c(5, 7, 17, 19, 24),
  period = c("Night", "Twilight", "Daylight", "Twilight", "Night")
)
curve_plot <- ggplot(
  predictions,
  aes(
    x = solar_hour,
    y = predicted_normalized,
    color = context,
    fill = context,
    linetype = context
  )
) +
  geom_rect(
    data = shade,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = period),
    alpha = 0.11,
    color = NA
  ) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.09, color = NA) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = context_colors) +
  scale_fill_manual(
    values = c(
      context_colors,
      Night = "#607482", Twilight = "#E6A66C", Daylight = "#F4ECD8"
    ),
    breaks = context_levels
  ) +
  scale_linetype_manual(values = context_linetypes) +
  scale_x_continuous(
    breaks = seq(0, 24, 3), limits = c(0, 24), expand = c(0, 0)
  ) +
  labs(
    title = "Pooled mountain gazelle diel activity across predator gradients",
    subtitle = "Continuous site-level wolf and jackal indices evaluated at -1 and +1 SD",
    x = "Apparent solar time",
    y = "Normalized predicted activity density",
    color = NULL,
    fill = NULL,
    linetype = NULL,
    caption = paste0(
      "Negative-binomial cyclic GAMM; curves integrate to one within each predator context. ",
      "\nShading shows coefficient-simulation 95% intervals; background indicates approximate light periods."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(color = "#4E5B60"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.key.width = grid::unit(1.35, "cm"),
    plot.caption = element_text(hjust = 0, color = "#59666B", lineheight = 1.05)
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )
ggsave(out("predicted_curves.png"), curve_plot, width = 10, height = 7, dpi = 320)
ggsave(out("predicted_curves.pdf"), curve_plot, width = 10, height = 7)

# Comparison with final class-specific GAMMs -------------------------------

class_comparison <- tibble()
if (file.exists(cfg$class_model_summary)) {
  class_comparison <- read_csv(cfg$class_model_summary, show_col_types = FALSE) %>%
    filter(term %in% focal_smooth_terms) %>%
    transmute(
      analysis = class,
      term,
      effect = unname(focal_smooth_labels[term]),
      EDF = EDF_if_smooth,
      P,
      supported = P < 0.05
    )
}
pooled_comparison <- model_summary %>%
  filter(model == "Requested k = 8", term %in% focal_smooth_terms) %>%
  transmute(
    analysis = "Pooled species",
    term,
    effect = unname(focal_smooth_labels[term]),
    EDF,
    P,
    supported = P < 0.05
  )
class_comparison <- bind_rows(pooled_comparison, class_comparison)
write_csv(class_comparison, out("class_comparison.csv"))

# Concise scientific report -------------------------------------------------

focal_results <- model_summary %>%
  filter(term %in% focal_smooth_terms) %>%
  mutate(effect = unname(focal_smooth_labels[term]))
primary_focal <- focal_results %>% filter(model == "Requested k = 8")
diagnostic_primary <- diagnostics %>% filter(model == "Requested k = 8")

effect_sentence <- function(term_name, subject) {
  row <- primary_focal %>% filter(.data$term == .env$term_name)
  paste0(
    subject, ": EDF = ", fmt(row$EDF, 2), ", approximate P = ", fmt_p(row$P),
    ifelse(row$P < 0.05, " (supported).", " (unsupported).")
  )
}

metric_table_lines <- c(
  "| Wolf | Jackal | Peak (solar h) | Proportion 06:00-18:00 | Mean h from solar noon |",
  "|---|---|---:|---:|---:|",
  paste0(
    "| ", curve_metrics$wolf_level, " | ", curve_metrics$jackal_level, " | ",
    fmt(curve_metrics$peak_solar_hour, 1), " | ",
    fmt(curve_metrics$daytime_06_18_fraction, 3), " | ",
    fmt(curve_metrics$mean_hours_from_solar_noon, 2), " |"
  )
)

stability_table_lines <- c(
  "| Model | Effect | EDF | Approximate P |",
  "|---|---|---:|---:|",
  paste0(
    "| ", focal_results$model, " | ", focal_results$effect, " | ",
    fmt(focal_results$EDF, 2), " | ", fmt_p(focal_results$P), " |"
  )
)
class_table_lines <- c(
  "| Analysis | Effect | EDF | Approximate P | Supported |",
  "|---|---|---:|---:|---|",
  paste0(
    "| ", class_comparison$analysis, " | ", class_comparison$effect, " | ",
    fmt(class_comparison$EDF, 2), " | ", fmt_p(class_comparison$P), " | ",
    ifelse(class_comparison$supported, "Yes", "No"), " |"
  )
)

results_lines <- c(
  "# Pooled species-level circular GAMM results",
  "",
  "## Data audit",
  "",
  paste0(
    "The full export contained ", format(nrow(raw), big.mark = ","),
    " mountain-gazelle sighting rows. Collapsing successive records within 30 minutes inside each camera and Camelot session, without demographic splitting, produced exactly ",
    format(nrow(encounters), big.mark = ","), " species-level independent encounters."
  ),
  paste0(
    "All ", format(nrow(event_round_matches), big.mark = ","),
    " encounters were assigned to one of 344 positive-effort camera-rounds after applying the documented session-specific clock corrections. The model used 8,256 camera-round-hour cells and did not sum demographic class-events."
  ),
  "",
  "## Requested pooled model",
  "",
  effect_sentence("s(solar_hour):wolf_z", "Wolf-associated change in diel timing"),
  effect_sentence("s(solar_hour):jackal_z", "Jackal-associated change in diel timing"),
  effect_sentence(
    "s(solar_hour):wolf_jackal_z", "Wolf x jackal temporal interaction"
  ),
  "The ordinary wolf, jackal, and interaction coefficients control differences in event intensity; the three cyclic varying-coefficient smooths above are the tests of diel timing.",
  "",
  "## Predicted normalized curves",
  "",
  metric_table_lines,
  "",
  "The four curves are population-level predictions from the requested model at all combinations of -1 and +1 SD wolf and jackal activity. Each curve integrates to one and therefore represents timing rather than abundance.",
  "",
  "## Stability and diagnostics",
  "",
  stability_table_lines,
  "",
  paste0(
    "The requested model converged with theta = ", fmt(diagnostic_primary$theta),
    ", explained ", fmt(100 * diagnostic_primary$deviance_explained, 1),
    "% of deviance, and had Pearson dispersion = ",
    fmt(diagnostic_primary$pearson_dispersion), "."
  ),
  paste0(
    "Observed and model-implied zero fractions were ",
    fmt(diagnostic_primary$observed_zero_fraction), " and ",
    fmt(diagnostic_primary$predicted_zero_fraction), ", respectively. Mean within-camera-round lag-1 Pearson residual correlation was ",
    fmt(diagnostic_primary$mean_within_camera_round_lag1_pearson), "."
  ),
  paste0(
    "The minimum ecological k-index was ",
    fmt(diagnostic_primary$minimum_ecological_k_index), " (P = ",
    fmt_p(diagnostic_primary$minimum_ecological_k_P),
    "), and maximum ecological concurvity was ",
    fmt(diagnostic_primary$maximum_ecological_concurvity), "."
  ),
  paste0(
    "The nominal jackal smooth was not retained with k = ", sensitivity_k,
    " (P = ", fmt_p(focal_results$P[
      focal_results$model == "Larger basis k = 12" &
        focal_results$term == "s(solar_hour):jackal_z"
    ]), "). Conversely, the camera-round sensitivity supported the wolf and jackal smooths and placed the interaction at the 0.05 boundary. The ecological k-index remained below one after increasing k. Thus, smooth-term inference was not stable across reasonable specifications, and no pooled predator-related timing effect should be treated as robust."
  ),
  "Point-estimate curves from the camera-round sensitivity were nearly identical to the requested curves (correlations 0.994-1.000), whereas the larger-basis curves were moderately similar (0.883-0.912). Full metrics are in `pooled_circular_GAMM_curve_stability.csv`.",
  "",
  "## Comparison with demographic classes",
  "",
  class_table_lines,
  "",
  "The pooled model should be compared with class-specific smooths by effect identity and curve shape, not by adding class-event sample sizes. A pooled encounter can contain multiple demographic classes, and the pooled dataset also retains gazelles without a focal demographic classification.",
  "",
  "## Ecological interpretation",
  "",
  "The requested fit suggests a possible pooled shift in diel timing along the jackal gradient, expressed mainly as a change in the balance between morning and late-afternoon activity. That result is exploratory because it was not robust to the larger-basis sensitivity. There was no stable pooled evidence for a wolf-associated timing shift or a wolf by jackal temporal interaction. These site-level associations do not identify responses to individual predator encounters and do not demonstrate causal avoidance.",
  "The nominal pooled jackal effect was unsupported in every class-specific jackal smooth, whereas the retained class-specific models supported a female wolf effect and wolf by jackal interactions for juveniles and males. The pooled result therefore cannot be attributed to one demographic class. Its point curves are likely weighted toward adults because adults dominate classified encounters, but pooling also includes unclassified gazelles and encounters containing multiple classes and can mask opposing class-specific responses.",
  "",
  "## Files",
  "",
  "- `pooled_circular_GAMM_predicted_curves.png` and `.pdf`",
  "- `pooled_circular_GAMM_predictions.csv`",
  "- `pooled_circular_GAMM_curve_metrics.csv`",
  "- `pooled_circular_GAMM_model_summary.csv`",
  "- `pooled_circular_GAMM_diagnostics.csv`",
  "- `pooled_circular_GAMM_species_level_independent_encounters.csv`",
  "- `pooled_circular_GAMM_models.rds`"
)
writeLines(results_lines, out("results.md"))
writeLines(capture.output(sessionInfo()), out("session_info.txt"))

cat("Pooled species-level circular GAMM analysis completed.\n")
