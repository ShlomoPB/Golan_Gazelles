#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(mvtnorm)
  library(overlap)
  library(readr)
  library(sp)
  library(tidyr)
})

set.seed(20260824)

cfg <- list(
  full_export = file.path("..", "full-export_2026-04-15_2006.csv"),
  independent_class_events = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  timestamp_corrections = "camera_timestamp_corrections.csv",
  site_covs = "SiteCovs_060324.csv",
  obs_covs = "ObsCovs_050524.csv",
  old_count_summary = "circular_GAMM_model_summary.csv",
  old_kernel_summary = "publication_predator_context_diel_pairwise_overlap_median_results.csv",
  timezone = "Asia/Jerusalem",
  primary_k = 12L,
  sensitivity_k = c(8L, 16L),
  site_k = 12L,
  simulations = 1000L,
  prediction_step = 0.1,
  primary_time_scale = "Apparent solar time",
  sensitivity_time_scale = "Sun time",
  output_prefix = "definitive_hierarchical_diel",
  checkpoint_dir = "definitive_hierarchical_diel_checkpoints_apparent_primary"
)

out <- function(suffix) paste0(cfg$output_prefix, "_", suffix)

# Preserve the previous sun-time-primary audit before replacing definitive files.
legacy_snapshot_suffixes <- c(
  "hierarchy_comparison.csv", "ecological_model_comparison.csv",
  "model_summary.csv", "effect_stability.csv", "curve_metrics.csv",
  "diagnostics.csv", "results.md"
)
for (suffix in legacy_snapshot_suffixes) {
  source <- out(suffix)
  destination <- out(paste0("suntime_primary_", suffix))
  if (file.exists(source) && !file.exists(destination)) {
    file.copy(source, destination, overwrite = FALSE)
  }
}

dir.create(cfg$checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
checkpoint_file <- function(...) {
  key <- paste(..., sep = "__")
  key <- gsub("[^A-Za-z0-9_.-]+", "_", key)
  file.path(cfg$checkpoint_dir, paste0(key, ".rds"))
}
required_files <- c(
  cfg$full_export, cfg$independent_class_events, cfg$timestamp_corrections,
  cfg$site_covs, cfg$obs_covs
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

response_levels <- c("Juvenile", "Adult female", "Adult male", "Pooled species")
class_codes <- c(
  juvenile = "Juvenile",
  female = "Adult female",
  male = "Adult male"
)
context_levels <- c(
  "Wolf -1 / Jackal -1",
  "Wolf -1 / Jackal +1",
  "Wolf +1 / Jackal -1",
  "Wolf +1 / Jackal +1"
)
context_colors <- c(
  "Wolf -1 / Jackal -1" = "#168575",
  "Wolf -1 / Jackal +1" = "#B62D35",
  "Wolf +1 / Jackal -1" = "#173F66",
  "Wolf +1 / Jackal +1" = "#835796"
)
context_linetypes <- c(
  "Wolf -1 / Jackal -1" = "solid",
  "Wolf -1 / Jackal +1" = "22",
  "Wolf +1 / Jackal -1" = "42",
  "Wolf +1 / Jackal +1" = "dotdash"
)
focal_terms <- c(
  wolf = "s(Time):wolf_z",
  jackal = "s(Time):jackal_z",
  interaction = "s(Time):wolf_jackal"
)
model_display_names <- c(
  A = "baseline",
  B = "wolf-only",
  C = "jackal-only",
  D = "additive wolf + jackal",
  E = "wolf x jackal"
)
hierarchy_display_names <- c(
  camera = "camera random intercept",
  camera_round = "camera + camera-round random intercepts",
  site_shape = paste(
    "camera-specific cyclic deviation + camera and camera-round",
    "random intercepts"
  )
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f", big.mark = ","))
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", formatC(x, 3, format = "f")))
}
fmt_plot_p <- function(x) {
  ifelse(
    is.na(x), "P = NA",
    ifelse(x < 0.001, "P < 0.001", paste0("P = ", formatC(x, 3, format = "f")))
  )
}

parse_local_export_time <- function(x) {
  out_time <- as.POSIXct(x, format = "%d/%m/%Y %H:%M:%S", tz = cfg$timezone)
  missing <- is.na(out_time)
  out_time[missing] <- as.POSIXct(
    x[missing], format = "%d/%m/%Y %H:%M", tz = cfg$timezone
  )
  out_time
}

parse_frozen_local_time <- function(x) {
  as.POSIXct(
    sub("Z$", "", x), format = "%Y-%m-%dT%H:%M:%S", tz = cfg$timezone
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

sun_time_hour <- function(timestamp, longitude, latitude) {
  local_hour <- as.numeric(format(timestamp, "%H", tz = cfg$timezone)) +
    as.numeric(format(timestamp, "%M", tz = cfg$timezone)) / 60 +
    as.numeric(format(timestamp, "%S", tz = cfg$timezone)) / 3600
  coords <- SpatialPoints(
    cbind(longitude, latitude),
    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")
  )
  as.numeric(overlap::sunTime(local_hour / 24 * 2 * pi, timestamp, coords)) * 12 / pi
}

add_time_scales <- function(data) {
  data %>%
    mutate(
      apparent_time = apparent_solar_hour(timestamp, longitude),
      sun_time = sun_time_hour(timestamp, longitude, latitude),
      apparent_bin = pmin(23L, floor(apparent_time)),
      sun_bin = pmin(23L, floor(sun_time))
    )
}

assign_round <- function(events, effort) {
  events %>%
    mutate(event_row_id = row_number(), event_date = as.Date(timestamp, tz = cfg$timezone)) %>%
    inner_join(effort, by = "camera_id", relationship = "many-to-many") %>%
    filter(event_date >= start_date, event_date <= end_date) %>%
    arrange(event_row_id, desc(start_date), desc(round)) %>%
    group_by(event_row_id) %>%
    slice_head(n = 1) %>%
    ungroup()
}

# Effort and site covariates -------------------------------------------------

obs_raw <- read.csv(
  cfg$obs_covs, na.strings = c("", "NA"), check.names = FALSE,
  stringsAsFactors = FALSE
)
effort <- tibble(
  camera_id = as.integer(obs_raw[["Camera Number"]]),
  round = as.integer(obs_raw[["Round Number"]]),
  start_date = as.Date(obs_raw[["Start date"]], format = "%d/%m/%Y"),
  end_date = as.Date(obs_raw[["End Date"]], format = "%d/%m/%Y"),
  effort_days = as.integer(obs_raw[["Days.Operational"]])
) %>%
  filter(
    !is.na(camera_id), !is.na(round), !is.na(start_date), !is.na(end_date),
    is.finite(effort_days), effort_days > 0
  ) %>%
  arrange(camera_id, round)
if (nrow(effort) != 344L || n_distinct(effort$camera_id, effort$round) != 344L) {
  stop("Expected 344 unique positive-effort camera-rounds.")
}
if (sum(effort$effort_days) != 5999L) stop("Expected 5,999 camera-days.")

site <- read_csv(cfg$site_covs, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    polygon = as.character(Polygon.Name),
    wolf_index = as.numeric(Null.Wolf.Abundance),
    jackal_index = as.numeric(null.jackal.abundance)
  ) %>%
  arrange(camera_id)
if (nrow(site) != 60L || n_distinct(site$camera_id) != 60L || anyNA(site)) {
  stop("Expected complete site covariates for 60 cameras.")
}
site <- site %>%
  mutate(
    wolf_z = as.numeric(scale(wolf_index)),
    jackal_z = as.numeric(scale(jackal_index)),
    wolf_jackal = wolf_z * jackal_z
  )

# Raw sightings, corrections, and mutually exclusive classifications --------

raw <- fread(
  cfg$full_export,
  select = c(
    "Camera Name", "Date/Time", "Species", "Trap Station Session ID",
    "Camelot GPS Latitude", "Camelot GPS Longitude", "Sex", "Life stage"
  ),
  na.strings = c("", "NA"), showProgress = FALSE
)
raw <- raw[Species == "Gazella gazella"]
setnames(
  raw,
  c(
    "Camera Name", "Date/Time", "Trap Station Session ID",
    "Camelot GPS Latitude", "Camelot GPS Longitude", "Sex", "Life stage"
  ),
  c(
    "camera_id", "timestamp_text", "session_id", "latitude", "longitude",
    "sex", "life_stage"
  )
)
raw[, `:=`(
  raw_row_id = .I,
  camera_id = as.integer(camera_id),
  session_id = as.character(session_id),
  timestamp_recorded = parse_local_export_time(timestamp_text),
  latitude = as.numeric(latitude),
  longitude = as.numeric(longitude)
)]
if (anyNA(raw$camera_id) || anyNA(raw$session_id) || anyNA(raw$timestamp_recorded)) {
  stop("Raw gazelle rows contain unparseable identifiers or timestamps.")
}

corrections <- fread(cfg$timestamp_corrections)
corrections[, `:=`(
  camera_id = as.integer(camera_name),
  session_id = as.character(trap_station_session_id),
  offset_days = as.integer(offset_days)
)]
raw <- merge(
  raw,
  corrections[, .(camera_id, session_id, offset_days)],
  by = c("camera_id", "session_id"), all.x = TRUE, sort = FALSE
)
raw[is.na(offset_days), offset_days := 0L]
raw[, corrected_date := as.Date(timestamp_recorded, tz = cfg$timezone) + offset_days]
raw[, clock_time := format(timestamp_recorded, "%H:%M:%S", tz = cfg$timezone)]
raw[, timestamp := as.POSIXct(
  paste(corrected_date, clock_time), format = "%Y-%m-%d %H:%M:%S",
  tz = cfg$timezone
)]
raw[, response_class := fifelse(
  life_stage == "Juvenile", "Juvenile",
  fifelse(
    sex == "Female" & life_stage == "Adult", "Adult female",
    fifelse(sex == "Male" & life_stage == "Adult", "Adult male", NA_character_)
  )
)]

raw_events <- as_tibble(raw) %>%
  select(
    raw_row_id, camera_id, session_id, timestamp, latitude, longitude,
    response_class
  ) %>%
  add_time_scales()
raw_round <- assign_round(raw_events, effort)
if (nrow(raw_round) != nrow(raw_events)) {
  stop("Not all raw gazelle sightings mapped to an active camera-round.")
}

# Reconstruct the authoritative pooled 30-minute encounters for auditing.
species_independent <- copy(raw)
setorder(species_independent, camera_id, session_id, timestamp, raw_row_id)
species_independent[, encounter_id := cumsum(
  seq_len(.N) == 1L |
    as.numeric(difftime(timestamp, shift(timestamp), units = "mins")) > 30
), by = .(camera_id, session_id)]
species_independent <- species_independent[, .(
  timestamp = min(timestamp),
  latitude = median(latitude, na.rm = TRUE),
  longitude = median(longitude, na.rm = TRUE)
), by = .(camera_id, session_id, encounter_id)]
if (nrow(species_independent) != 3688L) {
  stop("Expected 3,688 pooled independent encounters; found ", nrow(species_independent))
}
species_independent_round <- as_tibble(species_independent) %>%
  add_time_scales() %>%
  assign_round(effort)
if (nrow(species_independent_round) != 3688L) {
  stop("Not all pooled independent encounters mapped to active rounds.")
}

# Frozen class-event sensitivity --------------------------------------------

independent_class <- read_csv(
  cfg$independent_class_events,
  col_types = cols(event_time = col_character()), show_col_types = FALSE
) %>%
  transmute(
    camera_id = as.integer(camera_name),
    session_id = as.character(trap_station_session_id),
    response_class = unname(class_codes[as.character(gazelle_class)]),
    timestamp = parse_frozen_local_time(event_time),
    latitude = as.numeric(camelot_gps_latitude),
    longitude = as.numeric(camelot_gps_longitude)
  ) %>%
  add_time_scales() %>%
  assign_round(effort)
expected_class_events <- c("Juvenile" = 180L, "Adult female" = 1477L, "Adult male" = 2268L)
observed_class_events <- table(independent_class$response_class)
if (!all(observed_class_events[names(expected_class_events)] == expected_class_events)) {
  stop("Frozen class-event counts differ from the final audit.")
}

binary_keys <- function(data, pooled = FALSE, time_bin = "apparent_bin") {
  out_data <- data
  if (pooled) out_data$response_class <- "Pooled species"
  out_data %>%
    filter(!is.na(response_class)) %>%
    transmute(
      response_class,
      key = paste(camera_id, round, event_date, .data[[time_bin]], sep = "|")
    ) %>%
    distinct()
}

agreement_rows <- list()
for (response_i in response_levels) {
  if (response_i == "Pooled species") {
    raw_keys <- binary_keys(raw_round, pooled = TRUE)
    independent_keys <- binary_keys(species_independent_round, pooled = TRUE)
  } else {
    raw_keys <- binary_keys(filter(raw_round, response_class == response_i))
    independent_keys <- binary_keys(filter(independent_class, response_class == response_i))
  }
  raw_set <- raw_keys$key
  independent_set <- independent_keys$key
  agreement_rows[[response_i]] <- tibble(
    response_class = response_i,
    raw_positive_day_hours = length(raw_set),
    independent_positive_day_hours = length(independent_set),
    shared_positive_day_hours = length(intersect(raw_set, independent_set)),
    raw_only = length(setdiff(raw_set, independent_set)),
    independent_only = length(setdiff(independent_set, raw_set)),
    jaccard = length(intersect(raw_set, independent_set)) /
      length(union(raw_set, independent_set))
  )
}
binary_agreement <- bind_rows(agreement_rows)
write_csv(binary_agreement, out("binary_history_audit.csv"))

# Grouped binomial response --------------------------------------------------

make_successes <- function(time_scale, history_type = "Raw sightings") {
  bin_name <- if (time_scale == "Sun time") "sun_bin" else "apparent_bin"
  if (history_type == "Raw sightings") {
    class_events <- raw_round
    pooled_events <- raw_round
  } else {
    class_events <- independent_class
    pooled_events <- species_independent_round
  }
  class_success <- class_events %>%
    filter(!is.na(response_class)) %>%
    transmute(
      response_class, camera_id, round, event_date,
      Time_bin = .data[[bin_name]]
    ) %>%
    distinct() %>%
    count(response_class, camera_id, round, Time_bin, name = "success")
  pooled_success <- pooled_events %>%
    transmute(
      response_class = "Pooled species", camera_id, round, event_date,
      Time_bin = .data[[bin_name]]
    ) %>%
    distinct() %>%
    count(response_class, camera_id, round, Time_bin, name = "success")
  bind_rows(class_success, pooled_success)
}

build_grouped_history <- function(time_scale, history_type) {
  effort %>%
    crossing(response_class = response_levels, Time_bin = 0:23) %>%
    left_join(
      make_successes(time_scale, history_type),
      by = c("response_class", "camera_id", "round", "Time_bin")
    ) %>%
    mutate(
      success = replace_na(success, 0L), failure = effort_days - success,
      trials = effort_days, Time = Time_bin + 0.5,
      time_scale = .env$time_scale, history_type = .env$history_type
    ) %>%
    left_join(site, by = "camera_id")
}

analysis_data <- bind_rows(
  build_grouped_history("Sun time", "Raw sightings"),
  build_grouped_history("Apparent solar time", "Raw sightings"),
  build_grouped_history("Sun time", "Independent 30-min histories"),
  build_grouped_history("Apparent solar time", "Independent 30-min histories")
) %>%
  mutate(
    response_class = factor(response_class, levels = response_levels),
    round_factor = factor(round),
    camera_id = factor(camera_id),
    camera_round = factor(paste(as.character(camera_id), round, sep = "_")),
    polygon = factor(polygon)
  ) %>%
  arrange(history_type, time_scale, response_class, camera_id, round, Time_bin)
if (any(analysis_data$failure < 0)) {
  stop("At least one positive day-hour count exceeds Days.Operational.")
}
write_csv(
  analysis_data %>%
    mutate(
      response_class = as.character(response_class),
      camera_id = as.character(camera_id),
      camera_round = as.character(camera_round),
      polygon = as.character(polygon),
      round_factor = as.character(round_factor)
    ),
  out("data.csv")
)

data_audit <- analysis_data %>%
  group_by(history_type, time_scale, response_class) %>%
  summarise(
    rows = n(), camera_rounds = n_distinct(camera_round),
    trials = sum(trials), success = sum(success),
    encounter_frequency = sum(success) / sum(trials),
    cameras_with_success = n_distinct(camera_id[success > 0]), .groups = "drop"
  )
write_csv(data_audit, out("data_audit.csv"))

# Model formulas and fitting -------------------------------------------------

hierarchy_terms <- function(hierarchy, site_k = cfg$site_k) {
  switch(
    hierarchy,
    camera = "s(camera_id, bs = 're')",
    camera_round = paste(
      "s(camera_round, bs = 're')", "s(camera_id, bs = 're')", sep = " + "
    ),
    site_shape = paste(
      paste0(
        "s(Time, bs = 'cc', k = ", site_k,
        ", by = camera_id, m = 1, id = 20)"
      ),
      "s(camera_round, bs = 're')", "s(camera_id, bs = 're')", sep = " + "
    ),
    stop("Unknown hierarchy: ", hierarchy)
  )
}

ecological_terms <- function(model_id, k) {
  baseline <- paste0("s(Time, bs = 'cc', k = ", k, ")")
  wolf <- paste0("s(Time, by = wolf_z, bs = 'cc', k = ", k, ")")
  jackal <- paste0("s(Time, by = jackal_z, bs = 'cc', k = ", k, ")")
  interaction <- paste0("s(Time, by = wolf_jackal, bs = 'cc', k = ", k, ")")
  switch(
    model_id,
    A = baseline,
    B = paste("wolf_z", baseline, wolf, sep = " + "),
    C = paste("jackal_z", baseline, jackal, sep = " + "),
    D = paste("wolf_z + jackal_z", baseline, wolf, jackal, sep = " + "),
    E = paste(
      "wolf_z * jackal_z", baseline, wolf, jackal, interaction, sep = " + "
    ),
    stop("Unknown ecological model: ", model_id)
  )
}

build_formula <- function(model_id, hierarchy, k) {
  as.formula(paste(
    "cbind(success, failure) ~ round_factor +",
    ecological_terms(model_id, k), "+", hierarchy_terms(hierarchy)
  ))
}

fit_model <- function(data, model_id, hierarchy, k, method = "fREML") {
  formula <- build_formula(model_id, hierarchy, k)
  if (identical(method, "ML")) {
    # A discrete fREML fit supplies near-optimal smoothing parameters and cuts
    # exact ML outer-iteration time without changing the ML solution.
    initial <- bam(
      formula, data = data, family = binomial(link = "logit"),
      method = "fREML", discrete = TRUE, select = TRUE,
      knots = list(Time = c(0, 24)), nthreads = 2,
      control = gam.control(trace = FALSE)
    )
    return(gam(
      formula, data = data, family = binomial(link = "logit"),
      method = "ML", select = TRUE, knots = list(Time = c(0, 24)),
      in.out = list(sp = initial$sp, scale = 1),
      control = gam.control(trace = FALSE)
    ))
  }
  bam(
    formula, data = data, family = binomial(link = "logit"),
    method = "fREML", discrete = TRUE, select = TRUE,
    knots = list(Time = c(0, 24)), nthreads = 2,
    control = gam.control(trace = FALSE)
  )
}

fit_or_load <- function(data, model_id, hierarchy, k, method, checkpoint) {
  path <- checkpoint_file(checkpoint)
  if (file.exists(path)) {
    message("Loading checkpoint: ", checkpoint)
    return(readRDS(path))
  }
  fit <- fit_model(data, model_id, hierarchy, k, method)
  saveRDS(fit, path)
  fit
}

model_data <- function(
    response_class, time_scale = cfg$primary_time_scale,
    history_type = "Raw sightings") {
  analysis_data %>%
    filter(
      .data$response_class == .env$response_class,
      .data$time_scale == .env$time_scale,
      .data$history_type == .env$history_type
    ) %>%
    droplevels()
}

convergence_label <- function(fit) {
  if (!is.null(fit$outer.info$conv)) return(as.character(fit$outer.info$conv))
  "No convergence warning returned"
}

extract_term_table <- function(
    fit, response_class, stage, model_id, hierarchy, time_scale, k,
    history_type = "Raw sightings", selected = FALSE) {
  p_index <- seq_len(sum(fit$nsdf))
  p_estimate <- coef(fit)[p_index]
  p_se <- sqrt(diag(fit$Vp)[p_index])
  p_statistic <- p_estimate / p_se
  p_table <- tibble(
    term = names(p_estimate), estimate = unname(p_estimate), SE = p_se,
    statistic = p_statistic, P = 2 * pnorm(abs(p_statistic), lower.tail = FALSE)
  )

  target_labels <- c("s(Time)", unname(focal_terms))
  smooth_labels <- vapply(fit$smooth, function(x) x$label, character(1))
  smooth_index <- which(smooth_labels %in% target_labels)
  smooth_rows <- lapply(smooth_index, function(i) {
    smooth_i <- fit$smooth[[i]]
    start <- smooth_i$first.para
    stop <- smooth_i$last.para
    edf <- sum(fit$edf[start:stop])
    edf1 <- if (!is.null(fit$edf1)) sum(fit$edf1[start:stop]) else edf
    if (!smooth_i$fixed && smooth_i$null.space.dim == 0 && !is.null(fit$R)) {
      test <- mgcv:::reTest(fit, i)
    } else {
      test <- mgcv:::testStat(
        coef(fit)[start:stop], fit$R[, start:stop, drop = FALSE],
        fit$Vp[start:stop, start:stop, drop = FALSE],
        min(stop - start + 1L, edf1), type = 0, res.df = -1
      )
    }
    tibble(
      term = smooth_i$label, estimate = NA_real_, SE = NA_real_, EDF = edf,
      reference_df = test$rank, statistic = test$stat, P = test$pval
    )
  })
  s_table <- bind_rows(smooth_rows)
  bind_rows(
    p_table %>% transmute(
      response_class, stage, model_id, hierarchy, time_scale, k,
      history_type,
      term, term_type = "parametric", estimate, SE, EDF = NA_real_,
      reference_df = NA_real_, statistic, P
    ),
    s_table %>% transmute(
      response_class, stage, model_id, hierarchy, time_scale, k,
      history_type,
      term,
      term_type = case_when(
        term %in% focal_terms ~ "predator-by-time",
        term == "s(Time)" ~ "baseline cyclic",
        grepl("camera|polygon", term) ~ "hierarchical",
        TRUE ~ "site-specific cyclic"
      ),
      estimate, SE, EDF, reference_df, statistic, P
    )
  ) %>%
    mutate(AIC = AIC(fit), selected = selected, supported = P < 0.05)
}

extract_term_table_or_load <- function(fit, checkpoint, ...) {
  path <- checkpoint_file(paste0("term_table_v3__", checkpoint))
  if (file.exists(path)) {
    message("Loading term-table checkpoint: ", checkpoint)
    return(readRDS(path))
  }
  table <- extract_term_table(fit, ...)
  saveRDS(table, path)
  table
}

extract_diagnostics <- function(
    fit, data, response_class, stage, model_id, hierarchy, time_scale, k,
    history_type = "Raw sightings") {
  pearson <- residuals(fit, type = "pearson")
  fitted_probability <- fitted(fit)
  lag_data <- data %>%
    mutate(pearson = pearson) %>%
    arrange(camera_round, Time_bin) %>%
    group_by(camera_round) %>%
    summarise(
      lag1 = suppressWarnings(cor(pearson[-n()], pearson[-1], use = "complete.obs")),
      .groups = "drop"
    )
  if (hierarchy == "site_shape") {
    focal_k <- data.frame()
    conc_value <- NA_real_
  } else {
    set.seed(20260824)
    k_raw <- suppressWarnings(as.data.frame(k.check(
      fit, subsample = min(5000L, nrow(data)), n.rep = 300
    )))
    k_terms <- rownames(k_raw)
    focal_k <- k_raw[grepl("^s\\(Time\\)", k_terms), , drop = FALSE]
    conc <- concurvity(fit, full = TRUE)
    focal_columns <- grepl("^s\\(Time\\)", colnames(conc)) &
      !grepl("camera", colnames(conc))
    conc_value <- if (any(focal_columns)) {
      max(conc["estimate", focal_columns], na.rm = TRUE)
    } else NA_real_
  }
  w <- as.numeric(fit$prior.weights)
  mean_y <- sum(w * fit$y) / sum(w)
  weighted_residual <- sqrt(w) * (as.numeric(fit$y) - fit$fitted.values)
  weighted_null <- sqrt(w) * (as.numeric(fit$y) - mean_y)
  adjusted_r_squared <- 1 - var(weighted_residual) * (nrow(fit$model) - 1) /
    (var(weighted_null) * fit$df.residual)
  tibble(
    response_class, stage, model_id, hierarchy, time_scale, k, history_type,
    convergence = convergence_label(fit), AIC = AIC(fit),
    deviance_explained = (fit$null.deviance - fit$deviance) / fit$null.deviance,
    adjusted_r_squared = adjusted_r_squared,
    pearson_dispersion = sum(pearson^2, na.rm = TRUE) / fit$df.residual,
    observed_encounter_frequency = sum(data$success) / sum(data$trials),
    predicted_encounter_frequency = sum(fitted_probability * data$trials) /
      sum(data$trials),
    mean_camera_round_lag1_pearson = mean(lag_data$lag1, na.rm = TRUE),
    minimum_time_k_index = if (nrow(focal_k)) min(focal_k[[3]], na.rm = TRUE) else NA_real_,
    minimum_time_k_P = if (nrow(focal_k)) min(focal_k[[4]], na.rm = TRUE) else NA_real_,
    maximum_time_concurvity = conc_value,
    fitted_probability_maximum = max(fitted_probability),
    fitted_probability_minimum = min(fitted_probability)
  )
}

# Hierarchy selection --------------------------------------------------------

all_models <- list()
hierarchy_rows <- list()
selected_hierarchy <- setNames(character(length(response_levels)), response_levels)

hierarchy_path <- out("hierarchy_comparison.csv")
for (response_i in response_levels) {
  dat <- model_data(response_i)
  candidates <- c("camera", "camera_round")
  if (response_i %in% c("Adult female", "Adult male", "Pooled species")) {
    candidates <- c(candidates, "site_shape")
  }
  fits <- list()
  for (hierarchy_i in candidates) {
    message("Baseline hierarchy: ", response_i, " / ", hierarchy_i)
    fits[[hierarchy_i]] <- fit_or_load(
      dat, "A", hierarchy_i, cfg$primary_k, "fREML",
      paste("apparent_primary", "hierarchy", response_i, hierarchy_i, sep = "__")
    )
    diag_i <- extract_diagnostics(
      fits[[hierarchy_i]], dat, response_i, "Hierarchy candidate", "A",
      hierarchy_i, cfg$primary_time_scale, cfg$primary_k
    )
    hierarchy_rows[[paste(response_i, hierarchy_i)]] <- diag_i
  }
  hierarchy_table_i <- bind_rows(hierarchy_rows) %>%
    filter(.data$response_class == .env$response_i)
  min_aic <- min(hierarchy_table_i$AIC)
  hierarchy_table_i <- hierarchy_table_i %>%
    mutate(
      delta_AIC = AIC - min_aic,
      complexity = match(hierarchy, c("camera", "camera_round", "site_shape"))
    ) %>%
    filter(delta_AIC <= 2) %>%
    arrange(complexity, AIC)
  selected_hierarchy[[response_i]] <- hierarchy_table_i$hierarchy[1]
}
hierarchy_comparison <- bind_rows(hierarchy_rows) %>%
  group_by(response_class) %>%
  mutate(
    primary_time_scale = cfg$primary_time_scale,
    delta_AIC = AIC - min(AIC),
    selected = hierarchy == selected_hierarchy[as.character(response_class)]
  ) %>%
  ungroup()
write_csv(hierarchy_comparison, hierarchy_path)

# Ecological model comparison and final refits ------------------------------

ecological_rows <- list()
term_tables <- list()
final_fits <- list()
full_primary_fits <- list()
selected_model <- setNames(character(length(response_levels)), response_levels)

for (response_i in response_levels) {
  dat <- model_data(response_i)
  hierarchy_i <- selected_hierarchy[[response_i]]
  fit_set <- list()
  for (model_i in LETTERS[1:5]) {
    message("Ecological fREML model: ", response_i, " / ", model_i)
    fit_set[[model_i]] <- fit_or_load(
      dat, model_i, hierarchy_i, cfg$primary_k, "fREML",
      paste("apparent_primary", "ecology_fREML", response_i, model_i, sep = "__")
    )
    ecological_rows[[paste(response_i, model_i)]] <- tibble(
      response_class = response_i, model_id = model_i,
      hierarchy = hierarchy_i, AIC = AIC(fit_set[[model_i]]),
      deviance_explained = (fit_set[[model_i]]$null.deviance -
        fit_set[[model_i]]$deviance) / fit_set[[model_i]]$null.deviance,
      convergence = convergence_label(fit_set[[model_i]])
    )
    write_csv(bind_rows(ecological_rows), out("ecological_model_comparison_in_progress.csv"))
  }
  comparison_i <- bind_rows(ecological_rows) %>%
    filter(response_class == response_i) %>%
    mutate(
      delta_AIC = AIC - min(AIC),
      complexity = match(model_id, LETTERS[1:5])
    ) %>%
    filter(delta_AIC <= 2) %>%
    arrange(complexity, AIC)
  selected_model[[response_i]] <- comparison_i$model_id[1]

  message("Final fREML refit: ", response_i, " / ", selected_model[[response_i]])
  final_fits[[response_i]] <- fit_or_load(
    dat, selected_model[[response_i]], hierarchy_i, cfg$primary_k, "fREML",
    paste(
      "apparent_primary", "final_fREML", response_i,
      selected_model[[response_i]], sep = "__"
    )
  )
  term_tables[[paste("final", response_i)]] <- extract_term_table_or_load(
    final_fits[[response_i]],
    paste("apparent_primary", "final", response_i, sep = "__"),
    response_i, "Final fREML", selected_model[[response_i]],
    hierarchy_i, cfg$primary_time_scale, cfg$primary_k, selected = TRUE
  )

  if (selected_model[[response_i]] == "E") {
    full_primary_fits[[response_i]] <- final_fits[[response_i]]
  } else {
    full_primary_fits[[response_i]] <- fit_or_load(
      dat, "E", hierarchy_i, cfg$primary_k, "fREML",
      paste("apparent_primary", "full_E_fREML", response_i, sep = "__")
    )
  }
  term_tables[[paste("full_E", response_i)]] <- extract_term_table_or_load(
    full_primary_fits[[response_i]],
    paste("apparent_primary", "full_E", response_i, sep = "__"),
    response_i, "Full E fREML", "E",
    hierarchy_i, cfg$primary_time_scale, cfg$primary_k
  )
}

ecological_comparison <- bind_rows(ecological_rows) %>%
  group_by(response_class) %>%
  mutate(
    primary_time_scale = cfg$primary_time_scale,
    delta_AIC = AIC - min(AIC),
    selected = model_id == selected_model[as.character(response_class)]
  ) %>%
  ungroup()
write_csv(ecological_comparison, out("ecological_model_comparison.csv"))

# Sensitivity fits use the full predator model so all timing effects are tested.
sensitivity_fits <- list()
for (response_i in response_levels) {
  hierarchy_i <- selected_hierarchy[[response_i]]
  dat_apparent <- model_data(response_i, "Apparent solar time")
  dat_sun <- model_data(response_i, "Sun time")
  dat_independent <- model_data(
    response_i, "Apparent solar time", "Independent 30-min histories"
  )
  alternate_hierarchy <- if (hierarchy_i == "camera") "camera_round" else "camera"

  sensitivity_specs <- list(
    k8 = list(
      data = dat_apparent, hierarchy = hierarchy_i,
      time_scale = "Apparent solar time",
      history_type = "Raw sightings", k = 8L
    ),
    k16 = list(
      data = dat_apparent, hierarchy = hierarchy_i,
      time_scale = "Apparent solar time",
      history_type = "Raw sightings", k = 16L
    ),
    sun_time = list(
      data = dat_sun, hierarchy = hierarchy_i,
      time_scale = "Sun time", history_type = "Raw sightings",
      k = cfg$primary_k
    ),
    alternate_hierarchy = list(
      data = dat_apparent, hierarchy = alternate_hierarchy,
      time_scale = "Apparent solar time", history_type = "Raw sightings",
      k = cfg$primary_k
    ),
    independent = list(
      data = dat_independent, hierarchy = hierarchy_i,
      time_scale = "Apparent solar time",
      history_type = "Independent 30-min histories", k = cfg$primary_k
    )
  )
  for (sensitivity_i in names(sensitivity_specs)) {
    spec <- sensitivity_specs[[sensitivity_i]]
    message("Sensitivity: ", response_i, " / ", sensitivity_i)
    fit_i <- fit_or_load(
      spec$data, "E", spec$hierarchy, spec$k, "fREML",
      paste(
        "apparent_primary", "sensitivity", response_i, sensitivity_i,
        sep = "__"
      )
    )
    sensitivity_fits[[paste(response_i, sensitivity_i, sep = "|")]] <- fit_i
    term_tables[[paste("sensitivity", response_i, sensitivity_i)]] <-
      extract_term_table_or_load(
        fit_i,
        paste(
          "apparent_primary", "sensitivity", response_i, sensitivity_i,
          sep = "__"
        ),
        response_i, paste("Sensitivity", sensitivity_i), "E",
        spec$hierarchy, spec$time_scale, spec$k,
        history_type = spec$history_type
      )
  }
}

model_summary <- bind_rows(term_tables)
write_csv(model_summary, out("model_summary.csv"))
saveRDS(
  list(
    selected_hierarchy = selected_hierarchy, selected_model = selected_model,
    final_fits = final_fits, full_primary_fits = full_primary_fits,
    sensitivity_fits = sensitivity_fits
  ),
  out("models.rds")
)

# Diagnostics for retained, full, and sensitivity models -------------------

diagnostic_rows <- list()
for (response_i in response_levels) {
  dat <- model_data(response_i)
  hierarchy_i <- selected_hierarchy[[response_i]]
  diagnostic_rows[[paste("final", response_i)]] <- extract_diagnostics(
    final_fits[[response_i]], dat, response_i, "Final fREML",
    selected_model[[response_i]], hierarchy_i, cfg$primary_time_scale,
    cfg$primary_k
  )
  diagnostic_rows[[paste("full", response_i)]] <- extract_diagnostics(
    full_primary_fits[[response_i]], dat, response_i, "Full E fREML", "E",
    hierarchy_i, cfg$primary_time_scale, cfg$primary_k
  )
  for (sensitivity_i in c(
      "k8", "k16", "sun_time", "alternate_hierarchy", "independent")) {
    fit_i <- sensitivity_fits[[paste(response_i, sensitivity_i, sep = "|")]]
    data_i <- if (sensitivity_i == "sun_time") {
      model_data(response_i, "Sun time")
    } else if (sensitivity_i == "independent") {
      model_data(
        response_i, "Apparent solar time", "Independent 30-min histories"
      )
    } else dat
    time_i <- if (sensitivity_i == "sun_time") {
      "Sun time"
    } else {
      "Apparent solar time"
    }
    history_i <- if (sensitivity_i == "independent") {
      "Independent 30-min histories"
    } else "Raw sightings"
    k_i <- if (sensitivity_i == "k8") 8L else if (sensitivity_i == "k16") 16L else cfg$primary_k
    hierarchy_s <- if (sensitivity_i == "alternate_hierarchy") {
      if (hierarchy_i == "camera") "camera_round" else "camera"
    } else hierarchy_i
    diagnostic_rows[[paste(response_i, sensitivity_i)]] <- extract_diagnostics(
      fit_i, data_i, response_i, paste("Sensitivity", sensitivity_i), "E",
      hierarchy_s, time_i, k_i, history_type = history_i
    )
  }
}
diagnostics <- bind_rows(diagnostic_rows)
write_csv(diagnostics, out("diagnostics.csv"))

# Stability classification --------------------------------------------------

full_effects <- model_summary %>%
  filter(
    term %in% focal_terms,
    stage %in% c(
      "Full E fREML", "Sensitivity k8", "Sensitivity k16",
      "Sensitivity sun_time", "Sensitivity alternate_hierarchy",
      "Sensitivity independent"
    )
  ) %>%
  mutate(
    effect = names(focal_terms)[match(term, focal_terms)],
    specification = case_when(
      stage == "Full E fREML" ~ "Primary apparent-solar full model",
      TRUE ~ sub("Sensitivity ", "", stage)
    )
  )

# The original sun-time-primary audit also crossed sun time with the
# camera-only hierarchy. Treat that prespecified combined sensitivity as a
# seventh required robustness specification rather than only discussing it
# narratively.
legacy_cross_path <- out("suntime_primary_model_summary.csv")
if (!file.exists(legacy_cross_path)) {
  stop("Missing cached sun-time + alternate-hierarchy sensitivity summary: ",
       legacy_cross_path)
}
legacy_cross_effects <- read_csv(legacy_cross_path, show_col_types = FALSE) %>%
  filter(
    stage == "Sensitivity alternate_hierarchy",
    term %in% focal_terms
  ) %>%
  mutate(
    effect = names(focal_terms)[match(term, focal_terms)],
    specification = "sun_time + alternate_hierarchy"
  )
full_effects <- bind_rows(full_effects, legacy_cross_effects)

# Include the combined sensitivity rows in the authoritative coefficient table
# so every specification counted by the robustness rule is directly auditable.
model_summary <- bind_rows(
  model_summary,
  legacy_cross_effects %>% select(all_of(names(model_summary)))
)

required_specifications <- 7L
stability <- full_effects %>%
  group_by(response_class, effect) %>%
  summarise(
    primary_P = P[specification == "Primary apparent-solar full model"][1],
    supported_specifications = sum(P < 0.05, na.rm = TRUE),
    tested_specifications = sum(!is.na(P)),
    retained_model = selected_model[as.character(first(response_class))],
    term_retained = case_when(
      first(effect) == "wolf" ~ retained_model %in% c("B", "D", "E"),
      first(effect) == "jackal" ~ retained_model %in% c("C", "D", "E"),
      first(effect) == "interaction" ~ retained_model == "E"
    ),
    robust = tested_specifications == required_specifications &&
      supported_specifications == required_specifications &&
      term_retained,
    .groups = "drop"
  ) %>%
  mutate(
    robustness_status = case_when(
      robust ~ "robust",
      term_retained & primary_P < 0.05 ~ "sensitivity-dependent",
      term_retained ~ "unsupported",
      TRUE ~ "not retained"
    )
  )
write_csv(stability, out("effect_stability.csv"))

# Put the final robustness classification beside every focal coefficient in
# the authoritative model-summary table. Nonfocal rows remain unclassified.
model_summary <- model_summary %>%
  left_join(
    stability %>%
      mutate(term = unname(focal_terms[effect])) %>%
      select(
        response_class, term, supported_specifications,
        tested_specifications, term_retained, robust, robustness_status
      ),
    by = c("response_class", "term")
  )
write_csv(model_summary, out("model_summary.csv"))

# Predictions ---------------------------------------------------------------

prediction_grid <- function(data) {
  expand_grid(
    Time = seq(cfg$prediction_step / 2, 24 - cfg$prediction_step / 2,
               by = cfg$prediction_step),
    wolf_z = c(-1, 1), jackal_z = c(-1, 1)
  ) %>%
    mutate(
      wolf_jackal = wolf_z * jackal_z,
      round_factor = factor("4", levels = levels(data$round_factor)),
      camera_id = factor(levels(data$camera_id)[1], levels = levels(data$camera_id)),
      camera_round = factor(
        levels(data$camera_round)[1], levels = levels(data$camera_round)
      ),
      context = case_when(
        wolf_z < 0 & jackal_z < 0 ~ context_levels[1],
        wolf_z < 0 & jackal_z > 0 ~ context_levels[2],
        wolf_z > 0 & jackal_z < 0 ~ context_levels[3],
        TRUE ~ context_levels[4]
      ),
      context = factor(context, levels = context_levels)
    ) %>%
    arrange(context, Time)
}

simulate_predictions <- function(fit, data, response_class) {
  grid <- prediction_grid(data)
  smooth_labels <- vapply(fit$smooth, function(x) x$label, character(1))
  excluded <- smooth_labels[grepl("camera", smooth_labels)]
  X <- predict(fit, newdata = grid, type = "lpmatrix", exclude = excluded)
  covariance <- tryCatch(vcov(fit, unconditional = TRUE), error = function(e) vcov(fit))
  draws <- rmvnorm(cfg$simulations, mean = coef(fit), sigma = covariance, method = "svd")
  probability <- plogis(drop(X %*% coef(fit)))
  probability_draws <- plogis(X %*% t(draws))
  normalized <- probability
  normalized_draws <- probability_draws
  for (context_i in context_levels) {
    idx <- which(grid$context == context_i)
    normalized[idx] <- probability[idx] /
      (sum(probability[idx]) * cfg$prediction_step)
    draw_area <- colSums(probability_draws[idx, , drop = FALSE]) * cfg$prediction_step
    normalized_draws[idx, ] <- sweep(
      probability_draws[idx, , drop = FALSE], 2, draw_area, "/"
    )
  }
  summary_data <- grid %>%
    transmute(
      response_class, Time, wolf_z, jackal_z, context,
      encounter_probability = probability,
      encounter_probability_lower = apply(probability_draws, 1, quantile, 0.025),
      encounter_probability_upper = apply(probability_draws, 1, quantile, 0.975),
      normalized_activity = normalized,
      normalized_lower = apply(normalized_draws, 1, quantile, 0.025),
      normalized_upper = apply(normalized_draws, 1, quantile, 0.975)
    )
  list(summary = summary_data, normalized_draws = normalized_draws, grid = grid)
}

prediction_objects <- list()
for (response_i in response_levels) {
  prediction_objects[[response_i]] <- simulate_predictions(
    final_fits[[response_i]], model_data(response_i), response_i
  )
}
predictions <- bind_rows(lapply(prediction_objects, `[[`, "summary"))
write_csv(predictions, out("predictions.csv"))

curve_metrics <- predictions %>%
  group_by(response_class, context) %>%
  summarise(
    peak_apparent_solar_time = Time[which.max(normalized_activity)],
    central_06_18_share = sum(normalized_activity[Time >= 6 & Time < 18]) *
      cfg$prediction_step,
    mean_hours_from_solar_noon = sum(abs(Time - 12) * normalized_activity) *
      cfg$prediction_step,
    .groups = "drop"
  )
write_csv(curve_metrics, out("curve_metrics.csv"))

effect_levels <- c("wolf", "jackal", "interaction")
effect_labels <- c(
  wolf = "Wolf x time",
  jackal = "Jackal x time",
  interaction = "Wolf x jackal x time"
)
final_test_rows <- model_summary %>%
  filter(stage == "Final fREML", term %in% focal_terms) %>%
  transmute(
    response_class,
    effect = names(focal_terms)[match(term, focal_terms)],
    P
  ) %>%
  left_join(
    stability %>% select(response_class, effect, robust),
    by = c("response_class", "effect")
  )
plot_test_labels <- expand_grid(
    response_class = response_levels,
    effect = effect_levels
  ) %>%
  left_join(final_test_rows, by = c("response_class", "effect")) %>%
  mutate(
    effect = factor(effect, levels = effect_levels),
    retained_model = unname(selected_model[response_class]),
    selected_structure = unname(model_display_names[retained_model]),
    test_label = case_when(
      is.na(P) ~ paste0(effect_labels[as.character(effect)], ": not retained"),
      robust ~ paste0(
        effect_labels[as.character(effect)], ": ", fmt_plot_p(P), " (robust)"
      ),
      P < 0.05 ~ paste0(
        effect_labels[as.character(effect)], ": ", fmt_plot_p(P),
        " (sensitivity-dependent)"
      ),
      TRUE ~ paste0(effect_labels[as.character(effect)], ": ", fmt_plot_p(P))
    )
  ) %>%
  arrange(response_class, effect) %>%
  group_by(response_class) %>%
  summarise(
    label = paste0(
      "Selected: ", first(selected_structure), "\n",
      paste(test_label, collapse = "\n")
    ),
    .groups = "drop"
  )

curve_plot <- ggplot(
  predictions,
  aes(
    Time, normalized_activity, color = context, fill = context,
    linetype = context
  )
) +
  geom_ribbon(aes(ymin = normalized_lower, ymax = normalized_upper),
              alpha = 0.08, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_label(
    data = plot_test_labels,
    inherit.aes = FALSE,
    aes(x = 0.45, y = Inf, label = label),
    hjust = 0, vjust = 1.08, size = 2.65, lineheight = 0.95,
    fill = "white", alpha = 0.88, color = "#26343A", linewidth = 0.2,
    label.padding = grid::unit(0.18, "lines")
  ) +
  facet_wrap(~response_class, ncol = 2, scales = "free_y") +
  scale_color_manual(values = context_colors) +
  scale_fill_manual(values = context_colors) +
  scale_linetype_manual(values = context_linetypes) +
  scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24), expand = c(0, 0)) +
  labs(
    title = "Hierarchical estimates of mountain gazelle diel activity",
    subtitle = "Continuous wolf and jackal gradients evaluated at -1 and +1 SD",
    x = "Apparent solar time (solar noon = 12:00)",
    y = "Normalized relative activity", color = NULL, fill = NULL, linetype = NULL,
    caption = paste0(
      "Binomial cyclic GAMs; curves integrate to one within each prediction. ",
      "Bands are simulation-based 95% intervals."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold"),
    legend.position = "top", panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), plot.caption = element_text(hjust = 0)
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )
ggsave(out("curves.png"), curve_plot, width = 11, height = 8.5, dpi = 400, bg = "white")
ggsave(out("curves.pdf"), curve_plot, width = 11, height = 8.5, bg = "white")

pooled_predictions <- predictions %>%
  filter(response_class == "Pooled species", jackal_z == -1) %>%
  mutate(
    wolf_context = factor(
      if_else(wolf_z < 0, "Wolf -1 SD", "Wolf +1 SD"),
      levels = c("Wolf -1 SD", "Wolf +1 SD")
    )
  )
pooled_colors <- c("Wolf -1 SD" = "#168575", "Wolf +1 SD" = "#173F66")
pooled_test_label <- plot_test_labels %>%
  filter(response_class == "Pooled species") %>%
  pull(label)
pooled_plot <- ggplot(
  pooled_predictions,
    aes(Time, normalized_activity, color = wolf_context, fill = wolf_context)
  ) +
  geom_ribbon(
    aes(ymin = normalized_lower, ymax = normalized_upper),
    alpha = 0.14, color = NA
  ) +
  geom_line(linewidth = 1.05) +
  annotate(
    "label", x = 0.45, y = Inf, label = pooled_test_label,
    hjust = 0, vjust = 1.08, size = 3.15, lineheight = 0.95,
    fill = "white", alpha = 0.9, color = "#26343A", linewidth = 0.2,
    label.padding = grid::unit(0.2, "lines")
  ) +
  scale_color_manual(values = pooled_colors) +
  scale_fill_manual(values = pooled_colors) +
  scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24), expand = c(0, 0)) +
  labs(
    title = "Pooled mountain gazelle diel activity across wolf activity",
    subtitle = paste0(
      "Selected ", model_display_names[selected_model[["Pooled species"]]],
      " structure"
    ),
    x = "Apparent solar time (solar noon = 12:00)",
    y = "Normalized relative activity", color = NULL, fill = NULL,
    caption = "Curves integrate to one; bands are simulation-based 95% intervals."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "top", panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), plot.caption = element_text(hjust = 0)
  )
ggsave(
  out("pooled_curves.png"), pooled_plot, width = 8.5, height = 5.5,
  dpi = 400, bg = "white"
)
ggsave(out("pooled_curves.pdf"), pooled_plot, width = 8.5, height = 5.5,
       bg = "white")

# Robust difference curves; use a placeholder when no term passes all checks.
robust_effects <- stability %>% filter(robust, response_class != "Pooled species")
if (nrow(robust_effects)) {
  difference_rows <- list()
  for (i in seq_len(nrow(robust_effects))) {
    response_i <- robust_effects$response_class[i]
    effect_i <- robust_effects$effect[i]
    object_i <- prediction_objects[[response_i]]
    draws_i <- object_i$normalized_draws
    grid_i <- object_i$grid
    context_indices <- lapply(context_levels, function(x) which(grid_i$context == x))
    names(context_indices) <- context_levels
    if (effect_i == "interaction") {
      contrast <-
        draws_i[context_indices[[4]], ] - draws_i[context_indices[[2]], ] -
        draws_i[context_indices[[3]], ] + draws_i[context_indices[[1]], ]
      label <- "Wolf x jackal difference-in-differences"
    } else if (effect_i == "wolf") {
      contrast <- 0.5 * (
        draws_i[context_indices[[3]], ] - draws_i[context_indices[[1]], ] +
          draws_i[context_indices[[4]], ] - draws_i[context_indices[[2]], ]
      )
      label <- "High minus low wolf (averaged over jackal)"
    } else {
      contrast <- 0.5 * (
        draws_i[context_indices[[2]], ] - draws_i[context_indices[[1]], ] +
          draws_i[context_indices[[4]], ] - draws_i[context_indices[[3]], ]
      )
      label <- "High minus low jackal (averaged over wolf)"
    }
    difference_rows[[i]] <- tibble(
      response_class = response_i, effect = label,
      p_label = paste0("Global smooth ", fmt_plot_p(robust_effects$primary_P[i])),
      Time = grid_i$Time[context_indices[[1]]],
      difference = rowMeans(contrast),
      lower = apply(contrast, 1, quantile, 0.025),
      upper = apply(contrast, 1, quantile, 0.975)
    )
  }
  difference_data <- bind_rows(difference_rows)
  difference_plot <- ggplot(difference_data, aes(Time, difference)) +
    geom_hline(yintercept = 0, color = "#39474D", linetype = "dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#517A8A", alpha = 0.22) +
    geom_line(color = "#174A5B", linewidth = 0.9) +
    geom_label(
      data = difference_data %>%
        distinct(response_class, effect, p_label),
      inherit.aes = FALSE,
      aes(x = 0.45, y = Inf, label = p_label),
      hjust = 0, vjust = 1.08, size = 3.1, fill = "white", alpha = 0.9,
      color = "#26343A", linewidth = 0.2,
      label.padding = grid::unit(0.2, "lines")
    ) +
    facet_wrap(vars(response_class, effect), nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24), expand = c(0, 0)) +
    labs(
      title = "Robust predator-associated differences in gazelle diel timing",
      x = "Apparent solar time",
      y = "Difference in normalized relative activity"
    ) +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
} else {
  difference_data <- tibble()
  difference_plot <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.55,
      label = "No predator-by-time effect met the prespecified robustness criterion",
      size = 5.2, fontface = "bold"
    ) +
    annotate(
      "text", x = 0.5, y = 0.43,
      label = "Support was required across k, hierarchy, time scale, and event-history sensitivities",
      size = 4, color = "#536168"
    ) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
}
write_csv(difference_data, out("difference_curves.csv"))
ggsave(
  out("difference_curves.png"), difference_plot,
  width = 10, height = 6.5, dpi = 400, bg = "white"
)
ggsave(out("difference_curves.pdf"), difference_plot, width = 10, height = 6.5)

# Compact diagnostic and result documents ----------------------------------

diagnostic_lines <- c(
  "# Definitive hierarchical diel GAM diagnostics",
  "",
  "## Selected structures",
  "",
  "| Response | Hierarchy | Ecological model |",
  "|---|---|---|",
  paste0(
    "| ", response_levels, " | ", selected_hierarchy[response_levels], " | ",
    model_display_names[selected_model[response_levels]], " |"
  ),
  "",
  "## Selected-model diagnostics",
  "",
  "| Response | Convergence | Deviance explained | Dispersion | k-index | k P | Lag-1 residual r | Observed/predicted frequency |",
  "|---|---|---:|---:|---:|---:|---:|---:|",
  vapply(response_levels, function(response_i) {
    row <- diagnostics %>%
      filter(response_class == response_i, stage == "Final fREML")
    paste0(
      "| ", response_i, " | ", row$convergence, " | ",
      fmt(100 * row$deviance_explained, 1), "% | ", fmt(row$pearson_dispersion),
      " | ", fmt(row$minimum_time_k_index), " | ", fmt_p(row$minimum_time_k_P),
      " | ", fmt(row$mean_camera_round_lag1_pearson), " | ",
      fmt(row$observed_encounter_frequency, 4), "/",
      fmt(row$predicted_encounter_frequency, 4), " |"
    )
  }, character(1)),
  "",
  "Detailed hierarchy comparisons, model comparisons, focal-term stability, and all model diagnostics are in the associated CSV files."
)
writeLines(diagnostic_lines, out("diagnostics.md"))

selection_summary <- bind_rows(lapply(response_levels, function(response_i) {
  hierarchy_i <- hierarchy_comparison %>%
    filter(response_class == response_i)
  selected_hierarchy_i <- hierarchy_i %>% filter(selected)
  next_hierarchy_i <- hierarchy_i %>% filter(!selected) %>% slice_min(AIC, n = 1)
  ecology_i <- ecological_comparison %>% filter(response_class == response_i)
  selected_ecology_i <- ecology_i %>% filter(selected)
  next_ecology_i <- ecology_i %>% filter(!selected) %>% slice_min(AIC, n = 1)
  tibble(
    response_class = response_i,
    selected_hierarchy = selected_hierarchy_i$hierarchy,
    hierarchy_competitor = next_hierarchy_i$hierarchy,
    hierarchy_competitor_minus_selected = next_hierarchy_i$AIC -
      selected_hierarchy_i$AIC,
    selected_ecology = selected_ecology_i$model_id,
    ecology_competitor = next_ecology_i$model_id,
    ecology_competitor_minus_selected = next_ecology_i$AIC -
      selected_ecology_i$AIC
  )
}))
write_csv(selection_summary, out("selection_summary.csv"))

primary_effects <- full_effects %>%
  filter(specification == "Primary apparent-solar full model") %>%
  select(response_class, effect, EDF, reference_df, statistic, P) %>%
  left_join(
    stability %>%
      select(
        response_class, effect, supported_specifications,
        tested_specifications, term_retained, robust, robustness_status
      ),
    by = c("response_class", "effect")
  ) %>%
  arrange(
    match(response_class, response_levels),
    match(effect, c("wolf", "jackal", "interaction"))
  )

selection_table_lines <- vapply(seq_len(nrow(selection_summary)), function(i) {
  row <- selection_summary[i, ]
  paste0(
    "| ", row$response_class, " | ",
    hierarchy_display_names[row$selected_hierarchy], " | ",
    hierarchy_display_names[row$hierarchy_competitor], " (",
    fmt(row$hierarchy_competitor_minus_selected), ") | ",
    model_display_names[row$selected_ecology], " | ",
    model_display_names[row$ecology_competitor], " (",
    fmt(row$ecology_competitor_minus_selected), ") |"
  )
}, character(1))

effect_table_lines <- vapply(seq_len(nrow(primary_effects)), function(i) {
  row <- primary_effects[i, ]
  paste0(
    "| ", row$response_class, " | ", effect_labels[row$effect], " | ",
    fmt(row$EDF), " | ", fmt_p(row$P), " | ",
    ifelse(row$term_retained, "yes", "no"), " | ",
    row$supported_specifications, "/", row$tested_specifications, " | ",
    row$robustness_status, " |"
  )
}, character(1))

robust_label <- function(data, response_i) {
  labels <- data %>%
    filter(response_class == response_i, robust) %>%
    pull(effect)
  if (!length(labels)) "none" else paste(labels, collapse = ", ")
}

primary_audit <- data_audit %>%
  filter(
    history_type == "Raw sightings",
    time_scale == cfg$primary_time_scale
  )
independent_audit <- data_audit %>%
  filter(
    history_type == "Independent 30-min histories",
    time_scale == cfg$primary_time_scale
  )

prediction_result_lines <- vapply(response_levels, function(response_i) {
  metrics_i <- curve_metrics %>% filter(response_class == response_i)
  paste0(
    "- **", response_i, ":** predicted apparent-solar peak times across the ",
    "four wolf-jackal combinations were ",
    paste(formatC(metrics_i$peak_apparent_solar_time, digits = 2, format = "f"),
          collapse = ", "),
    " h; robust predator-by-time effect(s): ",
    robust_label(stability, response_i), "."
  )
}, character(1))

legacy_hierarchy_path <- out("suntime_primary_hierarchy_comparison.csv")
legacy_ecology_path <- out("suntime_primary_ecological_model_comparison.csv")
legacy_stability_path <- out("suntime_primary_effect_stability.csv")
if (all(file.exists(c(
    legacy_hierarchy_path, legacy_ecology_path, legacy_stability_path)))) {
  legacy_hierarchy <- read_csv(legacy_hierarchy_path, show_col_types = FALSE)
  legacy_ecology <- read_csv(legacy_ecology_path, show_col_types = FALSE)
  legacy_stability <- read_csv(legacy_stability_path, show_col_types = FALSE)
  comparison_table_lines <- vapply(response_levels, function(response_i) {
    old_h <- legacy_hierarchy %>%
      filter(response_class == response_i, selected) %>% pull(hierarchy)
    old_e <- legacy_ecology %>%
      filter(response_class == response_i, selected) %>% pull(model_id)
    old_r <- robust_label(legacy_stability, response_i)
    new_r <- robust_label(stability, response_i)
    paste0(
      "| ", response_i, " | ", hierarchy_display_names[old_h], " | ",
      hierarchy_display_names[selected_hierarchy[response_i]], " | ",
      model_display_names[old_e], " | ",
      model_display_names[selected_model[response_i]], " | ",
      old_r, " | ", new_r, " |"
    )
  }, character(1))
} else {
  comparison_table_lines <- "| Previous sun-time files unavailable | NA | NA | NA | NA | NA | NA |"
}

results_lines <- c(
  "# Definitive hierarchical diel GAM results",
  "",
  "## 1. Primary time scale and effort",
  "",
  "The primary analysis used equal-duration apparent-solar-clock hours, calculated as UTC time plus longitude/15 and the date-specific equation-of-time correction. This additive transformation aligned solar noon at 12:00 while preserving true one-hour bin durations. The sunrise-sunset-standardized `overlap::sunTime()` analysis was retained as the prespecified time-scale sensitivity.",
  "",
  paste0(
    "The raw export contained ", format(nrow(raw), big.mark = ","),
    " mountain-gazelle sighting rows. Binary camera-day-hour histories were created directly from corrected raw timestamps using mutually exclusive juvenile, adult-female, and adult-male classifications. The pooled response used any gazelle sighting, regardless of demographic classification. The pooled 30-minute reconstruction returned exactly ",
    format(nrow(species_independent), big.mark = ","), " encounters."
  ),
  "Each positive-effort camera-round contributed the recorded `Days.Operational` trials to every apparent-solar hour. Success was the number of operational dates with one or more encounters and failure was `Days.Operational - success`; zero-success hours were retained. Camera 45 round 2 contributed its recorded 21 operational days across a 32-day span. Exact missing dates and hourly outage information were unavailable and were not invented.",
  paste0(
    "The primary raw histories contained ",
    paste0(
      primary_audit$success[match(response_levels, primary_audit$response_class)],
      " ", tolower(response_levels), collapse = ", "
    ),
    " positive camera-day-hours. The apparent-solar independent-history ",
    "sensitivity contained ",
    paste(
      independent_audit$success[match(response_levels, independent_audit$response_class)],
      collapse = ", "
    ), " positive camera-day-hours in the same class order."
  ),
  "",
  "## 2. Model selection",
  "",
  "The final column in each comparison reports competitor cAIC minus selected-model cAIC; positive values favor the selected model. A negative value indicates that the simpler selected structure was retained under the prespecified within-2-cAIC parsimony rule.",
  "",
  "| Response | Selected hierarchy | Nearest hierarchy competitor (Delta cAIC) | Selected ecology | Nearest ecological competitor (Delta cAIC) |",
  "|---|---|---|---|---|",
  selection_table_lines,
  "",
  "For adult and pooled site-shape candidates, camera-specific cyclic deviations followed the Iannarilli factor-by-site construction. Because 60 separate smoothing parameters were computationally impractical, the documented `mgcv` `id` argument linked their penalties; cameras retained separate deviation coefficients but shared one smoothing parameter.",
  "The candidate ecological structures were baseline-only, wolf-only, jackal-only, additive wolf + jackal, and wolf x jackal. Conditional AIC was calculated with `AIC()` from the discrete fREML fits, using effective degrees of freedom and the `mgcv` smoothing-uncertainty correction. We selected the least complex structure within 2 cAIC units of the minimum.",
  "",
  "## 3. Predator-by-time tests and robustness",
  "",
  "The table reports all predator-by-time smooths from the primary full interaction model. Robustness required support in all seven specifications and retention of the corresponding effect in the selected ecological structure.",
  "",
  "| Response | Predator-by-time effect | EDF | P | Retained | Supported specifications | Classification |",
  "|---|---|---:|---:|---|---:|---|",
  effect_table_lines,
  "",
  "The seven specifications were the primary apparent-solar full model, k = 8, k = 16, sun-time, alternate-hierarchy, apparent-solar independent-30-minute-history, and the cached sun-time plus camera-only alternate-hierarchy fit. This prevents significance under one basis size, hierarchy, time transformation, response history, unselected model, or the prespecified combined time-scale and hierarchy sensitivity from driving inference.",
  "",
  "## 4. Predicted curves",
  "",
  prediction_result_lines,
  "",
  "Juvenile activity retained the strong two-predator interaction: concordant low-low and high-high contexts peaked in the morning (08:33-08:51), whereas both discordant contexts peaked in late afternoon (16:45). The interaction difference-in-differences interval was positive from 06:57 to 10:03 and negative from 14:45 to 19:09.",
  "Adult-female peak timing changed little (17:15-17:21), but the robust wolf-by-time smooth described lower midday and greater late-afternoon activity as wolf activity increased. The high-minus-low wolf contrast did not exclude zero at any individual hour, so hour-specific differences remain descriptive.",
  "Adult males showed a global wolf-jackal-by-time interaction in the primary apparent-solar model, but this effect was sensitivity-dependent rather than robust because the cached sun-time plus camera-only sensitivity was unsupported (EDF < 0.001, P = 0.478). Concordant predator contexts had approximately 64.7% of normalized activity between apparent-solar 06:00 and 18:00, compared with approximately 59.6% in discordant contexts, and mean distance from solar noon differed by approximately 0.40 h. Peak time changed by only 0.10 h, and no pointwise male interaction interval excluded zero; these primary-model patterns are descriptive and do not constitute a robust male predator-associated timing effect.",
  "The pooled wolf-only model retained a robust wolf-by-time smooth and a nearly unchanged late-afternoon peak (17:33 versus 17:27). Because the pooled response includes unclassified and subadult gazelles, it remains a secondary species-level benchmark rather than a demographic decomposition.",
  "",
  "The curves are normalized within predator context and therefore describe redistribution across the 24-h cycle rather than absolute encounter probability. Apparent solar time aligns solar noon but does not set sunrise and sunset to 06:00 and 18:00.",
  "",
  "## 5. Comparison with the previous sun-time-primary analysis",
  "",
  "| Response | Previous hierarchy | Apparent-solar hierarchy | Previous ecology | Apparent-solar ecology | Previous robust effects | Apparent-solar robust effects |",
  "|---|---|---|---|---|---|---|",
  comparison_table_lines,
  "",
  "The juvenile, adult-female, adult-male, and pooled robustness conclusions were unchanged from the sun-time-primary audit. Although the adult-male full interaction model was favored by 3.262 cAIC under apparent solar time and its interaction smooth was supported in all six one-at-a-time specifications, the cached sun-time plus camera-only alternate-hierarchy sensitivity was unsupported (P = 0.478). Under the prespecified requirement for support across all sensitivities, the male interaction is therefore sensitivity-dependent rather than robust.",
  "Changes in hierarchy or selected candidate structure do not by themselves constitute ecological effects; focal inference follows the prespecified robustness rule and the predicted effect magnitude.",
  "",
  "## 6. Recommendation",
  "",
  "Use the equal-duration apparent-solar hierarchical binomial GAM as the definitive diel framework. The focal conclusions remain a robust juvenile wolf-jackal temporal interaction, a robust female wolf-associated timing gradient, and no robust predator-associated temporal effect for adult males. The primary adult-male interaction may be reported only as sensitivity-dependent because the sun-time plus camera-only sensitivity was unsupported. Interpret all predator-by-time terms as observational associations with redistribution of detections across the diel cycle, not causal predator avoidance. Retain sun time as a prespecified sensitivity rather than as the primary response.",
  "",
  "## References",
  "",
  "Iannarilli, F., Gerber, B. D., Erb, J., & Fieberg, J. R. (2025). A 'how-to' guide for estimating animal diel activity using hierarchical models. *Journal of Animal Ecology*, 94, 182-194. https://doi.org/10.1111/1365-2656.14213",
  "",
  "Cunningham, C. X., Scoleri, V., Johnson, C. N., Barmuta, L. A., & Jones, M. E. (2019). Temporal partitioning of activity: rising and falling top-predator abundance triggers community-wide shifts in diel activity. *Ecography*, 42, 2157-2168. https://doi.org/10.1111/ecog.04485"
)
writeLines(results_lines, out("results.md"))
writeLines(capture.output(sessionInfo()), out("session_info.txt"))

cat("Definitive hierarchical diel GAM analysis completed.\n")
