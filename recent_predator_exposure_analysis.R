#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(survival)
  library(tidyr)
})

set.seed(20260825)

cfg <- list(
  full_export = file.path("..", "full-export_2026-04-15_2006.csv"),
  class_events = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  corrections = "camera_timestamp_corrections.csv",
  obs_covs = "ObsCovs_050524.csv",
  old_pamm_contrasts = "publication_gazelle_pressure_recurrent_events_pamm_window_contrasts.csv",
  old_pamm_permutations = "publication_gazelle_pressure_recurrent_events_permutation_tests.csv",
  timezone = "Asia/Jerusalem",
  interaction_min_dual_cases = 30L
)

output_summary <- "recent_predator_exposure_summary.csv"
output_results <- "recent_predator_exposure_results.md"
output_figure <- "recent_predator_exposure_figure.png"
output_data_audit <- "recent_predator_exposure_data_audit.csv"

required_files <- unlist(cfg[c(
  "full_export", "class_events", "corrections", "obs_covs",
  "old_pamm_contrasts", "old_pamm_permutations"
)])
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "))
}

class_levels <- c("Juvenile", "Adult female", "Adult male")
class_lookup <- c(
  juvenile = "Juvenile",
  female = "Adult female",
  male = "Adult male"
)
term_levels <- c("wolf_0_6h", "wolf_6_24h", "jackal_0_6h", "jackal_6_24h")
term_labels <- c(
  wolf_0_6h = "Wolf: 0-6 h",
  wolf_6_24h = "Wolf: >6-24 h",
  jackal_0_6h = "Jackal: 0-6 h",
  jackal_6_24h = "Jackal: >6-24 h"
)

parse_local_export_time <- function(x) {
  out <- as.POSIXct(x, format = "%d/%m/%Y %H:%M:%S", tz = cfg$timezone)
  missing <- is.na(out)
  out[missing] <- as.POSIXct(
    x[missing], format = "%d/%m/%Y %H:%M", tz = cfg$timezone
  )
  out
}

parse_frozen_local_time <- function(x) {
  as.POSIXct(sub("Z$", "", x), format = "%Y-%m-%dT%H:%M:%S", tz = cfg$timezone)
}

equation_of_time <- function(day_of_year) {
  b <- 2 * pi * (day_of_year - 81) / 364
  9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
}

fmt <- function(x, digits = 2) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f"))
  )
}

# Operational dates --------------------------------------------------------

obs_raw <- read.csv(
  cfg$obs_covs, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
effort <- data.table(
  camera_id = as.integer(obs_raw[["Camera Number"]]),
  round = as.integer(obs_raw[["Round Number"]]),
  start_date = as.Date(obs_raw[["Start date"]], format = "%d/%m/%Y"),
  end_date = as.Date(obs_raw[["End Date"]], format = "%d/%m/%Y"),
  days_operational = as.integer(obs_raw[["Days.Operational"]])
)[
  !is.na(camera_id) & !is.na(round) & !is.na(start_date) & !is.na(end_date) &
    is.finite(days_operational) & days_operational > 0
]
effort[, nominal_days := as.integer(end_date - start_date)]
effort[, day_discrepancy := nominal_days - days_operational]

ambiguous_effort <- effort[day_discrepancy != 0]
if (
  nrow(effort) != 344L || sum(effort$days_operational) != 5999L ||
    nrow(ambiguous_effort) != 1L ||
    ambiguous_effort$camera_id != 45L || ambiguous_effort$round != 2L ||
    ambiguous_effort$day_discrepancy != 11L
) {
  stop("Operational-date audit differs from the authoritative survey record.")
}

# Camera 45 round 2 has 21 known operational days but 11 unknown outage dates;
# excluding it avoids inventing date-level zero observations.
known_effort <- effort[!(camera_id == 45L & round == 2L)]
operational_dates <- known_effort[, .(
  solar_date = seq(start_date, end_date - 1, by = "day")
), by = .(camera_id, round, days_operational)]
if (nrow(operational_dates) != sum(known_effort$days_operational)) {
  stop("Known operational-date expansion does not match Days.Operational.")
}

# Corrected raw timestamps and independent predator events -----------------

required_export_columns <- c(
  "Camera Name", "Date/Time", "Species", "Trap Station Session ID",
  "Camelot GPS Latitude", "Camelot GPS Longitude"
)
raw <- fread(
  cfg$full_export,
  select = required_export_columns,
  na.strings = c("", "NA"),
  showProgress = FALSE
)
setnames(
  raw,
  required_export_columns,
  c("camera_id", "timestamp_text", "species", "session_id", "latitude", "longitude")
)
raw[, `:=`(
  raw_row_id = .I,
  camera_id = as.integer(camera_id),
  session_id = as.character(session_id),
  timestamp_recorded = parse_local_export_time(timestamp_text),
  latitude = as.numeric(latitude),
  longitude = as.numeric(longitude)
)]
raw <- raw[camera_id %in% 1:60 & !is.na(timestamp_recorded)]

corrections <- fread(cfg$corrections)
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
  paste(corrected_date, clock_time),
  format = "%Y-%m-%d %H:%M:%S", tz = cfg$timezone
)]

camera_locations <- raw[
  is.finite(longitude) & is.finite(latitude),
  .(longitude = median(longitude), latitude = median(latitude)),
  by = camera_id
][order(camera_id)]
if (
  nrow(camera_locations) != 60L ||
    !identical(camera_locations$camera_id, 1:60) || anyNA(camera_locations)
) {
  stop("Expected complete coordinates for cameras 1-60.")
}

predator_rows <- raw[
  species %chin% c("Canis lupus", "Canis aureus"),
  .(camera_id, session_id, timestamp, species, raw_row_id)
]
predator_rows[, event_date := as.Date(timestamp, tz = cfg$timezone)]
predator_rows <- merge(
  predator_rows,
  known_effort[, .(camera_id, round, start_date, end_date)],
  by = "camera_id", allow.cartesian = TRUE, sort = FALSE
)[event_date >= start_date & event_date < end_date]
predator_rows[, predator := fifelse(species == "Canis lupus", "Wolf", "Jackal")]
setorder(
  predator_rows, camera_id, round, session_id, predator, timestamp, raw_row_id
)
predator_rows[, gap_minutes := as.numeric(
  difftime(timestamp, shift(timestamp), units = "mins")
), by = .(camera_id, round, session_id, predator)]
predator_rows[, event_id := cumsum(
  is.na(gap_minutes) | gap_minutes > 30
), by = .(camera_id, round, session_id, predator)]
predator_events <- predator_rows[, .(
  event_time = min(timestamp),
  source_rows = .N
), by = .(camera_id, round, session_id, predator, event_id)]
predator_events[, event_time_num := as.numeric(event_time)]
setorder(predator_events, camera_id, predator, event_time_num)

# Equal-duration apparent-solar hourly intervals ---------------------------

hourly <- merge(
  operational_dates,
  camera_locations[, .(camera_id, longitude, latitude)],
  by = "camera_id", all.x = TRUE, sort = FALSE
)
hourly <- hourly[, .(solar_hour = 0:23), by = .(
  camera_id, round, solar_date, longitude, latitude
)]
hourly[, day_of_year := as.integer(format(solar_date, "%j"))]
hourly[, solar_offset_hours := longitude / 15 + equation_of_time(day_of_year) / 60]
hourly[, solar_midnight_utc := as.numeric(as.POSIXct(
  paste(solar_date, "00:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"
))]
hourly[, interval_start_num :=
  solar_midnight_utc - solar_offset_hours * 3600 + solar_hour * 3600]
hourly[, interval_end_num := interval_start_num + 3600]
setorder(hourly, camera_id, interval_start_num)
hourly[, hour_id := .I]
hourly[, stratum_id := paste(camera_id, round, solar_hour, sep = "_")]

if (nrow(hourly) != 24L * sum(known_effort$days_operational)) {
  stop("Hourly history does not equal 24 intervals per known operational day.")
}

# Independent demographic events use the frozen, corrected class-event file.
class_events <- fread(
  cfg$class_events,
  colClasses = list(character = "event_time")
)
class_events[, `:=`(
  camera_id = as.integer(camera_name),
  session_id = as.character(trap_station_session_id),
  gazelle_class = unname(class_lookup[as.character(gazelle_class)]),
  event_time = parse_frozen_local_time(event_time)
)]
class_events[, event_time_num := as.numeric(event_time)]
if (anyNA(class_events[, .(camera_id, session_id, gazelle_class, event_time_num)])) {
  stop("Frozen class-event file contains unparseable values.")
}
expected_class_events <- c(
  "Juvenile" = 180L, "Adult female" = 1477L, "Adult male" = 2268L
)
observed_class_events <- table(class_events$gazelle_class)
if (!all(observed_class_events[names(expected_class_events)] == expected_class_events)) {
  stop("Frozen class-event totals differ from the authoritative audit.")
}

map_events_to_hours <- function(events, intervals) {
  mapped <- vector("list", length = 60L)
  for (camera_i in 1:60) {
    event_i <- events[camera_id == camera_i]
    hour_i <- intervals[camera_id == camera_i]
    if (!nrow(event_i) || !nrow(hour_i)) next
    idx <- findInterval(event_i$event_time_num, hour_i$interval_start_num)
    valid <- idx > 0L
    valid[valid] <- event_i$event_time_num[valid] < hour_i$interval_end_num[idx[valid]]
    if (any(valid)) {
      mapped[[camera_i]] <- data.table(
        hour_id = hour_i$hour_id[idx[valid]],
        camera_id = camera_i,
        gazelle_class = event_i$gazelle_class[valid]
      )
    }
  }
  unique(rbindlist(mapped, fill = TRUE))
}

positive_hours <- map_events_to_hours(class_events, hourly)
mapped_event_audit <- class_events[, .(source_events = .N), by = gazelle_class]
mapped_event_audit <- merge(
  mapped_event_audit,
  positive_hours[, .(positive_hours = .N), by = gazelle_class],
  by = "gazelle_class", all.x = TRUE
)

# Recent exposures are evaluated at the start of each response hour. This
# prevents predator detections later in the gazelle hour from being treated as
# preceding exposure.
exposure_rows <- vector("list", length = 60L)
for (camera_i in 1:60) {
  starts <- hourly[camera_id == camera_i, interval_start_num]
  wolf_times <- predator_events[
    camera_id == camera_i & predator == "Wolf", event_time_num
  ]
  jackal_times <- predator_events[
    camera_id == camera_i & predator == "Jackal", event_time_num
  ]
  wolf_times <- sort(wolf_times)
  jackal_times <- sort(jackal_times)

  window_flags <- function(times) {
    if (!length(times)) {
      return(list(recent = integer(length(starts)), earlier = integer(length(starts))))
    }
    at_start <- findInterval(starts, times)
    before_6 <- findInterval(starts - 6 * 3600, times)
    before_24 <- findInterval(starts - 24 * 3600, times)
    list(
      recent = as.integer((at_start - before_6) > 0L),
      earlier = as.integer((before_6 - before_24) > 0L)
    )
  }

  wolf <- window_flags(wolf_times)
  jackal <- window_flags(jackal_times)
  exposure_rows[[camera_i]] <- data.table(
    camera_id = camera_i,
    interval_start_num = starts,
    wolf_0_6h = wolf$recent,
    wolf_6_24h = wolf$earlier,
    jackal_0_6h = jackal$recent,
    jackal_6_24h = jackal$earlier
  )
}
exposures <- rbindlist(exposure_rows)
hourly <- merge(
  hourly, exposures,
  by = c("camera_id", "interval_start_num"), all.x = TRUE, sort = FALSE
)
setorder(hourly, hour_id)
hourly[, `:=`(
  wolf_any_24h = as.integer(wolf_0_6h == 1L | wolf_6_24h == 1L),
  jackal_any_24h = as.integer(jackal_0_6h == 1L | jackal_6_24h == 1L)
)]
hourly[, exposure_category := fcase(
  wolf_any_24h == 0L & jackal_any_24h == 0L, "Neither",
  wolf_any_24h == 1L & jackal_any_24h == 0L, "Wolf only",
  wolf_any_24h == 0L & jackal_any_24h == 1L, "Jackal only",
  default = "Both"
)]
hourly[, exposure_category := factor(
  exposure_category, levels = c("Neither", "Wolf only", "Jackal only", "Both")
)]

class_hourly <- hourly[
  rep(seq_len(nrow(hourly)), each = length(class_levels))
]
class_hourly[, gazelle_class := factor(
  rep(class_levels, times = nrow(hourly)), levels = class_levels
)]
class_hourly[, response := 0L]
class_hourly[
  positive_hours,
  on = .(hour_id, camera_id, gazelle_class),
  response := 1L
]

# Conditional case-crossover models ---------------------------------------

fit_class_model <- function(class_name) {
  dat <- class_hourly[gazelle_class == class_name]
  informative <- dat[, .(
    cases = sum(response), controls = sum(response == 0L)
  ), by = stratum_id][cases > 0L & controls > 0L]
  dat <- dat[stratum_id %in% informative$stratum_id]
  dat[, camera_cluster := factor(camera_id)]

  fit <- clogit(
    response ~ wolf_0_6h + wolf_6_24h + jackal_0_6h + jackal_6_24h +
      strata(stratum_id) + cluster(camera_cluster),
    data = dat,
    method = "efron",
    na.action = na.fail,
    model = TRUE, x = TRUE, y = TRUE
  )

  coefficient_table <- as.data.frame(summary(fit)$coefficients)
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL
  se_column <- if ("robust se" %in% names(coefficient_table)) "robust se" else "se(coef)"
  p_column <- grep("Pr\\(", names(coefficient_table), value = TRUE)[1]
  coefficient_table <- as_tibble(coefficient_table) %>%
    transmute(
      gazelle_class = class_name,
      result_type = "Primary exposure effect",
      term,
      estimate = coef,
      std_error = .data[[se_column]],
      odds_ratio = exp(estimate),
      lower = exp(estimate - 1.96 * std_error),
      upper = exp(estimate + 1.96 * std_error),
      statistic = estimate / std_error,
      p_value = .data[[p_column]]
    )

  vc <- vcov(fit)
  beta <- coef(fit)
  make_contrast <- function(name, weights, type) {
    weights <- weights[names(beta)]
    est <- sum(weights * beta)
    se <- sqrt(as.numeric(t(weights) %*% vc %*% weights))
    tibble(
      gazelle_class = class_name,
      result_type = type,
      term = name,
      estimate = est,
      std_error = se,
      odds_ratio = exp(est),
      lower = exp(est - 1.96 * se),
      upper = exp(est + 1.96 * se),
      statistic = est / se,
      p_value = 2 * pnorm(-abs(est / se))
    )
  }
  zero_weights <- setNames(rep(0, length(beta)), names(beta))
  contrast_rows <- list()
  w <- zero_weights; w["wolf_0_6h"] <- 1; w["wolf_6_24h"] <- -1
  contrast_rows[[1]] <- make_contrast(
    "wolf_0_6h_vs_6_24h", w, "Within-predator lag contrast"
  )
  w <- zero_weights; w["jackal_0_6h"] <- 1; w["jackal_6_24h"] <- -1
  contrast_rows[[2]] <- make_contrast(
    "jackal_0_6h_vs_6_24h", w, "Within-predator lag contrast"
  )
  w <- zero_weights; w["wolf_0_6h"] <- 1; w["jackal_0_6h"] <- -1
  contrast_rows[[3]] <- make_contrast(
    "wolf_vs_jackal_0_6h", w, "Between-predator contrast"
  )
  w <- zero_weights; w["wolf_6_24h"] <- 1; w["jackal_6_24h"] <- -1
  contrast_rows[[4]] <- make_contrast(
    "wolf_vs_jackal_6_24h", w, "Between-predator contrast"
  )

  case_counts <- rbindlist(lapply(term_levels, function(term_i) {
    data.table(
      gazelle_class = class_name,
      term = term_i,
      exposed_gazelle_hours = dat[response == 1L & get(term_i) == 1L, .N],
      exposed_model_hours = dat[get(term_i) == 1L, .N]
    )
  }))
  coefficient_table <- left_join(
    coefficient_table, as_tibble(case_counts),
    by = c("gazelle_class", "term")
  )
  contrast_table <- bind_rows(contrast_rows) %>%
    mutate(
      exposed_gazelle_hours = NA_integer_,
      exposed_model_hours = NA_integer_
    )
  model_audit <- data.table(
    gazelle_class = class_name,
    source_positive_hours = class_hourly[gazelle_class == class_name, sum(response)],
    model_positive_hours = dat[, sum(response)],
    model_hours = nrow(dat),
    informative_strata = uniqueN(dat$stratum_id),
    cameras_with_cases = uniqueN(dat[response == 1L, camera_id]),
    dual_exposure_positive_hours = dat[
      response == 1L & wolf_any_24h == 1L & jackal_any_24h == 1L, .N
    ],
    converged = all(is.finite(coef(fit)))
  )
  list(
    fit = fit,
    results = bind_rows(coefficient_table, contrast_table),
    audit = model_audit
  )
}

model_objects <- setNames(lapply(class_levels, fit_class_model), class_levels)
model_results <- bind_rows(lapply(model_objects, `[[`, "results"))
model_audit <- rbindlist(lapply(model_objects, `[[`, "audit"))

interaction_adequate <- all(
  model_audit$dual_exposure_positive_hours >= cfg$interaction_min_dual_cases
)
interaction_fitted <- FALSE

model_results <- model_results %>%
  left_join(as_tibble(model_audit), by = "gazelle_class") %>%
  mutate(
    term_label = recode(
      term,
      !!!as.list(term_labels),
      wolf_0_6h_vs_6_24h = "Wolf: 0-6 h versus >6-24 h",
      jackal_0_6h_vs_6_24h = "Jackal: 0-6 h versus >6-24 h",
      wolf_vs_jackal_0_6h = "Wolf versus jackal: 0-6 h",
      wolf_vs_jackal_6_24h = "Wolf versus jackal: >6-24 h"
    ),
    inference_flag = case_when(
      gazelle_class == "Juvenile" & result_type == "Primary exposure effect" &
        exposed_gazelle_hours < 30 ~ "Sparse exposed juvenile outcomes",
      p_value < 0.05 ~ "Nominal support",
      TRUE ~ "Unsupported"
    )
  )
write_csv(model_results, output_summary)

exposure_category_counts <- hourly[, .N, by = exposure_category][
  order(exposure_category)
]
exposure_term_counts <- rbindlist(lapply(term_levels, function(term_i) {
  data.table(term = term_i, hourly_observations = hourly[get(term_i) == 1L, .N])
}))
case_exposure_counts <- class_hourly[response == 1L, .(
  positive_hours = .N,
  wolf_0_6h = sum(wolf_0_6h),
  wolf_6_24h = sum(wolf_6_24h),
  jackal_0_6h = sum(jackal_0_6h),
  jackal_6_24h = sum(jackal_6_24h),
  wolf_only_24h = sum(exposure_category == "Wolf only"),
  jackal_only_24h = sum(exposure_category == "Jackal only"),
  both_24h = sum(exposure_category == "Both"),
  neither_24h = sum(exposure_category == "Neither")
), by = gazelle_class]

data_audit <- rbindlist(list(
  data.table(
    audit_type = "Survey effort",
    item = c(
      "Recorded camera-days", "Known camera-days analyzed",
      "Known hourly observations", "Excluded ambiguous camera-round"
    ),
    value = c(
      sum(effort$days_operational), sum(known_effort$days_operational),
      nrow(hourly), 1
    ),
    detail = c(
      "All positive-effort rounds",
      "Camera 45 round 2 excluded",
      "24 true one-hour apparent-solar intervals per known day",
      "Camera 45 round 2: 21 operational days across 32 nominal days"
    )
  ),
  predator_events[, .(
    audit_type = "Independent predator events",
    item = predator,
    value = .N,
    detail = paste0(uniqueN(camera_id), " cameras")
  ), by = predator],
  exposure_category_counts[, .(
    audit_type = "Hourly 24-h exposure category",
    item = as.character(exposure_category),
    value = N,
    detail = paste0(fmt(100 * N / sum(N), 1), "% of analyzed hours")
  )],
  exposure_term_counts[, .(
    audit_type = "Hourly exposure indicator",
    item = term,
    value = hourly_observations,
    detail = paste0(fmt(100 * hourly_observations / nrow(hourly), 1), "% of analyzed hours")
  )]
), fill = TRUE)
fwrite(data_audit, output_data_audit)

# Figure -------------------------------------------------------------------

figure_data <- model_results %>%
  filter(result_type == "Primary exposure effect") %>%
  mutate(
    gazelle_class = factor(gazelle_class, levels = class_levels),
    term_label = factor(
      term_label,
      levels = rev(unname(term_labels[term_levels]))
    ),
    predator = if_else(grepl("^Wolf", term_label), "Wolf", "Jackal")
  )

exposure_figure <- ggplot(
  figure_data,
  aes(x = odds_ratio, y = term_label, color = predator)
) +
  geom_vline(xintercept = 1, color = "#4F5B5F", linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper), width = 0.13, linewidth = 0.65,
    orientation = "y"
  ) +
  geom_point(size = 2.4) +
  facet_wrap(~gazelle_class, nrow = 1) +
  scale_x_log10() +
  scale_color_manual(values = c(Wolf = "#315F8C", Jackal = "#C26B2E"), guide = "none") +
  labs(
    title = "Gazelle detections following recent local predator passage",
    subtitle = paste0(
      "Conditional within camera, survey round, and apparent-solar hour; ",
      "points are odds ratios with 95% cluster-robust confidence intervals"
    ),
    x = "Odds ratio for gazelle detection",
    y = NULL,
    caption = paste0(
      "Exposure was based on 30-minute-independent predator events before the ",
      "start of each response hour. Camera 45 round 2 was excluded because exact ",
      "outage dates were unavailable."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4D5A5E"),
    plot.caption = element_text(hjust = 0, color = "#5F696B"),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(8, 10, 8, 8)
  )
ggsave(output_figure, exposure_figure, width = 11, height = 5.6, dpi = 400, bg = "white")

# Existing PAMM audit and comparative report -------------------------------

old_contrasts <- read_csv(cfg$old_pamm_contrasts, show_col_types = FALSE) %>%
  filter(
    direction == "Pressure to gazelle",
    primary %in% c("Wolf", "Jackal"),
    window %in% c("0-6 hours", "6-24 hours")
  )
old_permutation <- read_csv(cfg$old_pamm_permutations, show_col_types = FALSE) %>%
  filter(
    direction == "Pressure to gazelle",
    primary %in% c("Wolf", "Jackal")
  )

pamm_rows <- old_contrasts %>%
  left_join(
    old_permutation %>%
      select(primary, secondary, permutation_p_value, p_bh_global),
    by = c("primary", "secondary")
  ) %>%
  arrange(match(primary, c("Wolf", "Jackal")), secondary, window)
pamm_table <- vapply(seq_len(nrow(pamm_rows)), function(i) {
  row <- pamm_rows[i, ]
  paste0(
    "| ", row$primary, " | ", row$secondary, " | ", row$window, " | ",
    fmt(row$hazard_ratio), " [", fmt(row$lwr), ", ", fmt(row$upr), "] | ",
    fmt_p(row$permutation_p_value), " | ", fmt_p(row$p_bh_global), " |"
  )
}, character(1))

exposure_table <- exposure_category_counts[, paste0(
  "| ", exposure_category, " | ", format(N, big.mark = ","), " | ",
  fmt(100 * N / sum(N), 1), "% |"
)]

primary_rows <- model_results %>%
  filter(result_type == "Primary exposure effect") %>%
  arrange(match(gazelle_class, class_levels), match(term, term_levels))
primary_table <- vapply(seq_len(nrow(primary_rows)), function(i) {
  row <- primary_rows[i, ]
  paste0(
    "| ", row$gazelle_class, " | ", row$term_label, " | ",
    fmt(row$odds_ratio), " [", fmt(row$lower), ", ", fmt(row$upper), "] | ",
    fmt_p(row$p_value), " | ", row$exposed_gazelle_hours, " | ",
    row$model_positive_hours, " | ", row$inference_flag, " |"
  )
}, character(1))

contrast_rows <- model_results %>%
  filter(result_type != "Primary exposure effect") %>%
  arrange(match(gazelle_class, class_levels), result_type, term)
contrast_table <- vapply(seq_len(nrow(contrast_rows)), function(i) {
  row <- contrast_rows[i, ]
  paste0(
    "| ", row$gazelle_class, " | ", row$term_label, " | ",
    fmt(row$odds_ratio), " [", fmt(row$lower), ", ", fmt(row$upper), "] | ",
    fmt_p(row$p_value), " |"
  )
}, character(1))

case_count_rows <- vapply(seq_len(nrow(case_exposure_counts)), function(i) {
  row <- case_exposure_counts[i]
  paste0(
    "| ", row$gazelle_class, " | ", row$positive_hours, " | ",
    row$wolf_0_6h, " | ", row$wolf_6_24h, " | ",
    row$jackal_0_6h, " | ", row$jackal_6_24h, " | ", row$both_24h, " |"
  )
}, character(1))

nominal_primary <- primary_rows %>% filter(p_value < 0.05, lower > 1 | upper < 1)
supported_text <- if (nrow(nominal_primary)) {
  paste0(
    "Nominal primary associations were: ",
    paste0(
      nominal_primary$gazelle_class, " ", nominal_primary$term_label,
      " (OR = ", fmt(nominal_primary$odds_ratio), ")",
      collapse = "; "
    ), ". These are unadjusted, prespecified coefficient-level tests."
  )
} else {
  "No primary recent-exposure coefficient had a 95% confidence interval excluding one."
}

interaction_text <- if (interaction_adequate) {
  paste0(
    "Dual-exposure outcome counts met the numerical threshold in every class, but ",
    "wolf-jackal interactions were not fitted because they were not part of the ",
    "prespecified primary estimand and would multiply four window-specific interaction ",
    "choices. The additive model is retained."
  )
} else {
  paste0(
    "Wolf-jackal recent-exposure interactions were not fitted. The prespecified ",
    "adequacy threshold was at least ", cfg$interaction_min_dual_cases,
    " dual-exposure positive hours in every class; observed class counts were ",
    paste0(model_audit$gazelle_class, " = ", model_audit$dual_exposure_positive_hours,
           collapse = ", "), "."
  )
}

report <- c(
  "# Recent local predator exposure analysis",
  "",
  "## 1. PAMM audit",
  "",
  "The existing wolf-jackal PAMM was `analyze_gazelle_pressure_recurrent_events.R`, with results in `publication_gazelle_pressure_recurrent_events_results.md`. It modeled the recurrent-event rate of a demographic gazelle class after each focal predator timestamp at the same camera and Camelot camera session. Exact duplicate export rows at a timestamp were collapsed, but the retained wolf-jackal run used every distinct timestamp rather than 30-minute-independent events. Follow-up began at a wolf or jackal timestamp, ended at the next same-predator event, camera-session end, or 14 days, and could cross calendar dates. Thus, the 0-6 h interval necessarily represented the same or immediately following date, but the complete risk interval was not restricted to one date.",
  "",
  "The response was the number of fawn, adult-female, or adult-male timestamp visits in piecewise follow-up intervals, with log interval duration as an offset. Elapsed time since the predator event, clock hour, and day of year were smooth terms. Camera and camera-session random-effect smooths were fitted; all focal wolf/jackal models retained both frailties. Wolf and jackal were analyzed in separate directional dyads, so the PAMM contained no joint wolf and jackal main effects or interaction. Its 0-6 h and 6-24 h summaries were contrasts against a late, taxon-specific data-supported reference interval, not against otherwise matched hours with no recent predator passage.",
  "",
  "| Predator | Gazelle class | PAMM window | Hazard ratio (95% CI) | Randomization P | Global BH q |",
  "|---|---|---|---:|---:|---:|",
  pamm_table,
  "",
  "Although several model-based PAMM intervals excluded one, no pressure-to-gazelle association passed the camera-session circular-shift test after global false-discovery-rate adjustment. The randomization was the stronger PAMM check because it preserved site use, event totals, clock time, and within-taxon clustering.",
  "",
  "## 2. Is this analysis redundant with the PAMM?",
  "",
  "No. The PAMM already addressed the broad question of directional temporal clustering after predator records, but it did not estimate the requested matched hourly odds under explicit 0-6 h and >6-24 h independent-passage exposures. The new analysis uses 30-minute-independent predator events, binary hourly gazelle encounter histories, simultaneous additive wolf and jackal terms, and exact matching by camera, survey round, and apparent-solar hour. It should be viewed as a focused alternative estimand and robustness analysis, not as an unrelated second test.",
  "",
  "## 3. Data structure and exposure counts",
  "",
  paste0(
    "We analyzed ", format(nrow(hourly), big.mark = ","),
    " camera-hours from ", format(nrow(operational_dates), big.mark = ","),
    " known operational camera-days. Camera 45 round 2 was excluded because its ",
    "21 operational days occurred within 32 nominal dates but the 11 outage dates ",
    "were unavailable. Excluding this round prevented unobserved outage hours from ",
    "being coded as gazelle absences. Other rounds matched their start-inclusive, ",
    "end-exclusive date spans exactly."
  ),
  paste0(
    "Predator exposures were derived from ",
    predator_events[predator == "Wolf", .N], " independent wolf events and ",
    predator_events[predator == "Jackal", .N],
    " independent jackal events. A predator record after the start of a response ",
    "hour could not classify that hour as exposed."
  ),
  "",
  "| Recent exposure category (preceding 24 h) | Camera-hours | Percent |",
  "|---|---:|---:|",
  exposure_table,
  "",
  "| Class | Positive hours | Wolf 0-6 h | Wolf >6-24 h | Jackal 0-6 h | Jackal >6-24 h | Both predators within 24 h |",
  "|---|---:|---:|---:|---:|---:|---:|",
  case_count_rows,
  "",
  interaction_text,
  "",
  "## 4. Model choice",
  "",
  "For each demographic class, we fitted a conditional logistic case-crossover model with recent wolf and jackal exposure indicators in the same additive model. Strata were camera x survey round x apparent-solar hour. Consequently, each coefficient compared exposed and unexposed operational dates at the same camera, during the same one-hour diel bin and survey period. Cluster-robust standard errors accounted for repeated strata within cameras. This design more directly removes stable camera-level habitat and predator-index differences than a conventional random-intercept GLMM and controls the diel cycle without requiring an additional estimated smooth.",
  "",
  "The response was at least one independent class encounter in the hour. Exposures were evaluated at the hour start, making the design conservative for predator detections occurring earlier within the same hour. Odds ratios are conditional associations and do not identify individual gazelles or demonstrate causal avoidance or attraction.",
  "Reported P values are nominal and were not adjusted across classes and lag windows; inference therefore emphasizes confidence intervals, exposed outcomes, and the prespecified within-predator lag contrasts.",
  "",
  "## 5. Results by demographic class",
  "",
  "| Class | Recent predator exposure | Odds ratio (95% CI) | P | Exposed positive hours | Model positive hours | Assessment |",
  "|---|---|---:|---:|---:|---:|---|",
  primary_table,
  "",
  supported_text,
  "None of the 0-6 h coefficients or within-predator lag contrasts had P < 0.05. The analysis therefore did not resolve an immediate reactive change in hourly gazelle detection following predator passage. The two adult-male >6-24 h estimates are modest nominal signals; because neither differed from its corresponding 0-6 h coefficient, they do not establish a delayed response.",
  "",
  "Lag and between-predator contrasts compare coefficients within the same fitted model:",
  "",
  "| Class | Contrast | Ratio of odds ratios (95% CI) | P |",
  "|---|---|---:|---:|",
  contrast_table,
  "",
  "Juvenile inference should be judged from exposed positive-hour counts and interval width, not from P values alone. An unsupported juvenile coefficient is not evidence that juveniles fail to respond when the relevant exposure-outcome combination is sparse.",
  "",
  "## 6. Comparison with circular GAM and previous dynamic analysis",
  "",
  "- **Main hierarchical circular GAM:** asks whether diel shape changes along continuous, site-level wolf and jackal activity gradients after accounting for hierarchical camera structure. Predator indices are spatially constant and do not represent a recent passage.",
  "- **Previous camera x round x hour dynamic GAMM:** asks whether gazelle timing covaries with recurring local predator activity patterns averaged within camera-round-hour combinations. It represents repeated local temporal niches, not the occurrence of one immediately preceding predator event.",
  "- **Existing PAMM:** asks how the recurrent gazelle timestamp-visit hazard changes smoothly after a focal predator timestamp in a same-camera session, with dyads fitted separately and a circular-shift null test.",
  "- **New recent-exposure analysis:** asks whether binary gazelle detection odds differ after an actual independent local predator passage in the previous 0-6 h or >6-24 h, relative to matched dates at the same camera, round, and apparent-solar hour.",
  "",
  "## 7. Manuscript recommendation",
  "",
  if (nrow(nominal_primary)) {
    paste0(
      "Treat this as a focused supplementary acute-response analysis. The modest ",
      "adult-male >6-24 h associations should be labeled nominal and should not be ",
      "interpreted as a delayed response because the direct lag contrasts were ",
      "unsupported. The analysis should not replace the circular GAM's site-level ",
      "diel inference or be described as causal avoidance or attraction."
    )
  } else {
    paste0(
      "Do not add a second null acute-response analysis to the main text. The result ",
      "is most useful as a supplementary robustness check alongside the earlier PAMM. ",
      "It does not prove absence of avoidance; it indicates that this design did not ",
      "resolve a consistent hourly response after recent local predator passage."
    )
  },
  "",
  "The manuscript itself was not edited.",
  "",
  "## Outputs",
  "",
  "- `recent_predator_exposure_summary.csv`",
  "- `recent_predator_exposure_results.md`",
  "- `recent_predator_exposure_figure.png`",
  "- `recent_predator_exposure_data_audit.csv`"
)
writeLines(report, output_results)

saveRDS(model_objects, "recent_predator_exposure_models.rds")
writeLines(capture.output(sessionInfo()), "recent_predator_exposure_session_info.txt")

cat("Recent predator exposure analysis completed.\n")
