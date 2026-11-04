#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
  library(patchwork)
})

# Directional camera-trap interaction analysis following the recurrent-event
# PAMM framework of Ferry et al. (2024; doi:10.1111/2041-210X.14338).

args <- commandArgs(trailingOnly = TRUE)
cfg <- list(
  full_export = file.path("..", "full-export_2026-04-15_2006.csv"),
  obs_cov_file = "ObsCovs_050524.csv",
  prefix = "publication_gazelle_pressure_recurrent_events",
  timezone = "Asia/Jerusalem",
  event_gap_minutes = "30",
  max_follow_days = "14",
  permutation_reps = "499",
  seed = "20260804",
  split_date = "auto",
  save_models = "true",
  audit_only = "false",
  censor_competing_pressures = "false",
  sampling_frame = "camelot_sessions",
  event_definition = "timestamp",
  variance_family = "auto",
  pressures = "Wolf,Jackal,Boar,Human"
)

for (arg in args) {
  pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(pieces) == 2) {
    cfg[[sub("^--", "", pieces[1])]] <- pieces[2]
  }
}

event_gap_minutes <- as.numeric(cfg$event_gap_minutes)
max_follow_days <- as.numeric(cfg$max_follow_days)
permutation_reps <- as.integer(cfg$permutation_reps)
seed <- as.integer(cfg$seed)
save_models <- tolower(cfg$save_models) %in% c("true", "1", "yes")
audit_only <- tolower(cfg$audit_only) %in% c("true", "1", "yes")
censor_competing_pressures <- tolower(cfg$censor_competing_pressures) %in% c("true", "1", "yes")
split_date <- if (tolower(cfg$split_date) == "auto") as.Date(NA) else as.Date(cfg$split_date)
sampling_frame <- tolower(cfg$sampling_frame)
event_definition <- tolower(cfg$event_definition)
variance_family <- tolower(cfg$variance_family)
if (!sampling_frame %in% c("camelot_sessions", "obs_rounds")) {
  stop("sampling_frame must be 'camelot_sessions' or 'obs_rounds'.")
}
if (!event_definition %in% c("timestamp", "gap")) {
  stop("event_definition must be 'timestamp' or 'gap'.")
}
if (variance_family == "auto") {
  variance_family <- if (event_definition == "timestamp") "quasipoisson" else "poisson"
}
if (!variance_family %in% c("poisson", "quasipoisson")) {
  stop("variance_family must be 'auto', 'poisson', or 'quasipoisson'.")
}
set.seed(seed)

class_levels <- c("Fawn", "Female", "Male")
available_pressure_levels <- c("Wolf", "Jackal", "Boar", "Human", "Dog")
pressure_levels <- trimws(strsplit(cfg$pressures, ",", fixed = TRUE)[[1]])
pressure_levels <- pressure_levels[nzchar(pressure_levels)]
unknown_pressures <- setdiff(pressure_levels, available_pressure_levels)
if (!length(pressure_levels) || length(unknown_pressures)) {
  stop(
    "pressures must contain one or more of: ",
    paste(available_pressure_levels, collapse = ", "),
    if (length(unknown_pressures)) paste0(". Unknown: ", paste(unknown_pressures, collapse = ", ")) else ""
  )
}
taxon_levels <- c(class_levels, pressure_levels)
direction_levels <- c("Pressure to gazelle", "Gazelle to pressure")
all_pressure_colors <- c(
  Wolf = "#315f8c",
  Jackal = "#c26b2e",
  Boar = "#5b6f45",
  Human = "#7b4f3b",
  Dog = "#a05a2c"
)
pressure_colors <- all_pressure_colors[pressure_levels]
class_colors <- c(Fawn = "#168a79", Female = "#bd3438", Male = "#183b5b")

fmt_num <- function(x, digits = 2) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f"))
  )
}

parse_local_datetime <- function(x, timezone) {
  as.POSIXct(x, format = "%d/%m/%Y %H:%M", tz = timezone)
}

safe_quantile <- function(x, prob) {
  if (!length(x) || all(!is.finite(x))) return(NA_real_)
  unname(stats::quantile(x[is.finite(x)], prob, names = FALSE, type = 8))
}

theme_publication <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      legend.position = "bottom"
    )
}

required_columns <- c(
  "Camera Name", "Date/Time", "Sighting Quantity", "Species",
  "Trap Station Session ID", "Session Start Date", "Session End Date",
  "Sex", "Life stage"
)

message("Reading full Camelot export: ", normalizePath(cfg$full_export))
raw <- fread(
  cfg$full_export,
  select = required_columns,
  na.strings = c("", "NA"),
  showProgress = FALSE
)

missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns)) {
  stop("Full export is missing columns: ", paste(missing_columns, collapse = ", "))
}

setnames(
  raw,
  required_columns,
  c(
    "camera", "event_time", "quantity", "species", "session_id",
    "session_start", "session_end", "sex", "life_stage"
  )
)

raw[, `:=`(
  camera = as.integer(camera),
  session_id = as.character(session_id),
  event_time = parse_local_datetime(event_time, cfg$timezone),
  session_start = parse_local_datetime(session_start, cfg$timezone),
  session_end = parse_local_datetime(session_end, cfg$timezone),
  quantity = suppressWarnings(as.numeric(quantity)),
  sex = tolower(trimws(as.character(sex))),
  life_stage = tolower(trimws(as.character(life_stage)))
)]
raw[!is.finite(quantity), quantity := 1]

human_labels <- c(
  "Homo sapiens", "Homo sapiens-solider", "Bicycle Rider",
  "Vehicle military", "Vehicle civilian", "Vehicle construction"
)
direct_human_labels <- c("Homo sapiens", "Homo sapiens-solider")

pressure_species <- list(
  Wolf = "Canis lupus",
  Jackal = "Canis aureus",
  Boar = "Sus scrofa",
  Human = human_labels,
  Dog = "Canis familiaris"
)
selected_pressure_species <- unname(unlist(pressure_species[pressure_levels], use.names = FALSE))
selected_species <- c("Gazella gazella", selected_pressure_species)

focal_raw <- raw[
  species %chin% selected_species
]
focal_raw <- focal_raw[!is.na(camera) & !is.na(session_id) & !is.na(event_time)]

focal_raw[, taxon := fcase(
  species == "Gazella gazella" & life_stage == "juvenile", "Fawn",
  species == "Gazella gazella" & sex == "female" & !life_stage %chin% c("juvenile", "sub-adult"), "Female",
  species == "Gazella gazella" & sex == "male" & !life_stage %chin% c("juvenile", "sub-adult"), "Male",
  species == "Canis lupus", "Wolf",
  species == "Canis aureus", "Jackal",
  species == "Sus scrofa", "Boar",
  species %chin% human_labels, "Human",
  species == "Canis familiaris", "Dog",
  default = NA_character_
)]

if (sampling_frame == "camelot_sessions") {
  # Camelot session boundaries define the observation-at-risk intervals.
  focal_raw <- focal_raw[
    !is.na(session_start) & !is.na(session_end) &
      session_end > session_start & event_time >= session_start & event_time <= session_end
  ]
  focal_raw[, camera_session := paste(camera, session_id, sep = "_")]
  sessions <- focal_raw[, .(
    session_start = min(session_start),
    session_end = max(session_end),
    inconsistent_start = uniqueN(session_start) > 1,
    inconsistent_end = uniqueN(session_end) > 1
  ), by = .(camera, session_id, camera_session)]
  survey_start_date <- as.Date(min(sessions$session_start), tz = cfg$timezone)
  survey_end_date <- as.Date(max(sessions$session_end), tz = cfg$timezone)
} else {
  # Optional compatibility mode for the abundance-model sampling frame.
  obs_raw <- fread(cfg$obs_cov_file, na.strings = c("", "NA"), showProgress = FALSE)
  required_obs_columns <- c(
    "Camera Number", "Round Number", "Start date", "End Date", "Days.Operational"
  )
  missing_obs_columns <- setdiff(required_obs_columns, names(obs_raw))
  if (length(missing_obs_columns)) {
    stop("Observation-covariate file is missing columns: ", paste(missing_obs_columns, collapse = ", "))
  }
  active_obs <- obs_raw[, .(
    camera = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    round_start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    round_end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    days_operational = suppressWarnings(as.numeric(Days.Operational))
  )][
    is.finite(days_operational) & days_operational > 0 &
      !is.na(round_start_date) & !is.na(round_end_date)
  ]
  focal_raw[, event_row_id := .I]
  focal_raw[, event_date := as.Date(event_time, tz = cfg$timezone)]
  round_candidates <- merge(
    focal_raw[, .(event_row_id, camera, event_date)],
    active_obs,
    by = "camera",
    allow.cartesian = TRUE
  )[event_date >= round_start_date & event_date <= round_end_date]
  setorder(round_candidates, event_row_id, -round_start_date)
  round_map <- round_candidates[, .SD[1], by = event_row_id]
  focal_raw <- merge(
    focal_raw,
    round_map[, .(event_row_id, round, round_start_date, round_end_date)],
    by = "event_row_id",
    all = FALSE
  )
  focal_raw[, original_session_id := session_id]
  focal_raw[, `:=`(
    session_start = as.POSIXct(round_start_date, tz = cfg$timezone),
    session_end = as.POSIXct(round_end_date + 1, tz = cfg$timezone),
    session_id = paste0("round_", round, "__camelot_", original_session_id),
    camera_session = paste(camera, paste0("round_", round), sep = "_")
  )]
  sessions <- focal_raw[, .(
    session_id = paste0("round_", round[1]),
    session_start = min(session_start),
    session_end = max(session_end),
    inconsistent_start = FALSE,
    inconsistent_end = FALSE
  ), by = .(camera, camera_session)]
  survey_start_date <- min(active_obs$round_start_date)
  survey_end_date <- max(active_obs$round_end_date)
}
survey_midpoint_date <- survey_start_date + floor(as.numeric(survey_end_date - survey_start_date) / 2)
if (is.na(split_date)) split_date <- survey_midpoint_date

classification_audit <- rbindlist(list(
  raw[
    species %chin% selected_species,
    .(stage = "Raw target rows", n = .N),
    by = .(source_label = species)
  ],
  focal_raw[
    !is.na(taxon),
    .(stage = paste("Valid", sampling_frame, "rows"), n = .N),
    by = .(source_label = taxon)
  ]
), use.names = TRUE, fill = TRUE)

unclassified_gazelle <- focal_raw[species == "Gazella gazella" & is.na(taxon), .N]
focal_raw <- focal_raw[!is.na(taxon) & taxon %chin% taxon_levels]

# Exact duplicate export rows are repeated annotations/images, not additional
# temporal detections. Distinct timestamps remain separate in timestamp mode.
deduplicated <- focal_raw[
  , .(
    quantity = max(quantity, na.rm = TRUE),
    raw_rows = .N,
    direct_human = any(species %chin% direct_human_labels)
  ),
  by = .(camera, session_id, camera_session, session_start, session_end, taxon, event_time)
]
setorder(deduplicated, camera, session_id, taxon, event_time)
if (event_definition == "timestamp") {
  deduplicated[, analysis_event_id := seq_len(.N), by = .(camera, session_id, taxon)]
} else {
  deduplicated[, gap_minutes := as.numeric(difftime(event_time, shift(event_time), units = "mins")),
    by = .(camera, session_id, taxon)]
  deduplicated[, analysis_event_id := cumsum(is.na(gap_minutes) | gap_minutes > event_gap_minutes),
    by = .(camera, session_id, taxon)]
}

events <- deduplicated[
  , .(
    event_time = min(event_time),
    quantity = max(quantity, na.rm = TRUE),
    exact_records = sum(raw_rows),
    direct_human = any(direct_human)
  ),
  by = .(
    camera, session_id, camera_session, session_start, session_end,
    taxon, analysis_event_id
  )
]
setorder(events, camera, session_id, taxon, event_time)

event_audit <- events[, .(
  analysis_events = .N,
  cameras = uniqueN(camera),
  sessions = uniqueN(camera_session),
  first_event = min(event_time),
  last_event = max(event_time),
  median_raw_rows_per_event = as.numeric(stats::median(exact_records)),
  max_raw_rows_per_event = as.numeric(max(exact_records)),
  direct_person_events = if (taxon[1] == "Human") as.integer(sum(direct_human)) else NA_integer_
), by = taxon]
event_audit <- merge(
  event_audit,
  focal_raw[, .(raw_annotated_rows = .N), by = taxon],
  by = "taxon",
  all.x = TRUE
)
event_audit <- merge(
  event_audit,
  deduplicated[, .(distinct_timestamps = .N), by = taxon],
  by = "taxon",
  all.x = TRUE
)
event_audit[, exact_duplicate_rows_collapsed := raw_annotated_rows - distinct_timestamps]
valid_session_count <- nrow(sessions)
inconsistent_session_count <- sum(sessions$inconsistent_start | sessions$inconsistent_end)
event_audit[, `:=`(
  raw_rows_total = raw[,.N],
  valid_sessions = valid_session_count,
  inconsistent_session_dates = inconsistent_session_count,
  unclassified_gazelle_rows = unclassified_gazelle,
  sampling_frame = sampling_frame,
  event_definition = event_definition,
  event_gap_minutes = if (event_definition == "gap") event_gap_minutes else NA_real_
)]
event_audit[, taxon_order := match(taxon, taxon_levels)]
setorder(event_audit, taxon_order)
event_audit[, taxon_order := NULL]

fwrite(classification_audit, paste0(cfg$prefix, "_classification_audit.csv"))
fwrite(event_audit, paste0(cfg$prefix, "_event_audit.csv"))
fwrite(
  events[, .(
    camera, session_id, taxon, event_time, quantity, exact_records, direct_human
  )],
  paste0(
    cfg$prefix,
    if (event_definition == "gap") "_independent_events.csv" else "_timestamp_visit_events.csv"
  )
)

message("Analysis event counts:")
print(event_audit[, .(taxon, analysis_events, cameras, sessions)])
if (audit_only) quit(save = "no", status = 0)

make_recurrent_surveys <- function(event_data, primary, secondary, tertiary, max_days) {
  primary_events <- event_data[taxon == primary]
  secondary_events <- event_data[taxon == secondary]
  tertiary_events <- event_data[taxon %chin% tertiary]
  if (!nrow(primary_events)) return(list(surveys = data.table(), lags = data.table()))

  grouped_primary <- split(primary_events, primary_events$camera_session)
  grouped_secondary <- split(secondary_events$event_time, secondary_events$camera_session)
  grouped_tertiary <- split(tertiary_events$event_time, tertiary_events$camera_session)
  survey_out <- vector("list", length(grouped_primary))
  lag_out <- vector("list", length(grouped_primary))
  out_index <- 1L

  for (session_key in names(grouped_primary)) {
    p_dat <- grouped_primary[[session_key]][order(event_time)]
    p_times <- p_dat$event_time
    s_times <- sort(grouped_secondary[[session_key]])
    t_times <- sort(grouped_tertiary[[session_key]])
    if (is.null(s_times)) s_times <- as.POSIXct(character(), tz = cfg$timezone)
    if (is.null(t_times)) t_times <- as.POSIXct(character(), tz = cfg$timezone)

    local_surveys <- vector("list", length(p_times))
    local_lags <- vector("list", length(p_times))
    for (j in seq_along(p_times)) {
      start_time <- p_times[j]
      next_primary <- if (j < length(p_times)) p_times[j + 1L] else as.POSIXct(NA)
      future_tertiary <- t_times[t_times > start_time]
      next_tertiary <- if (length(future_tertiary)) future_tertiary[1] else as.POSIXct(NA)
      administrative_end <- min(
        p_dat$session_end[j],
        start_time + max_days * 86400,
        na.rm = TRUE
      )
      candidate_ends <- c(administrative_end, next_primary, next_tertiary)
      candidate_ends <- candidate_ends[!is.na(candidate_ends)]
      end_time <- min(candidate_ends)
      duration_days <- as.numeric(difftime(end_time, start_time, units = "days"))
      if (!is.finite(duration_days) || duration_days <= 1 / 1440) next

      survey_id <- paste(primary, secondary, session_key, j, sep = "__")
      censor_reason <- if (!is.na(next_primary) && end_time == next_primary) {
        "next_primary"
      } else if (!is.na(next_tertiary) && end_time == next_tertiary) {
        "tertiary"
      } else if (end_time == p_dat$session_end[j]) {
        "session_end"
      } else {
        "maximum_follow_up"
      }
      local_surveys[[j]] <- data.table(
        survey_id = survey_id,
        camera = p_dat$camera[j],
        session_id = p_dat$session_id[j],
        camera_session = session_key,
        primary = primary,
        secondary = secondary,
        survey_start = start_time,
        survey_end = end_time,
        duration_days = duration_days,
        censor_reason = censor_reason,
        period = ifelse(as.Date(start_time, tz = cfg$timezone) <= split_date, "Early", "Late")
      )

      observed_secondary <- s_times[s_times >= start_time & s_times <= end_time]
      if (length(observed_secondary)) {
        local_lags[[j]] <- data.table(
          survey_id = survey_id,
          target_time = observed_secondary,
          lag_days = as.numeric(difftime(observed_secondary, start_time, units = "days"))
        )
      }
    }
    survey_out[[out_index]] <- rbindlist(local_surveys, fill = TRUE)
    lag_out[[out_index]] <- rbindlist(local_lags, fill = TRUE)
    out_index <- out_index + 1L
  }

  list(
    surveys = rbindlist(survey_out, fill = TRUE),
    lags = rbindlist(lag_out, fill = TRUE)
  )
}

cut_hours <- unique(c(
  seq(0, 24, by = 2),
  seq(30, 72, by = 6),
  seq(84, 168, by = 12),
  seq(192, max_follow_days * 24, by = 24),
  max_follow_days * 24
))
cut_hours <- sort(cut_hours[cut_hours >= 0 & cut_hours <= max_follow_days * 24])
if (tail(cut_hours, 1) < max_follow_days * 24) {
  cut_hours <- c(cut_hours, max_follow_days * 24)
}
cut_days <- cut_hours / 24
intervals <- data.table(
  bin_id = seq_len(length(cut_days) - 1L),
  cut_start = head(cut_days, -1L),
  cut_end = tail(cut_days, -1L)
)

make_ped <- function(surveys, lags) {
  if (!nrow(surveys)) return(data.table())
  n_intervals <- nrow(intervals)
  ped <- surveys[rep(seq_len(.N), each = n_intervals)]
  ped[, bin_id := rep(seq_len(n_intervals), times = nrow(surveys))]
  ped <- intervals[ped, on = "bin_id"]
  ped <- ped[cut_start < duration_days]
  ped[, `:=`(
    tstart = cut_start,
    tend = pmin(cut_end, duration_days)
  )]
  ped[, exposure_days := tend - tstart]
  ped <- ped[exposure_days > 0]
  ped[, event_count := 0L]

  if (nrow(lags)) {
    lags <- copy(lags)
    lags[, bin_id := pmin(
      findInterval(lag_days, cut_days, all.inside = TRUE),
      n_intervals
    )]
    lag_counts <- lags[, .(event_count = .N), by = .(survey_id, bin_id)]
    ped[lag_counts, on = .(survey_id, bin_id), event_count := i.event_count]
  }

  ped[, midpoint_days := (tstart + tend) / 2]
  ped[, midpoint_time := survey_start + midpoint_days * 86400]
  ped[, clock_hour := as.numeric(format(midpoint_time, "%H", tz = cfg$timezone)) +
    as.numeric(format(midpoint_time, "%M", tz = cfg$timezone)) / 60]
  ped[, day_of_year := as.numeric(format(midpoint_time, "%j", tz = cfg$timezone))]
  ped[, log_elapsed := log1p(midpoint_days * 24)]
  ped[, `:=`(
    camera = factor(camera),
    camera_session = factor(camera_session),
    period = factor(period, levels = c("Early", "Late"))
  )]
  ped[]
}

fit_pamm <- function(ped) {
  if (nrow(ped) < 100 || sum(ped$event_count) < 20 || uniqueN(ped$camera) < 3) {
    return(NULL)
  }
  k_elapsed <- min(8L, max(4L, uniqueN(round(ped$log_elapsed, 5)) - 1L))
  analysis_family <- if (variance_family == "quasipoisson") {
    quasipoisson(link = "log")
  } else {
    poisson(link = "log")
  }
  tryCatch(
    bam(
      event_count ~ s(log_elapsed, k = k_elapsed, bs = "tp") +
        s(clock_hour, k = 8, bs = "cc") +
        s(day_of_year, k = 10, bs = "cc") +
        s(camera, bs = "re") +
        s(camera_session, bs = "re"),
      data = ped,
      family = analysis_family,
      offset = log(exposure_days),
      knots = list(clock_hour = c(0, 24), day_of_year = c(0.5, 366.5)),
      method = "fREML",
      discrete = TRUE,
      nthreads = 2
    ),
    error = function(e) {
      warning("Full PAMM failed; fitting camera-only frailty: ", conditionMessage(e))
      tryCatch(
        bam(
          event_count ~ s(log_elapsed, k = k_elapsed, bs = "tp") +
            s(clock_hour, k = 8, bs = "cc") +
            s(day_of_year, k = 10, bs = "cc") +
            s(camera, bs = "re"),
          data = ped,
          family = analysis_family,
          offset = log(exposure_days),
          knots = list(clock_hour = c(0, 24), day_of_year = c(0.5, 366.5)),
          method = "fREML",
          discrete = TRUE,
          nthreads = 2
        ),
        error = function(e2) {
          warning("Simplified PAMM failed: ", conditionMessage(e2))
          NULL
        }
      )
    }
  )
}

prediction_matrix <- function(model, elapsed_hours, template) {
  newdata <- data.frame(
    period = factor("Early", levels = levels(template$period)),
    log_elapsed = log1p(elapsed_hours),
    clock_hour = 12,
    day_of_year = 183,
    camera = factor(levels(template$camera)[1], levels = levels(template$camera)),
    camera_session = factor(
      levels(template$camera_session)[1],
      levels = levels(template$camera_session)
    )
  )
  random_terms <- intersect(c("s(camera)", "s(camera_session)"), model$smooth |> vapply(
    function(x) x$label,
    character(1)
  ))
  predict(model, newdata = newdata, type = "lpmatrix", exclude = random_terms)
}

make_curve_and_contrasts <- function(model, ped, primary, secondary, direction) {
  if (is.null(model)) return(list(curve = data.table(), contrasts = data.table(), smooth = data.table()))
  support <- ped[, .(
    exposure_days = sum(exposure_days),
    risk_surveys = uniqueN(survey_id),
    risk_cameras = uniqueN(camera),
    supported_end_days = max(tend)
  ), by = bin_id]
  support <- support[
    risk_surveys >= max(30, ceiling(0.05 * uniqueN(ped$survey_id))) &
      risk_cameras >= max(5, ceiling(0.20 * uniqueN(ped$camera))) &
      exposure_days >= 1
  ]
  max_supported_days <- if (nrow(support)) max(support$supported_end_days) else {
    safe_quantile(ped$tend, 0.90)
  }
  max_supported_days <- min(max_follow_days, max(max_supported_days, 0.5))
  ref_end_hours <- max_supported_days * 24
  ref_start_hours <- max(6, 0.75 * ref_end_hours)
  curve_hours <- unique(c(
    seq(0.25, 24, length.out = 80),
    seq(25, ref_end_hours, length.out = 120)
  ))
  curve_hours <- sort(curve_hours[curve_hours <= ref_end_hours])
  ref_hours <- seq(ref_start_hours, ref_end_hours, length.out = 60)
  x_curve <- prediction_matrix(model, curve_hours, ped)
  x_ref <- prediction_matrix(model, ref_hours, ped)
  ref_vector <- colMeans(x_ref)
  contrast_matrix <- sweep(x_curve, 2, ref_vector, "-")
  beta <- coef(model)
  vc <- vcov(model, unconditional = TRUE)
  estimate <- as.numeric(contrast_matrix %*% beta)
  variance <- rowSums((contrast_matrix %*% vc) * contrast_matrix)
  se <- sqrt(pmax(variance, 0))
  curve <- data.table(
    primary = primary,
    secondary = secondary,
    direction = direction,
    elapsed_hours = curve_hours,
    hazard_ratio = exp(estimate),
    lwr = exp(estimate - 1.96 * se),
    upr = exp(estimate + 1.96 * se),
    reference_start_hours = ref_start_hours,
    reference_end_hours = ref_end_hours,
    reference_window = paste0(
      formatC(ref_start_hours / 24, digits = 1, format = "f"), "-",
      formatC(ref_end_hours / 24, digits = 1, format = "f"), " days"
    )
  )

  windows <- list(
    "0-6 hours" = seq(0.25, 6, length.out = 30),
    "6-24 hours" = seq(6.25, 24, length.out = 40),
    "1-3 days" = seq(25, 72, length.out = 50),
    "3-7 days" = seq(73, min(168, max_follow_days * 24), length.out = 50)
  )
  contrasts <- rbindlist(lapply(names(windows), function(window_name) {
    hours <- windows[[window_name]]
    hours <- hours[hours <= ref_end_hours]
    if (!length(hours)) return(NULL)
    x_window <- colMeans(prediction_matrix(model, hours, ped))
    contrast <- x_window - ref_vector
    est <- as.numeric(contrast %*% beta)
    est_se <- sqrt(max(as.numeric(contrast %*% vc %*% contrast), 0))
    data.table(
      primary = primary,
      secondary = secondary,
      direction = direction,
      window = window_name,
      log_hazard_ratio = est,
      std_error = est_se,
      hazard_ratio = exp(est),
      lwr = exp(est - 1.96 * est_se),
      upr = exp(est + 1.96 * est_se),
      reference_start_hours = ref_start_hours,
      reference_end_hours = ref_end_hours
    )
  }), fill = TRUE)

  smry <- summary(model)
  smooth_row <- grep("^s\\(log_elapsed\\)", rownames(smry$s.table))
  smooth <- if (length(smooth_row)) {
    data.table(
      primary = primary,
      secondary = secondary,
      direction = direction,
      edf = smry$s.table[smooth_row[1], "edf"],
      statistic = smry$s.table[smooth_row[1], ncol(smry$s.table) - 1L],
      p_value = pmax(0, pmin(1, smry$s.table[smooth_row[1], ncol(smry$s.table)]))
    )
  } else {
    data.table(
      primary = primary, secondary = secondary, direction = direction,
      edf = NA_real_, statistic = NA_real_, p_value = NA_real_
    )
  }
  list(curve = curve, contrasts = contrasts, smooth = smooth)
}

lag_summary <- function(surveys, lags, primary, secondary, direction) {
  first_lags <- if (nrow(lags)) lags[, .(first_lag_days = min(lag_days)), by = survey_id] else data.table()
  survey_first <- merge(surveys, first_lags, by = "survey_id", all.x = TRUE)
  observed <- survey_first[!is.na(first_lag_days), first_lag_days]
  data.table(
    primary = primary,
    secondary = secondary,
    direction = direction,
    n_primary_surveys = nrow(surveys),
    n_cameras = uniqueN(surveys$camera),
    recurrent_secondary_events = nrow(lags),
    surveys_with_secondary = sum(!is.na(survey_first$first_lag_days)),
    proportion_within_6h = mean(!is.na(survey_first$first_lag_days) & survey_first$first_lag_days <= 0.25),
    proportion_within_24h = mean(!is.na(survey_first$first_lag_days) & survey_first$first_lag_days <= 1),
    median_first_lag_hours_observed = if (length(observed)) stats::median(observed) * 24 else NA_real_,
    q25_first_lag_hours_observed = safe_quantile(observed, 0.25) * 24,
    q75_first_lag_hours_observed = safe_quantile(observed, 0.75) * 24,
    median_follow_up_days = stats::median(surveys$duration_days),
    censored_by_next_primary = sum(surveys$censor_reason == "next_primary"),
    censored_by_tertiary = sum(surveys$censor_reason == "tertiary")
  )
}

prepare_shift_groups <- function(target_events, session_table) {
  target <- merge(
    target_events[, .(camera_session, event_time)],
    session_table[, .(camera_session, session_start, session_end)],
    by = "camera_session",
    all.x = TRUE
  )
  target[, duration_days_floor := pmax(1L, floor(as.numeric(difftime(
    session_end, session_start, units = "days"
  ))))]
  target[, elapsed_seconds := as.numeric(difftime(event_time, session_start, units = "secs"))]
  target[, day_index := pmin(floor(elapsed_seconds / 86400), duration_days_floor - 1L)]
  target[, time_of_day_seconds := elapsed_seconds - floor(elapsed_seconds / 86400) * 86400]
  split(target, target$camera_session)
}

assign_lags_to_surveys <- function(time_groups, survey_groups) {
  all_lags <- numeric()
  for (session_key in intersect(names(time_groups), names(survey_groups))) {
    times <- sort(time_groups[[session_key]])
    surveys <- survey_groups[[session_key]][order(survey_start)]
    starts <- as.numeric(surveys$survey_start)
    ends <- as.numeric(surveys$survey_end)
    indices <- findInterval(as.numeric(times), starts)
    valid <- indices > 0L
    indices_valid <- indices[valid]
    times_valid <- as.numeric(times[valid])
    valid_end <- times_valid <= ends[indices_valid]
    if (any(valid_end)) {
      all_lags <- c(
        all_lags,
        (times_valid[valid_end] - starts[indices_valid[valid_end]]) / 86400
      )
    }
  }
  all_lags
}

permutation_test <- function(surveys, target_events, observed_lags, reps) {
  if (!nrow(surveys) || !nrow(target_events) || reps < 1) {
    return(list(test = data.table(), profile = data.table()))
  }
  survey_groups <- split(surveys, surveys$camera_session)
  target_groups <- prepare_shift_groups(target_events, sessions)
  observed_lag_values <- observed_lags$lag_days
  early_exposure <- sum(pmax(0, pmin(surveys$duration_days, 1)))
  late_exposure <- sum(pmax(0, pmin(surveys$duration_days, 7) - 1))
  if (early_exposure <= 0 || late_exposure <= 0) {
    return(list(test = data.table(), profile = data.table()))
  }

  calc_stat <- function(lag_values) {
    early_count <- sum(lag_values >= 0 & lag_values <= 1)
    late_count <- sum(lag_values > 1 & lag_values <= 7)
    log(((early_count + 0.5) / early_exposure) / ((late_count + 0.5) / late_exposure))
  }

  bin_exposure <- surveys[, {
    pmax(0, pmin(duration_days, intervals$cut_end) - intervals$cut_start)
  }, by = survey_id]$V1
  bin_exposure <- colSums(matrix(bin_exposure, ncol = nrow(intervals), byrow = TRUE))
  bin_exposure <- pmax(bin_exposure, 1e-10)

  profile_ratio <- function(lag_values) {
    bins <- pmin(findInterval(lag_values, cut_days, all.inside = TRUE), nrow(intervals))
    counts <- tabulate(bins, nbins = nrow(intervals))
    rates <- (counts + 0.5) / bin_exposure
    ref_bins <- intervals$cut_start >= min(7, max_follow_days / 2)
    rates / exp(mean(log(rates[ref_bins])))
  }

  observed_stat <- calc_stat(observed_lag_values)
  observed_profile <- profile_ratio(observed_lag_values)
  permutation_stats <- numeric(reps)
  permutation_profiles <- matrix(NA_real_, nrow = reps, ncol = nrow(intervals))

  for (b in seq_len(reps)) {
    shifted_times <- lapply(target_groups, function(group_data) {
      n_days <- unique(group_data$duration_days_floor)[1]
      offset <- if (n_days > 1) sample.int(n_days - 1L, 1L) else 0L
      shifted_day <- (group_data$day_index + offset) %% n_days
      group_data$session_start + shifted_day * 86400 + group_data$time_of_day_seconds
    })
    shifted_lags <- assign_lags_to_surveys(shifted_times, survey_groups)
    permutation_stats[b] <- calc_stat(shifted_lags)
    permutation_profiles[b, ] <- profile_ratio(shifted_lags)
  }

  null_center <- stats::median(permutation_stats, na.rm = TRUE)
  deviation <- observed_stat - null_center
  p_value <- (1 + sum(abs(permutation_stats - null_center) >= abs(deviation), na.rm = TRUE)) /
    (1 + sum(is.finite(permutation_stats)))
  test <- data.table(
    observed_log_short_long_ratio = observed_stat,
    null_median_log_short_long_ratio = null_center,
    interaction_log_ratio = deviation,
    interaction_ratio = exp(deviation),
    null_adjusted_lwr = exp(observed_stat - safe_quantile(permutation_stats, 0.975)),
    null_adjusted_upr = exp(observed_stat - safe_quantile(permutation_stats, 0.025)),
    permutation_p_value = p_value,
    observed_events_0_1d = sum(observed_lag_values >= 0 & observed_lag_values <= 1),
    observed_events_1_7d = sum(observed_lag_values > 1 & observed_lag_values <= 7),
    exposure_0_1d = early_exposure,
    exposure_1_7d = late_exposure,
    permutation_reps = reps
  )
  profile <- data.table(
    bin_id = intervals$bin_id,
    elapsed_hours = (intervals$cut_start + intervals$cut_end) / 2 * 24,
    observed_relative_rate = observed_profile,
    null_median = apply(permutation_profiles, 2, stats::median, na.rm = TRUE),
    null_lwr = apply(permutation_profiles, 2, safe_quantile, prob = 0.025),
    null_upr = apply(permutation_profiles, 2, safe_quantile, prob = 0.975),
    exposure_days = bin_exposure
  )
  list(test = test, profile = profile)
}

all_lag_summaries <- list()
all_smooth_tests <- list()
all_contrasts <- list()
all_curves <- list()
all_permutation_tests <- list()
all_null_profiles <- list()
all_diagnostics <- list()
model_fits <- list()
result_index <- 1L

dyads <- rbindlist(list(
  CJ(pressure = pressure_levels, class = class_levels)[
    , .(primary = pressure, secondary = class, direction = "Pressure to gazelle")
  ],
  CJ(pressure = pressure_levels, class = class_levels)[
    , .(primary = class, secondary = pressure, direction = "Gazelle to pressure")
  ]
))

for (i in seq_len(nrow(dyads))) {
  primary <- dyads$primary[i]
  secondary <- dyads$secondary[i]
  direction <- dyads$direction[i]
  focal_pressure <- if (primary %chin% pressure_levels) primary else secondary
  tertiary <- if (censor_competing_pressures) setdiff(pressure_levels, focal_pressure) else character()
  message(sprintf("[%d/%d] %s -> %s", i, nrow(dyads), primary, secondary))

  recurrent <- make_recurrent_surveys(
    events,
    primary = primary,
    secondary = secondary,
    tertiary = tertiary,
    max_days = max_follow_days
  )
  surveys_i <- recurrent$surveys
  lags_i <- recurrent$lags
  if (!nrow(surveys_i)) next
  ped_i <- make_ped(surveys_i, lags_i)
  model_i <- fit_pamm(ped_i)
  curve_output <- make_curve_and_contrasts(
    model_i, ped_i, primary, secondary, direction
  )
  permutation_output <- permutation_test(
    surveys_i,
    events[taxon == secondary],
    lags_i,
    permutation_reps
  )

  id_columns <- list(
    primary = primary,
    secondary = secondary,
    direction = direction
  )
  all_lag_summaries[[result_index]] <- lag_summary(
    surveys_i, lags_i, primary, secondary, direction
  )
  all_smooth_tests[[result_index]] <- curve_output$smooth
  all_contrasts[[result_index]] <- curve_output$contrasts
  all_curves[[result_index]] <- curve_output$curve
  all_permutation_tests[[result_index]] <- cbind(
    as.data.table(id_columns),
    permutation_output$test
  )
  all_null_profiles[[result_index]] <- cbind(
    as.data.table(id_columns),
    permutation_output$profile
  )

  if (!is.null(model_i)) {
    pearson_dispersion <- sum(residuals(model_i, type = "pearson")^2, na.rm = TRUE) /
      model_i$df.residual
    all_diagnostics[[result_index]] <- data.table(
      primary = primary,
      secondary = secondary,
      direction = direction,
      ped_rows = nrow(ped_i),
      events = sum(ped_i$event_count),
      cameras = uniqueN(ped_i$camera),
      deviance_explained = summary(model_i)$dev.expl,
      pearson_dispersion = pearson_dispersion,
      variance_family = variance_family,
      converged = isTRUE(model_i$converged),
      full_session_frailty = "s(camera_session)" %chin% vapply(model_i$smooth, function(x) x$label, character(1))
    )
    if (save_models) {
      model_fits[[paste(primary, secondary, sep = "__")]] <- model_i
    }
  }
  result_index <- result_index + 1L
}

lag_summaries <- rbindlist(all_lag_summaries, fill = TRUE)
smooth_tests <- rbindlist(all_smooth_tests, fill = TRUE)
window_contrasts <- rbindlist(all_contrasts, fill = TRUE)
hazard_curves <- rbindlist(all_curves, fill = TRUE)
permutation_tests <- rbindlist(all_permutation_tests, fill = TRUE)
null_profiles <- rbindlist(all_null_profiles, fill = TRUE)
diagnostics <- rbindlist(all_diagnostics, fill = TRUE)

if (nrow(permutation_tests)) {
  permutation_tests[, p_bh_global := stats::p.adjust(permutation_p_value, method = "BH")]
  permutation_tests[, p_bh_direction := stats::p.adjust(permutation_p_value, method = "BH"), by = direction]
}

output_tables <- c(
  lag_summaries = "lag_summary",
  smooth_tests = "pamm_smooth_tests",
  window_contrasts = "pamm_window_contrasts",
  hazard_curves = "pamm_hazard_curves",
  permutation_tests = "permutation_tests",
  null_profiles = "permutation_null_profiles",
  diagnostics = "model_diagnostics"
)
for (object_name in names(output_tables)) {
  fwrite(get(object_name), paste0(cfg$prefix, "_", output_tables[[object_name]], ".csv"))
}
if (save_models) saveRDS(model_fits, paste0(cfg$prefix, "_model_fits.rds"))

plot_curves <- function(direction_value, color_variable, palette, title, subtitle, filename_suffix) {
  plot_data <- hazard_curves[direction == direction_value]
  if (!nrow(plot_data)) return(NULL)
  if (direction_value == "Pressure to gazelle") {
    plot_data[, `:=`(
      pressure = factor(primary, levels = pressure_levels),
      gazelle_class = factor(secondary, levels = class_levels),
      color_group = secondary
    )]
    facet_formula <- pressure ~ gazelle_class
  } else {
    plot_data[, `:=`(
      pressure = factor(secondary, levels = pressure_levels),
      gazelle_class = factor(primary, levels = class_levels),
      color_group = secondary
    )]
    facet_formula <- pressure ~ gazelle_class
  }
  p <- ggplot(plot_data, aes(elapsed_hours / 24, hazard_ratio, color = color_group, fill = color_group)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.16, color = NA) +
    geom_line(linewidth = 0.8) +
    facet_grid(facet_formula, scales = "free_y") +
    scale_x_continuous(breaks = c(0, 1, 3, 7, 14)) +
    scale_color_manual(values = palette, guide = "none") +
    scale_fill_manual(values = palette, guide = "none") +
    coord_cartesian(ylim = c(0, quantile(plot_data$upr[is.finite(plot_data$upr)], 0.98, na.rm = TRUE))) +
    labs(
      x = "Days since primary event",
      y = "Adjusted visitation hazard ratio",
      title = title,
      subtitle = subtitle
    ) +
    theme_publication()
  ggsave(paste0(cfg$prefix, filename_suffix, ".png"), p, width = 11.5, height = 8.2, dpi = 320)
  ggsave(paste0(cfg$prefix, filename_suffix, ".pdf"), p, width = 11.5, height = 8.2)
  p
}

plot_curves(
  "Pressure to gazelle",
  "gazelle_class",
  class_colors,
  paste0("Gazelle visitation after ", paste(tolower(pressure_levels), collapse = ", "), " events"),
  "PAMM hazard ratios use each dyad's late data-supported interval; bands are 95% confidence intervals",
  "_gazelle_response_curves"
)
plot_curves(
  "Gazelle to pressure",
  "pressure",
  pressure_colors,
  paste0(paste(pressure_levels, collapse = ", "), " visitation after gazelle events"),
  "Reverse-direction PAMM hazard ratios distinguish temporal order from symmetric co-use",
  "_pressure_response_curves"
)

if (nrow(permutation_tests)) {
  permutation_plot_data <- copy(permutation_tests)
  permutation_plot_data[, pressure := ifelse(primary %chin% pressure_levels, primary, secondary)]
  permutation_plot_data[, gazelle_class := ifelse(primary %chin% class_levels, primary, secondary)]
  permutation_plot_data[, pair := factor(
    paste(pressure, gazelle_class, sep = " / "),
    levels = rev(as.vector(outer(pressure_levels, class_levels, paste, sep = " / ")))
  )]
  permutation_plot_data[, direction := factor(direction, levels = direction_levels)]
  p_permutation <- ggplot(
    permutation_plot_data,
    aes(interaction_ratio, pair, xmin = null_adjusted_lwr, xmax = null_adjusted_upr, color = direction)
  ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey35") +
    geom_errorbar(orientation = "y", width = 0.16, position = position_dodge(width = 0.46)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.46)) +
    scale_x_log10() +
    scale_color_manual(values = c("Pressure to gazelle" = "#b33a3a", "Gazelle to pressure" = "#315f8c")) +
    labs(
      x = "Observed 0-1 day vs 1-7 day rate ratio, relative to randomized timing",
      y = NULL,
      color = "Direction",
      title = "Directional short-lag associations after preserving camera and time of day",
      subtitle = "Intervals are 95% randomization envelopes, not model-based confidence intervals"
    ) +
    theme_publication(10.5)
  ggsave(paste0(cfg$prefix, "_permutation_effects.png"), p_permutation, width = 10.5, height = 8.2, dpi = 320)
  ggsave(paste0(cfg$prefix, "_permutation_effects.pdf"), p_permutation, width = 10.5, height = 8.2)
}

main_results <- merge(
  window_contrasts[window == "0-6 hours"],
  permutation_tests,
  by = c("primary", "secondary", "direction"),
  all = TRUE
)
main_results[, `:=`(
  direction_order = match(direction, direction_levels),
  primary_order = match(primary, taxon_levels),
  secondary_order = match(secondary, taxon_levels)
)]
setorder(main_results, direction_order, primary_order, secondary_order)
main_results[, c("direction_order", "primary_order", "secondary_order") := NULL]

event_description <- if (event_definition == "gap") {
  paste0(
    "Exact duplicate export rows were collapsed, then consecutive records within ",
    event_gap_minutes,
    " min at one camera-session-taxon were grouped as one independent event."
  )
} else {
  "Exact duplicate export rows at one camera-session-taxon timestamp were collapsed, but every distinct timestamp was retained; no independence interval was imposed."
}
event_table_label <- if (event_definition == "gap") "Independent events" else "Distinct timestamp visits"
variance_description <- if (variance_family == "quasipoisson") {
  "quasi-Poisson uncertainty accounts for extra-Poisson variation."
} else {
  "Poisson piecewise-event rates were modeled with camera and camera-session frailties."
}
pressure_notes <- c(
  paste0("- Focal pressure set: `", paste(pressure_levels, collapse = ", "), "`."),
  if ("Human" %chin% pressure_levels) "- Human-associated activity includes people, soldiers, bicycles, and civilian, military, and construction vehicles, matching the earlier human-pressure covariate." else NULL,
  if ("Dog" %chin% pressure_levels) "- Dog activity is the export category `Canis familiaris`; the available metadata do not distinguish livestock-guarding, accompanied, roaming, or feral dogs." else NULL
)

results_lines <- c(
  paste0("# Directional gazelle-", paste(tolower(pressure_levels), collapse = "-"), " recurrent-event analysis"),
  "",
  "## Analysis",
  "",
  paste0("- Input: `", cfg$full_export, "`."),
  paste0("- Sampling frame: `", sampling_frame, "`, from ", survey_start_date, " to ", survey_end_date, "."),
  paste0("- Event definition: `", event_definition, "`. ", event_description),
  paste0("- Variance family: `", variance_family, "`; ", variance_description),
  paste0("- Maximum follow-up after each primary event: ", max_follow_days, " days, additionally censored by the next primary event or Camelot camera-session end."),
  paste0("- Competing-pressure censoring: `", censor_competing_pressures, "`; the default dyadic analysis does not discard follow-up when another pressure taxon is detected."),
  pressure_notes,
  "- Primary model: directional recurrent-event piecewise exponential additive mixed model (PAMM), with smooth effects of elapsed time, local time of day, and day of year, plus camera and camera-session frailties.",
  paste0("- Randomization test: ", permutation_reps, " integer-day circular shifts of secondary events within camera sessions, preserving camera use, event totals, and clock time."),
  "- The randomization effect compares the 0-1 day rate with the 1-7 day rate relative to neutral shifted schedules. Ratios below 1 indicate fewer short lags than expected; ratios above 1 indicate more short lags.",
  paste0("- Randomization p-values were adjusted by the Benjamini-Hochberg method across all ", nrow(dyads), " directional tests and, secondarily, within each ", nrow(dyads) / length(direction_levels), "-test directional family."),
  "",
  "## Event audit",
  "",
  paste0("| Taxon | ", event_table_label, " | Cameras | Sessions |"),
  "|---|---:|---:|---:|"
)
for (i in seq_len(nrow(event_audit))) {
  results_lines <- c(results_lines, sprintf(
    "| %s | %d | %d | %d |",
    event_audit$taxon[i], event_audit$analysis_events[i],
    event_audit$cameras[i], event_audit$sessions[i]
  ))
}

results_lines <- c(
  results_lines,
  "",
  "## Short-lag results",
  "",
  "| Direction | Primary | Secondary | PAMM HR, 0-6 h vs late supported interval | 95% CI | Reference days | Randomization ratio | p | BH q, global | BH q, direction |",
  "|---|---|---|---:|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(main_results))) {
  results_lines <- c(results_lines, sprintf(
    "| %s | %s | %s | %s | [%s, %s] | %s-%s | %s | %s | %s | %s |",
    main_results$direction[i], main_results$primary[i], main_results$secondary[i],
    fmt_num(main_results$hazard_ratio[i]), fmt_num(main_results$lwr[i]), fmt_num(main_results$upr[i]),
    fmt_num(main_results$reference_start_hours[i] / 24, 1),
    fmt_num(main_results$reference_end_hours[i] / 24, 1),
    fmt_num(main_results$interaction_ratio[i]), fmt_p(main_results$permutation_p_value[i]),
    fmt_p(main_results$p_bh_global[i]), fmt_p(main_results$p_bh_direction[i])
  ))
}

global_discoveries <- main_results[p_bh_global < 0.05]
direction_discoveries <- main_results[p_bh_direction < 0.05]
global_labels <- if (nrow(global_discoveries)) {
  paste(paste0(global_discoveries$primary, " to ", global_discoveries$secondary), collapse = ", ")
} else {
  "none"
}
direction_labels <- if (nrow(direction_discoveries)) {
  paste(paste0(direction_discoveries$primary, " to ", direction_discoveries$secondary), collapse = ", ")
} else {
  "none"
}
global_inference <- if (nrow(global_discoveries)) {
  paste0(
    "At a 5% false-discovery-rate threshold across all ", nrow(main_results),
    " tests, the retained associations were: ", global_labels, "."
  )
} else {
  paste0(
    "No association met a 5% false-discovery-rate threshold across all ",
    nrow(main_results), " tests."
  )
}
interaction_guardrail <- if (event_definition == "gap") {
  paste0(
    "The PAMM curves describe directional changes in rates of ", event_gap_minutes,
    "-min independent camera events after another taxon was recorded. The camera-session randomization retains spatial intensity, diel timing, and within-taxon temporal clustering. Neither analysis demonstrates individual-level following, predation, protection, or causal behavioral avoidance."
  )
} else {
  "The PAMM curves describe directional changes in camera trigger rates after another taxon was recorded. Because consecutive timestamps may represent the same individual or group, this timestamp-level analysis estimates detection-sequence association rather than independent animal encounters. The camera-session randomization is the stronger test of fine-scale temporal association because it retains spatial intensity, diel timing, and within-taxon temporal clustering. Neither analysis demonstrates predation, individual-level pursuit, or causal behavioral avoidance."
}

results_lines <- c(
  results_lines,
  "",
  "## Multiplicity-aware inference",
  "",
  paste0(global_inference, " Within the two prespecified directional families, the 5% threshold retained: ", direction_labels, ". Direction-family adjustment is secondary and should not replace the global result."),
  "",
  "## Interpretation guardrails",
  "",
  interaction_guardrail,
  "",
  "The simple median observed lag is conditional on a secondary event occurring before censoring and must not be interpreted as an unbiased waiting time. It is included only as a descriptive quantity; the recurrent-event PAMM is the inferential analysis.",
  "",
  "## Key outputs",
  "",
  paste0("- `", cfg$prefix, "_gazelle_response_curves.png`"),
  paste0("- `", cfg$prefix, "_pressure_response_curves.png`"),
  paste0("- `", cfg$prefix, "_permutation_effects.png`"),
  paste0("- `", cfg$prefix, "_lag_summary.csv`"),
  paste0("- `", cfg$prefix, "_pamm_window_contrasts.csv`"),
  paste0("- `", cfg$prefix, "_permutation_tests.csv`"),
  "",
  "## Method references",
  "",
  "Ferry, N., Dupont, P., Bender, A., & Heurich, M. (2024). Introducing recurrent event analyses to assess species interactions based on camera-trap data: A comparison with time-to-first-event approaches. *Methods in Ecology and Evolution*, 15, 1233-1246. https://doi.org/10.1111/2041-210X.14338",
  "",
  "Niedballa, J., et al. (2019). Assessing analytical methods for detecting spatiotemporal interactions between species from camera trapping data. *Remote Sensing in Ecology and Conservation*. https://doi.org/10.1002/rse2.107"
)
writeLines(results_lines, paste0(cfg$prefix, "_results.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = paste0(cfg$prefix, "_session_info.txt"))
message("Analysis complete. Prefix: ", cfg$prefix)
