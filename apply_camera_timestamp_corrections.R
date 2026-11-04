#!/usr/bin/env Rscript

# Apply documented camera-clock date corrections to the prepared independent
# events while preserving the original export and local clock time of day.

options(stringsAsFactors = FALSE)

source_file <- "Gazella gazella_independent_events_30min.csv"
correction_file <- "camera_timestamp_corrections.csv"
obs_file <- "ObsCovs_050524.csv"
output_file <- "Gazella gazella_independent_events_30min_timestamp_corrected.csv"
audit_file <- "publication_camera_timestamp_correction_audit.csv"

required <- c(source_file, correction_file, obs_file)
missing_required <- required[!file.exists(required)]
if (length(missing_required)) {
  stop("Missing correction inputs: ", paste(missing_required, collapse = ", "))
}

events <- read.csv(source_file, check.names = FALSE)
corrections <- read.csv(correction_file, check.names = FALSE)
obs <- read.csv(obs_file, na.strings = c("NA", ""), check.names = FALSE)

event_key <- paste(events$camera_name, events$trap_station_session_id, sep = ":")
correction_key <- paste(
  corrections$camera_name,
  corrections$trap_station_session_id,
  sep = ":"
)
if (anyDuplicated(correction_key)) stop("Duplicate timestamp-correction keys.")

correction_index <- match(event_key, correction_key)
affected <- !is.na(correction_index)
missing_sessions <- setdiff(correction_key, unique(event_key))
if (length(missing_sessions)) {
  stop("Correction sessions absent from event file: ", paste(missing_sessions, collapse = ", "))
}

original_time <- events$event_time
recorded_date <- as.Date(substr(original_time, 1, 10), format = "%Y-%m-%d")
if (anyNA(recorded_date)) stop("One or more event dates could not be parsed.")

offset_days <- integer(nrow(events))
offset_days[affected] <- corrections$offset_days[correction_index[affected]]
corrected_date <- recorded_date + offset_days

# Camera clocks had the wrong date but valid local time of day. Reconstructing
# the timestamp string avoids a daylight-saving shift in the recorded hour.
events$event_time[affected] <- paste0(
  format(corrected_date[affected], "%Y-%m-%d"),
  substr(original_time[affected], 11, nchar(original_time[affected]))
)

obs_check <- data.frame(
  camera_name = as.integer(obs[["Camera Number"]]),
  round = as.integer(obs[["Round Number"]]),
  start_date = as.Date(obs[["Start date"]], format = "%d/%m/%Y"),
  end_date = as.Date(obs[["End Date"]], format = "%d/%m/%Y"),
  stringsAsFactors = FALSE
)

audit <- data.frame(
  event_row = which(affected),
  camera_name = events$camera_name[affected],
  trap_station_session_id = events$trap_station_session_id[affected],
  gazelle_class = events$gazelle_class[affected],
  event_id = events$event_id[affected],
  event_count = events$event_count[affected],
  original_event_time = original_time[affected],
  corrected_event_time = events$event_time[affected],
  offset_days = offset_days[affected],
  expected_round = corrections$expected_round[correction_index[affected]],
  reason = corrections$reason[correction_index[affected]],
  stringsAsFactors = FALSE
)

audit$within_expected_round <- mapply(
  function(camera, round, date) {
    idx <- obs_check$camera_name == camera & obs_check$round == round
    if (sum(idx) != 1L) return(FALSE)
    !is.na(obs_check$start_date[idx]) &&
      !is.na(obs_check$end_date[idx]) &&
      date >= obs_check$start_date[idx] &&
      date <= obs_check$end_date[idx]
  },
  camera = audit$camera_name,
  round = audit$expected_round,
  date = corrected_date[affected]
)
if (!all(audit$within_expected_round)) {
  bad <- audit[!audit$within_expected_round, ]
  stop(
    "Corrected events outside expected rounds: ",
    paste(paste(bad$camera_name, bad$trap_station_session_id, sep = ":"), collapse = ", ")
  )
}

write.csv(events, output_file, row.names = FALSE, na = "")
write.csv(audit, audit_file, row.names = FALSE, na = "")

message(
  "Corrected ", sum(affected), " independent-event rows across ",
  length(unique(event_key[affected])), " camera sessions."
)
