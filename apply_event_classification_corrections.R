#!/usr/bin/env Rscript

# Remove documented juvenile records that were duplicated into sex-based adult
# classes. Juvenile copies are retained; only the duplicate adult-class rows
# are excluded from the analysis-ready independent-event file.

options(stringsAsFactors = FALSE)

source_file <- "Gazella gazella_independent_events_30min_timestamp_corrected.csv"
exclusion_file <- "event_classification_exclusions.csv"
output_file <- "Gazella gazella_independent_events_30min_analysis_ready.csv"
audit_file <- "publication_event_classification_correction_audit.csv"
count_files <- c(
  juvenile = "gazelle_juvenile_abundance_60x7.csv",
  female = "gazelle_female_abundance_60x7.csv",
  male = "gazelle_male_abundance_60x7.csv"
)

required <- c(source_file, exclusion_file, count_files)
missing_required <- required[!file.exists(required)]
if (length(missing_required)) {
  stop("Missing classification-correction inputs: ", paste(missing_required, collapse = ", "))
}

events <- read.csv(source_file, check.names = FALSE)
exclusions <- read.csv(exclusion_file, check.names = FALSE)

make_key <- function(camera, session, class, time) {
  paste(camera, session, class, time, sep = ":")
}
event_key <- make_key(
  events$camera_name,
  events$trap_station_session_id,
  events$gazelle_class,
  events$event_time
)
exclusion_key <- make_key(
  exclusions$camera_name,
  exclusions$trap_station_session_id,
  exclusions$gazelle_class,
  exclusions$event_time
)
if (anyDuplicated(exclusion_key)) stop("Duplicate event-classification exclusion keys.")

match_count <- vapply(exclusion_key, function(key) sum(event_key == key), integer(1))
if (any(match_count != 1L)) {
  bad <- exclusion_key[match_count != 1L]
  stop("Each classification exclusion must match exactly one event: ", paste(bad, collapse = ", "))
}

excluded <- event_key %in% exclusion_key
audit <- merge(
  events[excluded, ],
  exclusions,
  by = c("camera_name", "trap_station_session_id", "gazelle_class", "event_time"),
  all.x = TRUE,
  sort = FALSE
)
audit$action <- "Excluded duplicate adult-class row; juvenile row retained"
events <- events[!excluded, ]

read_count_total <- function(path) {
  x <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  sum(as.matrix(x[paste0("X", 1:7)]), na.rm = TRUE)
}
expected_totals <- vapply(count_files, read_count_total, numeric(1))
observed_totals <- tapply(events$event_count, events$gazelle_class, sum, na.rm = TRUE)
observed_totals <- observed_totals[names(expected_totals)]
if (anyNA(observed_totals) || !identical(as.numeric(observed_totals), as.numeric(expected_totals))) {
  stop(
    "Analysis-ready event totals do not match count-history totals. Expected ",
    paste(names(expected_totals), expected_totals, sep = "=", collapse = ", "),
    "; observed ",
    paste(names(observed_totals), observed_totals, sep = "=", collapse = ", ")
  )
}

write.csv(events, output_file, row.names = FALSE, na = "")
write.csv(audit, audit_file, row.names = FALSE, na = "")

message(
  "Excluded ", sum(excluded), " duplicated juvenile rows from adult classes. ",
  "Analysis-ready animal totals match all three count histories."
)
