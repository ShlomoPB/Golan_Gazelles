#!/usr/bin/env Rscript

# Prepare an internal, student-ready copy of the all-species Camelot export.
# Every source row is retained. Documented camera-clock corrections are added
# transparently, and no source file is modified.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source_file <- "../full-export_2026-04-15_2006.csv"
gazelle_flash_file <- "../full-export_2026-04-15_1717.csv"
gazelle_correction_file <- "camera_timestamp_corrections.csv"
deployment_file <- "ObsCovs_050524.csv"
camera_site_file <- "publication_camera_array_map_sites.csv"
output_dir <- "student_data_all_species"

required_files <- c(
  source_file,
  gazelle_flash_file,
  gazelle_correction_file,
  deployment_file,
  camera_site_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input(s): ", paste(missing_files, collapse = ", "))
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

parse_dmy_date <- function(x) {
  as.Date(substr(x, 1, 10), format = "%d/%m/%Y")
}

to_iso_local <- function(x) {
  parsed_date <- parse_dmy_date(x)
  time_text <- substr(x, 12, nchar(x))
  ifelse(
    is.na(parsed_date) | !nzchar(time_text),
    NA_character_,
    paste(format(parsed_date, "%Y-%m-%d"), time_text)
  )
}

replace_date_keep_clock <- function(x, offset_days) {
  parsed_date <- parse_dmy_date(x)
  time_text <- substr(x, 12, nchar(x))
  corrected_date <- parsed_date + offset_days
  ifelse(
    is.na(corrected_date) | !nzchar(time_text),
    NA_character_,
    paste(format(corrected_date, "%Y-%m-%d"), time_text)
  )
}

raw <- read_csv(
  source_file,
  col_types = cols(.default = col_character()),
  na = c("", "NA"),
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "Camera ID", "Camera Name", "Attention Needed Flag", "Camera Check Flag",
  "Date/Time", "Sighting Created", "Sighting Quantity", "Species",
  "Trap Station ID", "Camelot GPS Latitude", "Camelot GPS Longitude",
  "Trap Station Name", "Trap Station Session Camera ID", "Session End Date",
  "Trap Station Session ID", "Session Start Date", "Sex", "Life stage"
)
if (!setequal(names(raw), required_columns)) {
  stop("The all-species export schema differs from the expected 18 columns.")
}
if (anyNA(parse_dmy_date(raw[["Date/Time"]]))) {
  stop("One or more source sighting dates could not be parsed.")
}

# Preserve the source label and add transparent, student-facing name fields.
# Common names harmonize duplicate scientific/common labels but do not alter
# the recorded classification or imply a new taxonomic identification.
species_name_map <- tribble(
  ~species_recorded_label, ~species_label_corrected, ~species_common_name,
  "Alectoris chukar",      "Alectoris chukar",       "Chukar",
  "Bicycle Rider",        "Bicycle rider",         "Bicycle rider",
  "Bird unknown",         "Unidentified bird",     "Unidentified bird",
  "Bos taurus",           "Bos taurus",            "Domestic cattle",
  "Bubulcus ibis",        "Bubulcus ibis",         "Cattle egret",
  "Canid unknown",        "Unidentified canid",    "Unidentified canid",
  "Canis aureus",         "Canis aureus",          "Golden jackal",
  "Canis familiaris",     "Canis familiaris",      "Domestic dog",
  "Canis lupus",          "Canis lupus",           "Gray wolf",
  "Capra aegagrus_hircus", "Capra aegagrus hircus", "Domestic goat",
  "Chukar partridge",     "Chukar partridge",      "Chukar",
  "Corvus cornix",        "Corvus cornix",         "Hooded crow",
  "Equus africanus",      "Equus africanus",       "Domestic donkey",
  "Equus ferus caballus", "Equus ferus caballus", "Domestic horse",
  "Erinaceus concolor",   "Erinaceus concolor",    "Southern white-breasted hedgehog",
  "Felis silvesteris",    "Felis silvestris",      "Wildcat",
  "Gazella gazella",      "Gazella gazella",       "Mountain gazelle",
  "Herpestes ichneumon",  "Herpestes ichneumon",  "Egyptian mongoose",
  "Homo sapiens",         "Homo sapiens",          "Human",
  "Homo sapiens-solider", "Homo sapiens-soldier",  "Soldier",
  "Hyaena hyaena",        "Hyaena hyaena",         "Striped hyena",
  "Hystrix indica",       "Hystrix indica",        "Indian crested porcupine",
  "Lanius collurio",      "Lanius collurio",       "Red-backed shrike",
  "Lepus capensis",       "Lepus capensis",        "Cape hare",
  "Martes Foina",         "Martes foina",          "Beech marten",
  "Meles meles",          "Meles meles",           "European badger",
  "Ovis aries",           "Ovis aries",            "Domestic sheep",
  "Procavia capensis",    "Procavia capensis",     "Rock hyrax",
  "Red-backed shrike",    "Red-backed shrike",     "Red-backed shrike",
  "Rodent Small",         "Small rodent",          "Small rodent",
  "Stellagama stellio",   "Stellagama stellio",    "Starred agama",
  "Sus scrofa",           "Sus scrofa",            "Wild boar",
  "Testudo graeca",       "Testudo graeca",        "Spur-thighed tortoise",
  "Unidentified Animal",  "Unidentified animal",   "Unidentified animal",
  "Vehicle civilian",     "Civilian vehicle",      "Civilian vehicle",
  "Vehicle construction", "Construction vehicle",  "Construction vehicle",
  "Vehicle military",     "Military vehicle",      "Military vehicle",
  "Vulpes vulpes",        "Vulpes vulpes",         "Red fox"
)

recorded_species <- sort(unique(na.omit(raw[["Species"]])))
unmapped_species <- setdiff(recorded_species, species_name_map$species_recorded_label)
unused_mappings <- setdiff(species_name_map$species_recorded_label, recorded_species)
if (length(unmapped_species) || length(unused_mappings)) {
  stop(
    "Species-name mapping does not exactly match the source labels. Unmapped: ",
    paste(unmapped_species, collapse = ", "), "; unused: ",
    paste(unused_mappings, collapse = ", ")
  )
}

base_corrections <- read_csv(
  gazelle_correction_file,
  col_types = cols(
    camera_name = col_integer(),
    trap_station_session_id = col_integer(),
    offset_days = col_integer(),
    expected_round = col_integer(),
    reason = col_character()
  ),
  show_col_types = FALSE
)

# The prior table was scoped to sessions containing gazelle events. Camera 24
# session 25 contains no gazelle event but is the first deployment from the
# same clock that was set to 1 January. The authoritative deployment began on
# 26 August, exactly 238 days later, so the same documented correction applies.
all_species_extension <- tibble(
  camera_name = 24L,
  trap_station_session_id = 25L,
  offset_days = 238L,
  expected_round = 1L,
  reason = paste(
    "All-species extension: first camera-24 session from the clock set to",
    "01 January; authoritative deployment began 26 August"
  )
)

corrections <- bind_rows(base_corrections, all_species_extension) %>%
  arrange(camera_name, expected_round, trap_station_session_id)

correction_keys <- paste(
  corrections$camera_name,
  corrections$trap_station_session_id,
  sep = ":"
)
if (anyDuplicated(correction_keys)) stop("Duplicate correction keys detected.")
if (!all(corrections$camera_name %in% c(22L, 24L))) {
  stop("Unexpected camera in the timestamp-correction table.")
}

station_camera_id <- as.integer(raw[["Camera Name"]])
session_id <- as.integer(raw[["Trap Station Session ID"]])
source_keys <- paste(station_camera_id, session_id, sep = ":")
correction_index <- match(source_keys, correction_keys)
correction_applied <- !is.na(correction_index)

missing_correction_sessions <- setdiff(correction_keys, unique(source_keys))
if (length(missing_correction_sessions)) {
  stop(
    "Correction session(s) missing from all-species export: ",
    paste(missing_correction_sessions, collapse = ", ")
  )
}

offset_days <- integer(nrow(raw))
expected_round <- rep(NA_integer_, nrow(raw))
correction_reason <- rep(NA_character_, nrow(raw))
offset_days[correction_applied] <-
  corrections$offset_days[correction_index[correction_applied]]
expected_round[correction_applied] <-
  corrections$expected_round[correction_index[correction_applied]]
correction_reason[correction_applied] <-
  corrections$reason[correction_index[correction_applied]]

original_date <- parse_dmy_date(raw[["Date/Time"]])
corrected_date <- original_date + offset_days
datetime_original_local <- to_iso_local(raw[["Date/Time"]])
datetime_corrected_local <- replace_date_keep_clock(
  raw[["Date/Time"]],
  offset_days
)

if (any(substr(datetime_original_local, 12, 16) !=
        substr(datetime_corrected_local, 12, 16))) {
  stop("A timestamp correction changed the recorded local clock time.")
}
if (any(corrected_date[!correction_applied] != original_date[!correction_applied])) {
  stop("An unflagged row changed date.")
}

deployment_raw <- read_csv(
  deployment_file,
  col_types = cols(.default = col_character()),
  na = c("", "NA"),
  show_col_types = FALSE
)
deployments <- deployment_raw %>%
  transmute(
    station_camera_id = as.integer(.data[["Camera Number"]]),
    survey_round = as.integer(.data[["Round Number"]]),
    deployment_start_date = as.Date(.data[["Start date"]], "%d/%m/%Y"),
    deployment_end_date = as.Date(.data[["End Date"]], "%d/%m/%Y"),
    days_operational = suppressWarnings(as.numeric(.data[["Days.Operational"]])),
    camera_model = .data[["Camera.Model"]],
    camera_sensitivity = .data[["Camera.Sensitivity"]],
    capture_delay = .data[["Capture.Delay"]],
    vegetation_structure = .data[["Vegetation.Structure"]],
    path_type = .data[["Path.Type"]],
    area_behind_camera = .data[["Area.Behind.Camera"]],
    bottleneck = .data[["Bottleneck."]],
    livestock_presence = .data[["Livestock.Present."]],
    agricultural_presence = .data[["Agricultural.Presence"]]
  ) %>%
  filter(!is.na(deployment_start_date), !is.na(deployment_end_date)) %>%
  arrange(station_camera_id, deployment_start_date, survey_round)

# Assign by corrected sighting date. Iterating deployments in chronological
# order makes a boundary date belong to the later deployment, matching the
# session-start convention. Documented correction sessions are then assigned
# their verified expected round explicitly.
survey_round <- rep(NA_integer_, nrow(raw))
assignment_method <- rep("outside documented deployment dates", nrow(raw))
for (camera in sort(unique(deployments$station_camera_id))) {
  camera_rows <- which(station_camera_id == camera)
  camera_deployments <- deployments %>%
    filter(station_camera_id == camera) %>%
    arrange(deployment_start_date, survey_round)
  for (j in seq_len(nrow(camera_deployments))) {
    in_deployment <-
      corrected_date[camera_rows] >= camera_deployments$deployment_start_date[j] &
      corrected_date[camera_rows] <= camera_deployments$deployment_end_date[j]
    matched_rows <- camera_rows[which(in_deployment)]
    survey_round[matched_rows] <- camera_deployments$survey_round[j]
    assignment_method[matched_rows] <- "corrected date within deployment"
  }
}
survey_round[correction_applied] <- expected_round[correction_applied]
assignment_method[correction_applied] <- "documented corrected session"

if (!all(survey_round[correction_applied] == expected_round[correction_applied])) {
  stop("A corrected session was not assigned to its documented survey round.")
}

deployment_details <- deployments %>%
  rename(expected_round = survey_round)

sites <- read_csv(
  camera_site_file,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  transmute(
    station_camera_id = as.integer(.data[["Camera.ID"]]),
    sampling_stratum = sampling_stratum,
    stratum_type = stratum_type,
    latitude_authoritative = as.numeric(latitude),
    longitude_authoritative = as.numeric(longitude)
  )
if (n_distinct(sites$station_camera_id) != nrow(sites)) {
  stop("Authoritative camera-site table has duplicate camera IDs.")
}

student_data <- tibble(
  source_row_id = seq_len(nrow(raw)),
  station_camera_id = station_camera_id,
  hardware_camera_id_exported = as.integer(raw[["Camera ID"]]),
  trap_station_id_exported = as.integer(raw[["Trap Station ID"]]),
  trap_station_name_exported = raw[["Trap Station Name"]],
  trap_station_session_id = session_id,
  trap_station_session_camera_id =
    as.integer(raw[["Trap Station Session Camera ID"]]),
  survey_round = survey_round,
  round_assignment_method = assignment_method,
  datetime_original_local = datetime_original_local,
  datetime_corrected_local = datetime_corrected_local,
  date_correction_applied = correction_applied,
  date_correction_days = offset_days,
  date_correction_reason = correction_reason,
  sighting_created_exported = to_iso_local(raw[["Sighting Created"]]),
  species = raw[["Species"]],
  species_recorded = !is.na(raw[["Species"]]),
  sighting_quantity = suppressWarnings(as.numeric(raw[["Sighting Quantity"]])),
  sex = raw[["Sex"]],
  life_stage = raw[["Life stage"]],
  attention_needed_flag = raw[["Attention Needed Flag"]],
  camera_check_flag = raw[["Camera Check Flag"]],
  latitude_exported = suppressWarnings(as.numeric(raw[["Camelot GPS Latitude"]])),
  longitude_exported = suppressWarnings(as.numeric(raw[["Camelot GPS Longitude"]])),
  session_start_exported_original = to_iso_local(raw[["Session Start Date"]]),
  session_end_exported_original = to_iso_local(raw[["Session End Date"]])
) %>%
  left_join(
    species_name_map,
    by = c("species" = "species_recorded_label")
  ) %>%
  relocate(species_label_corrected, species_common_name, .after = species) %>%
  left_join(
    deployment_details,
    by = c(
      "station_camera_id",
      "survey_round" = "expected_round"
    )
  ) %>%
  left_join(sites, by = "station_camera_id") %>%
  relocate(
    deployment_start_date,
    deployment_end_date,
    days_operational,
    .after = round_assignment_method
  ) %>%
  relocate(
    latitude_authoritative,
    longitude_authoritative,
    sampling_stratum,
    stratum_type,
    .after = longitude_exported
  )

if (nrow(student_data) != nrow(raw)) stop("Row preservation check failed.")
if (anyNA(student_data$station_camera_id)) stop("Missing station camera IDs.")
if (anyNA(student_data$datetime_corrected_local)) {
  stop("Missing corrected sighting timestamps.")
}

corrected_cameras <- sort(unique(student_data$station_camera_id[
  student_data$date_correction_applied
]))
if (!identical(corrected_cameras, c(22L, 24L))) {
  stop("Corrections were applied to unexpected cameras.")
}

# Flash was omitted from the all-species export but retained in a gazelle-only
# export with the same camera, session, and original sighting timestamps. Join
# a flash value only where that exact timestamp group has one unambiguous flash
# status. This recovers all but 30 gazelle rows and a small number of exact-time
# co-detections; all other species rows remain explicitly unavailable.
flash_raw <- read_csv(
  gazelle_flash_file,
  col_types = cols(.default = col_character()),
  na = c("", "NA"),
  show_col_types = FALSE,
  progress = FALSE
)
if (!all(c(
  "Camera Name", "Trap Station Session ID", "Date/Time", "Flash", "Species"
) %in% names(flash_raw))) {
  stop("Gazelle flash export lacks one or more required columns.")
}
if (n_distinct(flash_raw$Species, na.rm = TRUE) != 1L ||
    unique(na.omit(flash_raw$Species)) != "Gazella gazella") {
  stop("The flash source is not the expected gazelle-only export.")
}

flash_map <- flash_raw %>%
  transmute(
    station_camera_id = as.integer(.data[["Camera Name"]]),
    trap_station_session_id = as.integer(.data[["Trap Station Session ID"]]),
    datetime_original_local = to_iso_local(.data[["Date/Time"]]),
    flash_export_value = .data[["Flash"]]
  ) %>%
  group_by(
    station_camera_id,
    trap_station_session_id,
    datetime_original_local
  ) %>%
  summarise(
    flash_value_count = n_distinct(flash_export_value, na.rm = TRUE),
    flash_status = if (flash_value_count == 1L) {
      sort(unique(na.omit(flash_export_value)))[1]
    } else {
      NA_character_
    },
    .groups = "drop"
  )

student_data <- student_data %>%
  left_join(
    flash_map,
    by = c(
      "station_camera_id",
      "trap_station_session_id",
      "datetime_original_local"
    )
  ) %>%
  mutate(
    flash_used = case_when(
      grepl("^Flash fired", flash_status) ~ TRUE,
      grepl("^Flash did not fire", flash_status) ~ FALSE,
      TRUE ~ NA
    ),
    flash_data_source = case_when(
      flash_value_count == 1L ~ "gazelle_export_exact_camera_session_time",
      flash_value_count > 1L ~ "ambiguous_gazelle_export_time",
      TRUE ~ "unavailable_in_all_species_export"
    )
  ) %>%
  select(-flash_value_count)

gazelle_flash_check <- student_data %>%
  filter(species == "Gazella gazella") %>%
  count(flash_data_source)
if (sum(gazelle_flash_check$n) !=
    sum(student_data$species == "Gazella gazella", na.rm = TRUE)) {
  stop("Gazelle flash reconciliation changed the gazelle row total.")
}
if (sum(
  gazelle_flash_check$n[
    gazelle_flash_check$flash_data_source == "unavailable_in_all_species_export"
  ]
) != 0L) {
  stop("One or more gazelle rows failed to match the flash export.")
}

# Concise final file for the student. The authoritative coordinates are the
# analysis-ready locations; exported coordinates are retained to expose any
# source discrepancy rather than silently replacing it.
student_final <- student_data %>%
  transmute(
    source_row_id,
    camera_id = station_camera_id,
    camera_session_id = trap_station_session_id,
    camera_round = survey_round,
    latitude_authoritative,
    longitude_authoritative,
    latitude_exported,
    longitude_exported,
    deployment_start_date,
    deployment_end_date,
    days_operational,
    date_time = datetime_corrected_local,
    date_correction_applied,
    date_correction_days,
    quantity = sighting_quantity,
    species_recorded = species,
    species = species_label_corrected,
    species_common_name,
    sex,
    life_stage,
    flash_used,
    flash_status,
    flash_data_source,
    polygon_name = sampling_stratum,
    polygon_type = stratum_type,
    camera_model,
    camera_sensitivity,
    capture_delay,
    vegetation_structure,
    path_type,
    area_behind_camera,
    bottleneck,
    livestock_presence,
    agricultural_presence
  )

if (nrow(student_final) != nrow(raw)) {
  stop("Final student CSV failed the source-row preservation check.")
}
write_csv(
  student_final,
  file.path(output_dir, "golan_all_species_student_final.csv"),
  na = ""
)

full_output <- file.path(
  output_dir,
  "golan_all_species_raw_rows_date_corrected.csv"
)
write_csv(student_data, full_output, na = "")
compressed_output <- paste0(full_output, ".gz")
write_csv(student_data, compressed_output, na = "")

# Exact duplicate source rows are retained above. This compact copy records
# their multiplicity but must not be mistaken for a 30-minute independent-
# event dataset.
collapsed_data <- student_data %>%
  group_by(across(-source_row_id)) %>%
  summarise(
    first_source_row_id = min(source_row_id),
    raw_row_multiplicity = n(),
    .groups = "drop"
  ) %>%
  relocate(first_source_row_id, raw_row_multiplicity)

collapsed_output <- file.path(
  output_dir,
  "golan_all_species_exact_rows_collapsed_date_corrected.csv"
)
write_csv(collapsed_data, collapsed_output, na = "")

correction_audit <- student_data %>%
  filter(date_correction_applied) %>%
  mutate(
    original_date = as.Date(substr(datetime_original_local, 1, 10)),
    corrected_date = as.Date(substr(datetime_corrected_local, 1, 10)),
    within_authoritative_deployment =
      corrected_date >= deployment_start_date &
      corrected_date <= deployment_end_date
  ) %>%
  group_by(
    station_camera_id,
    trap_station_session_id,
    date_correction_days,
    survey_round,
    date_correction_reason,
    deployment_start_date,
    deployment_end_date
  ) %>%
  summarise(
    source_rows = n(),
    original_date_min = min(original_date),
    original_date_max = max(original_date),
    corrected_date_min = min(corrected_date),
    corrected_date_max = max(corrected_date),
    rows_within_authoritative_deployment =
      sum(within_authoritative_deployment, na.rm = TRUE),
    rows_outside_authoritative_deployment =
      sum(!within_authoritative_deployment, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(station_camera_id, survey_round)
write_csv(
  correction_audit,
  file.path(output_dir, "golan_all_species_date_correction_audit.csv"),
  na = ""
)
write_csv(
  corrections,
  file.path(output_dir, "golan_all_species_timestamp_corrections_applied.csv"),
  na = ""
)

unassigned_round_summary <- student_data %>%
  filter(is.na(survey_round)) %>%
  mutate(
    corrected_date = as.Date(substr(datetime_corrected_local, 1, 10)),
    species_for_summary = coalesce(species, "[missing species label]")
  ) %>%
  group_by(
    station_camera_id,
    trap_station_session_id,
    corrected_date,
    species_for_summary
  ) %>%
  summarise(
    source_rows = n(),
    first_source_row_id = min(source_row_id),
    .groups = "drop"
  ) %>%
  rename(species = species_for_summary) %>%
  arrange(station_camera_id, trap_station_session_id, corrected_date, species)
write_csv(
  unassigned_round_summary,
  file.path(output_dir, "golan_all_species_unassigned_round_rows_summary.csv"),
  na = ""
)

species_summary <- student_data %>%
  mutate(
    species_for_summary = coalesce(species, "[missing species label]"),
    species_label_corrected = coalesce(
      species_label_corrected,
      "[missing species label]"
    ),
    species_common_name = coalesce(
      species_common_name,
      "[missing species label]"
    )
  ) %>%
  group_by(
    species_for_summary,
    species_label_corrected,
    species_common_name
  ) %>%
  summarise(
    raw_rows = n(),
    cameras = n_distinct(station_camera_id),
    summed_sighting_quantity = sum(sighting_quantity, na.rm = TRUE),
    corrected_timestamp_rows = sum(date_correction_applied),
    .groups = "drop"
  ) %>%
  rename(species = species_for_summary) %>%
  arrange(desc(raw_rows), species)
write_csv(
  species_summary,
  file.path(output_dir, "golan_all_species_species_summary.csv"),
  na = ""
)
write_csv(
  species_name_map,
  file.path(output_dir, "golan_all_species_species_name_dictionary.csv"),
  na = ""
)

data_dictionary <- tribble(
  ~variable, ~description,
  "source_row_id", "Stable row number in the 2026-04-15 20:06 source export.",
  "station_camera_id", "Camera-trap station number from Camera Name; use this as the site ID.",
  "hardware_camera_id_exported", "Camelot Camera ID retained from the export; this is not the station ID.",
  "trap_station_id_exported", "Camelot trap-station database ID retained from the export.",
  "trap_station_name_exported", "Camelot trap-station name retained from the export.",
  "trap_station_session_id", "Camelot deployment/session identifier used for timestamp corrections.",
  "trap_station_session_camera_id", "Camelot session-camera database identifier.",
  "survey_round", "Survey round assigned from corrected date and the authoritative deployment table; documented correction sessions use their verified expected round.",
  "round_assignment_method", "Whether round came from the corrected date or a documented corrected session; rows outside deployment dates remain unassigned.",
  "deployment_start_date", "Authoritative deployment start date from ObsCovs_050524.csv.",
  "deployment_end_date", "Authoritative deployment end date from ObsCovs_050524.csv.",
  "days_operational", "Recorded operational days for the camera-round.",
  "datetime_original_local", "Original camera-clock sighting timestamp, standardized to ISO order without timezone conversion.",
  "datetime_corrected_local", "Date-corrected camera-clock timestamp. Clock time of day is preserved exactly.",
  "date_correction_applied", "TRUE only for documented camera 22 or 24 sessions requiring correction.",
  "date_correction_days", "Whole days added to the original camera date: 148 for affected camera-22 sessions and 238 for camera 24.",
  "date_correction_reason", "Documented reason for the date correction.",
  "sighting_created_exported", "Camelot database creation timestamp, standardized only; it is not the camera sighting time.",
  "species", "Species or non-wildlife category exactly as recorded in the source export.",
  "species_label_corrected", "Recorded label with evident spelling, capitalization, or word-order errors corrected; the source label remains in species.",
  "species_common_name", "Consistent English common name or non-wildlife category derived from the recorded label; duplicate scientific/common labels are harmonized here.",
  "species_recorded", "Whether the export contains a nonmissing species/category label.",
  "sighting_quantity", "Recorded number of individuals for the sighting row.",
  "sex", "Recorded sex category.",
  "life_stage", "Recorded life-stage category.",
  "attention_needed_flag", "Camelot attention-needed flag retained from the export.",
  "camera_check_flag", "Camelot camera-check flag retained from the export.",
  "latitude_exported", "Latitude in the Camelot export; retained for provenance.",
  "longitude_exported", "Longitude in the Camelot export; retained for provenance.",
  "latitude_authoritative", "Authoritative camera-map latitude; preferred for analysis.",
  "longitude_authoritative", "Authoritative camera-map longitude; preferred for analysis.",
  "sampling_stratum", "Named sampling stratum from the authoritative camera-site table.",
  "stratum_type", "Sampling-stratum category from the authoritative camera-site table.",
  "session_start_exported_original", "Original Camelot session start; known to inherit clock-date errors for affected cameras and retained only for provenance.",
  "session_end_exported_original", "Original Camelot session end; known to inherit clock-date errors for affected cameras and retained only for provenance.",
  "camera_model", "Camera model recorded for the assigned camera-round.",
  "camera_sensitivity", "Camera sensitivity setting recorded for the assigned camera-round.",
  "capture_delay", "Capture-delay setting recorded for the assigned camera-round.",
  "vegetation_structure", "Vegetation structure recorded for the assigned camera-round.",
  "path_type", "Path type recorded for the assigned camera-round.",
  "area_behind_camera", "Area behind camera recorded for the assigned camera-round.",
  "bottleneck", "Bottleneck indicator recorded for the assigned camera-round.",
  "livestock_presence", "Livestock-presence category recorded for the assigned camera-round.",
  "agricultural_presence", "Agricultural-presence category recorded for the assigned camera-round.",
  "first_source_row_id", "Collapsed-file only: first source row represented by the exact-record group.",
  "raw_row_multiplicity", "Collapsed-file only: number of exactly repeated raw rows represented. This is not an independent-event count."
)
write_csv(
  data_dictionary,
  file.path(output_dir, "golan_all_species_data_dictionary.csv"),
  na = ""
)

source_exact_duplicate_rows <- sum(duplicated(raw))
unassigned_rows <- sum(is.na(student_data$survey_round))
qc_summary <- tibble(
  metric = c(
    "Source file",
    "Source MD5",
    "Source rows",
    "Output rows",
    "Source columns",
    "Output columns",
    "Output CSV MD5",
    "Compressed output MD5",
    "Station cameras",
    "Distinct nonmissing species/category labels",
    "Rows with missing species/category label",
    "Exactly duplicated source rows beyond first copy",
    "Exact-record groups in collapsed file",
    "Corrected rows total",
    "Corrected camera 22 rows",
    "Corrected camera 24 rows",
    "Corrected rows outside nominal deployment dates",
    "Rows with unchanged local clock time",
    "Rows without assigned survey round",
    "Corrected timestamp minimum",
    "Corrected timestamp maximum",
    "Authoritative coordinate rows missing"
  ),
  value = c(
    normalizePath(source_file, winslash = "/"),
    unname(tools::md5sum(source_file)),
    nrow(raw),
    nrow(student_data),
    ncol(raw),
    ncol(student_data),
    unname(tools::md5sum(full_output)),
    unname(tools::md5sum(compressed_output)),
    n_distinct(student_data$station_camera_id),
    n_distinct(student_data$species, na.rm = TRUE),
    sum(!student_data$species_recorded),
    source_exact_duplicate_rows,
    nrow(collapsed_data),
    sum(student_data$date_correction_applied),
    sum(student_data$date_correction_applied & student_data$station_camera_id == 22L),
    sum(student_data$date_correction_applied & student_data$station_camera_id == 24L),
    sum(correction_audit$rows_outside_authoritative_deployment),
    sum(
      substr(student_data$datetime_original_local, 12, 16) ==
        substr(student_data$datetime_corrected_local, 12, 16)
    ),
    unassigned_rows,
    min(student_data$datetime_corrected_local),
    max(student_data$datetime_corrected_local),
    sum(
      is.na(student_data$latitude_authoritative) |
        is.na(student_data$longitude_authoritative)
    )
  )
)
write_csv(
  qc_summary,
  file.path(output_dir, "golan_all_species_qc_summary.csv"),
  na = ""
)

readme <- c(
  "# Golan all-species camera-trap export: student data package",
  "",
  paste0(
    "Source: `", basename(source_file), "` (", format(nrow(raw), big.mark = ","),
    " raw rows; source MD5 `", unname(tools::md5sum(source_file)), "`)."
  ),
  "The source export was not modified.",
  "",
  "## Primary file",
  "",
  "`golan_all_species_raw_rows_date_corrected.csv` preserves every source row.",
  "The identical `.csv.gz` copy is recommended for transfer because it is smaller.",
  "Use `datetime_corrected_local` for sighting time and `station_camera_id` for site.",
  "The original timestamp remains in `datetime_original_local`, and every date",
  "change is documented by the correction flag, day offset, reason, and audit table.",
  "",
  "Camera 22 affected sessions received +148 days. Camera 24 received +238 days",
  "for all six sessions, including first-session ID 25, which was absent from the",
  "gazelle-only correction table because that session had no gazelle event.",
  paste0(
    "In total, ", format(sum(student_data$date_correction_applied), big.mark = ","),
    " all-species rows were date-corrected. Recorded clock time of day was unchanged."
  ),
  paste0(
    format(sum(correction_audit$rows_outside_authoritative_deployment), big.mark = ","),
    " corrected camera-24 rows (", nrow(correction_audit %>%
      filter(rows_outside_authoritative_deployment > 0)),
    " sessions) fall on the evening immediately before the nominal deployment",
    " start date. They retain the established +238-day correction and are assigned",
    " by documented session; no additional one-day shift was invented."
  ),
  "",
  "## Important limitations",
  "",
  paste0(
    "The export contains ", format(source_exact_duplicate_rows, big.mark = ","),
    " exact repeated rows beyond the first copy. These may represent repeated image",
    " frames or export duplication; the export has no media identifier with which to",
    " distinguish them. They are retained in the primary file."
  ),
  "`golan_all_species_exact_rows_collapsed_date_corrected.csv` is a convenience copy",
  "with `raw_row_multiplicity`, but it is not an independent-detection dataset.",
  "A species-specific independence rule must be applied for event-rate analyses.",
  "",
  "Species labels, including human and vehicle categories and missing labels, are",
  "preserved exactly as exported in the full file's `species` column. In the concise",
  "student file, `species` is corrected, `species_recorded` preserves the source label,",
  "and `species_common_name` provides the harmonized common name. No biological",
  "filtering was performed. Exact",
  "coordinates are included for internal project use and should",
  "not be redistributed without considering sensitive-species and site confidentiality.",
  "",
  "The exported session start/end fields are retained only for provenance because",
  "they inherit camera-clock date errors. Prefer the authoritative deployment dates",
  "and `survey_round` columns.",
  paste0(
    format(unassigned_rows, big.mark = ","),
    " rows fall outside all documented deployment intervals and therefore retain",
    " a missing survey round; they are itemized in the unassigned-row summary."
  ),
  "",
  "## Supporting files",
  "",
  "- `golan_all_species_date_correction_audit.csv`: correction totals and date ranges.",
  "- `golan_all_species_timestamp_corrections_applied.csv`: correction rules used.",
  "- `golan_all_species_species_summary.csv`: row and quantity totals by recorded label.",
  "- `golan_all_species_species_name_dictionary.csv`: source, corrected, and common-name mapping.",
  "- `golan_all_species_unassigned_round_rows_summary.csv`: rows outside documented deployments.",
  "- `golan_all_species_data_dictionary.csv`: column definitions.",
  "- `golan_all_species_qc_summary.csv`: reproducibility and integrity checks.",
  "",
  "Generated by `wrangle_golan_all_species_student_data.R`."
)
writeLines(readme, file.path(output_dir, "README.md"), useBytes = TRUE)

message(
  "Wrote ", format(nrow(student_data), big.mark = ","),
  " row-preserving records and ", format(nrow(collapsed_data), big.mark = ","),
  " exact-record groups to ", output_dir, "."
)
