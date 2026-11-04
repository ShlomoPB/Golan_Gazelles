#!/usr/bin/env Rscript

# Non-destructive publication audit. This script reads existing inputs and
# outputs, but does not refit models, move files, or modify source datasets.

options(stringsAsFactors = FALSE)

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

read_count_matrix <- function(path) {
  x <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  y <- as.matrix(x[paste0("X", 1:7)])
  storage.mode(y) <- "numeric"
  rownames(y) <- x$camera_name
  y
}

assign_active_round <- function(event_camera, event_date, active_obs) {
  idx <- which(
    active_obs$camera_id == event_camera &
      event_date >= active_obs$start_date &
      event_date <= active_obs$end_date
  )
  if (!length(idx)) return(NA_integer_)
  idx <- idx[which.max(active_obs$start_date[idx])]
  active_obs$round[idx]
}

class_label <- c(juvenile = "Fawn", female = "Adult female", male = "Adult male")

site_file <- "SiteCovs_060324.csv"
obs_file <- "ObsCovs_050524.csv"
event_file <- "Gazella gazella_independent_events_30min_analysis_ready.csv"
count_files <- c(
  Fawn = "gazelle_juvenile_abundance_60x7.csv",
  `Adult female` = "gazelle_female_abundance_60x7.csv",
  `Adult male` = "gazelle_male_abundance_60x7.csv"
)
full_export_file <- file.path("..", "full-export_2026-04-15_2006.csv")
curated_export_file <- file.path("..", "full-export_2026-04-15_1717.csv")

required <- c(site_file, obs_file, event_file, count_files)
missing_required <- required[!file.exists(required)]
if (length(missing_required)) {
  stop("Missing required audit inputs: ", paste(missing_required, collapse = ", "))
}

site <- read.csv(site_file, na.strings = c("NA", ""), check.names = FALSE)
obs <- read.csv(obs_file, na.strings = c("NA", ""), check.names = FALSE)
events <- read.csv(event_file, na.strings = c("NA", ""), check.names = FALSE)
counts <- lapply(count_files, read_count_matrix)

obs_audit <- data.frame(
  camera_id = as.integer(obs[["Camera Number"]]),
  round = as.integer(obs[["Round Number"]]),
  start_date = as.Date(obs[["Start date"]], format = "%d/%m/%Y"),
  end_date = as.Date(obs[["End Date"]], format = "%d/%m/%Y"),
  days_operational = safe_num(obs[["Days.Operational"]])
)
active_obs <- obs_audit[
  is.finite(obs_audit$days_operational) &
    obs_audit$days_operational > 0 &
    !is.na(obs_audit$start_date) &
    !is.na(obs_audit$end_date),
]

event_timestamp <- as.POSIXct(
  sub("Z$", "", events$event_time),
  format = "%Y-%m-%dT%H:%M:%S",
  tz = "Asia/Jerusalem"
)
event_date <- as.Date(format(event_timestamp, tz = "Asia/Jerusalem"))

solar_elevation <- function(ts, latitude, longitude) {
  doy <- as.integer(format(ts, "%j", tz = "UTC"))
  hour <- as.numeric(format(ts, "%H", tz = "UTC")) +
    as.numeric(format(ts, "%M", tz = "UTC")) / 60 +
    as.numeric(format(ts, "%S", tz = "UTC")) / 3600
  gamma <- 2 * pi / 365 * (doy - 1 + (hour - 12) / 24)
  eq_time <- 229.18 * (
    0.000075 +
      0.001868 * cos(gamma) -
      0.032077 * sin(gamma) -
      0.014615 * cos(2 * gamma) -
      0.040849 * sin(2 * gamma)
  )
  decl <- 0.006918 -
    0.399912 * cos(gamma) +
    0.070257 * sin(gamma) -
    0.006758 * cos(2 * gamma) +
    0.000907 * sin(2 * gamma) -
    0.002697 * cos(3 * gamma) +
    0.00148 * sin(3 * gamma)
  true_solar_minutes <- (hour * 60 + eq_time + 4 * longitude) %% 1440
  hour_angle <- true_solar_minutes / 4 - 180
  lat_rad <- latitude * pi / 180
  ha_rad <- hour_angle * pi / 180
  cos_zenith <- sin(lat_rad) * sin(decl) +
    cos(lat_rad) * cos(decl) * cos(ha_rad)
  cos_zenith <- pmin(1, pmax(-1, cos_zenith))
  90 - acos(cos_zenith) * 180 / pi
}

event_round <- mapply(
  assign_active_round,
  event_camera = as.integer(events$camera_name),
  event_date = event_date,
  MoreArgs = list(active_obs = active_obs)
)
events$active_round <- as.integer(event_round)
events$valid_active_round <- !is.na(events$active_round)
events$class <- unname(class_label[as.character(events$gazelle_class)])
events$event_count <- safe_num(events$event_count)
count_summary <- do.call(
  rbind,
  lapply(names(counts), function(class_i) {
    y <- counts[[class_i]]
    camera_totals <- rowSums(y, na.rm = TRUE)
    data.frame(
      dataset = basename(count_files[[class_i]]),
      analysis_population = "Final 60 x 7 count history",
      class = class_i,
      unit = "Summed classified animals in retained camera-round cells",
      count = sum(y, na.rm = TRUE),
      event_rows = NA_integer_,
      cameras = sum(camera_totals > 0),
      nonmissing_camera_rounds = sum(!is.na(y)),
      missing_camera_rounds = sum(is.na(y)),
      filter = "NA camera-rounds excluded; zeros retained",
      stringsAsFactors = FALSE
    )
  })
)

event_summary <- do.call(
  rbind,
  lapply(names(class_label), function(code) {
    idx <- events$gazelle_class == code
    data.frame(
      dataset = basename(event_file),
      analysis_population = "All prepared 30-min independent events",
      class = unname(class_label[code]),
      unit = "Independent event rows and summed classified animals",
      count = sum(events$event_count[idx], na.rm = TRUE),
      event_rows = sum(idx),
      cameras = length(unique(events$camera_name[idx])),
      nonmissing_camera_rounds = NA_integer_,
      missing_camera_rounds = NA_integer_,
      filter = "No active-camera-round filter",
      stringsAsFactors = FALSE
    )
  })
)

valid_event_summary <- do.call(
  rbind,
  lapply(names(class_label), function(code) {
    idx <- events$gazelle_class == code & events$valid_active_round
    data.frame(
      dataset = basename(event_file),
      analysis_population = "Diel: events within documented active camera-rounds",
      class = unname(class_label[code]),
      unit = "Independent event rows and summed classified animals",
      count = sum(events$event_count[idx], na.rm = TRUE),
      event_rows = sum(idx),
      cameras = length(unique(events$camera_name[idx])),
      nonmissing_camera_rounds = NA_integer_,
      missing_camera_rounds = NA_integer_,
      filter = "Local event date within positive-effort camera-round; boundary assigned to later round",
      stringsAsFactors = FALSE
    )
  })
)

fawn_y <- counts[["Fawn"]]
female_y <- counts[["Adult female"]]
informative <- is.finite(fawn_y) & is.finite(female_y) & (fawn_y + female_y > 0)
informative_cameras <- rowSums(informative) > 0

beta_summary <- data.frame(
  dataset = paste(basename(count_files[c("Fawn", "Adult female")]), collapse = " + "),
  analysis_population = "Primary beta-binomial informative camera-rounds",
  class = c("Fawn", "Adult female"),
  unit = "Summed classified animals in informative camera-rounds",
  count = c(sum(fawn_y[informative]), sum(female_y[informative])),
  event_rows = NA_integer_,
  cameras = sum(informative_cameras),
  nonmissing_camera_rounds = sum(informative),
  missing_camera_rounds = NA_integer_,
  filter = "Both class counts finite, combined count > 0, complete model covariates",
  stringsAsFactors = FALSE
)

count_provenance <- rbind(event_summary, valid_event_summary, count_summary, beta_summary)
write_csv(count_provenance, "publication_count_provenance.csv")

descriptive_ratio <- data.frame(
  analysis_population = c(
    "Final 60 x 7 count histories",
    "All prepared 30-min independent events",
    "Diel events within date-resolved active camera-rounds"
  ),
  fawn_count = c(
    sum(fawn_y, na.rm = TRUE),
    sum(events$event_count[events$gazelle_class == "juvenile"], na.rm = TRUE),
    sum(events$event_count[
      events$gazelle_class == "juvenile" & events$valid_active_round
    ], na.rm = TRUE)
  ),
  adult_female_count = c(
    sum(female_y, na.rm = TRUE),
    sum(events$event_count[events$gazelle_class == "female"], na.rm = TRUE),
    sum(events$event_count[
      events$gazelle_class == "female" & events$valid_active_round
    ], na.rm = TRUE)
  ),
  status = c(
    "AUTHORITATIVE descriptive composition index and beta-binomial counts",
    "Source-event sensitivity; includes events outside documented active rounds",
    "Diel-population sensitivity; excludes events not mapped to date-resolved active rounds"
  ),
  stringsAsFactors = FALSE
)
descriptive_ratio$fawn_per_adult_female <-
  descriptive_ratio$fawn_count / descriptive_ratio$adult_female_count
descriptive_ratio$fawns_per_100_adult_females <-
  100 * descriptive_ratio$fawn_per_adult_female
write_csv(
  descriptive_ratio,
  "publication_descriptive_fawn_female_ratio.csv"
)

excluded_events <- events[!events$valid_active_round, c(
  "camera_name", "gazelle_class", "class", "event_time", "event_count"
)]
write_csv(excluded_events, "publication_events_excluded_from_active_rounds.csv")

reconciliation_rows <- list()
reconciliation_i <- 1L
for (class_i in names(counts)) {
  code_i <- names(class_label)[class_label == class_i]
  y <- counts[[class_i]]
  for (camera_i in rownames(y)) {
    for (round_i in seq_len(ncol(y))) {
      event_idx <-
        events$gazelle_class == code_i &
        events$valid_active_round &
        as.integer(events$camera_name) == as.integer(camera_i) &
        events$active_round == round_i
      event_sum <- sum(events$event_count[event_idx], na.rm = TRUE)
      event_rows_i <- sum(event_idx)
      count_value <- y[camera_i, round_i]
      reconciliation_rows[[reconciliation_i]] <- data.frame(
        class = class_i,
        camera_id = as.integer(camera_i),
        round = round_i,
        count_history_value = count_value,
        mapped_event_rows = event_rows_i,
        mapped_event_count_sum = event_sum,
        difference_event_minus_history = if (
          is.finite(count_value)
        ) event_sum - count_value else NA_real_,
        agreement = if (
          is.na(count_value)
        ) "Count history is NA" else if (event_sum == count_value) {
          "Exact"
        } else {
          "Different"
        },
        stringsAsFactors = FALSE
      )
      reconciliation_i <- reconciliation_i + 1L
    }
  }
}
count_reconciliation <- do.call(rbind, reconciliation_rows)
write_csv(count_reconciliation, "publication_count_history_event_reconciliation.csv")

data_audit <- data.frame(
  metric = c(
    "site_covariate_rows",
    "unique_site_ids",
    "survey_rounds",
    "expected_camera_rounds",
    "observation_covariate_rows",
    "duplicate_camera_round_rows",
    "positive_effort_camera_rounds",
    "missing_effort_camera_rounds",
    "date_resolved_positive_effort_camera_rounds",
    "positive_effort_rounds_missing_dates",
    "total_positive_effort_camera_days",
    "date_resolved_positive_effort_camera_days",
    "survey_start",
    "survey_end",
    "survey_midpoint",
    "prepared_independent_event_rows",
    "active_round_independent_event_rows",
    "events_excluded_outside_active_rounds",
    "beta_binomial_informative_camera_rounds",
    "beta_binomial_contributing_cameras",
    "beta_binomial_contributing_polygons"
  ),
  value = c(
    nrow(site),
    length(unique(site$Camera.ID)),
    length(unique(obs_audit$round)),
    length(unique(obs_audit$camera_id)) * length(unique(obs_audit$round)),
    nrow(obs_audit),
    sum(duplicated(obs_audit[c("camera_id", "round")])),
    sum(is.finite(obs_audit$days_operational) & obs_audit$days_operational > 0),
    sum(!is.finite(obs_audit$days_operational) | obs_audit$days_operational <= 0),
    nrow(active_obs),
    sum(
      is.finite(obs_audit$days_operational) &
        obs_audit$days_operational > 0 &
        (is.na(obs_audit$start_date) | is.na(obs_audit$end_date))
    ),
    sum(obs_audit$days_operational[is.finite(obs_audit$days_operational) & obs_audit$days_operational > 0]),
    sum(active_obs$days_operational),
    as.character(min(active_obs$start_date)),
    as.character(max(active_obs$end_date)),
    as.character(min(active_obs$start_date) + floor(as.numeric(max(active_obs$end_date) - min(active_obs$start_date)) / 2)),
    nrow(events),
    sum(events$valid_active_round),
    sum(!events$valid_active_round),
    sum(informative),
    sum(informative_cameras),
    length(unique(site$Polygon.Name[match(as.integer(rownames(fawn_y))[informative_cameras], site$Camera.ID)]))
  ),
  stringsAsFactors = FALSE
)
write_csv(data_audit, "publication_data_audit_summary.csv")

# Raw classification inventory. These are raw sighting rows, not independent
# events or unique animals. This table documents how legacy labels entered the
# project and prevents them from being silently pooled with focal classes.
classification_audit <- data.frame()
if (file.exists(full_export_file)) {
  collapse_30min_stats <- function(timestamp, group_key, quantity) {
    pieces <- split(seq_along(timestamp), group_key)
    stats <- vapply(pieces, function(piece_idx) {
      ord <- order(timestamp[piece_idx])
      t <- timestamp[piece_idx][ord]
      q <- quantity[piece_idx][ord]
      keep <- !is.na(t)
      t <- t[keep]
      q <- q[keep]
      if (!length(t)) return(c(events = 0, animals = 0))
      event_id <- cumsum(c(TRUE, as.numeric(diff(t), units = "mins") > 30))
      event_max <- tapply(q, event_id, max, na.rm = TRUE)
      event_max[!is.finite(event_max)] <- 0
      c(events = length(unique(event_id)), animals = sum(event_max))
    }, numeric(2))
    rowSums(stats)
  }

  if (requireNamespace("data.table", quietly = TRUE)) {
    raw <- data.table::fread(
      full_export_file,
      select = c(
        "Camera Name", "Date/Time", "Sighting Quantity", "Species",
        "Trap Station Session ID", "Sex", "Life stage"
      ),
      showProgress = FALSE
    )
  } else {
    raw <- read.csv(full_export_file, check.names = FALSE)
    raw <- raw[c(
      "Camera Name", "Date/Time", "Sighting Quantity", "Species",
      "Trap Station Session ID", "Sex", "Life stage"
    )]
  }
  raw <- raw[raw$Species == "Gazella gazella",]
  raw$Sex[is.na(raw$Sex) | raw$Sex == ""] <- "Missing"
  raw[["Life stage"]][is.na(raw[["Life stage"]]) | raw[["Life stage"]] == ""] <- "Missing"
  key <- paste(raw$Sex, raw[["Life stage"]], sep = " | ")
  split_rows <- split(seq_len(nrow(raw)), key)
  classification_audit <- do.call(
    rbind,
    lapply(split_rows, function(idx) {
      group_key <- paste(
        raw[["Camera Name"]][idx], raw[["Trap Station Session ID"]][idx],
        sep = "_"
      )
      timestamp <- as.POSIXct(
        raw[["Date/Time"]][idx],
        format = "%d/%m/%Y %H:%M",
        tz = "Asia/Jerusalem"
      )
      collapsed <- collapse_30min_stats(
        timestamp,
        group_key,
        safe_num(raw[["Sighting Quantity"]][idx])
      )
      data.frame(
        sex = raw$Sex[idx[1]],
        life_stage = raw[["Life stage"]][idx[1]],
        raw_sighting_rows = length(idx),
        summed_sighting_quantity = sum(safe_num(raw[["Sighting Quantity"]][idx]), na.rm = TRUE),
        collapsed_30min_events = unname(collapsed["events"]),
        collapsed_30min_summed_event_max = unname(collapsed["animals"]),
        cameras = length(unique(raw[["Camera Name"]][idx])),
        current_focal_class = if (
          raw$Sex[idx[1]] == "Female" && raw[["Life stage"]][idx[1]] == "Adult"
        ) "Potential adult-female source" else if (
          raw$Sex[idx[1]] == "Male" && raw[["Life stage"]][idx[1]] == "Adult"
        ) "Potential adult-male source" else if (
          raw[["Life stage"]][idx[1]] == "Juvenile"
        ) "Potential juvenile source; exact curation must be verified" else if (
          raw$Sex[idx[1]] == "Not Adult Male" || raw[["Life stage"]][idx[1]] == "Sub-adult"
        ) "Excluded from focal adult-male class unless explicitly documented" else {
          "Unidentifiable/nonfocal"
        },
        stringsAsFactors = FALSE
      )
    })
  )
  classification_audit <- classification_audit[
    order(-classification_audit$raw_sighting_rows),
  ]
}
write_csv(classification_audit, "publication_raw_gazelle_classification_audit.csv")

subadult_provenance <- data.frame()
export_files <- c(full_export_file, curated_export_file)
export_files <- export_files[file.exists(export_files)]
if (length(export_files) && exists("collapse_30min_stats")) {
  subadult_rows <- list()
  subadult_i <- 1L
  for (export_i in export_files) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      export_data <- data.table::fread(
        export_i,
        select = c(
          "Camera Name", "Date/Time", "Sighting Quantity", "Species",
          "Trap Station Session ID", "Sex", "Life stage"
        ),
        showProgress = FALSE
      )
    } else {
      export_data <- read.csv(export_i, check.names = FALSE)
    }
    export_data <- export_data[export_data$Species == "Gazella gazella",]
    definitions <- list(
      `Sex recorded as Not Adult Male` = export_data$Sex == "Not Adult Male",
      `Life stage recorded as Sub-adult` = export_data[["Life stage"]] == "Sub-adult",
      `Either subadult-related label` =
        export_data$Sex == "Not Adult Male" |
        export_data[["Life stage"]] == "Sub-adult"
    )
    for (definition_i in names(definitions)) {
      idx <- which(definitions[[definition_i]] %in% TRUE)
      collapsed <- collapse_30min_stats(
        timestamp_i <- as.POSIXct(
          export_data[["Date/Time"]][idx],
          format = "%d/%m/%Y %H:%M",
          tz = "Asia/Jerusalem"
        ),
        paste(
          export_data[["Camera Name"]][idx],
          export_data[["Trap Station Session ID"]][idx],
          sep = "_"
        ),
        safe_num(export_data[["Sighting Quantity"]][idx])
      )
      active_i <- !is.na(mapply(
        assign_active_round,
        event_camera = as.integer(export_data[["Camera Name"]][idx]),
        event_date = as.Date(format(timestamp_i, tz = "Asia/Jerusalem")),
        MoreArgs = list(active_obs = active_obs)
      ))
      collapsed_active <- collapse_30min_stats(
        timestamp_i[active_i],
        paste(
          export_data[["Camera Name"]][idx][active_i],
          export_data[["Trap Station Session ID"]][idx][active_i],
          sep = "_"
        ),
        safe_num(export_data[["Sighting Quantity"]][idx][active_i])
      )
      subadult_rows[[subadult_i]] <- data.frame(
        source_dataset = basename(export_i),
        definition = definition_i,
        raw_sighting_rows = length(idx),
        raw_sighting_quantity = sum(safe_num(export_data[["Sighting Quantity"]][idx]), na.rm = TRUE),
        collapsed_30min_events = unname(collapsed["events"]),
        collapsed_30min_summed_event_max = unname(collapsed["animals"]),
        active_round_collapsed_events = unname(collapsed_active["events"]),
        active_round_summed_event_max = unname(collapsed_active["animals"]),
        treatment_in_current_focal_analyses = "Excluded",
        stringsAsFactors = FALSE
      )
      subadult_i <- subadult_i + 1L
    }
  }
  subadult_provenance <- do.call(rbind, subadult_rows)
}
write_csv(subadult_provenance, "publication_subadult_male_provenance.csv")

variable_crosswalk <- data.frame(
  project_variable = c(
    "Null.Wolf.Abundance", "null.jackal.abundance", "human.abundance",
    "Open.Spaces", "Built.and.Roads", "Nature.Reserves",
    "Cattle.Sightings", "Days.Operational", "gazelle_class",
    "juvenile", "female", "male", "Early", "Late"
  ),
  standardized_variable = c(
    "Null.Wolf.Abundance.z / wolf_z", "null.jackal.abundance.z / jackal_z", "human_z",
    "Open.Spaces.z / open_z", "Built.and.Roads.z / built_z", "Nature.Reserves.z / reserve_z",
    "cattle_sightings.z", "log_days_operational.z / log_days", NA,
    NA, NA, NA, NA, NA
  ),
  recommended_manuscript_term = c(
    "site-level wolf activity/relative-abundance index",
    "site-level jackal activity/relative-abundance index",
    "site-level human activity index",
    "open-space cover", "built-and-road cover", "nature-reserve cover",
    "cattle sightings", "operational camera days", "gazelle demographic class",
    "fawn", "adult female", "adult male", "Early survey period", "Late survey period"
  ),
  caution = c(
    "Generated covariate; derivation and uncertainty propagation require documentation",
    "Generated covariate; derivation and uncertainty propagation require documentation",
    "Do not call abundance until derivation is verified",
    "Continuous site covariate", "Continuous site covariate", "Continuous site covariate",
    "Detection covariate", "Effort/detection covariate; log transformed in NMix",
    "Do not pool classes", "Legacy files may call this juvenile", "Classified adults only",
    "Not Adult Male and sub-adult records are not part of this class",
    "Survey period, not measured temperature", "Survey period, not measured temperature"
  ),
  stringsAsFactors = FALSE
)
write_csv(variable_crosswalk, "publication_variable_crosswalk.csv")

inventory <- file.info(list.files(".", full.names = FALSE, all.files = FALSE))
inventory$file <- rownames(inventory)
inventory$extension <- tools::file_ext(inventory$file)
inventory$category <- ifelse(
  grepl("gazelle_predator.*overlap|predator_overlap", inventory$file, ignore.case = TRUE),
  "LEGACY/EXCLUDED: gazelle-predator temporal overlap",
  ifelse(
    grepl("smoke", inventory$file, ignore.case = TRUE),
    "EXPLORATORY/SMOKE",
    ifelse(
      grepl("spPGOcc|spoccupancy|occupancy", inventory$file, ignore.case = TRUE),
      "EXPLORATORY/SUPERSEDED OCCUPANCY",
      ifelse(
        grepl("PAMM|recurrent|pressure", inventory$file, ignore.case = TRUE),
        "EXPLORATORY/NOT IN CURRENT MANUSCRIPT",
        "UNCLASSIFIED: review required"
      )
    )
  )
)
file_inventory <- inventory[c("file", "extension", "size", "mtime", "category")]
rownames(file_inventory) <- NULL
write_csv(file_inventory, "publication_file_inventory.csv")

input_files <- c(
  site_file, obs_file, event_file, unname(count_files),
  "Camera Coordinates.xlsx", full_export_file, curated_export_file,
  "publication_wolf_jackal_interaction_nmix_spabundance_fits.rds",
  "gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds",
  "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
)
input_files <- input_files[file.exists(input_files)]
input_info <- file.info(input_files)
checksums <- data.frame(
  file = input_files,
  bytes = input_info$size,
  modified = format(input_info$mtime, "%Y-%m-%d %H:%M:%S %z"),
  md5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)
write_csv(checksums, "publication_input_checksums.csv")

package_names <- c(
  "spAbundance", "unmarked", "glmmTMB", "TMB", "mgcv", "activity",
  "overlap", "data.table", "ggplot2"
)
package_versions <- data.frame(
  package = package_names,
  installed_version = vapply(package_names, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(utils::packageVersion(pkg))
    } else {
      NA_character_
    }
  }, character(1)),
  built = vapply(package_names, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      built_i <- utils::packageDescription(pkg)$Built
      if (is.null(built_i)) NA_character_ else as.character(built_i)
    } else {
      NA_character_
    }
  }, character(1)),
  library_path = vapply(package_names, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      dirname(find.package(pkg))
    } else {
      NA_character_
    }
  }, character(1)),
  audit_date = as.character(Sys.Date()),
  note = "",
  stringsAsFactors = FALSE
)
tmb_version <- package_versions$installed_version[package_versions$package == "TMB"]
glmm_tmb_built_tmb <- "1.9.18"
package_versions$note[package_versions$package == "glmmTMB"] <- if (
  length(tmb_version) == 1L && !is.na(tmb_version) &&
    identical(tmb_version, glmm_tmb_built_tmb)
) {
  paste0("Compatible with installed TMB ", tmb_version)
} else {
  paste0(
    "Version mismatch: glmmTMB was built with TMB ", glmm_tmb_built_tmb,
    "; installed TMB is ", tmb_version
  )
}
write_csv(package_versions, "publication_package_versions.csv")

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, "publication_audit_session_info.txt")

cat("Publication phase-1 audit complete.\n")
cat("  active camera-rounds:", nrow(active_obs), "\n")
cat("  camera-days:", sum(active_obs$days_operational), "\n")
cat("  valid diel events:", sum(events$valid_active_round), "\n")
cat("  beta-binomial informative rounds:", sum(informative), "\n")
