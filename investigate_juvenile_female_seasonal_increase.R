#!/usr/bin/env Rscript

# Exploratory audit of the Early/Late increase in the camera-based
# juvenile:female composition index. The published model object is read but
# not overwritten or refitted. All new outputs are written to a separate
# directory.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tidyr)
})

output_dir <- "juvenile_female_seasonal_audit"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

juvenile_file <- "gazelle_juvenile_abundance_60x7.csv"
female_file <- "gazelle_female_abundance_60x7.csv"
site_file <- "SiteCovs_060324.csv"
obs_file <- "ObsCovs_050524.csv"
event_file <- "Gazella gazella_independent_events_30min_analysis_ready.csv"
fit_file <- "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
raw_export_file <- file.path("..", "full-export_2026-04-15_2006.csv")

required_files <- c(
  juvenile_file, female_file, site_file, obs_file, event_file, fit_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "))
}

z_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
}

safe_ratio <- function(numerator, denominator) {
  ifelse(is.finite(denominator) & denominator > 0, numerator / denominator, NA_real_)
}

read_count_long <- function(path, value_name) {
  read.csv(path, na.strings = c("NA", ""), check.names = FALSE) %>%
    as_tibble() %>%
    mutate(camera_id = as.integer(camera_name)) %>%
    pivot_longer(
      cols = all_of(paste0("X", 1:7)),
      names_to = "round_label",
      values_to = value_name
    ) %>%
    mutate(round = as.integer(sub("^X", "", round_label))) %>%
    select(camera_id, round, all_of(value_name))
}

build_season_lookup <- function(obs) {
  obs <- obs %>%
    mutate(
      camera_id = as.integer(`Camera Number`),
      round = as.integer(`Round Number`),
      start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
      end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
      mid_date_observed = start_date + floor(as.numeric(end_date - start_date) / 2)
    )
  survey_start <- min(obs$start_date, na.rm = TRUE)
  survey_end <- max(obs$end_date, na.rm = TRUE)
  split_date <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)
  round_mid <- obs %>%
    filter(!is.na(mid_date_observed)) %>%
    group_by(round) %>%
    summarise(round_mid = median(as.numeric(mid_date_observed)), .groups = "drop")
  obs <- obs %>%
    left_join(round_mid, by = "round") %>%
    mutate(
      mid_date = as.Date(
        ifelse(
          is.na(mid_date_observed),
          round_mid,
          as.numeric(mid_date_observed)
        ),
        origin = "1970-01-01"
      ),
      season = factor(
        if_else(mid_date <= split_date, "Early (hotter)", "Late (mild-cold)"),
        levels = c("Early (hotter)", "Late (mild-cold)")
      ),
      period = factor(
        if_else(season == "Early (hotter)", "Early", "Late"),
        levels = c("Early", "Late")
      )
    )
  attr(obs, "split_date") <- split_date
  attr(obs, "survey_start") <- survey_start
  attr(obs, "survey_end") <- survey_end
  obs
}

juvenile <- read_count_long(juvenile_file, "juvenile")
female <- read_count_long(female_file, "female")
counts <- full_join(juvenile, female, by = c("camera_id", "round"))
if (any(xor(is.na(counts$juvenile), is.na(counts$female)))) {
  stop("Juvenile and adult-female histories have different missing-data patterns.")
}

obs_raw <- read.csv(
  obs_file,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
obs_lookup <- build_season_lookup(obs_raw) %>%
  transmute(
    camera_id,
    round,
    start_date,
    end_date,
    mid_date,
    season,
    period,
    days_operational = as.numeric(Days.Operational),
    path_type = factor(Path.Type),
    camera_model = factor(Camera.Model),
    vegetation_structure = factor(Vegetation.Structure),
    area_behind_camera = factor(Area.Behind.Camera)
  )

site_covs <- read_csv(site_file, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    Polygon.Name = factor(Polygon.Name),
    wolf_z = z_scale(Null.Wolf.Abundance),
    jackal_z = z_scale(null.jackal.abundance),
    open_z = z_scale(Open.Spaces),
    built_z = z_scale(Built.and.Roads),
    reserve_z = z_scale(Nature.Reserves),
    human_z = z_scale(human.abundance)
  )

camera_round <- counts %>%
  left_join(obs_lookup, by = c("camera_id", "round")) %>%
  left_join(site_covs, by = "camera_id") %>%
  mutate(
    valid_counts = is.finite(juvenile) & is.finite(female),
    positive_effort = is.finite(days_operational) & days_operational > 0,
    informative = valid_counts & juvenile + female > 0
  )

valid_data <- camera_round %>%
  filter(valid_counts, positive_effort)

model_data <- camera_round %>%
  filter(
    is.finite(juvenile),
    is.finite(female),
    juvenile + female > 0,
    complete.cases(
      wolf_z,
      jackal_z,
      open_z,
      built_z,
      reserve_z,
      human_z,
      season,
      path_type,
      camera_model,
      vegetation_structure,
      area_behind_camera,
      Polygon.Name
    )
  ) %>%
  arrange(camera_id, round) %>%
  mutate(
    row_id = row_number(),
    fawn = juvenile,
    camera_id = factor(camera_id),
    Polygon.Name = droplevels(Polygon.Name),
    season = droplevels(season),
    period = droplevels(period),
    path_type = droplevels(path_type),
    camera_model = droplevels(camera_model),
    vegetation_structure = droplevels(vegetation_structure),
    area_behind_camera = droplevels(area_behind_camera),
    round_factor = factor(round),
    round_centered = round - mean(round)
  )

published_fits <- readRDS(fit_file)
primary_fit <- published_fits[["seasonal_additive_stable"]]
if (!inherits(primary_fit, "glmmTMB")) {
  stop("The published season-adjusted additive model was not found.")
}
primary_frame <- model.frame(primary_fit)
primary_response <- model.response(primary_frame)
if (
  nrow(model_data) != nobs(primary_fit) ||
    !all(as.numeric(primary_response[, 1]) == as.numeric(model_data$juvenile)) ||
    !all(as.numeric(primary_response[, 2]) == as.numeric(model_data$female)) ||
    !identical(as.character(primary_frame$camera_id), as.character(model_data$camera_id)) ||
    !identical(as.character(primary_frame$season), as.character(model_data$season))
) {
  stop("The reconstructed exploratory frame does not match the published model frame.")
}

model_keys <- model_data %>%
  transmute(camera_id = as.integer(as.character(camera_id)), round, model_included = TRUE)
valid_data <- valid_data %>%
  left_join(model_keys, by = c("camera_id", "round")) %>%
  mutate(model_included = replace_na(model_included, FALSE))

summarise_camera_rounds <- function(data, grouping_variables) {
  data %>%
    group_by(across(all_of(grouping_variables))) %>%
    summarise(
      operational_cameras = n_distinct(camera_id),
      camera_nights = sum(days_operational),
      valid_camera_rounds = n(),
      informative_camera_rounds = sum(juvenile + female > 0),
      juvenile_count = sum(juvenile),
      adult_female_count = sum(female),
      juvenile_rate_per_100_camera_nights = 100 * juvenile_count / camera_nights,
      adult_female_rate_per_100_camera_nights = 100 * adult_female_count / camera_nights,
      juvenile_female_ratio = safe_ratio(juvenile_count, adult_female_count),
      juveniles_per_100_females = 100 * juvenile_female_ratio,
      model_camera_rounds = sum(model_included),
      model_cameras = n_distinct(camera_id[model_included]),
      model_camera_nights = sum(days_operational[model_included]),
      model_juvenile_count = sum(juvenile[model_included]),
      model_adult_female_count = sum(female[model_included]),
      .groups = "drop"
    )
}

round_summary <- summarise_camera_rounds(valid_data, c("round")) %>%
  left_join(
    valid_data %>%
      group_by(round) %>%
      summarise(
        period = paste(sort(unique(as.character(period))), collapse = "/"),
        round_start = min(start_date, na.rm = TRUE),
        round_end = max(end_date, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "round"
  ) %>%
  relocate(round, period, round_start, round_end)

period_summary <- summarise_camera_rounds(valid_data, c("period")) %>%
  mutate(period = as.character(period))

write_csv(round_summary, file.path(output_dir, "juvenile_female_summary_by_round.csv"))
write_csv(period_summary, file.path(output_dir, "juvenile_female_summary_by_period.csv"))

# Camera composition, concentration, and paired-period summaries.
camera_period <- valid_data %>%
  group_by(camera_id, period) %>%
  summarise(
    camera_rounds = n(),
    camera_nights = sum(days_operational),
    juvenile_count = sum(juvenile),
    adult_female_count = sum(female),
    juvenile_rate_per_100_camera_nights = 100 * juvenile_count / camera_nights,
    adult_female_rate_per_100_camera_nights = 100 * adult_female_count / camera_nights,
    juveniles_per_100_females = 100 * safe_ratio(juvenile_count, adult_female_count),
    .groups = "drop"
  )
write_csv(camera_period, file.path(output_dir, "camera_period_composition.csv"))

camera_totals <- valid_data %>%
  group_by(camera_id) %>%
  summarise(
    camera_nights = sum(days_operational),
    juvenile_count = sum(juvenile),
    adult_female_count = sum(female),
    .groups = "drop"
  )

concentration_for <- function(data, label) {
  ranked <- data %>%
    group_by(camera_id) %>%
    summarise(juvenile_count = sum(juvenile), .groups = "drop") %>%
    arrange(desc(juvenile_count), camera_id) %>%
    mutate(
      rank = row_number(),
      share = juvenile_count / sum(juvenile_count),
      cumulative_share = cumsum(share)
    )
  tibble(
    population = label,
    juvenile_count = sum(ranked$juvenile_count),
    cameras_with_juveniles = sum(ranked$juvenile_count > 0),
    top_camera_id = ranked$camera_id[1],
    top_camera_count = ranked$juvenile_count[1],
    top_camera_share = ranked$share[1],
    top_3_camera_share = sum(head(ranked$share, 3)),
    top_5_camera_share = sum(head(ranked$share, 5)),
    cameras_for_50_percent = which(ranked$cumulative_share >= 0.50)[1],
    cameras_for_80_percent = which(ranked$cumulative_share >= 0.80)[1]
  )
}

camera_concentration <- bind_rows(
  concentration_for(valid_data, "Overall"),
  concentration_for(filter(valid_data, period == "Early"), "Early"),
  concentration_for(filter(valid_data, period == "Late"), "Late")
)
write_csv(camera_concentration, file.path(output_dir, "juvenile_camera_concentration.csv"))

paired_camera_ids <- camera_period %>%
  count(camera_id, name = "periods") %>%
  filter(periods == 2) %>%
  pull(camera_id)

paired_camera_wide <- camera_period %>%
  filter(camera_id %in% paired_camera_ids) %>%
  mutate(period = as.character(period)) %>%
  pivot_wider(
    names_from = period,
    values_from = c(
      camera_rounds,
      camera_nights,
      juvenile_count,
      adult_female_count,
      juvenile_rate_per_100_camera_nights,
      adult_female_rate_per_100_camera_nights,
      juveniles_per_100_females
    ),
    names_sep = "_"
  ) %>%
  mutate(
    juvenile_rate_change =
      juvenile_rate_per_100_camera_nights_Late -
      juvenile_rate_per_100_camera_nights_Early,
    adult_female_rate_change =
      adult_female_rate_per_100_camera_nights_Late -
      adult_female_rate_per_100_camera_nights_Early,
    ratio_change = juveniles_per_100_females_Late - juveniles_per_100_females_Early
  )
write_csv(paired_camera_wide, file.path(output_dir, "paired_camera_period_changes.csv"))

period_composition_ratio <- function(data) {
  totals <- data %>%
    group_by(period) %>%
    summarise(
      juvenile = sum(juvenile),
      female = sum(female),
      ratio = juvenile / female,
      .groups = "drop"
    )
  late_ratio <- totals$ratio[totals$period == "Late"]
  early_ratio <- totals$ratio[totals$period == "Early"]
  late_ratio / early_ratio
}

observed_leave_one_camera_out <- bind_rows(
  tibble(
    omitted_camera = NA_integer_,
    observed_late_early_ratio_of_ratios = period_composition_ratio(valid_data)
  ),
  bind_rows(lapply(sort(unique(valid_data$camera_id)), function(camera_i) {
    tibble(
      omitted_camera = camera_i,
      observed_late_early_ratio_of_ratios = period_composition_ratio(
        filter(valid_data, camera_id != camera_i)
      )
    )
  }))
)
write_csv(
  observed_leave_one_camera_out,
  file.path(output_dir, "observed_leave_one_camera_out.csv")
)

# Exploratory model sensitivities. The first row is extracted from the saved
# published fit; all subsequent fits are new and are explicitly diagnostic.
control <- glmmTMBControl(
  optCtrl = list(iter.max = 10000, eval.max = 10000),
  profile = TRUE,
  rank_check = "adjust"
)

primary_formula <- formula(primary_fit)

extract_term <- function(fit, term_pattern) {
  coefficients <- as.data.frame(summary(fit)$coefficients$cond)
  term <- grep(term_pattern, rownames(coefficients), value = TRUE)
  if (length(term) != 1L) {
    return(tibble(
      term = NA_character_, estimate = NA_real_, std_error = NA_real_,
      odds_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
      p_value = NA_real_
    ))
  }
  estimate <- coefficients[term, "Estimate"]
  std_error <- coefficients[term, "Std. Error"]
  tibble(
    term = term,
    estimate = estimate,
    std_error = std_error,
    odds_ratio = exp(estimate),
    conf_low = exp(estimate - 1.96 * std_error),
    conf_high = exp(estimate + 1.96 * std_error),
    p_value = coefficients[term, "Pr(>|z|)"]
  )
}

diagnostic_row <- function(fit, label, status, term_pattern, data = NULL) {
  term <- extract_term(fit, term_pattern)
  bind_cols(
    tibble(
      analysis = label,
      status = status,
      n_observations = stats::nobs(fit),
      n_cameras = if (is.null(data)) {
        n_distinct(model_data$camera_id)
      } else {
        n_distinct(data$camera_id)
      },
      AIC = stats::AIC(fit),
      convergence_code = fit$fit$convergence,
      positive_definite_hessian = isTRUE(fit$sdr$pdHess),
      warning_or_error = NA_character_
    ),
    term
  )
}

fit_diagnostic <- function(data, fit_formula, label, term_pattern) {
  data <- droplevels(data)
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      glmmTMB(
        fit_formula,
        family = betabinomial(link = "logit"),
        data = data,
        control = control
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      fit = NULL,
      row = tibble(
        analysis = label,
        status = "Exploratory sensitivity",
        n_observations = nrow(data),
        n_cameras = n_distinct(data$camera_id),
        AIC = NA_real_,
        convergence_code = NA_integer_,
        positive_definite_hessian = FALSE,
        warning_or_error = conditionMessage(fit),
        term = NA_character_, estimate = NA_real_, std_error = NA_real_,
        odds_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
        p_value = NA_real_
      )
    ))
  }
  row <- diagnostic_row(
    fit,
    label,
    "Exploratory sensitivity",
    term_pattern,
    data
  )
  row$warning_or_error <- if (length(warnings)) {
    paste(unique(warnings), collapse = " | ")
  } else {
    NA_character_
  }
  list(fit = fit, row = row)
}

primary_row <- diagnostic_row(
  primary_fit,
  "Published primary model",
  "Published primary",
  "^seasonLate",
  model_data
)

top_cameras <- camera_totals %>%
  arrange(desc(juvenile_count), camera_id) %>%
  pull(camera_id)
top_camera_rounds <- model_data %>%
  arrange(desc(juvenile), camera_id, round) %>%
  transmute(camera_id = as.character(camera_id), round) %>%
  slice_head(n = 3)
paired_model_cameras <- model_data %>%
  distinct(camera_id, season) %>%
  count(camera_id, name = "periods") %>%
  filter(periods == 2) %>%
  pull(camera_id)

subset_definitions <- list(
  `Exclude highest-juvenile camera` = filter(
    model_data,
    as.integer(as.character(camera_id)) != top_cameras[1]
  ),
  `Exclude three highest-juvenile cameras` = filter(
    model_data,
    !as.integer(as.character(camera_id)) %in% top_cameras[1:3]
  ),
  `Exclude highest-juvenile camera-round` = anti_join(
    model_data,
    top_camera_rounds[1, ],
    by = c("camera_id", "round")
  ),
  `Exclude three highest-juvenile camera-rounds` = anti_join(
    model_data,
    top_camera_rounds,
    by = c("camera_id", "round")
  ),
  `Exclude survey round 7` = filter(model_data, round != 7),
  `Exclude survey rounds 6 and 7` = filter(model_data, !round %in% c(6, 7)),
  `Cameras represented in both periods` = filter(
    model_data,
    camera_id %in% paired_model_cameras
  )
)

subset_sensitivity <- lapply(names(subset_definitions), function(label) {
  fit_diagnostic(
    subset_definitions[[label]],
    primary_formula,
    label,
    "^seasonLate"
  )
})

message("Fitting leave-one-camera-out exploratory sensitivities...")
leave_one_model <- lapply(sort(unique(as.character(model_data$camera_id))), function(camera_i) {
  fit_diagnostic(
    filter(model_data, as.character(camera_id) != camera_i),
    primary_formula,
    paste0("Leave out camera ", camera_i),
    "^seasonLate"
  )
})

# The full adjustment structure becomes singular when camera 6 (or the three
# highest-juvenile cameras) is removed. A reduced, explicitly labeled backup
# retains the focal ecological terms and nesting while omitting sparse camera-
# design nuisance factors, allowing the direction of those exclusions to be
# inspected without presenting them as replacements for the primary model.
focused_formula <- cbind(fawn, female) ~
  wolf_z + jackal_z + open_z + season +
  (1 | Polygon.Name / camera_id)
focused_definitions <- list(
  `Focused adjustment baseline` = model_data,
  `Focused adjustment: exclude highest-juvenile camera` = filter(
    model_data,
    as.integer(as.character(camera_id)) != top_cameras[1]
  ),
  `Focused adjustment: exclude three highest-juvenile cameras` = filter(
    model_data,
    !as.integer(as.character(camera_id)) %in% top_cameras[1:3]
  )
)
focused_sensitivity <- lapply(names(focused_definitions), function(label) {
  fit_diagnostic(
    focused_definitions[[label]],
    focused_formula,
    label,
    "^seasonLate"
  )
})

model_sensitivity <- bind_rows(
  primary_row,
  bind_rows(lapply(subset_sensitivity, `[[`, "row")),
  bind_rows(lapply(focused_sensitivity, `[[`, "row")),
  bind_rows(lapply(leave_one_model, `[[`, "row"))
)
write_csv(model_sensitivity, file.path(output_dir, "season_effect_model_sensitivity.csv"))

# Continuous and factor-round sensitivity models.
no_time_formula <- update(primary_formula, . ~ . - season)
continuous_round_formula <- update(no_time_formula, . ~ . + round_centered)
factor_round_formula <- update(no_time_formula, . ~ . + round_factor)

no_time_result <- fit_diagnostic(
  model_data,
  no_time_formula,
  "No temporal term",
  "a^"
)
continuous_round_result <- fit_diagnostic(
  model_data,
  continuous_round_formula,
  "Continuous survey round",
  "^round_centered$"
)
factor_round_result <- fit_diagnostic(
  model_data,
  factor_round_formula,
  "Survey round as factor",
  "^round_factor"
)

lrt_row <- function(smaller, larger, comparison) {
  if (is.null(smaller) || is.null(larger)) {
    return(tibble(
      comparison, df_difference = NA_real_, likelihood_ratio = NA_real_,
      p_value = NA_real_, smaller_AIC = NA_real_, larger_AIC = NA_real_
    ))
  }
  smaller_log_likelihood <- logLik(smaller)
  larger_log_likelihood <- logLik(larger)
  df_difference <- attr(larger_log_likelihood, "df") -
    attr(smaller_log_likelihood, "df")
  likelihood_ratio <- 2 * (
    as.numeric(larger_log_likelihood) - as.numeric(smaller_log_likelihood)
  )
  tibble(
    comparison,
    df_difference = df_difference,
    likelihood_ratio = likelihood_ratio,
    p_value = stats::pchisq(
      likelihood_ratio,
      df = df_difference,
      lower.tail = FALSE
    ),
    smaller_AIC = AIC(smaller),
    larger_AIC = AIC(larger)
  )
}

temporal_model_tests <- bind_rows(
  lrt_row(
    no_time_result$fit,
    primary_fit,
    "Published binary Early/Late term"
  ),
  lrt_row(
    no_time_result$fit,
    continuous_round_result$fit,
    "Continuous survey-round trend"
  ),
  lrt_row(
    no_time_result$fit,
    factor_round_result$fit,
    "Any survey-round differences"
  )
)
write_csv(temporal_model_tests, file.path(output_dir, "temporal_model_tests.csv"))

round_model_coefficients <- bind_rows(
  primary_row %>% mutate(model = "Published binary Early/Late"),
  continuous_round_result$row %>% mutate(model = "Continuous survey round")
)
if (!is.null(factor_round_result$fit)) {
  factor_coefficients <- as.data.frame(
    summary(factor_round_result$fit)$coefficients$cond
  ) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    filter(grepl("^round_factor", term)) %>%
    transmute(
      model = "Survey round as factor",
      term,
      estimate = Estimate,
      std_error = `Std. Error`,
      odds_ratio = exp(estimate),
      conf_low = exp(estimate - 1.96 * std_error),
      conf_high = exp(estimate + 1.96 * std_error),
      p_value = `Pr(>|z|)`
    )
} else {
  factor_coefficients <- tibble()
}
write_csv(
  bind_rows(
    round_model_coefficients %>%
      select(model, term, estimate, std_error, odds_ratio, conf_low, conf_high, p_value),
    factor_coefficients
  ),
  file.path(output_dir, "temporal_model_coefficients.csv")
)

# Date-cutoff sensitivity around the 17 October split.
cutoff_dates <- as.Date(c(
  "2020-10-03", "2020-10-10", "2020-10-17", "2020-10-24", "2020-10-31"
))
cutoff_formula <- update(no_time_formula, . ~ . + period_alt)
cutoff_results <- lapply(cutoff_dates, function(cutoff_date) {
  cutoff_data <- model_data %>%
    mutate(
      period_alt = factor(
        if_else(mid_date <= cutoff_date, "Early", "Late"),
        levels = c("Early", "Late")
      )
    )
  result <- fit_diagnostic(
    cutoff_data,
    cutoff_formula,
    paste0("Cutoff ", cutoff_date),
    "^period_altLate$"
  )
  result$row %>%
    mutate(
      cutoff_date = cutoff_date,
      early_rows = sum(cutoff_data$period_alt == "Early"),
      late_rows = sum(cutoff_data$period_alt == "Late")
    )
})
cutoff_sensitivity <- bind_rows(cutoff_results)
write_csv(cutoff_sensitivity, file.path(output_dir, "period_cutoff_sensitivity.csv"))

# Cross-class 30-minute episode and group-size summaries.
assign_active_round <- function(event_camera, event_date, active_obs) {
  index <- which(
    active_obs$camera_id == event_camera &
      event_date >= active_obs$start_date &
      event_date <= active_obs$end_date
  )
  if (!length(index)) return(NA_integer_)
  index <- index[which.max(active_obs$start_date[index])]
  active_obs$round[index]
}

events <- read.csv(
  event_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(
    camera_id = as.integer(camera_name),
    session_id = as.integer(trap_station_session_id),
    timestamp = as.POSIXct(
      sub("Z$", "", event_time),
      format = "%Y-%m-%dT%H:%M:%S",
      tz = "Asia/Jerusalem"
    ),
    event_date = as.Date(format(timestamp, tz = "Asia/Jerusalem")),
    event_count = as.numeric(event_count),
    event_row_id = row_number()
  )

active_obs <- obs_lookup %>%
  filter(
    is.finite(days_operational),
    days_operational > 0,
    !is.na(start_date),
    !is.na(end_date)
  )
events$round <- mapply(
  assign_active_round,
  event_camera = events$camera_id,
  event_date = events$event_date,
  MoreArgs = list(active_obs = active_obs)
)
events <- events %>%
  left_join(
    obs_lookup %>% select(camera_id, round, period),
    by = c("camera_id", "round")
  )

focal_events <- events %>%
  filter(
    gazelle_class %in% c("juvenile", "female"),
    !is.na(timestamp),
    is.finite(event_count),
    !is.na(round),
    !is.na(period)
  ) %>%
  arrange(camera_id, session_id, timestamp) %>%
  group_by(camera_id, session_id) %>%
  mutate(
    gap_minutes = as.numeric(difftime(timestamp, lag(timestamp), units = "mins")),
    cross_class_episode = cumsum(is.na(gap_minutes) | gap_minutes > 30)
  ) %>%
  ungroup()

episodes <- focal_events %>%
  group_by(camera_id, session_id, cross_class_episode) %>%
  summarise(
    round = first(round),
    period = first(period),
    episode_start = min(timestamp, na.rm = TRUE),
    episode_end = max(timestamp, na.rm = TRUE),
    juvenile_present = any(gazelle_class == "juvenile"),
    adult_female_present = any(gazelle_class == "female"),
    juvenile_count = if (juvenile_present) {
      max(event_count[gazelle_class == "juvenile"], na.rm = TRUE)
    } else 0,
    adult_female_count = if (adult_female_present) {
      max(event_count[gazelle_class == "female"], na.rm = TRUE)
    } else 0,
    .groups = "drop"
  ) %>%
  mutate(
    episode_type = case_when(
      juvenile_present & adult_female_present ~ "Juvenile and adult female",
      juvenile_present ~ "Juvenile only",
      TRUE ~ "Adult female only"
    )
  )

episode_summary_by_round <- episodes %>%
  group_by(round) %>%
  summarise(
    juvenile_episodes = sum(juvenile_present),
    juvenile_with_adult_female = sum(juvenile_present & adult_female_present),
    juvenile_without_adult_female = sum(juvenile_present & !adult_female_present),
    adult_female_only = sum(!juvenile_present & adult_female_present),
    proportion_juvenile_episodes_with_adult_female =
      juvenile_with_adult_female / juvenile_episodes,
    .groups = "drop"
  )
episode_summary_by_period <- episodes %>%
  group_by(period) %>%
  summarise(
    juvenile_episodes = sum(juvenile_present),
    juvenile_with_adult_female = sum(juvenile_present & adult_female_present),
    juvenile_without_adult_female = sum(juvenile_present & !adult_female_present),
    adult_female_only = sum(!juvenile_present & adult_female_present),
    proportion_juvenile_episodes_with_adult_female =
      juvenile_with_adult_female / juvenile_episodes,
    .groups = "drop"
  )

exact_timestamp_overlap <- focal_events %>%
  group_by(camera_id, session_id, timestamp) %>%
  summarise(
    round = first(round),
    period = first(period),
    juvenile_present = any(gazelle_class == "juvenile"),
    adult_female_present = any(gazelle_class == "female"),
    .groups = "drop"
  )
exact_timestamp_overlap_summary <- exact_timestamp_overlap %>%
  group_by(period) %>%
  summarise(
    exact_shared_timestamps = sum(juvenile_present & adult_female_present),
    .groups = "drop"
  )

event_group_size <- focal_events %>%
  mutate(
    demographic_class = recode(
      gazelle_class,
      juvenile = "Juvenile",
      female = "Adult female"
    )
  ) %>%
  group_by(round, period, demographic_class) %>%
  summarise(
    independent_event_rows = n(),
    summed_animals = sum(event_count),
    mean_group_size = mean(event_count),
    median_group_size = median(event_count),
    maximum_group_size = max(event_count, na.rm = TRUE),
    .groups = "drop"
  )
event_group_size_period <- focal_events %>%
  mutate(
    demographic_class = recode(
      gazelle_class,
      juvenile = "Juvenile",
      female = "Adult female"
    )
  ) %>%
  group_by(period, demographic_class) %>%
  summarise(
    independent_event_rows = n(),
    summed_animals = sum(event_count),
    mean_group_size = mean(event_count),
    median_group_size = median(event_count),
    maximum_group_size = max(event_count, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(episode_summary_by_round, file.path(output_dir, "cross_class_episodes_by_round.csv"))
write_csv(episode_summary_by_period, file.path(output_dir, "cross_class_episodes_by_period.csv"))
write_csv(event_group_size, file.path(output_dir, "event_group_sizes_by_round.csv"))
write_csv(event_group_size_period, file.path(output_dir, "event_group_sizes_by_period.csv"))
write_csv(
  exact_timestamp_overlap_summary,
  file.path(output_dir, "exact_timestamp_overlap_by_period.csv")
)

# Existing metadata contain broad life-stage labels but no finer juvenile age
# category, body-size field, or image path. Summarize those labels by period.
metadata_availability <- tibble(
  field = c(
    "Sex", "Life stage", "Sighting Quantity", "Date/Time",
    "Body-size or juvenile-age notes", "Image filename or path"
  ),
  available = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  implication = c(
    "Available as broad categorical metadata",
    "Adult, Juvenile, Sub-adult, or missing; no finer juvenile ages",
    "Available but not an age estimate",
    "Available for seasonal assignment",
    "Cannot test juvenile maturation or identify neonates",
    "Targeted photograph review cannot be linked reproducibly from this export"
  )
)

classification_metadata <- tibble()
if (file.exists(raw_export_file)) {
  raw <- fread(
    raw_export_file,
    select = c(
      "Camera Name", "Date/Time", "Sighting Quantity", "Species",
      "Trap Station Session ID", "Sex", "Life stage"
    ),
    showProgress = FALSE
  )
  raw <- raw[Species == "Gazella gazella"]
  raw[, timestamp := as.POSIXct(
    `Date/Time`,
    format = "%d/%m/%Y %H:%M",
    tz = "Asia/Jerusalem"
  )]
  raw[, event_date := as.Date(format(timestamp, tz = "Asia/Jerusalem"))]
  raw[, round := mapply(
    assign_active_round,
    event_camera = as.integer(`Camera Name`),
    event_date = event_date,
    MoreArgs = list(active_obs = active_obs)
  )]
  raw <- merge(
    raw,
    as.data.table(obs_lookup %>% select(camera_id, round, period)),
    by.x = c("Camera Name", "round"),
    by.y = c("camera_id", "round"),
    all.x = TRUE,
    sort = FALSE
  )
  raw[is.na(Sex) | Sex == "", Sex := "Missing"]
  raw[is.na(`Life stage`) | `Life stage` == "", `Life stage` := "Missing"]
  classification_metadata <- raw[
    !is.na(period),
    .(
      raw_sighting_rows = .N,
      summed_sighting_quantity = sum(as.numeric(`Sighting Quantity`), na.rm = TRUE),
      cameras = uniqueN(`Camera Name`)
    ),
    by = .(period, sex = Sex, life_stage = `Life stage`)
  ] %>%
    as_tibble() %>%
    arrange(period, desc(raw_sighting_rows))
}
write_csv(
  classification_metadata,
  file.path(output_dir, "classification_metadata_by_period.csv")
)
write_csv(
  metadata_availability,
  file.path(output_dir, "juvenile_age_metadata_availability.csv")
)

# Diagnostic figure 1: numerator, denominator, and ratio through time.
round_plot_data <- bind_rows(
  round_summary %>%
    transmute(
      round,
      metric = "Juvenile encounter rate\n(per 100 camera-nights)",
      value = juvenile_rate_per_100_camera_nights
    ),
  round_summary %>%
    transmute(
      round,
      metric = "Adult-female encounter rate\n(per 100 camera-nights)",
      value = adult_female_rate_per_100_camera_nights
    ),
  round_summary %>%
    transmute(
      round,
      metric = "Observed juveniles per\n100 adult females",
      value = juveniles_per_100_females
    )
) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "Juvenile encounter rate\n(per 100 camera-nights)",
        "Adult-female encounter rate\n(per 100 camera-nights)",
        "Observed juveniles per\n100 adult females"
      )
    )
  )

round_plot <- ggplot(round_plot_data, aes(round, value)) +
  geom_line(linewidth = 0.8, color = "#1F4E5F") +
  geom_point(size = 2.3, shape = 21, fill = "white", color = "#1F4E5F") +
  geom_text(
    aes(label = sprintf("%.1f", value)),
    vjust = -0.65,
    size = 2.8,
    color = "#263238"
  ) +
  facet_grid(metric ~ ., scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = 1:7, limits = c(0.8, 7.2)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.20))) +
  labs(
    x = "Survey round",
    y = NULL,
    title = "Juvenile and adult-female detections across survey rounds",
    subtitle = paste0(
      "Early/Late was assigned from each camera-round midpoint; ",
      "rounds 3 and 4 contain both periods"
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "juvenile_female_round_diagnostics.png"),
  round_plot,
  width = 8.2,
  height = 7.2,
  dpi = 360
)
ggsave(
  file.path(output_dir, "juvenile_female_round_diagnostics.pdf"),
  round_plot,
  width = 8.2,
  height = 7.2
)

# Diagnostic figure 2: paired operational cameras.
paired_plot_data <- camera_period %>%
  filter(camera_id %in% paired_camera_ids) %>%
  select(
    camera_id,
    period,
    juvenile_rate_per_100_camera_nights,
    adult_female_rate_per_100_camera_nights,
    juveniles_per_100_females
  ) %>%
  pivot_longer(
    cols = c(
      juvenile_rate_per_100_camera_nights,
      adult_female_rate_per_100_camera_nights,
      juveniles_per_100_females
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      juvenile_rate_per_100_camera_nights = "Juvenile rate per 100 camera-nights",
      adult_female_rate_per_100_camera_nights = "Adult-female rate per 100 camera-nights",
      juveniles_per_100_females = "Juveniles per 100 adult females"
    )
  ) %>%
  filter(is.finite(value))

paired_medians <- paired_plot_data %>%
  group_by(metric, period) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")

paired_plot <- ggplot(
  paired_plot_data,
  aes(period, value, group = camera_id)
) +
  geom_line(color = "#7B8588", alpha = 0.28, linewidth = 0.35) +
  geom_point(color = "#7B8588", alpha = 0.45, size = 1.25) +
  geom_line(
    data = paired_medians,
    aes(x = period, y = value, group = 1),
    inherit.aes = FALSE,
    color = "#B34135",
    linewidth = 1.15
  ) +
  geom_point(
    data = paired_medians,
    aes(x = period, y = value),
    inherit.aes = FALSE,
    color = "#B34135",
    fill = "white",
    shape = 21,
    size = 2.7,
    stroke = 0.8
  ) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_y_sqrt() +
  labs(
    x = NULL,
    y = NULL,
    title = "Within-camera Early-to-Late changes",
    subtitle = paste0(
      "Thin lines are cameras operating in both periods; red lines connect period medians; ",
      "square-root y-scales"
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "paired_camera_period_diagnostics.png"),
  paired_plot,
  width = 10.4,
  height = 4.5,
  dpi = 360
)
ggsave(
  file.path(output_dir, "paired_camera_period_diagnostics.pdf"),
  paired_plot,
  width = 10.4,
  height = 4.5
)

# Machine-readable audit summary for report generation and independent review.
audit_metrics <- tibble(
  metric = c(
    "Survey midpoint used by source code",
    "All positive-effort camera-rounds",
    "All positive-effort camera-nights",
    "All informative camera-rounds",
    "Primary model camera-rounds",
    "Primary model cameras",
    "Primary model juvenile count",
    "Primary model adult-female count",
    "Cameras operating in both periods",
    "Mapped juvenile/female class-event rows",
    "Class-event rows outside documented active dates",
    "Exact timestamps containing both classes"
  ),
  value = c(
    as.character(attr(build_season_lookup(obs_raw), "split_date")),
    sum(valid_data$positive_effort),
    sum(valid_data$days_operational),
    sum(valid_data$juvenile + valid_data$female > 0),
    nrow(model_data),
    n_distinct(model_data$camera_id),
    sum(model_data$juvenile),
    sum(model_data$female),
    length(paired_camera_ids),
    nrow(focal_events),
    sum(
      events$gazelle_class %in% c("juvenile", "female") & is.na(events$round)
    ),
    sum(
      exact_timestamp_overlap$juvenile_present &
        exact_timestamp_overlap$adult_female_present
    )
  )
)
write_csv(audit_metrics, file.path(output_dir, "audit_metrics.csv"))

message("Exploratory seasonal audit outputs written to: ", output_dir)
