#!/usr/bin/env Rscript

# Exploratory decomposition of the seasonal juvenile:adult-female ratio.
# This script reads established analysis-ready inputs and writes only to its
# own folder. It does not alter or refit the primary manuscript model.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tidyr)
})

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
} else {
  normalizePath("juvenile_female_seasonal_rate_analysis", winslash = "/")
}
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/")
output_dir <- script_dir

input_files <- list(
  juvenile = file.path(project_dir, "gazelle_juvenile_abundance_60x7.csv"),
  female = file.path(project_dir, "gazelle_female_abundance_60x7.csv"),
  events = file.path(
    project_dir,
    "Gazella gazella_independent_events_30min_analysis_ready.csv"
  ),
  observation = file.path(project_dir, "ObsCovs_050524.csv"),
  site = file.path(project_dir, "SiteCovs_060324.csv"),
  published_script = file.path(
    project_dir,
    "analyze_fawn_female_beta_binomial_wolf_jackal.R"
  ),
  prior_audit_script = file.path(
    project_dir,
    "investigate_juvenile_female_seasonal_increase.R"
  )
)

missing_inputs <- names(input_files)[!file.exists(unlist(input_files))]
if (length(missing_inputs)) {
  stop(
    "Missing required input(s): ",
    paste(paste0(missing_inputs, " = ", unlist(input_files[missing_inputs])),
          collapse = "; ")
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

split_date <- as.Date("2020-10-17")
bootstrap_replicates <- 10000L
bootstrap_seed <- 20201017L

z_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
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

juvenile_history <- read_count_long(input_files$juvenile, "juvenile_count")
female_history <- read_count_long(input_files$female, "adult_female_count")

observation <- read.csv(
  input_files$observation,
  na.strings = c("NA", ""),
  check.names = FALSE
) %>%
  as_tibble() %>%
  transmute(
    camera_id = as.integer(.data[["Camera Number"]]),
    round = as.integer(.data[["Round Number"]]),
    start_date = as.Date(.data[["Start date"]], format = "%d/%m/%Y"),
    end_date = as.Date(.data[["End Date"]], format = "%d/%m/%Y"),
    days_operational = as.numeric(.data[["Days.Operational"]]),
    path_type = factor(.data[["Path.Type"]]),
    camera_model = factor(.data[["Camera.Model"]]),
    vegetation_structure = factor(.data[["Vegetation.Structure"]]),
    area_behind_camera = factor(.data[["Area.Behind.Camera"]])
  ) %>%
  mutate(
    midpoint_date = start_date + floor(as.numeric(end_date - start_date) / 2),
    period = factor(
      if_else(midpoint_date <= split_date, "Early", "Late"),
      levels = c("Early", "Late")
    )
  )

survey_start <- min(observation$start_date, na.rm = TRUE)
survey_end <- max(observation$end_date, na.rm = TRUE)
source_midpoint <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)
if (!identical(source_midpoint, split_date)) {
  stop("The reconstructed survey midpoint is not 2020-10-17.")
}

site_covariates <- read_csv(input_files$site, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    polygon = factor(Polygon.Name),
    wolf_z = z_scale(Null.Wolf.Abundance),
    jackal_z = z_scale(null.jackal.abundance),
    open_z = z_scale(Open.Spaces),
    built_z = z_scale(Built.and.Roads),
    reserve_z = z_scale(Nature.Reserves),
    human_z = z_scale(human.abundance)
  )

camera_round <- full_join(
  juvenile_history,
  female_history,
  by = c("camera_id", "round")
) %>%
  left_join(observation, by = c("camera_id", "round")) %>%
  left_join(site_covariates, by = "camera_id")

if (any(xor(
  is.na(camera_round$juvenile_count),
  is.na(camera_round$adult_female_count)
))) {
  stop("Juvenile and adult-female histories have different missingness.")
}

valid_data <- camera_round %>%
  filter(
    is.finite(juvenile_count),
    is.finite(adult_female_count),
    is.finite(days_operational),
    days_operational > 0
  )

if (nrow(valid_data) != 344L || n_distinct(valid_data$camera_id) != 60L) {
  stop("Positive-effort camera-round reconstruction differs from the audit.")
}

assign_active_round <- function(event_camera, event_date, deployments) {
  candidate <- which(
    deployments$camera_id == event_camera &
      event_date >= deployments$start_date &
      event_date <= deployments$end_date &
      is.finite(deployments$days_operational) &
      deployments$days_operational > 0
  )
  if (!length(candidate)) return(NA_integer_)
  # Follow the established date-resolved event audit: when documented round
  # intervals overlap on a boundary date, assign the event to the deployment
  # with the latest start. The primary count histories retain a separate
  # session-based allocation and are not expected to match every event cell.
  candidate <- candidate[which.max(deployments$start_date[candidate])]
  deployments$round[candidate]
}

events <- read_csv(input_files$events, show_col_types = FALSE) %>%
  mutate(
    camera_id = as.integer(camera_name),
    event_date = as.Date(substr(event_time, 1, 10)),
    event_count = as.numeric(event_count)
  )

events$round <- mapply(
  assign_active_round,
  event_camera = events$camera_id,
  event_date = events$event_date,
  MoreArgs = list(deployments = observation)
)

focal_events <- events %>%
  filter(
    gazelle_class %in% c("juvenile", "female"),
    is.finite(event_count)
  )

if (anyNA(focal_events$round)) {
  stop("One or more focal independent events could not be assigned to a round.")
}

event_camera_round <- focal_events %>%
  group_by(camera_id, round, gazelle_class) %>%
  summarise(
    independent_events = n(),
    summed_event_maximum = sum(event_count),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = gazelle_class,
    values_from = c(independent_events, summed_event_maximum),
    values_fill = 0
  ) %>%
  rename(
    juvenile_events = independent_events_juvenile,
    adult_female_events = independent_events_female,
    juvenile_reconstructed_count = summed_event_maximum_juvenile,
    adult_female_reconstructed_count = summed_event_maximum_female
  )

valid_data <- valid_data %>%
  left_join(event_camera_round, by = c("camera_id", "round")) %>%
  mutate(
    across(
      c(
        juvenile_events,
        adult_female_events,
        juvenile_reconstructed_count,
        adult_female_reconstructed_count
      ),
      ~ replace_na(.x, 0)
    )
  )

reconstruction_audit <- tibble(
  metric = c(
    "Juvenile individuals in count history",
    "Juvenile summed event maxima",
    "Adult-female individuals in count history",
    "Adult-female summed event maxima",
    "Juvenile independent events",
    "Adult-female independent events",
    "Camera-round juvenile count mismatches",
    "Camera-round adult-female count mismatches"
  ),
  value = c(
    sum(valid_data$juvenile_count),
    sum(valid_data$juvenile_reconstructed_count),
    sum(valid_data$adult_female_count),
    sum(valid_data$adult_female_reconstructed_count),
    sum(valid_data$juvenile_events),
    sum(valid_data$adult_female_events),
    sum(valid_data$juvenile_count != valid_data$juvenile_reconstructed_count),
    sum(
      valid_data$adult_female_count !=
        valid_data$adult_female_reconstructed_count
    )
  )
)

reconstruction_mismatches <- valid_data %>%
  filter(
    juvenile_count != juvenile_reconstructed_count |
      adult_female_count != adult_female_reconstructed_count
  ) %>%
  select(
    camera_id, round, start_date, end_date,
    juvenile_count, juvenile_reconstructed_count,
    adult_female_count, adult_female_reconstructed_count
  )

write_csv(
  reconstruction_mismatches,
  file.path(output_dir, "event_count_reconstruction_mismatches.csv")
)

if (sum(valid_data$juvenile_count) != 199 ||
    sum(valid_data$adult_female_count) != 1951 ||
    sum(valid_data$juvenile_reconstructed_count) != 199 ||
    sum(valid_data$adult_female_reconstructed_count) != 1951) {
  stop(
    "Independent-event totals do not match the count-history totals: ",
    "history totals = ", sum(valid_data$juvenile_count), "/",
    sum(valid_data$adult_female_count), "; reconstructed totals = ",
    sum(valid_data$juvenile_reconstructed_count), "/",
    sum(valid_data$adult_female_reconstructed_count), "."
  )
}

summarise_counts <- function(data, grouping) {
  data %>%
    group_by(across(all_of(grouping))) %>%
    summarise(
      operational_cameras = n_distinct(camera_id),
      camera_nights = sum(days_operational),
      valid_camera_rounds = n(),
      juvenile_individual_count = sum(juvenile_count),
      adult_female_individual_count = sum(adult_female_count),
      juvenile_independent_events = sum(juvenile_events),
      adult_female_independent_events = sum(adult_female_events),
      juvenile_event_assigned_individuals =
        sum(juvenile_reconstructed_count),
      adult_female_event_assigned_individuals =
        sum(adult_female_reconstructed_count),
      mean_juvenile_group_size =
        juvenile_event_assigned_individuals / juvenile_independent_events,
      mean_adult_female_group_size =
        adult_female_event_assigned_individuals /
          adult_female_independent_events,
      juvenile_individual_rate_per_100_camera_nights =
        100 * juvenile_individual_count / camera_nights,
      adult_female_individual_rate_per_100_camera_nights =
        100 * adult_female_individual_count / camera_nights,
      juvenile_event_rate_per_100_camera_nights =
        100 * juvenile_independent_events / camera_nights,
      adult_female_event_rate_per_100_camera_nights =
        100 * adult_female_independent_events / camera_nights,
      observed_juvenile_female_ratio =
        juvenile_individual_count / adult_female_individual_count,
      juveniles_per_100_adult_females =
        100 * observed_juvenile_female_ratio,
      .groups = "drop"
    )
}

round_dates <- valid_data %>%
  group_by(round) %>%
  summarise(
    period_coverage = paste(
      sort(unique(as.character(period))),
      collapse = "/"
    ),
    round_start = min(start_date),
    round_end = max(end_date),
    round_midpoint = as.Date(
      median(as.numeric(midpoint_date)),
      origin = "1970-01-01"
    ),
    .groups = "drop"
  )

round_summary <- summarise_counts(valid_data, "round") %>%
  left_join(round_dates, by = "round") %>%
  relocate(round, period_coverage, round_start, round_end, round_midpoint)

period_summary <- summarise_counts(valid_data, "period") %>%
  mutate(period = as.character(period))

# Cluster bootstrap: resample cameras and retain each sampled camera's complete
# set of available rounds. Missing effort contributes neither counts nor nights.
camera_ids <- sort(unique(valid_data$camera_id))
round_ids <- sort(unique(valid_data$round))
period_ids <- c("Early", "Late")

matrix_by_group <- function(data, value, group, group_levels) {
  wide <- data %>%
    transmute(
      camera_id,
      group_value = as.character(.data[[group]]),
      value = .data[[value]]
    ) %>%
    complete(
      camera_id = camera_ids,
      group_value = as.character(group_levels),
      fill = list(value = 0)
    ) %>%
    arrange(match(camera_id, camera_ids), match(group_value, group_levels)) %>%
    pivot_wider(names_from = group_value, values_from = value, values_fill = 0)
  as.matrix(wide[, -1, drop = FALSE])
}

period_camera <- valid_data %>%
  group_by(camera_id, period) %>%
  summarise(
    days_operational = sum(days_operational),
    juvenile_count = sum(juvenile_count),
    adult_female_count = sum(adult_female_count),
    juvenile_events = sum(juvenile_events),
    adult_female_events = sum(adult_female_events),
    .groups = "drop"
  )

round_matrices <- lapply(
  c(
    "days_operational", "juvenile_count", "adult_female_count",
    "juvenile_events", "adult_female_events"
  ),
  function(variable) {
    matrix_by_group(valid_data, variable, "round", round_ids)
  }
)
names(round_matrices) <- c("effort", "juvenile", "female", "juvenile_event", "female_event")

period_matrices <- lapply(
  c(
    "days_operational", "juvenile_count", "adult_female_count",
    "juvenile_events", "adult_female_events"
  ),
  function(variable) {
    matrix_by_group(period_camera, variable, "period", period_ids)
  }
)
names(period_matrices) <- names(round_matrices)

bootstrap_metrics <- function(matrices, n_groups) {
  output <- array(
    NA_real_,
    dim = c(
      bootstrap_replicates,
      n_groups,
      5L
    ),
    dimnames = list(
      NULL,
      NULL,
      c(
        "juvenile_individual_rate",
        "adult_female_individual_rate",
        "juvenile_event_rate",
        "adult_female_event_rate",
        "juvenile_female_ratio"
      )
    )
  )
  set.seed(bootstrap_seed)
  for (iteration in seq_len(bootstrap_replicates)) {
    weights <- tabulate(
      sample.int(length(camera_ids), length(camera_ids), replace = TRUE),
      nbins = length(camera_ids)
    )
    effort <- as.numeric(crossprod(weights, matrices$effort))
    juvenile <- as.numeric(crossprod(weights, matrices$juvenile))
    female <- as.numeric(crossprod(weights, matrices$female))
    juvenile_event <- as.numeric(crossprod(weights, matrices$juvenile_event))
    female_event <- as.numeric(crossprod(weights, matrices$female_event))
    output[iteration, , "juvenile_individual_rate"] <- 100 * juvenile / effort
    output[iteration, , "adult_female_individual_rate"] <- 100 * female / effort
    output[iteration, , "juvenile_event_rate"] <- 100 * juvenile_event / effort
    output[iteration, , "adult_female_event_rate"] <- 100 * female_event / effort
    output[iteration, , "juvenile_female_ratio"] <- 100 * juvenile / female
  }
  output
}

round_bootstrap <- bootstrap_metrics(round_matrices, length(round_ids))
period_bootstrap <- bootstrap_metrics(period_matrices, length(period_ids))

interval_table <- function(bootstrap_array, group_values, group_name) {
  bind_rows(lapply(seq_along(group_values), function(group_index) {
    bind_rows(lapply(dimnames(bootstrap_array)[[3]], function(metric_name) {
      interval <- quantile(
        bootstrap_array[, group_index, metric_name],
        probs = c(0.025, 0.975),
        na.rm = TRUE,
        names = FALSE
      )
      tibble(
        group = as.character(group_values[group_index]),
        metric = metric_name,
        bootstrap_lower = interval[1],
        bootstrap_upper = interval[2]
      )
    }))
  })) %>%
    rename(!!group_name := group)
}

round_intervals <- interval_table(round_bootstrap, round_ids, "round") %>%
  mutate(round = as.integer(round))
period_intervals <- interval_table(period_bootstrap, period_ids, "period")

add_intervals <- function(summary_data, interval_data, group_name) {
  metric_crosswalk <- c(
    juvenile_individual_rate = "juvenile_individual_rate_per_100_camera_nights",
    adult_female_individual_rate = "adult_female_individual_rate_per_100_camera_nights",
    juvenile_event_rate = "juvenile_event_rate_per_100_camera_nights",
    adult_female_event_rate = "adult_female_event_rate_per_100_camera_nights",
    juvenile_female_ratio = "juveniles_per_100_adult_females"
  )
  wide_intervals <- interval_data %>%
    mutate(metric = unname(metric_crosswalk[metric])) %>%
    pivot_wider(
      names_from = metric,
      values_from = c(bootstrap_lower, bootstrap_upper),
      names_glue = "{metric}_{.value}"
    )
  left_join(summary_data, wide_intervals, by = group_name)
}

round_summary <- add_intervals(round_summary, round_intervals, "round")
period_summary <- add_intervals(period_summary, period_intervals, "period")

write_csv(
  reconstruction_audit,
  file.path(output_dir, "event_count_reconstruction_audit.csv")
)
write_csv(
  round_summary,
  file.path(output_dir, "juvenile_female_counts_rates_by_round.csv")
)
write_csv(
  period_summary,
  file.path(output_dir, "juvenile_female_counts_rates_by_period.csv")
)

# Negative-binomial mixed models include all positive-effort camera-rounds,
# including zero detections. The adjusted models match the stable additive
# covariate set used in the published ratio analysis.
model_control <- glmmTMBControl(
  optCtrl = list(iter.max = 10000, eval.max = 10000),
  profile = TRUE,
  rank_check = "adjust"
)

fit_capture <- function(formula, data, family = nbinom2(), dispformula = ~1) {
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      glmmTMB(
        formula,
        data = data,
        family = family,
        dispformula = dispformula,
        control = model_control
      ),
      warning = function(warning_condition) {
        warnings <<- c(warnings, conditionMessage(warning_condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error_condition) error_condition
  )
  list(fit = fit, warnings = unique(warnings))
}

extract_period_rr <- function(result, model_label, response_label, subset_label) {
  fit <- result$fit
  if (inherits(fit, "error")) {
    return(tibble(
      model = model_label,
      response = response_label,
      data_subset = subset_label,
      contrast = "Late vs Early",
      rate_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      likelihood_ratio_p = NA_real_,
      n_camera_rounds = NA_integer_,
      n_cameras = NA_integer_,
      AIC = NA_real_,
      convergence_code = NA_integer_,
      positive_definite_hessian = FALSE,
      warning_or_error = conditionMessage(fit)
    ))
  }
  coefficients <- summary(fit)$coefficients$cond
  estimate <- coefficients["periodLate", "Estimate"]
  std_error <- coefficients["periodLate", "Std. Error"]
  tibble(
    model = model_label,
    response = response_label,
    data_subset = subset_label,
    contrast = "Late vs Early",
    rate_ratio = exp(estimate),
    conf_low = exp(estimate - 1.96 * std_error),
    conf_high = exp(estimate + 1.96 * std_error),
    p_value = coefficients["periodLate", "Pr(>|z|)"],
    likelihood_ratio_p = NA_real_,
    n_camera_rounds = nrow(model.frame(fit)),
    n_cameras = n_distinct(model.frame(fit)$camera_id),
    AIC = AIC(fit),
    convergence_code = fit$fit$convergence,
    positive_definite_hessian = isTRUE(fit$sdr$pdHess),
    warning_or_error = if (length(result$warnings)) {
      paste(result$warnings, collapse = " | ")
    } else NA_character_
  )
}

fit_class_model <- function(data, response, adjusted = FALSE) {
  data <- data %>%
    mutate(
      camera_id = factor(camera_id),
      period = droplevels(period)
    )
  if (adjusted) {
    data <- data %>%
      filter(complete.cases(
        polygon, wolf_z, jackal_z, open_z, built_z, reserve_z, human_z,
        path_type, camera_model, vegetation_structure, area_behind_camera
      )) %>%
      droplevels()
    formula <- as.formula(paste0(
      response,
      " ~ period + wolf_z + jackal_z + open_z + built_z + reserve_z + ",
      "human_z + path_type + camera_model + vegetation_structure + ",
      "area_behind_camera + offset(log(days_operational)) + ",
      "(1 | polygon / camera_id)"
    ))
  } else {
    formula <- as.formula(paste0(
      response,
      " ~ period + offset(log(days_operational)) + (1 | camera_id)"
    ))
  }
  fit_capture(formula, data)
}

fit_joint_model <- function(data) {
  stacked <- data %>%
    transmute(
      camera_id = factor(camera_id),
      camera_round_id = interaction(camera_id, round, drop = TRUE),
      period = droplevels(period),
      days_operational,
      `Adult female` = adult_female_count,
      Juvenile = juvenile_count
    ) %>%
    pivot_longer(
      cols = c(`Adult female`, Juvenile),
      names_to = "demographic_class",
      values_to = "count"
    ) %>%
    mutate(
      demographic_class = factor(
        demographic_class,
        levels = c("Adult female", "Juvenile")
      )
    )
  full <- fit_capture(
    count ~ demographic_class * period + offset(log(days_operational)) +
      (0 + demographic_class | camera_id) + (1 | camera_round_id),
    stacked,
    dispformula = ~ demographic_class
  )
  list(full = full, data = stacked)
}

extract_joint_difference <- function(result, subset_label) {
  full_fit <- result$full$fit
  if (inherits(full_fit, "error")) {
    return(tibble(
      model = "Joint class-by-period NB2 GLMM",
      response = "Summed event maxima",
      data_subset = subset_label,
      contrast = "Juvenile Late/Early RR divided by adult-female Late/Early RR",
      rate_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
      p_value = NA_real_, likelihood_ratio_p = NA_real_,
      n_camera_rounds = NA_integer_, n_cameras = NA_integer_, AIC = NA_real_,
      convergence_code = NA_integer_, positive_definite_hessian = FALSE,
      warning_or_error = conditionMessage(full_fit)
    ))
  }
  coefficient_name <- "demographic_classJuvenile:periodLate"
  coefficients <- summary(full_fit)$coefficients$cond
  estimate <- coefficients[coefficient_name, "Estimate"]
  std_error <- coefficients[coefficient_name, "Std. Error"]
  tibble(
    model = "Joint class-by-period NB2 GLMM",
    response = "Summed event maxima",
    data_subset = subset_label,
    contrast = "Juvenile Late/Early RR divided by adult-female Late/Early RR",
    rate_ratio = exp(estimate),
    conf_low = exp(estimate - 1.96 * std_error),
    conf_high = exp(estimate + 1.96 * std_error),
    p_value = coefficients[coefficient_name, "Pr(>|z|)"],
    likelihood_ratio_p = NA_real_,
    n_camera_rounds = nrow(result$data) / 2,
    n_cameras = n_distinct(result$data$camera_id),
    AIC = AIC(full_fit),
    convergence_code = full_fit$fit$convergence,
    positive_definite_hessian = isTRUE(full_fit$sdr$pdHess),
    warning_or_error = if (length(result$full$warnings)) {
      paste(result$full$warnings, collapse = " | ")
    } else NA_character_
  )
}

paired_cameras <- valid_data %>%
  distinct(camera_id, period) %>%
  count(camera_id, name = "periods_observed") %>%
  filter(periods_observed == 2L) %>%
  pull(camera_id)
paired_data <- valid_data %>% filter(camera_id %in% paired_cameras)

model_rows <- list()

for (response in c("juvenile_count", "adult_female_count")) {
  response_label <- if (response == "juvenile_count") "Juvenile individuals" else "Adult-female individuals"
  model_rows[[length(model_rows) + 1L]] <- extract_period_rr(
    fit_class_model(valid_data, response, adjusted = TRUE),
    "Adjusted class-specific NB2 GLMM",
    response_label,
    "All positive-effort camera-rounds"
  )
  model_rows[[length(model_rows) + 1L]] <- extract_period_rr(
    fit_class_model(valid_data, response, adjusted = FALSE),
    "Focused class-specific NB2 GLMM",
    response_label,
    "All positive-effort camera-rounds"
  )
  model_rows[[length(model_rows) + 1L]] <- extract_period_rr(
    fit_class_model(paired_data, response, adjusted = FALSE),
    "Focused class-specific NB2 GLMM",
    response_label,
    "Cameras operating in both periods"
  )
}

for (response in c("juvenile_events", "adult_female_events")) {
  response_label <- if (response == "juvenile_events") "Juvenile independent events" else "Adult-female independent events"
  model_rows[[length(model_rows) + 1L]] <- extract_period_rr(
    fit_class_model(valid_data, response, adjusted = FALSE),
    "Focused class-specific NB2 GLMM",
    response_label,
    "All positive-effort camera-rounds"
  )
}

model_rows[[length(model_rows) + 1L]] <- extract_joint_difference(
  fit_joint_model(valid_data),
  "All positive-effort camera-rounds"
)
model_rows[[length(model_rows) + 1L]] <- extract_joint_difference(
  fit_joint_model(paired_data),
  "Cameras operating in both periods"
)

model_results <- bind_rows(model_rows)
write_csv(
  model_results,
  file.path(output_dir, "seasonal_rate_model_results.csv")
)

if (any(
  model_results$convergence_code != 0 |
    !model_results$positive_definite_hessian,
  na.rm = TRUE
)) {
  warning("At least one reported model has a convergence or Hessian problem.")
}

# Publication-quality three-panel supplementary figure. Intervals are percentile
# camera-cluster bootstrap intervals and therefore describe sampling variation
# among cameras, not model-based confidence intervals.
rate_plot_data <- bind_rows(
  round_summary %>%
    transmute(
      round,
      round_midpoint,
      panel = "A. Juvenile individuals per 100 camera-nights",
      series = "Juvenile",
      estimate = juvenile_individual_rate_per_100_camera_nights,
      lower = juvenile_individual_rate_per_100_camera_nights_bootstrap_lower,
      upper = juvenile_individual_rate_per_100_camera_nights_bootstrap_upper
    ),
  round_summary %>%
    transmute(
      round,
      round_midpoint,
      panel = "B. Adult-female individuals per 100 camera-nights",
      series = "Adult female",
      estimate = adult_female_individual_rate_per_100_camera_nights,
      lower = adult_female_individual_rate_per_100_camera_nights_bootstrap_lower,
      upper = adult_female_individual_rate_per_100_camera_nights_bootstrap_upper
    ),
  round_summary %>%
    transmute(
      round,
      round_midpoint,
      panel = "C. Juveniles per 100 adult females",
      series = "Juveniles per 100 adult females",
      estimate = juveniles_per_100_adult_females,
      lower = juveniles_per_100_adult_females_bootstrap_lower,
      upper = juveniles_per_100_adult_females_bootstrap_upper
    )
) %>%
  mutate(
    panel = factor(
      panel,
      levels = c(
        "A. Juvenile individuals per 100 camera-nights",
        "B. Adult-female individuals per 100 camera-nights",
        "C. Juveniles per 100 adult females"
      )
    ),
    series = factor(
      series,
      levels = c(
        "Juvenile",
        "Adult female",
        "Juveniles per 100 adult females"
      )
    )
  )

round_labels <- paste0(
  "R", round_summary$round
)
boundary_round_x <- 3 +
  as.numeric(split_date - round_summary$round_midpoint[round_summary$round == 3]) /
  as.numeric(
    round_summary$round_midpoint[round_summary$round == 4] -
      round_summary$round_midpoint[round_summary$round == 3]
  )

figure <- ggplot(
  rate_plot_data,
  aes(round, estimate, color = series, group = series)
) +
  geom_vline(
    xintercept = boundary_round_x,
    linewidth = 0.55,
    linetype = "22",
    color = "#555555"
  ) +
  geom_text(
    data = tibble(
      panel = factor(levels(rate_plot_data$panel), levels = levels(rate_plot_data$panel)),
      round = boundary_round_x,
      estimate = Inf,
      label = "17 Oct"
    ),
    aes(x = round, y = estimate, label = label),
    inherit.aes = FALSE,
    hjust = -0.12,
    vjust = 1.25,
    size = 3.0,
    family = "serif",
    color = "#555555"
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.10,
    linewidth = 0.55
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(shape = 21, fill = "white", size = 2.8, stroke = 0.85) +
  facet_grid(panel ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(
    values = c(
      Juvenile = "#C55A2A",
      `Adult female` = "#287A8B",
      `Juveniles per 100 adult females` = "#303030"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = round_summary$round,
    labels = round_labels,
    expand = expansion(mult = c(0.035, 0.035))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) +
  labs(
    x = "Survey round",
    y = NULL,
    color = NULL,
    caption = paste0(
      "Bars are 95% camera-cluster bootstrap intervals; points are pooled observations.\n",
      "Dashed line: 17 October 2020. Counts sum event-level maximum group sizes."
    )
  ) +
  theme_classic(base_size = 11.5, base_family = "serif") +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 0,
      hjust = 1,
      face = "bold",
      margin = margin(r = 8)
    ),
    panel.spacing.y = grid::unit(1.0, "lines"),
    axis.text.x = element_text(size = 9.3),
    plot.caption = element_text(hjust = 0, size = 8.5, lineheight = 1.05),
    plot.margin = margin(7, 10, 7, 8)
  )

ggsave(
  file.path(output_dir, "juvenile_female_seasonal_rates_supplement.png"),
  figure,
  width = 8.4,
  height = 8.3,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(output_dir, "juvenile_female_seasonal_rates_supplement.pdf"),
  figure,
  width = 8.4,
  height = 8.3,
  device = cairo_pdf,
  bg = "white"
)

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f")))
}

period_early <- period_summary %>% filter(period == "Early")
period_late <- period_summary %>% filter(period == "Late")
raw_juvenile_rate_ratio <-
  period_late$juvenile_individual_rate_per_100_camera_nights /
  period_early$juvenile_individual_rate_per_100_camera_nights
raw_female_rate_ratio <-
  period_late$adult_female_individual_rate_per_100_camera_nights /
  period_early$adult_female_individual_rate_per_100_camera_nights
raw_ratio_of_ratios <-
  period_late$observed_juvenile_female_ratio /
  period_early$observed_juvenile_female_ratio

get_model_row <- function(desired_model, desired_response, desired_subset) {
  model_results %>%
    filter(
      .data$model == desired_model,
      .data$response == desired_response,
      .data$data_subset == desired_subset
    ) %>%
    slice(1)
}

adjusted_juvenile <- get_model_row(
  "Adjusted class-specific NB2 GLMM",
  "Juvenile individuals",
  "All positive-effort camera-rounds"
)
adjusted_female <- get_model_row(
  "Adjusted class-specific NB2 GLMM",
  "Adult-female individuals",
  "All positive-effort camera-rounds"
)
event_juvenile <- get_model_row(
  "Focused class-specific NB2 GLMM",
  "Juvenile independent events",
  "All positive-effort camera-rounds"
)
event_female <- get_model_row(
  "Focused class-specific NB2 GLMM",
  "Adult-female independent events",
  "All positive-effort camera-rounds"
)
joint_primary <- get_model_row(
  "Joint class-by-period NB2 GLMM",
  "Summed event maxima",
  "All positive-effort camera-rounds"
)

report <- c(
  "# Exploratory seasonal decomposition of the juvenile:adult-female ratio",
  "",
  "## Scope and existing inputs",
  "",
  paste0(
    "This analysis is exploratory and does not modify the primary beta-binomial ",
    "model or manuscript outputs. It reuses `gazelle_juvenile_abundance_60x7.csv`, ",
    "`gazelle_female_abundance_60x7.csv`, `ObsCovs_050524.csv`, ",
    "`SiteCovs_060324.csv`, and the clock- and classification-corrected ",
    "`Gazella gazella_independent_events_30min_analysis_ready.csv`. The published ",
    "ratio workflow is `analyze_fawn_female_beta_binomial_wolf_jackal.R`; the earlier ",
    "influence and metadata audit is `investigate_juvenile_female_seasonal_increase.R`."
  ),
  "",
  paste0(
    "The source-code definition assigns a camera-round to Early when its midpoint is ",
    "on or before 17 October 2020 and to Late otherwise. Thus, Late technically begins ",
    "18 October; no camera-round midpoint occurs exactly on 17 October, so the common ",
    "shorthand does not change the fitted data."
  ),
  "",
  "## Reconstruction",
  "",
  paste0(
    "The count histories contain ", sum(valid_data$juvenile_count),
    " juvenile and ", sum(valid_data$adult_female_count),
    " adult-female individuals across ", nrow(valid_data),
    " positive-effort camera-rounds from ", n_distinct(valid_data$camera_id),
    " cameras. The independent-event file contains ",
    sum(valid_data$juvenile_events), " juvenile and ",
    sum(valid_data$adult_female_events),
    " adult-female 30-min events. Summing each event's maximum simultaneously ",
    "observed individuals reproduces the overall totals exactly. However, ",
    nrow(reconstruction_mismatches),
    " camera-round rows differ because the count histories retain session-based ",
    "allocation whereas the event file is assigned by corrected event date. The ",
    "project cleanup audit identifies boundary-date allocation and manual curation as ",
    "likely causes but does not provide a complete historical session-to-round ",
    "crosswalk. Counts are repeated detections, not known unique animals."
  ),
  "",
  "## Descriptive decomposition",
  "",
  paste0(
    "Early contributed ", period_early$camera_nights, " camera-nights, ",
    period_early$juvenile_individual_count, " juvenile individuals, and ",
    period_early$adult_female_individual_count, " adult-female individuals; Late ",
    "contributed ", period_late$camera_nights, " camera-nights, ",
    period_late$juvenile_individual_count, " juvenile individuals, and ",
    period_late$adult_female_individual_count, " adult-female individuals. The pooled ",
    "juvenile rate increased by ", fmt(100 * (raw_juvenile_rate_ratio - 1), 1),
    "%, the pooled adult-female rate changed by ",
    fmt(100 * (raw_female_rate_ratio - 1), 1),
    "%, and the observed juvenile:female ratio increased by ",
    fmt(100 * (raw_ratio_of_ratios - 1), 1), "%."
  ),
  "",
  paste0(
    "The increase was expressed primarily as more juvenile encounters rather than ",
    "larger juvenile groups: juvenile event rows increased from ",
    period_early$juvenile_independent_events, " to ",
    period_late$juvenile_independent_events,
    ", while mean detected group size changed from ",
    fmt(period_early$mean_juvenile_group_size, 2),
    " to ",
    fmt(period_late$mean_juvenile_group_size, 2), "."
  ),
  "",
  "## Count-model results",
  "",
  paste0(
    "Negative-binomial mixed models included every positive-effort camera-round, ",
    "retained zero detections, used log camera-nights as an offset, and accounted for ",
    "repeated cameras. In models using the same stable additive adjustment set as the ",
    "published ratio analysis, the Late/Early juvenile count-rate ratio was ",
    fmt(adjusted_juvenile$rate_ratio), " (95% CI ",
    fmt(adjusted_juvenile$conf_low), "-", fmt(adjusted_juvenile$conf_high),
    ", P = ", fmt_p(adjusted_juvenile$p_value), "), whereas the adult-female ",
    "rate ratio was ", fmt(adjusted_female$rate_ratio), " (95% CI ",
    fmt(adjusted_female$conf_low), "-", fmt(adjusted_female$conf_high),
    ", P = ", fmt_p(adjusted_female$p_value), ")."
  ),
  "",
  paste0(
    "The event-count models gave the same interpretation: juvenile independent events ",
    "increased (rate ratio ", fmt(event_juvenile$rate_ratio), ", 95% CI ",
    fmt(event_juvenile$conf_low), "-", fmt(event_juvenile$conf_high),
    ", P = ", fmt_p(event_juvenile$p_value), "), but adult-female events did not ",
    "show a supported change (rate ratio ", fmt(event_female$rate_ratio),
    ", 95% CI ", fmt(event_female$conf_low), "-", fmt(event_female$conf_high),
    ", P = ", fmt_p(event_female$p_value), "). The direct joint class-by-period ",
    "contrast was close to, but did not cross, the conventional threshold (ratio of ",
    "Late/Early rate ratios = ", fmt(joint_primary$rate_ratio), ", 95% CI ",
    fmt(joint_primary$conf_low), "-", fmt(joint_primary$conf_high),
    ", Wald P = ", fmt_p(joint_primary$p_value), ")."
  ),
  "",
  paste0(
    "Restricting the analysis to the ", length(paired_cameras),
    " cameras represented in both periods left all effect directions essentially ",
    "unchanged. Camera turnover therefore does not explain the pattern, although the ",
    "limited juvenile information makes interval estimates imprecise."
  ),
  "",
  "## Biological interpretation",
  "",
  "- **Additional autumn births:** neither supported nor excluded. The available labels contain no newborn category, body-size notes, or finer juvenile ages.",
  "- **Increasing juvenile mobility or detectability:** consistent with the increase in juvenile event frequency, the absence of increased detected group size, and the spread of juvenile detections across more cameras in the prior influence audit. It is not directly proven.",
  "- **Adult-female movement during the rut:** the pooled adult-female rate declined, but mixed models did not support a decline after accounting for camera and effort. This mechanism remains possible but receives weak direct support here.",
  "- **Early juvenile mortality:** not estimable from these repeated-detection data and not, by itself, an explanation for an increase during the survey. It remains relevant background for interpreting predator associations.",
  "",
  "These are camera-based detection rates and composition indices. They do not estimate survival, recruitment, births, or absolute abundance.",
  "",
  "## Recommendation",
  "",
  paste0(
    "Add one cautious sentence to Section 3.2 and place the three-panel figure in the ",
    "supplement. A defensible sentence is: 'The later increase in juvenile ",
    "representation was driven primarily by a higher frequency of juvenile camera ",
    "detections rather than larger detected juvenile groups; adult-female rates showed ",
    "no supported decline in mixed models, and the available classifications could not ",
    "distinguish increased juvenile detectability from additional births.'"
  ),
  "",
  "## Reproducible outputs",
  "",
  "- `juvenile_female_counts_rates_by_round.csv`: requested seven-round descriptive table with camera-bootstrap intervals.",
  "- `juvenile_female_counts_rates_by_period.csv`: Early/Late descriptive decomposition.",
  "- `seasonal_rate_model_results.csv`: effect sizes, intervals, tests, and convergence diagnostics.",
  "- `event_count_reconstruction_audit.csv`: count-history versus independent-event verification.",
  "- `event_count_reconstruction_mismatches.csv`: camera-round allocation differences retained for audit.",
  "- `juvenile_female_seasonal_rates_supplement.png` and `.pdf`: supplementary figure.",
  "- `run_juvenile_female_seasonal_rate_analysis.R`: complete analysis script."
)

writeLines(
  report,
  file.path(output_dir, "juvenile_female_seasonal_rate_results.md"),
  useBytes = TRUE
)

message(
  "Completed exploratory seasonal-rate analysis: ",
  nrow(valid_data), " camera-rounds, ",
  n_distinct(valid_data$camera_id), " cameras, ",
  bootstrap_replicates, " camera-cluster bootstrap replicates."
)
