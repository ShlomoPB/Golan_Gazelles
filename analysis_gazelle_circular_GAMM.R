#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(mvtnorm)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_map <- list(
  event_file = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  site_cov_file = "SiteCovs_060324.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  output_dir = ".",
  local_timezone = "Asia/Jerusalem",
  primary_k = "8",
  juvenile_fallback_k = "6",
  juvenile_diagnostic_k = "10",
  prediction_step = "0.1",
  coefficient_simulations = "1500",
  site_shape_sensitivity = "true",
  seed = "20260820"
)
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    key <- sub("^--", "", parts[1])
    if (!key %in% names(arg_map)) stop("Unknown argument: ", key)
    arg_map[[key]] <- parts[2]
  }
}

event_file <- arg_map$event_file
site_cov_file <- arg_map$site_cov_file
obs_cov_file <- arg_map$obs_cov_file
output_dir <- arg_map$output_dir
local_timezone <- arg_map$local_timezone
primary_k <- as.integer(arg_map$primary_k)
juvenile_fallback_k <- as.integer(arg_map$juvenile_fallback_k)
juvenile_diagnostic_k <- as.integer(arg_map$juvenile_diagnostic_k)
prediction_step <- as.numeric(arg_map$prediction_step)
coefficient_simulations <- as.integer(arg_map$coefficient_simulations)
run_site_shape <- tolower(arg_map$site_shape_sensitivity) %in% c("true", "1", "yes")
seed <- as.integer(arg_map$seed)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
out_file <- function(name) file.path(output_dir, name)
set.seed(seed)

required_files <- c(event_file, site_cov_file, obs_cov_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

class_keys <- c("juvenile", "female", "male")
class_labels <- c(
  juvenile = "Juvenile",
  female = "Adult female",
  male = "Adult male"
)
expected_events <- c(juvenile = 180L, female = 1477L, male = 2268L)
response_columns <- c(
  juvenile = "juvenile_events",
  female = "female_events",
  male = "male_events"
)
context_levels <- c(
  "Low wolf / Low jackal",
  "Low wolf / High jackal",
  "High wolf / Low jackal",
  "High wolf / High jackal"
)
context_colors <- c(
  "Low wolf / Low jackal" = "#188977",
  "Low wolf / High jackal" = "#B12B32",
  "High wolf / Low jackal" = "#173B63",
  "High wolf / High jackal" = "#8A5A9E"
)
context_linetypes <- c(
  "Low wolf / Low jackal" = "solid",
  "Low wolf / High jackal" = "22",
  "High wolf / Low jackal" = "42",
  "High wolf / High jackal" = "dotdash"
)

fmt <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "NA",
    formatC(x, digits = digits, format = "f", big.mark = ",")
  )
}

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f"))
  )
}

equation_of_time <- function(day_of_year) {
  b <- 2 * pi * (day_of_year - 81) / 364
  9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
}

parse_local_camera_time <- function(x) {
  as.POSIXct(
    sub("Z$", "", x),
    format = "%Y-%m-%dT%H:%M:%S",
    tz = local_timezone
  )
}

apparent_solar_hour <- function(timestamp, longitude) {
  day_of_year <- as.integer(format(timestamp, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(timestamp, "%H", tz = "UTC")) +
    as.numeric(format(timestamp, "%M", tz = "UTC")) / 60 +
    as.numeric(format(timestamp, "%S", tz = "UTC")) / 3600
  (utc_hour + longitude / 15 + equation_of_time(day_of_year) / 60) %% 24
}

read_effort <- function(path) {
  read.csv(
    path,
    na.strings = c("NA", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    transmute(
      camera_id = as.integer(`Camera Number`),
      round = as.integer(`Round Number`),
      start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
      end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
      effort_days = as.numeric(`Days.Operational`)
    ) %>%
    filter(
      !is.na(camera_id),
      !is.na(round),
      !is.na(start_date),
      !is.na(end_date),
      is.finite(effort_days),
      effort_days > 0
    ) %>%
    arrange(camera_id, round)
}

effort <- read_effort(obs_cov_file)
if (
  nrow(effort) != 344L ||
    n_distinct(effort$camera_id, effort$round) != 344L
) {
  stop(
    "Expected 344 unique positive-effort camera-rounds; found ",
    nrow(effort),
    "."
  )
}

survey_start <- min(effort$start_date)
survey_end <- max(effort$end_date)
survey_split <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)

site_data <- read_csv(site_cov_file, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    polygon = as.character(Polygon.Name),
    jackal_index = as.numeric(null.jackal.abundance),
    wolf_index = as.numeric(Null.Wolf.Abundance)
  ) %>%
  arrange(camera_id)

if (
  nrow(site_data) != 60L ||
    n_distinct(site_data$camera_id) != 60L ||
    any(!complete.cases(site_data))
) {
  stop("Expected complete site and predator covariates for all 60 cameras.")
}

wolf_center <- mean(site_data$wolf_index)
wolf_scale <- sd(site_data$wolf_index)
jackal_center <- mean(site_data$jackal_index)
jackal_scale <- sd(site_data$jackal_index)
site_data <- site_data %>%
  mutate(
    wolf_z = (wolf_index - wolf_center) / wolf_scale,
    jackal_z = (jackal_index - jackal_center) / jackal_scale,
    wolf_jackal_z = wolf_z * jackal_z
  )

raw_events <- read_csv(
  event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  mutate(
    event_row_id = row_number(),
    camera_id = as.integer(camera_name),
    gazelle_class = as.character(gazelle_class),
    timestamp = parse_local_camera_time(event_time),
    event_date = as.Date(format(timestamp, tz = local_timezone)),
    solar_time = apparent_solar_hour(timestamp, camelot_gps_longitude),
    solar_bin = floor(solar_time),
    solar_hour = solar_bin + 0.5
  )

if (any(is.na(raw_events$timestamp)) || any(!is.finite(raw_events$solar_time))) {
  stop("The frozen event file contains unparseable timestamps or solar times.")
}

event_round_matches <- raw_events %>%
  select(
    event_row_id,
    camera_id,
    gazelle_class,
    timestamp,
    event_date,
    solar_time,
    solar_bin,
    solar_hour
  ) %>%
  inner_join(effort, by = "camera_id", relationship = "many-to-many") %>%
  filter(event_date >= start_date, event_date <= end_date) %>%
  arrange(event_row_id, desc(start_date), desc(round)) %>%
  group_by(event_row_id) %>%
  slice_head(n = 1) %>%
  ungroup()

unmatched_events <- anti_join(
  raw_events %>% select(event_row_id),
  event_round_matches %>% select(event_row_id),
  by = "event_row_id"
)
if (nrow(unmatched_events) > 0) {
  stop(
    nrow(unmatched_events),
    " frozen events could not be assigned to a positive-effort camera-round."
  )
}

event_totals <- event_round_matches %>%
  count(gazelle_class, name = "n_events")
observed_events <- setNames(event_totals$n_events, event_totals$gazelle_class)
if (
  nrow(event_round_matches) != sum(expected_events) ||
    !all(observed_events[names(expected_events)] == expected_events)
) {
  stop(
    "Validated class-event counts did not equal juvenile = 180, ",
    "female = 1477, and male = 2268."
  )
}

event_counts <- event_round_matches %>%
  count(camera_id, round, solar_bin, gazelle_class, name = "n_events") %>%
  pivot_wider(
    names_from = gazelle_class,
    values_from = n_events,
    values_fill = 0,
    names_glue = "{gazelle_class}_events"
  )

hourly_data <- effort %>%
  select(camera_id, round, start_date, end_date, effort_days) %>%
  crossing(solar_bin = 0:23) %>%
  mutate(solar_hour = solar_bin + 0.5) %>%
  left_join(site_data, by = "camera_id") %>%
  left_join(event_counts, by = c("camera_id", "round", "solar_bin")) %>%
  mutate(
    across(all_of(unname(response_columns)), ~ replace_na(.x, 0L)),
    camera_round = paste(camera_id, round, sep = "_"),
    midpoint_date = start_date + floor(as.numeric(end_date - start_date) / 2),
    period = ifelse(midpoint_date <= survey_split, "Early", "Late"),
    late_indicator = ifelse(period == "Late", 1, 0),
    camera_factor = factor(camera_id),
    camera_round_factor = factor(camera_round),
    polygon_factor = factor(polygon),
    round_factor = factor(round)
  ) %>%
  arrange(camera_id, round, solar_bin)

if (nrow(hourly_data) != 344L * 24L) {
  stop("Expected 8,256 camera-round-hour rows; found ", nrow(hourly_data), ".")
}
hourly_totals <- vapply(
  response_columns,
  function(column) sum(hourly_data[[column]]),
  numeric(1)
)
if (!all(hourly_totals == expected_events[names(hourly_totals)])) {
  stop("Hourly aggregation did not preserve the authoritative event totals.")
}

write_csv(
  hourly_data %>%
    select(
      camera_id,
      polygon,
      round,
      start_date,
      end_date,
      midpoint_date,
      period,
      solar_bin,
      solar_hour,
      effort_days,
      juvenile_events,
      female_events,
      male_events,
      wolf_index,
      jackal_index,
      wolf_z,
      jackal_z,
      wolf_jackal_z,
      camera_round
    ),
  out_file("circular_GAMM_hourly_data.csv")
)

base_formula <- function(k) {
  as.formula(paste0(
    "n_events ~ wolf_z * jackal_z + round_factor + ",
    "s(solar_hour, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = jackal_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_jackal_z, bs = 'cc', k = ", k, ") + ",
    "s(camera_round_factor, bs = 're') + ",
    "s(camera_factor, bs = 're') + ",
    "s(polygon_factor, bs = 're') + ",
    "offset(log(effort_days))"
  ))
}

site_shape_formula <- function(k) {
  as.formula(paste0(
    "n_events ~ wolf_z * jackal_z + round_factor + ",
    "s(solar_hour, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = jackal_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_jackal_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = camera_factor, bs = 'cc', k = 4, m = 1, id = 20) + ",
    "s(camera_round_factor, bs = 're') + ",
    "s(camera_factor, bs = 're') + ",
    "s(polygon_factor, bs = 're') + ",
    "offset(log(effort_days))"
  ))
}

season_shape_formula <- function(k) {
  as.formula(paste0(
    "n_events ~ wolf_z * jackal_z + round_factor + ",
    "s(solar_hour, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = late_indicator, bs = 'cc', k = 6) + ",
    "s(solar_hour, by = wolf_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = jackal_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_jackal_z, bs = 'cc', k = ", k, ") + ",
    "s(camera_round_factor, bs = 're') + ",
    "s(camera_factor, bs = 're') + ",
    "s(polygon_factor, bs = 're') + ",
    "offset(log(effort_days))"
  ))
}

season_site_shape_formula <- function(k) {
  as.formula(paste0(
    "n_events ~ wolf_z * jackal_z + round_factor + ",
    "s(solar_hour, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = late_indicator, bs = 'cc', k = 6) + ",
    "s(solar_hour, by = wolf_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = jackal_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = wolf_jackal_z, bs = 'cc', k = ", k, ") + ",
    "s(solar_hour, by = camera_factor, bs = 'cc', k = 4, m = 1, id = 20) + ",
    "s(camera_round_factor, bs = 're') + ",
    "s(camera_factor, bs = 're') + ",
    "s(polygon_factor, bs = 're') + ",
    "offset(log(effort_days))"
  ))
}

prepare_model_data <- function(class_key) {
  hourly_data %>%
    transmute(
      n_events = .data[[response_columns[[class_key]]]],
      solar_hour,
      effort_days,
      wolf_z,
      jackal_z,
      wolf_jackal_z,
      period,
      late_indicator,
      round_factor,
      camera_round_factor,
      camera_factor,
      polygon_factor,
      camera_id,
      round,
      solar_bin
    )
}

fit_circular_model <- function(model_data, k, formula_builder = base_formula) {
  bam(
    formula_builder(k),
    data = model_data,
    family = nb(link = "log"),
    method = "fREML",
    discrete = TRUE,
    select = TRUE,
    knots = list(solar_hour = c(0, 24)),
    control = gam.control(trace = FALSE)
  )
}

extract_convergence <- function(fit) {
  if (!is.null(fit$outer.info$conv)) return(as.character(fit$outer.info$conv))
  if (!is.null(fit$mgcv.conv)) {
    if (is.list(fit$mgcv.conv) && !is.null(fit$mgcv.conv$fully.converged)) {
      return(ifelse(
        fit$mgcv.conv$fully.converged,
        "full convergence",
        "not converged"
      ))
    }
    if (is.logical(fit$mgcv.conv) && length(fit$mgcv.conv) == 1L) {
      return(ifelse(fit$mgcv.conv, "full convergence", "not converged"))
    }
  }
  "No convergence warning returned"
}

extract_summary_table <- function(fit, class_key) {
  fit_summary <- summary(fit)
  parametric_raw <- as.data.frame(fit_summary$p.table)
  parametric <- tibble(
    class = class_labels[[class_key]],
    term = rownames(fit_summary$p.table),
    term_type = "parametric",
    estimate_if_parametric = parametric_raw[[1]],
    SE_if_parametric = parametric_raw[[2]],
    EDF_if_smooth = NA_real_,
    reference_df = NA_real_,
    test_statistic = parametric_raw[[3]],
    P = parametric_raw[[4]],
    supported = parametric_raw[[4]] < 0.05
  )

  smooth_raw <- as.data.frame(fit_summary$s.table)
  smooth_terms <- rownames(fit_summary$s.table)
  smooth <- tibble(
    class = class_labels[[class_key]],
    term = smooth_terms,
    term_type = case_when(
      grepl("camera|polygon", smooth_terms) ~ "random smooth",
      smooth_terms == "s(solar_hour)" ~ "baseline cyclic smooth",
      TRUE ~ "predator-by-time cyclic smooth"
    ),
    estimate_if_parametric = NA_real_,
    SE_if_parametric = NA_real_,
    EDF_if_smooth = smooth_raw[[1]],
    reference_df = smooth_raw[[2]],
    test_statistic = smooth_raw[[3]],
    P = smooth_raw[[4]],
    supported = smooth_raw[[4]] < 0.05
  )
  bind_rows(parametric, smooth)
}

extract_k_table <- function(fit, class_key) {
  k_result <- k.check(fit)
  k_raw <- as.data.frame(k_result)
  tibble(
    class = class_labels[[class_key]],
    term = rownames(k_result),
    k_prime = k_raw[[1]],
    edf = k_raw[[2]],
    k_index = k_raw[[3]],
    P = k_raw[[4]]
  )
}

extract_concurvity_table <- function(fit, class_key) {
  concurvity_matrix <- concurvity(fit, full = TRUE)
  as.data.frame(as.table(concurvity_matrix), responseName = "value") %>%
    transmute(
      class = class_labels[[class_key]],
      measure = as.character(Var1),
      term = as.character(Var2),
      value
    )
}

diagnostic_summary <- function(fit, model_data, class_key, k_used) {
  fit_summary <- summary(fit)
  fitted_counts <- fitted(fit)
  theta <- fit$family$getTheta(trans = TRUE)
  pearson_residuals <- residuals(fit, type = "pearson")
  predicted_zero <- (theta / (theta + fitted_counts))^theta

  lag_data <- model_data %>%
    mutate(pearson_residual = pearson_residuals) %>%
    arrange(camera_round_factor, solar_bin) %>%
    group_by(camera_round_factor) %>%
    summarise(
      lag1 = suppressWarnings(cor(
        pearson_residual[-n()],
        pearson_residual[-1],
        use = "complete.obs"
      )),
      .groups = "drop"
    )

  tibble(
    class = class_labels[[class_key]],
    rows = nrow(model_data),
    events = sum(model_data$n_events),
    positive_cells = sum(model_data$n_events > 0),
    k_used,
    theta,
    deviance_explained = fit_summary$dev.expl,
    adjusted_r_squared = fit_summary$r.sq,
    convergence = extract_convergence(fit),
    pearson_dispersion = sum(pearson_residuals^2, na.rm = TRUE) /
      fit$df.residual,
    observed_zero_fraction = mean(model_data$n_events == 0),
    mean_predicted_zero_probability = mean(predicted_zero),
    mean_within_camera_round_lag1_residual_correlation = mean(
      lag_data$lag1,
      na.rm = TRUE
    ),
    aic = AIC(fit)
  )
}

primary_models <- list()
model_data_list <- list()
model_summary_list <- list()
k_table_list <- list()
concurvity_list <- list()
diagnostic_list <- list()
juvenile_basis_audit <- tibble()
k_used_by_class <- setNames(rep(primary_k, length(class_keys)), class_keys)

for (class_key in class_keys) {
  message("Fitting primary circular GAMM for ", class_labels[[class_key]], "...")
  model_data <- prepare_model_data(class_key)
  fit <- fit_circular_model(model_data, primary_k)
  k_table <- extract_k_table(fit, class_key)

  ecological_k <- k_table %>%
    filter(grepl("solar_hour", term), is.finite(P))
  needs_juvenile_fallback <- class_key == "juvenile" &&
    extract_convergence(fit) == "not converged"

  if (needs_juvenile_fallback) {
    message("Refitting juvenile model with k = ", juvenile_fallback_k, "...")
    fit_fallback <- fit_circular_model(model_data, juvenile_fallback_k)
    fallback_k_table <- extract_k_table(fit_fallback, class_key)
    if (
      is.finite(AIC(fit_fallback)) &&
        (
          extract_convergence(fit_fallback) == "full convergence" ||
            AIC(fit_fallback) <= AIC(fit) + 2
        )
    ) {
      fit <- fit_fallback
      k_table <- fallback_k_table
      k_used <- juvenile_fallback_k
    } else {
      k_used <- primary_k
    }
  } else if (
    class_key == "juvenile" &&
      any(ecological_k$P < 0.025, na.rm = TRUE)
  ) {
    message(
      "Auditing juvenile basis adequacy with k = ",
      juvenile_diagnostic_k,
      "..."
    )
    fit_diagnostic <- fit_circular_model(model_data, juvenile_diagnostic_k)
    diagnostic_k_table <- extract_k_table(fit_diagnostic, class_key)
    ecological_k_diagnostic <- diagnostic_k_table %>%
      filter(grepl("solar_hour", term), is.finite(P))

    juvenile_basis_audit <- bind_rows(
      juvenile_basis_audit,
      tibble(
        candidate_k = c(primary_k, juvenile_diagnostic_k),
        aic = c(AIC(fit), AIC(fit_diagnostic)),
        minimum_temporal_k_check_p = c(
          min(ecological_k$P, na.rm = TRUE),
          min(ecological_k_diagnostic$P, na.rm = TRUE)
        ),
        deviance_explained = c(
          summary(fit)$dev.expl,
          summary(fit_diagnostic)$dev.expl
        ),
        convergence = c(
          extract_convergence(fit),
          extract_convergence(fit_diagnostic)
        )
      )
    )

    use_diagnostic_basis <-
      extract_convergence(fit_diagnostic) == "full convergence" &&
      min(ecological_k_diagnostic$P, na.rm = TRUE) >= 0.025 &&
      AIC(fit_diagnostic) <= AIC(fit) + 2
    if (use_diagnostic_basis) {
      fit <- fit_diagnostic
      k_table <- diagnostic_k_table
      k_used <- juvenile_diagnostic_k
    } else {
      k_used <- primary_k
    }
  } else {
    k_used <- primary_k
  }

  primary_models[[class_key]] <- fit
  k_used_by_class[[class_key]] <- k_used
  model_data_list[[class_key]] <- model_data
  model_summary_list[[class_key]] <- extract_summary_table(fit, class_key)
  k_table_list[[class_key]] <- k_table
  concurvity_list[[class_key]] <- extract_concurvity_table(fit, class_key)
  diagnostic_list[[class_key]] <- diagnostic_summary(
    fit,
    model_data,
    class_key,
    k_used
  )
}

if (nrow(juvenile_basis_audit) > 0) {
  write_csv(
    juvenile_basis_audit,
    out_file("circular_GAMM_juvenile_basis_audit.csv")
  )
}

model_summary <- bind_rows(model_summary_list)
k_diagnostics <- bind_rows(k_table_list)
concurvity_diagnostics <- bind_rows(concurvity_list)
diagnostics <- bind_rows(diagnostic_list)

write_csv(model_summary, out_file("circular_GAMM_model_summary.csv"))
write_csv(k_diagnostics, out_file("circular_GAMM_k_diagnostics.csv"))
write_csv(concurvity_diagnostics, out_file("circular_GAMM_concurvity.csv"))
write_csv(diagnostics, out_file("circular_GAMM_fit_diagnostics.csv"))
saveRDS(primary_models, out_file("circular_GAMM_models.rds"))

pdf(out_file("circular_GAMM_diagnostic_plots.pdf"), width = 9, height = 8)
gam_check_text <- list()
for (class_key in class_keys) {
  plot.new()
  title(main = paste("mgcv diagnostics:", class_labels[[class_key]]))
  gam_check_text[[class_key]] <- capture.output(
    gam.check(primary_models[[class_key]], rep = 0)
  )
}
dev.off()

site_shape_results <- tibble()
site_shape_models <- list()
random_intercept_models <- primary_models
if (run_site_shape) {
  for (class_key in class_keys) {
    message("Fitting site-shape sensitivity for ", class_labels[[class_key]], "...")
    sensitivity_fit <- tryCatch(
      fit_circular_model(
        model_data_list[[class_key]],
        k_used_by_class[[class_key]],
        formula_builder = site_shape_formula
      ),
      error = function(e) e
    )
    if (inherits(sensitivity_fit, "error")) {
      site_shape_results <- bind_rows(
        site_shape_results,
        tibble(
          class = class_labels[[class_key]],
          fitted = FALSE,
          convergence = conditionMessage(sensitivity_fit),
          primary_aic = AIC(primary_models[[class_key]]),
          site_shape_aic = NA_real_,
          delta_aic_site_shape_minus_primary = NA_real_,
          site_shape_total_edf = NA_real_,
          deviance_explained = NA_real_,
          selected_as_final = FALSE
        )
      )
    } else {
      site_shape_models[[class_key]] <- sensitivity_fit
      smooth_table <- as.data.frame(summary(sensitivity_fit)$s.table) %>%
        mutate(term = rownames(summary(sensitivity_fit)$s.table))
      site_shape_edf <- sum(
        smooth_table[[1]][grepl("camera_factor", smooth_table$term) &
          grepl("solar_hour", smooth_table$term)],
        na.rm = TRUE
      )
      site_shape_results <- bind_rows(
        site_shape_results,
        tibble(
          class = class_labels[[class_key]],
          fitted = TRUE,
          convergence = extract_convergence(sensitivity_fit),
          primary_aic = AIC(primary_models[[class_key]]),
          site_shape_aic = AIC(sensitivity_fit),
          delta_aic_site_shape_minus_primary = AIC(sensitivity_fit) -
            AIC(primary_models[[class_key]]),
          site_shape_total_edf = site_shape_edf,
          deviance_explained = summary(sensitivity_fit)$dev.expl,
          selected_as_final = class_key != "juvenile" &&
            extract_convergence(sensitivity_fit) == "full convergence" &&
            AIC(sensitivity_fit) - AIC(primary_models[[class_key]]) < -10
        )
      )
    }
  }
  saveRDS(site_shape_models, out_file("circular_GAMM_site_shape_sensitivity_models.rds"))
  write_csv(
    site_shape_results,
    out_file("circular_GAMM_site_shape_sensitivity.csv")
  )
}

# Adult camera-specific cyclic deviations were promoted only when they clearly
# improved the prespecified random-intercept model and converged. The juvenile
# site-shape model remains a sensitivity because only 180 events inform 60 sites.
if (nrow(site_shape_results) > 0) {
  for (class_key in class_keys) {
    class_label <- class_labels[[class_key]]
    selected <- site_shape_results %>%
      filter(class == class_label) %>%
      pull(selected_as_final)
    if (length(selected) == 1L && isTRUE(selected) && !is.null(site_shape_models[[class_key]])) {
      primary_models[[class_key]] <- site_shape_models[[class_key]]
    }
  }
}

# Rebuild all authoritative summaries from the selected final models.
model_summary_list <- list()
k_table_list <- list()
concurvity_list <- list()
diagnostic_list <- list()
for (class_key in class_keys) {
  model_summary_list[[class_key]] <- extract_summary_table(
    primary_models[[class_key]],
    class_key
  )
  k_table_list[[class_key]] <- extract_k_table(primary_models[[class_key]], class_key)
  concurvity_list[[class_key]] <- extract_concurvity_table(
    primary_models[[class_key]],
    class_key
  )
  diagnostic_list[[class_key]] <- diagnostic_summary(
    primary_models[[class_key]],
    model_data_list[[class_key]],
    class_key,
    k_used_by_class[[class_key]]
  ) %>%
    mutate(
      model_structure = ifelse(
        class_key != "juvenile" &&
          !identical(primary_models[[class_key]], random_intercept_models[[class_key]]),
        "camera-specific cyclic deviations",
        "random intercepts only"
      )
    )
}
model_summary <- bind_rows(model_summary_list)
k_diagnostics <- bind_rows(k_table_list)
concurvity_diagnostics <- bind_rows(concurvity_list)
diagnostics <- bind_rows(diagnostic_list)

write_csv(model_summary, out_file("circular_GAMM_model_summary.csv"))
write_csv(k_diagnostics, out_file("circular_GAMM_k_diagnostics.csv"))
write_csv(concurvity_diagnostics, out_file("circular_GAMM_concurvity.csv"))
write_csv(diagnostics, out_file("circular_GAMM_fit_diagnostics.csv"))
write_csv(site_shape_results, out_file("circular_GAMM_site_shape_sensitivity.csv"))
saveRDS(primary_models, out_file("circular_GAMM_models.rds"))
saveRDS(
  random_intercept_models,
  out_file("circular_GAMM_random_intercept_candidate_models.rds")
)

pdf(out_file("circular_GAMM_diagnostic_plots.pdf"), width = 9, height = 8)
gam_check_text <- list()
for (class_key in class_keys) {
  plot.new()
  title(main = paste("mgcv diagnostics:", class_labels[[class_key]]))
  gam_check_text[[class_key]] <- capture.output(
    gam.check(primary_models[[class_key]], rep = 0)
  )
}
dev.off()

# Control for the frozen Early/Late survey-period contrast as a difference in
# diel curve shape. This sensitivity does not estimate predator-by-period
# interactions; it asks whether the pooled predator smooths survive adjustment
# for broad seasonal timing changes beyond the round-specific intercepts.
season_shape_models <- list()
season_shape_rows <- list()
for (class_key in class_keys) {
  class_label <- class_labels[[class_key]]
  has_site_shape <- any(grepl(
    "solar_hour.*camera_factor",
    vapply(primary_models[[class_key]]$smooth, function(x) x$label, character(1))
  ))
  formula_builder <- if (has_site_shape) {
    season_site_shape_formula
  } else {
    season_shape_formula
  }
  message("Fitting Early/Late curve-shape sensitivity for ", class_label, "...")
  season_fit <- fit_circular_model(
    model_data_list[[class_key]],
    k_used_by_class[[class_key]],
    formula_builder = formula_builder
  )
  season_shape_models[[class_key]] <- season_fit

  primary_smooth <- as.data.frame(summary(primary_models[[class_key]])$s.table) %>%
    mutate(term = rownames(summary(primary_models[[class_key]])$s.table))
  season_smooth <- as.data.frame(summary(season_fit)$s.table) %>%
    mutate(term = rownames(summary(season_fit)$s.table))
  ecological_terms <- c(
    "s(solar_hour):wolf_z" = "Wolf time-varying effect",
    "s(solar_hour):jackal_z" = "Jackal time-varying effect",
    "s(solar_hour):wolf_jackal_z" = "Wolf x jackal time-varying effect",
    "s(solar_hour):late_indicator" = "Late versus Early curve-shape control"
  )
  for (term_i in names(ecological_terms)) {
    primary_index <- match(term_i, primary_smooth$term)
    season_index <- match(term_i, season_smooth$term)
    season_shape_rows[[length(season_shape_rows) + 1L]] <- tibble(
      class = class_label,
      effect = unname(ecological_terms[[term_i]]),
      term = term_i,
      primary_edf = ifelse(is.na(primary_index), NA_real_, primary_smooth[[1]][primary_index]),
      primary_P = ifelse(is.na(primary_index), NA_real_, primary_smooth[[4]][primary_index]),
      season_adjusted_edf = ifelse(is.na(season_index), NA_real_, season_smooth[[1]][season_index]),
      season_adjusted_P = ifelse(is.na(season_index), NA_real_, season_smooth[[4]][season_index]),
      convergence = extract_convergence(season_fit),
      delta_aic_season_adjusted_minus_primary = AIC(season_fit) -
        AIC(primary_models[[class_key]]),
      season_adjusted_deviance_explained = summary(season_fit)$dev.expl
    )
  }
}
season_shape_sensitivity <- bind_rows(season_shape_rows)
write_csv(
  season_shape_sensitivity,
  out_file("circular_GAMM_season_shape_sensitivity.csv")
)
saveRDS(
  season_shape_models,
  out_file("circular_GAMM_season_shape_sensitivity_models.rds")
)

simulate_predictions <- function(fit, model_data, class_key) {
  solar_grid <- seq(
    prediction_step / 2,
    24 - prediction_step / 2,
    by = prediction_step
  )
  new_data <- expand_grid(
    solar_hour = solar_grid,
    wolf_z = c(-1, 1),
    jackal_z = c(-1, 1)
  ) %>%
    mutate(
      wolf_jackal_z = wolf_z * jackal_z,
      round_factor = factor("4", levels = levels(model_data$round_factor)),
      camera_round_factor = factor(
        levels(model_data$camera_round_factor)[1],
        levels = levels(model_data$camera_round_factor)
      ),
      camera_factor = factor(
        levels(model_data$camera_factor)[1],
        levels = levels(model_data$camera_factor)
      ),
      polygon_factor = factor(
        levels(model_data$polygon_factor)[1],
        levels = levels(model_data$polygon_factor)
      ),
      effort_days = 1,
      context = case_when(
        wolf_z == -1 & jackal_z == -1 ~ "Low wolf / Low jackal",
        wolf_z == -1 & jackal_z == 1 ~ "Low wolf / High jackal",
        wolf_z == 1 & jackal_z == -1 ~ "High wolf / Low jackal",
        TRUE ~ "High wolf / High jackal"
      ),
      context = factor(context, levels = context_levels)
    ) %>%
    arrange(context, solar_hour)

  random_terms <- vapply(fit$smooth, function(x) x$label, character(1))
  random_terms <- random_terms[grepl("camera|polygon", random_terms)]
  linear_predictor_matrix <- predict(
    fit,
    newdata = new_data,
    type = "lpmatrix",
    exclude = random_terms
  )
  covariance <- tryCatch(
    vcov(fit, unconditional = TRUE),
    error = function(e) vcov(fit)
  )
  coefficient_draws <- rmvnorm(
    coefficient_simulations,
    mean = coef(fit),
    sigma = covariance,
    method = "svd"
  )

  point_raw <- exp(drop(linear_predictor_matrix %*% coef(fit)))
  draw_eta <- linear_predictor_matrix %*% t(coefficient_draws)
  draw_raw <- exp(pmin(draw_eta, 700))
  point_normalized <- point_raw
  draw_normalized <- draw_raw

  for (context_i in context_levels) {
    index <- which(new_data$context == context_i)
    point_area <- sum(point_raw[index]) * prediction_step
    draw_area <- colSums(draw_raw[index, , drop = FALSE]) * prediction_step
    point_normalized[index] <- point_raw[index] / point_area
    draw_normalized[index, ] <- sweep(
      draw_raw[index, , drop = FALSE],
      2,
      draw_area,
      "/"
    )
  }

  prediction_table <- new_data %>%
    transmute(
      class = class_labels[[class_key]],
      solar_hour,
      wolf_level = ifelse(wolf_z < 0, "Low (-1 SD)", "High (+1 SD)"),
      jackal_level = ifelse(jackal_z < 0, "Low (-1 SD)", "High (+1 SD)"),
      wolf_z,
      jackal_z,
      context,
      predicted_raw = point_raw,
      raw_CI_lower = apply(draw_raw, 1, quantile, probs = 0.025),
      raw_CI_upper = apply(draw_raw, 1, quantile, probs = 0.975),
      predicted_normalized = point_normalized,
      CI_lower = apply(draw_normalized, 1, quantile, probs = 0.025),
      CI_upper = apply(draw_normalized, 1, quantile, probs = 0.975)
    )

  context_indices <- setNames(
    lapply(context_levels, function(x) which(new_data$context == x)),
    context_levels
  )
  contrast_definitions <- list(
    list(
      label = "High - low wolf | jackal low",
      high = "High wolf / Low jackal",
      low = "Low wolf / Low jackal",
      focal = "Wolf"
    ),
    list(
      label = "High - low wolf | jackal high",
      high = "High wolf / High jackal",
      low = "Low wolf / High jackal",
      focal = "Wolf"
    ),
    list(
      label = "High - low jackal | wolf low",
      high = "Low wolf / High jackal",
      low = "Low wolf / Low jackal",
      focal = "Jackal"
    ),
    list(
      label = "High - low jackal | wolf high",
      high = "High wolf / High jackal",
      low = "High wolf / Low jackal",
      focal = "Jackal"
    )
  )

  difference_tables <- lapply(contrast_definitions, function(definition) {
    high_index <- context_indices[[definition$high]]
    low_index <- context_indices[[definition$low]]
    draw_difference <- draw_normalized[high_index, , drop = FALSE] -
      draw_normalized[low_index, , drop = FALSE]
    point_difference <- point_normalized[high_index] - point_normalized[low_index]
    lower <- apply(draw_difference, 1, quantile, probs = 0.025)
    upper <- apply(draw_difference, 1, quantile, probs = 0.975)
    tibble(
      class = class_labels[[class_key]],
      contrast = definition$label,
      focal_predator = definition$focal,
      solar_hour = solar_grid,
      difference_normalized = point_difference,
      CI_lower = lower,
      CI_upper = upper,
      interval_excludes_zero = lower > 0 | upper < 0
    )
  }) %>%
    bind_rows()

  metric_tables <- lapply(context_levels, function(context_i) {
    index <- context_indices[[context_i]]
    normalized_i <- draw_normalized[index, , drop = FALSE]
    daytime_index <- solar_grid >= 6 & solar_grid < 18
    daytime_draws <- colSums(
      normalized_i[daytime_index, , drop = FALSE]
    ) * prediction_step
    peak_draw_indices <- max.col(t(normalized_i), ties.method = "first")
    peak_draws <- solar_grid[peak_draw_indices]
    point_peak <- solar_grid[which.max(point_normalized[index])]
    tibble(
      class = class_labels[[class_key]],
      context = context_i,
      predicted_06_18_proportion = sum(
        point_normalized[index][daytime_index]
      ) * prediction_step,
      proportion_CI_lower = quantile(daytime_draws, 0.025),
      proportion_CI_upper = quantile(daytime_draws, 0.975),
      principal_peak_hour = point_peak,
      peak_hour_simulation_median = median(peak_draws),
      peak_hour_CI_lower = quantile(peak_draws, 0.025),
      peak_hour_CI_upper = quantile(peak_draws, 0.975)
    )
  }) %>%
    bind_rows()

  list(
    predictions = prediction_table,
    differences = difference_tables,
    metrics = metric_tables
  )
}

simulation_results <- lapply(class_keys, function(class_key) {
  message("Simulating population-level predictions for ", class_labels[[class_key]], "...")
  simulate_predictions(
    primary_models[[class_key]],
    model_data_list[[class_key]],
    class_key
  )
})
names(simulation_results) <- class_keys

predictions <- bind_rows(lapply(simulation_results, `[[`, "predictions"))
differences <- bind_rows(lapply(simulation_results, `[[`, "differences"))
derived_metrics <- bind_rows(lapply(simulation_results, `[[`, "metrics"))

write_csv(predictions, out_file("circular_GAMM_predictions.csv"))
write_csv(differences, out_file("circular_GAMM_difference_curves.csv"))
write_csv(derived_metrics, out_file("circular_GAMM_derived_metrics.csv"))

extract_supported_windows <- function(data) {
  data <- arrange(data, solar_hour)
  supported <- data$interval_excludes_zero
  if (!any(supported)) return(tibble())
  runs <- rle(supported)
  run_end <- cumsum(runs$lengths)
  run_start <- run_end - runs$lengths + 1L
  supported_runs <- which(runs$values)
  index_sets <- lapply(
    supported_runs,
    function(i) seq.int(run_start[i], run_end[i])
  )

  crosses_midnight <- rep(FALSE, length(index_sets))
  if (
    length(index_sets) > 1L &&
      supported[1] && supported[length(supported)]
  ) {
    index_sets[[1]] <- c(index_sets[[length(index_sets)]], index_sets[[1]])
    crosses_midnight[1] <- TRUE
    index_sets <- index_sets[-length(index_sets)]
    crosses_midnight <- crosses_midnight[-length(crosses_midnight)]
  }

  bind_rows(lapply(seq_along(index_sets), function(i) {
    indices <- index_sets[[i]]
    if (crosses_midnight[i]) {
      late_indices <- indices[data$solar_hour[indices] > 12]
      early_indices <- indices[data$solar_hour[indices] <= 12]
      start_hour <- max(0, min(data$solar_hour[late_indices]) - prediction_step / 2)
      end_hour <- min(24, max(data$solar_hour[early_indices]) + prediction_step / 2)
    } else {
      start_hour <- max(0, min(data$solar_hour[indices]) - prediction_step / 2)
      end_hour <- min(24, max(data$solar_hour[indices]) + prediction_step / 2)
    }
    tibble(
      start_hour = start_hour,
      end_hour = end_hour,
      crosses_midnight = crosses_midnight[i],
      duration_hours = length(indices) * prediction_step,
      direction = ifelse(
        mean(data$difference_normalized[indices]) > 0,
        "higher under high focal-predator activity",
        "lower under high focal-predator activity"
      )
    )
  }))
}

supported_windows <- differences %>%
  group_by(class, contrast, focal_predator) %>%
  group_modify(~ extract_supported_windows(.x)) %>%
  ungroup()
write_csv(
  supported_windows,
  out_file("circular_GAMM_supported_time_windows.csv")
)

validation_tables <- lapply(class_keys, function(class_key) {
  model_data <- model_data_list[[class_key]]
  fit <- primary_models[[class_key]]
  model_data %>%
    mutate(fitted_count = fitted(fit)) %>%
    group_by(solar_hour) %>%
    summarise(
      class = class_labels[[class_key]],
      observed_events = sum(n_events),
      fitted_events = sum(fitted_count),
      total_effort_days = sum(effort_days),
      observed_rate = observed_events / total_effort_days,
      fitted_rate = fitted_events / total_effort_days,
      .groups = "drop"
    )
}) %>%
  bind_rows()
write_csv(validation_tables, out_file("circular_GAMM_raw_validation.csv"))

class_factor_levels <- unname(class_labels[class_keys])
predictions <- predictions %>%
  mutate(
    class = factor(class, levels = class_factor_levels),
    context = factor(context, levels = context_levels)
  )
raw_rug <- event_round_matches %>%
  transmute(
    class = factor(unname(class_labels[gazelle_class]), levels = class_factor_levels),
    solar_hour = solar_time
  )

prediction_plot <- ggplot(
  predictions,
  aes(
    x = solar_hour,
    y = predicted_normalized,
    color = context,
    fill = context,
    linetype = context
  )
) +
  geom_ribbon(
    aes(ymin = CI_lower, ymax = CI_upper),
    alpha = 0.07,
    color = NA
  ) +
  geom_line(linewidth = 0.95) +
  geom_rug(
    data = raw_rug,
    inherit.aes = FALSE,
    aes(x = solar_hour),
    sides = "b",
    alpha = 0.055,
    color = "#263238",
    linewidth = 0.18,
    length = grid::unit(0.025, "npc")
  ) +
  facet_wrap(vars(class), ncol = 1, scales = "free_y") +
  scale_color_manual(values = context_colors, drop = FALSE, name = NULL) +
  scale_fill_manual(values = context_colors, drop = FALSE, name = NULL) +
  scale_linetype_manual(values = context_linetypes, drop = FALSE, name = NULL) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 3),
    expand = c(0, 0)
  ) +
  labs(
    x = "Apparent solar time",
    y = "Normalized predicted activity density",
    title = "Gazelle diel activity across continuous wolf-jackal gradients",
    subtitle = paste0(
      "Curves show population-level predictions at -1 and +1 SD; ",
      "rugs show observed independent class-events"
    ),
    caption = paste0(
      "Shading shows simulation-based 95% intervals. Curves are normalized ",
      "within predator combinations and describe timing, not abundance."
    )
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E3E6E6", linewidth = 0.35),
    panel.grid.major.x = element_line(color = "#ECEEEE", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 11.5),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "#4D5A5E"),
    plot.caption = element_text(color = "#5F696B", hjust = 0),
    legend.position = "top",
    legend.key.width = grid::unit(1.25, "cm")
  )

ggsave(
  out_file("circular_GAMM_predicted_curves.png"),
  prediction_plot,
  width = 10.5,
  height = 9.2,
  dpi = 320
)
ggsave(
  out_file("circular_GAMM_predicted_curves.pdf"),
  prediction_plot,
  width = 10.5,
  height = 9.2
)

differences <- differences %>%
  mutate(
    class = factor(class, levels = class_factor_levels),
    contrast = factor(
      contrast,
      levels = c(
        "High - low wolf | jackal low",
        "High - low wolf | jackal high",
        "High - low jackal | wolf low",
        "High - low jackal | wolf high"
      )
    )
  )
difference_plot <- ggplot(
  differences,
  aes(x = solar_hour, y = difference_normalized, color = focal_predator)
) +
  geom_hline(yintercept = 0, color = "#5D676A", linewidth = 0.45) +
  geom_ribbon(
    aes(ymin = CI_lower, ymax = CI_upper, fill = focal_predator),
    alpha = 0.12,
    color = NA
  ) +
  geom_line(linewidth = 0.75) +
  facet_grid(rows = vars(class), cols = vars(contrast), scales = "free_y") +
  scale_color_manual(values = c(Wolf = "#173B63", Jackal = "#B12B32"), name = NULL) +
  scale_fill_manual(values = c(Wolf = "#173B63", Jackal = "#B12B32"), name = NULL) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 6),
    expand = c(0, 0)
  ) +
  labs(
    x = "Apparent solar time",
    y = "Difference in normalized activity density",
    title = "Predator-associated differences in modeled gazelle diel activity",
    subtitle = "Positive values indicate greater relative activity under the high focal-predator condition",
    caption = "Shading shows simulation-based 95% intervals."
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E8E8", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 9),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "top"
  )

ggsave(
  out_file("circular_GAMM_difference_curves.png"),
  difference_plot,
  width = 14,
  height = 8.5,
  dpi = 320
)
ggsave(
  out_file("circular_GAMM_difference_curves.pdf"),
  difference_plot,
  width = 14,
  height = 8.5
)

validation_tables <- validation_tables %>%
  mutate(class = factor(class, levels = class_factor_levels))
validation_plot <- ggplot(validation_tables, aes(x = solar_hour)) +
  geom_point(
    aes(y = observed_rate),
    color = "#2F3E43",
    size = 1.7,
    alpha = 0.8
  ) +
  geom_line(
    aes(y = fitted_rate),
    color = "#B56A2D",
    linewidth = 0.9
  ) +
  facet_wrap(vars(class), ncol = 1, scales = "free_y") +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 3),
    expand = c(0, 0)
  ) +
  labs(
    x = "Apparent solar time",
    y = "Events per operational camera-day",
    title = "Raw hourly rates and fitted circular GAMM values",
    subtitle = "Points are observed aggregate rates; orange lines aggregate fitted values"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15)
  )
ggsave(
  out_file("circular_GAMM_raw_validation.png"),
  validation_plot,
  width = 9,
  height = 8,
  dpi = 320
)
ggsave(
  out_file("circular_GAMM_raw_validation.pdf"),
  validation_plot,
  width = 9,
  height = 8
)

ecological_smooths <- model_summary %>%
  filter(term_type == "predator-by-time cyclic smooth") %>%
  mutate(
    effect = case_when(
      grepl("wolf_jackal", term) ~ "Wolf x jackal time-varying effect",
      grepl("wolf_z", term) ~ "Wolf time-varying effect",
      grepl("jackal_z", term) ~ "Jackal time-varying effect",
      TRUE ~ term
    )
  )

max_ecological_concurvity <- concurvity_diagnostics %>%
  filter(
    measure == "estimate",
    grepl("solar_hour", term),
    !grepl("camera_factor", term)
  ) %>%
  group_by(class) %>%
  summarise(max_concurvity = max(value, na.rm = TRUE), .groups = "drop")

diagnostic_lines <- c(
  "# Circular GAMM diagnostics",
  "",
  "## Data and effort",
  "",
  paste0(
    "The frozen independent-event file contributed exactly ",
    sum(expected_events),
    " class-events: ",
    expected_events[["juvenile"]],
    " juvenile, ",
    format(expected_events[["female"]], big.mark = ","),
    " adult-female, and ",
    format(expected_events[["male"]], big.mark = ","),
    " adult-male events."
  ),
  paste0(
    "The analysis table contains ",
    format(nrow(hourly_data), big.mark = ","),
    " rows (344 positive-effort camera-rounds x 24 apparent-solar-hour bins)."
  ),
  paste0(
    "Each hourly cell used `Days.Operational` as its exposure because exact ",
    "within-day malfunction times were unavailable. Thus, the offset estimates ",
    "events in a given solar-hour bin per operational camera-day."
  ),
  paste0(
    "Wolf and jackal indices were standardized across the 60 cameras ",
    "(wolf mean = ",
    fmt(wolf_center),
    ", SD = ",
    fmt(wolf_scale),
    "; jackal mean = ",
    fmt(jackal_center),
    ", SD = ",
    fmt(jackal_scale),
    ")."
  ),
  "",
  "## Primary specification",
  "",
  paste0(
    "Separate negative-binomial `mgcv::bam` models used a log link, fREML, ",
    "cyclic cubic bases over 0-24 h, shrinkage (`select = TRUE`), factor(round), ",
    "and random intercept smooths for camera-round, camera, and polygon."
  ),
  paste0(
    "Because supported camera-specific timing heterogeneity materially improved ",
    "adult fit, the final female and male models also contained shared-penalty ",
    "camera-specific cyclic deviations (k = 4 per camera). The juvenile model ",
    "retained random intercepts only because 180 events were insufficient to ",
    "support 60 camera-level activity curves."
  ),
  "The ordinary wolf, jackal, and interaction terms describe overall event intensity; their time-varying cyclic smooths describe changes in diel shape.",
  "",
  "## Class-specific diagnostics",
  ""
)

for (class_key in class_keys) {
  class_label <- class_labels[[class_key]]
  diag_i <- diagnostics %>% filter(class == class_label)
  smooth_i <- ecological_smooths %>% filter(class == class_label)
  k_i <- k_diagnostics %>%
    filter(
      class == class_label,
      grepl("solar_hour", term),
      !grepl("camera_factor", term),
      is.finite(k_index)
    )
  conc_i <- max_ecological_concurvity %>% filter(class == class_label)
  diagnostic_lines <- c(
    diagnostic_lines,
    paste0("### ", class_label),
    "",
    paste0("Formula: `", paste(deparse(formula(primary_models[[class_key]])), collapse = " "), "`"),
    "",
    paste0(
      "The model used k = ",
      diag_i$k_used,
      " with `",
      diag_i$model_structure,
      "`",
      ", theta = ",
      fmt(diag_i$theta),
      ", and explained ",
      fmt(100 * diag_i$deviance_explained, 1),
      "% of deviance. Convergence status: `",
      diag_i$convergence,
      "`."
    ),
    paste0(
      "Observed zeros = ",
      fmt(100 * diag_i$observed_zero_fraction, 1),
      "%; mean model-implied zero probability = ",
      fmt(100 * diag_i$mean_predicted_zero_probability, 1),
      "%; Pearson dispersion = ",
      fmt(diag_i$pearson_dispersion),
      "."
    ),
    paste0(
      "Mean within-camera-round lag-1 Pearson residual correlation = ",
      fmt(diag_i$mean_within_camera_round_lag1_residual_correlation),
      ". Maximum estimated ecological-smooth concurvity = ",
      ifelse(nrow(conc_i), fmt(conc_i$max_concurvity), "NA"),
      "."
    ),
    paste0(
      "Temporal-smooth k-index range = ",
      fmt(min(k_i$k_index, na.rm = TRUE)),
      "-",
      fmt(max(k_i$k_index, na.rm = TRUE)),
      "; minimum k-check P = ",
      fmt_p(min(k_i$P, na.rm = TRUE)),
      "."
    ),
    "",
    "| Time-varying term | EDF | Approximate P |",
    "|---|---:|---:|",
    paste0(
      "| ",
      smooth_i$effect,
      " | ",
      fmt(smooth_i$EDF_if_smooth, 2),
      " | ",
      fmt_p(smooth_i$P),
      " |"
    ),
    "",
    "Native `mgcv` residual, k-index, and concurvity diagnostics were used; `gratia` and `DHARMa` were not required or installed.",
    ""
  )
}

if (nrow(site_shape_results) > 0) {
  diagnostic_lines <- c(
    diagnostic_lines,
    "## Camera-specific curve sensitivity",
    "",
    "A shared-penalty camera-specific cyclic deviation (k = 4 per camera) was fitted following the hierarchical timing logic of Iannarilli et al. (2025). It was promoted for adult classes only when it converged and improved AIC by more than 10; the juvenile version remained a sensitivity regardless of fit because only 180 events informed 60 cameras.",
    "",
    "| Class | Convergence | Delta AIC (site shape - random-intercept candidate) | Site-shape EDF | Deviance explained | Selected final |",
    "|---|---|---:|---:|---:|---|",
    paste0(
      "| ",
      site_shape_results$class,
      " | ",
      site_shape_results$convergence,
      " | ",
      fmt(site_shape_results$delta_aic_site_shape_minus_primary, 2),
      " | ",
      fmt(site_shape_results$site_shape_total_edf, 2),
      " | ",
      fmt(100 * site_shape_results$deviance_explained, 1),
      "% | ",
      site_shape_results$selected_as_final,
      " |"
    ),
    "",
    "AIC values from penalized fREML fits are treated descriptively rather than as formal likelihood-ratio tests. Primary inference remains based on the prespecified population-level smooths, shrinkage, diagnostics, and prediction simulations.",
    ""
  )
}
diagnostic_lines <- c(
  diagnostic_lines,
  "## Early/Late curve-shape sensitivity",
  "",
  paste0(
    "The frozen survey midpoint was `",
    survey_split,
    "`. A numeric Late-period cyclic difference smooth (k = 6) controlled broad ",
    "period-specific changes in diel shape while factor(round) continued to ",
    "control round-specific event intensity. This sensitivity did not fit ",
    "predator-by-period interactions."
  ),
  "",
  "| Class | Effect | Primary EDF | Primary P | Season-adjusted EDF | Season-adjusted P | Delta AIC |",
  "|---|---|---:|---:|---:|---:|---:|",
  paste0(
    "| ",
    season_shape_sensitivity$class,
    " | ",
    season_shape_sensitivity$effect,
    " | ",
    fmt(season_shape_sensitivity$primary_edf, 2),
    " | ",
    fmt_p(season_shape_sensitivity$primary_P),
    " | ",
    fmt(season_shape_sensitivity$season_adjusted_edf, 2),
    " | ",
    fmt_p(season_shape_sensitivity$season_adjusted_P),
    " | ",
    fmt(season_shape_sensitivity$delta_aic_season_adjusted_minus_primary, 2),
    " |"
  ),
  "",
  "The three supported predator-timing smooths retained the same inference after this adjustment. AIC differences are descriptive because the compared penalized fREML specifications differ in fixed-effect structure.",
  ""
)
writeLines(diagnostic_lines, out_file("circular_GAMM_diagnostics.md"))

support_counts <- ecological_smooths %>%
  group_by(class) %>%
  summarise(n_supported = sum(P < 0.05), .groups = "drop")
broad_difference_summary <- differences %>%
  group_by(class, contrast) %>%
  summarise(
    supported_hours = sum(interval_excludes_zero) * prediction_step,
    max_absolute_difference = max(abs(difference_normalized)),
    .groups = "drop"
  )

all_converged <- all(grepl("converg", tolower(diagnostics$convergence)))
adequate_k <- all(
  k_diagnostics$P[
    grepl("solar_hour", k_diagnostics$term) &
      !grepl("camera_factor", k_diagnostics$term) &
      is.finite(k_diagnostics$P)
  ] >= 0.01
)
acceptable_concurvity <- all(max_ecological_concurvity$max_concurvity < 0.9)
coherent_effects <- any(broad_difference_summary$supported_hours >= 2)
season_effects_stable <- all(c(
  season_shape_sensitivity %>%
    filter(
      class == "Juvenile",
      effect == "Wolf x jackal time-varying effect"
    ) %>%
    pull(season_adjusted_P) < 0.05,
  season_shape_sensitivity %>%
    filter(
      class == "Adult female",
      effect == "Wolf time-varying effect"
    ) %>%
    pull(season_adjusted_P) < 0.05,
  season_shape_sensitivity %>%
    filter(
      class == "Adult male",
      effect == "Wolf x jackal time-varying effect"
    ) %>%
    pull(season_adjusted_P) < 0.05
))
primary_recommendation <- if (
  all_converged && adequate_k && acceptable_concurvity &&
    coherent_effects && season_effects_stable
) {
  "GAMM"
} else {
  "KERNEL"
}

results_lines <- c(
  "# Circular GAMM results",
  "",
  "## Analysis summary",
  "",
  paste0(
    "We modeled independent class-event counts in ",
    format(nrow(hourly_data), big.mark = ","),
    " camera-round-hour cells with negative-binomial hierarchical cyclic GAMMs. ",
    "Predator covariates were continuous standardized site-level activity indices, not predator density."
  ),
  "Population-level predictions excluded all random intercepts, used round 4 and one operational camera-day, and were normalized within each predator combination for timing comparisons.",
  "",
  "## Time-varying predator effects",
  ""
)

for (class_key in class_keys) {
  class_label <- class_labels[[class_key]]
  smooth_i <- ecological_smooths %>% filter(class == class_label)
  metric_i <- derived_metrics %>% filter(class == class_label)
  diff_i <- broad_difference_summary %>% filter(class == class_label)
  get_effect <- function(pattern) smooth_i %>% filter(grepl(pattern, effect)) %>% slice(1)
  wolf_i <- get_effect("^Wolf time")
  jackal_i <- get_effect("^Jackal time")
  interaction_i <- get_effect("Wolf x jackal")
  results_lines <- c(
    results_lines,
    paste0("### ", class_label),
    "",
    paste0(
      "The wolf-by-time smooth had EDF = ",
      fmt(wolf_i$EDF_if_smooth, 2),
      " and approximate P = ",
      fmt_p(wolf_i$P),
      "; the jackal-by-time smooth had EDF = ",
      fmt(jackal_i$EDF_if_smooth, 2),
      " and P = ",
      fmt_p(jackal_i$P),
      "; and the wolf-by-jackal-by-time smooth had EDF = ",
      fmt(interaction_i$EDF_if_smooth, 2),
      " and P = ",
      fmt_p(interaction_i$P),
      "."
    ),
    paste0(
      "Across the four -1/+1 SD combinations, the modeled proportion of ",
      "normalized activity between 06:00 and 18:00 apparent solar time ranged ",
      "from ",
      fmt(min(metric_i$predicted_06_18_proportion), 2),
      " to ",
      fmt(max(metric_i$predicted_06_18_proportion), 2),
      ". Principal point-estimate peaks ranged from ",
      fmt(min(metric_i$principal_peak_hour), 1),
      " to ",
      fmt(max(metric_i$principal_peak_hour), 1),
      " h."
    ),
    paste0(
      "Simulation intervals excluded zero for ",
      fmt(max(diff_i$supported_hours), 1),
      " h at most among the four conditional predator difference curves."
    ),
    ""
  )
}

season_control <- season_shape_sensitivity %>%
  filter(effect == "Late versus Early curve-shape control")
season_predator <- season_shape_sensitivity %>%
  filter(effect != "Late versus Early curve-shape control")
results_lines <- c(
  results_lines,
  "## Early/Late curve-shape sensitivity",
  "",
  paste0(
    "We added a cyclic Late-versus-Early difference smooth using the frozen ",
    "survey midpoint (`",
    survey_split,
    "`) while retaining round-specific intensity effects. The period smooth ",
    "was unsupported for juveniles (EDF = ",
    fmt(season_control$season_adjusted_edf[season_control$class == "Juvenile"], 2),
    ", P = ",
    fmt_p(season_control$season_adjusted_P[season_control$class == "Juvenile"]),
    ") but supported for adult females (EDF = ",
    fmt(season_control$season_adjusted_edf[season_control$class == "Adult female"], 2),
    ", P ",
    fmt_p(season_control$season_adjusted_P[season_control$class == "Adult female"]),
    ") and adult males (EDF = ",
    fmt(season_control$season_adjusted_edf[season_control$class == "Adult male"], 2),
    ", P ",
    fmt_p(season_control$season_adjusted_P[season_control$class == "Adult male"]),
    ")."
  ),
  paste0(
    "The supported predator effects were stable: juvenile wolf-by-jackal-by-time ",
    "P ",
    fmt_p(season_predator$season_adjusted_P[
      season_predator$class == "Juvenile" &
        grepl("Wolf x jackal", season_predator$effect)
    ]),
    ", adult-female wolf-by-time P ",
    fmt_p(season_predator$season_adjusted_P[
      season_predator$class == "Adult female" &
        grepl("^Wolf time", season_predator$effect)
    ]),
    ", and adult-male wolf-by-jackal-by-time P = ",
    fmt_p(season_predator$season_adjusted_P[
      season_predator$class == "Adult male" &
        grepl("Wolf x jackal", season_predator$effect)
    ]),
    ". Thus, broad seasonal changes in diel shape did not account for the ",
    "predator-associated smooths."
  ),
  ""
)

results_lines <- c(
  results_lines,
  "## Relationship to the categorical kernel analysis",
  "",
  paste0(
    "The existing kernel analysis supported 0 of 3 eligible context comparisons ",
    "for juveniles, 4 of 6 for adult females, and 4 of 6 for adult males. The ",
    "GAMM should be judged by stability, continuous-gradient inference, and the ",
    "coherence of predicted temporal changes rather than exact agreement with ",
    "those categorical tests."
  ),
  "",
  "## Limitations",
  "",
  "`Days.Operational` was the best available exposure measure but assumes approximately equal exposure across the 24 solar-hour bins within a camera-round.",
  "Predator indices are camera-level relative activity gradients, so associations do not demonstrate direct predator avoidance or causal behavioral displacement.",
  "The juvenile sample is much smaller than the adult samples; unsupported juvenile effects should not be interpreted as proof of no response.",
  paste0(
    "The juvenile model retained substantial within-camera-round residual temporal ",
    "correlation (mean lag-1 Pearson correlation = ",
    fmt(diagnostics$mean_within_camera_round_lag1_residual_correlation[
      diagnostics$class == "Juvenile"
    ]),
    "), so its interaction should be treated as exploratory despite its smooth-term test."
  ),
  "Adult global smooth tests were supported, but their normalized pointwise difference intervals did not exclude zero at individual times; inference is therefore distribution-wide rather than evidence for sharply bounded hourly windows.",
  "Approximate smooth-term P values are interpreted with shrinkage, diagnostics, and simulation-based effect curves rather than as stand-alone evidence.",
  "",
  "## Recommendation for manuscript",
  "",
  paste0("**PRIMARY DIEL ANALYSIS RECOMMENDATION: ", primary_recommendation, "**"),
  "",
  if (primary_recommendation == "GAMM") {
    paste0(
      "The circular GAMMs converged, passed the prespecified basis and ",
      "concurvity checks, and yielded temporally coherent predator-associated ",
      "differences. They provide stronger primary inference by retaining the ",
      "continuous predator gradients, repeated camera-round structure, and effort."
      , " The supported predator smooths were also stable after controlling for ",
      "Early/Late changes in diel curve shape."
    )
  } else {
    paste0(
      "The circular GAMMs did not meet all prespecified stability and ecological-",
      "coherence criteria. The categorical kernel analysis should remain primary, ",
      "with these models treated as a sensitivity analysis rather than promoted ",
      "because of isolated approximate P values."
    )
  },
  "",
  "### Manuscript-ready Methods",
  "",
  paste0(
    "To test whether gazelle diel timing varied continuously across predator-activity ",
    "gradients, we fitted separate hierarchical cyclic generalized additive models ",
    "for juvenile, adult-female, and adult-male independent events following ",
    "Iannarilli et al. (2025). We aggregated 30-min independent class-events into ",
    "24 one-hour apparent-solar-time bins for each positive-effort camera-round and ",
    "retained zero-event cells. Models used a negative-binomial distribution with ",
    "a log link and the logarithm of operational camera-days as an offset. Each ",
    "model included ordinary standardized wolf and jackal activity terms and their ",
    "interaction, factor(round), a baseline cyclic cubic smooth of apparent solar ",
    "time, and cyclic varying-coefficient smooths for wolf activity, jackal activity, ",
    "and their interaction. Random-intercept smooths accounted for camera-round, ",
    "camera, and sampling-polygon clustering. Supported shared-penalty, camera-specific ",
    "cyclic deviations accounted for site-to-site variation in adult curve shape; ",
    "we omitted these deviations from the juvenile model because 180 events were ",
    "insufficient to estimate 60 camera-specific curves reliably. We fitted models ",
    "with `mgcv::bam` using fREML, shrinkage selection, and conservative basis ",
    "dimensions (k = 8 for adults). A diagnostic expansion from k = 8 to k = 10 ",
    "for juveniles resolved a marginal basis-adequacy warning without changing the ",
    "ecological inference. We evaluated convergence, basis adequacy, ",
    "concurvity, residual patterns, zero frequencies, and agreement between fitted ",
    "and observed hourly rates. As a seasonal sensitivity, we added a cyclic ",
    "Late-versus-Early difference smooth based on the frozen 17 October survey ",
    "midpoint; this controlled period-specific curve shape without fitting ",
    "predator-by-period interactions. We generated population-level curves at all ",
    "combinations of -1 and +1 SD wolf and jackal activity and obtained 95% ",
    "intervals from ",
    format(coefficient_simulations, big.mark = ","),
    " simulations of the model coefficients and covariance matrix."
  ),
  "",
  "### Manuscript-ready Results",
  "",
  paste0(
    "The hierarchical circular analysis included ",
    format(nrow(hourly_data), big.mark = ","),
    " camera-round-hour cells and ",
    format(sum(expected_events), big.mark = ","),
    " independent demographic-class events. Juvenile diel shape varied with the ",
    "wolf by jackal gradient (EDF = ",
    fmt(ecological_smooths$EDF_if_smooth[
      ecological_smooths$class == "Juvenile" &
        ecological_smooths$effect == "Wolf x jackal time-varying effect"
    ], 2),
    ", P ",
    fmt_p(ecological_smooths$P[
      ecological_smooths$class == "Juvenile" &
        ecological_smooths$effect == "Wolf x jackal time-varying effect"
    ]),
    "), whereas its wolf-only and jackal-only timing smooths were unsupported. ",
    "Predicted juvenile activity peaked in the morning (8.6 h) when predator ",
    "indices were concordant (both low or both high) and in late afternoon ",
    "(16.1 h) when indices were discordant. Given only 180 juvenile events and ",
    "remaining residual temporal dependence, we treated this interaction as exploratory."
  ),
  paste0(
    "Adult-female diel shape varied with wolf activity (EDF = ",
    fmt(ecological_smooths$EDF_if_smooth[
      ecological_smooths$class == "Adult female" &
        ecological_smooths$effect == "Wolf time-varying effect"
    ], 2),
    ", P ",
    fmt_p(ecological_smooths$P[
      ecological_smooths$class == "Adult female" &
        ecological_smooths$effect == "Wolf time-varying effect"
    ]),
    "), but not with jackal activity or their interaction. The point-estimate ",
    "principal peak shifted from 9.7 h at low wolf activity to 16.7 h at high ",
    "wolf activity, while the predicted 06:00-18:00 proportion remained nearly ",
    "constant (0.73), indicating redistribution within the diurnal-crepuscular ",
    "period rather than increased nocturnality."
  ),
  paste0(
    "Adult-male diel shape varied with the wolf by jackal gradient (EDF = ",
    fmt(ecological_smooths$EDF_if_smooth[
      ecological_smooths$class == "Adult male" &
        ecological_smooths$effect == "Wolf x jackal time-varying effect"
    ], 2),
    ", P = ",
    fmt_p(ecological_smooths$P[
      ecological_smooths$class == "Adult male" &
        ecological_smooths$effect == "Wolf x jackal time-varying effect"
    ]),
    "); wolf-only and jackal-only timing smooths were unsupported. Modeled ",
    "daytime activity was higher when predator indices were concordant (0.64) ",
    "than discordant (0.59), although the principal peak remained near 16.6-16.9 h. ",
    "The juvenile interaction, female wolf effect, and male interaction remained ",
    "supported after controlling for Early/Late curve shape. For both adult ",
    "classes, fully propagated pointwise intervals overlapped zero throughout ",
    "the day, so these results indicate global changes in curve shape rather than ",
    "precisely delimited hourly shifts."
  ),
  "",
  "### Proposed figure caption",
  "",
  paste0(
    "**Figure X. Modeled diel activity of mountain gazelle demographic classes ",
    "across continuous wolf and jackal activity gradients.** Curves show ",
    "population-level normalized predictions from class-specific negative-binomial ",
    "hierarchical cyclic GAMMs at -1 and +1 SD of each predator index; shading ",
    "shows simulation-based 95% intervals and rugs show observed independent ",
    "class-events. Curves integrate to one within each predator combination and ",
    "therefore describe event timing rather than event abundance."
  ),
  "",
  "### Kernel sensitivity placement",
  "",
  if (primary_recommendation == "GAMM") {
    paste0(
      "Move the four median-split kernel curves, pairwise overlap coefficients, ",
      "camera-cluster bootstrap intervals, and `compareCkern` tests to Supporting ",
      "Information. In the main text, state: 'Median-split kernel analyses provided ",
      "a categorical sensitivity check on the continuous hierarchical GAMM results.'"
    )
  } else {
    "Retain the current kernel Methods, Results, and principal figure; report the continuous GAMMs only as a model-based sensitivity analysis if space permits."
  }
)
writeLines(results_lines, out_file("circular_GAMM_results.md"))

writeLines(
  capture.output(sessionInfo()),
  out_file("circular_GAMM_session_info.txt")
)

cat("Circular GAMM analysis completed.\n")
cat("Primary recommendation:", primary_recommendation, "\n")
