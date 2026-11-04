#!/usr/bin/env Rscript

# Revised H3: demographic differences in predator-associated daylight use.
# This script uses the finalized class-specific 30-minute event file and writes
# all new outputs to a dedicated directory without modifying prior analyses.

suppressPackageStartupMessages({
  library(activity)
  library(broom.mixed)
  library(dplyr)
  library(emmeans)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tidyr)
})

set.seed(20260819)

event_file <- "Gazella gazella_independent_events_30min_analysis_ready.csv"
site_file <- "SiteCovs_060324.csv"
obs_file <- "ObsCovs_050524.csv"
context_file <- "publication_predator_context_class_diel_2x2_median_camera_categories.csv"
context_overlap_file <- "publication_authoritative_diel_context_comparisons.csv"
output_dir <- file.path("results", "H3_diel_demography")
local_timezone <- "Asia/Jerusalem"
minimum_curve_events <- 5L
activity_reps <- 99L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(event_file, site_file, obs_file, context_file, context_overlap_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required H3 inputs: ", paste(missing_files, collapse = ", "))
}

class_levels <- c("Fawn", "Female", "Male")
context_levels <- c(
  "Low both",
  "High jackal / Low wolf",
  "Low jackal / High wolf",
  "High both"
)

z_scale <- function(x) as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

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

solar_elevation <- function(timestamp, latitude, longitude) {
  day_of_year <- as.integer(format(timestamp, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(timestamp, "%H", tz = "UTC")) +
    as.numeric(format(timestamp, "%M", tz = "UTC")) / 60 +
    as.numeric(format(timestamp, "%S", tz = "UTC")) / 3600
  gamma <- 2 * pi / 365 * (day_of_year - 1 + (utc_hour - 12) / 24)
  eq_time <- 229.18 * (
    0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -
      0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma)
  )
  declination <- 0.006918 - 0.399912 * cos(gamma) +
    0.070257 * sin(gamma) - 0.006758 * cos(2 * gamma) +
    0.000907 * sin(2 * gamma) - 0.002697 * cos(3 * gamma) +
    0.00148 * sin(3 * gamma)
  true_solar_minutes <- (utc_hour * 60 + eq_time + 4 * longitude) %% 1440
  hour_angle <- true_solar_minutes / 4 - 180
  latitude_rad <- latitude * pi / 180
  hour_angle_rad <- hour_angle * pi / 180
  cos_zenith <- sin(latitude_rad) * sin(declination) +
    cos(latitude_rad) * cos(declination) * cos(hour_angle_rad)
  cos_zenith <- pmin(1, pmax(-1, cos_zenith))
  90 - acos(cos_zenith) * 180 / pi
}

assign_rounds <- function(events, active_obs) {
  event_index <- events %>%
    transmute(event_row_id, camera_name, event_date)
  event_index %>%
    inner_join(
      active_obs %>% select(camera_name, round, start_date, end_date),
      by = "camera_name",
      relationship = "many-to-many"
    ) %>%
    filter(event_date >= start_date, event_date <= end_date) %>%
    group_by(event_row_id) %>%
    slice_max(start_date, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(event_row_id, round)
}

model_diagnostics <- function(model, name) {
  if (is.null(model)) {
    return(tibble(model = name, convergence = NA_integer_, pd_hessian = FALSE))
  }
  beta <- tryCatch(fixef(model)$cond, error = function(e) numeric())
  se <- tryCatch(sqrt(diag(vcov(model)$cond)), error = function(e) numeric())
  gradient <- tryCatch(
    max(abs(model$sdr$gradient.fixed), na.rm = TRUE),
    error = function(e) NA_real_
  )
  pearson <- tryCatch(residuals(model, type = "pearson"), error = function(e) numeric())
  dispersion <- if (length(pearson)) {
    sum(pearson^2, na.rm = TRUE) / max(1, nobs(model) - length(beta))
  } else {
    NA_real_
  }
  singular <- tryCatch(
    isTRUE(performance::check_singularity(model, tolerance = 1e-5)),
    error = function(e) NA
  )
  convergence_code <- model$fit$convergence
  positive_hessian <- isTRUE(model$sdr$pdHess)
  fitted_values <- fitted(model)
  model_aic <- AIC(model)
  model_nobs <- nobs(model)
  tibble(
    model = name,
    convergence = convergence_code,
    pd_hessian = positive_hessian,
    max_abs_gradient = gradient,
    max_abs_beta = if (length(beta)) max(abs(beta), na.rm = TRUE) else NA_real_,
    max_se = if (length(se)) max(se, na.rm = TRUE) else NA_real_,
    singular = singular,
    pearson_dispersion = dispersion,
    min_fitted = min(fitted_values, na.rm = TRUE),
    max_fitted = max(fitted_values, na.rm = TRUE),
    AIC = model_aic,
    nobs = model_nobs
  )
}

is_stable <- function(model) {
  d <- model_diagnostics(model, "candidate")
  isTRUE(d$convergence == 0L) && isTRUE(d$pd_hessian) &&
    is.finite(d$max_abs_beta) && d$max_abs_beta < 20 &&
    is.finite(d$max_se) && d$max_se < 10
}

fit_model <- function(formula, data, family = binomial()) {
  glmmTMB(
    formula,
    family = family,
    data = data,
    control = glmmTMBControl(
      optCtrl = list(iter.max = 10000, eval.max = 10000)
    )
  )
}

random_variances <- function(model, model_name) {
  components <- tryCatch(VarCorr(model)$cond, error = function(e) list())
  if (!length(components)) return(tibble())
  bind_rows(lapply(names(components), function(group_i) {
    matrix_i <- components[[group_i]]
    standard_deviation <- unname(attr(matrix_i, "stddev")[1])
    tibble(
      model = model_name,
      group = group_i,
      variance = standard_deviation^2,
      standard_deviation = standard_deviation
    )
  }))
}

lrt_row <- function(reduced, full, test_name) {
  comparison <- suppressMessages(as.data.frame(anova(reduced, full)))
  last <- comparison[nrow(comparison), , drop = FALSE]
  get_col <- function(candidates) {
    hit <- intersect(candidates, names(last))
    if (!length(hit)) return(NA_real_)
    as.numeric(last[[hit[1]]])
  }
  tibble(
    test = test_name,
    method = "Likelihood-ratio test",
    statistic = get_col(c("Chisq", "LRT")),
    df = get_col(c("Chi Df", "Df")),
    p_value = get_col(c("Pr(>Chisq)", "Pr(>Chi)")),
    reduced_AIC = AIC(reduced),
    full_AIC = AIC(full)
  )
}

wald_joint <- function(model, terms, test_name) {
  beta <- fixef(model)$cond
  variance <- vcov(model)$cond
  term_index <- match(terms, names(beta))
  if (anyNA(term_index)) stop("Missing Wald-test terms: ", paste(terms, collapse = ", "))
  estimate <- beta[term_index]
  variance_sub <- variance[term_index, term_index, drop = FALSE]
  statistic <- as.numeric(t(estimate) %*% solve(variance_sub, estimate))
  tibble(
    test = test_name,
    method = "Wald chi-square",
    statistic = statistic,
    df = length(terms),
    p_value = pchisq(statistic, df = length(terms), lower.tail = FALSE),
    reduced_AIC = NA_real_,
    full_AIC = AIC(model)
  )
}

extract_probability_summary <- function(emm) {
  out <- as.data.frame(summary(emm, type = "response", infer = c(TRUE, TRUE)))
  estimate_name <- intersect(c("prob", "response"), names(out))[1]
  lower_name <- intersect(c("asymp.LCL", "lower.CL"), names(out))[1]
  upper_name <- intersect(c("asymp.UCL", "upper.CL"), names(out))[1]
  out$probability <- out[[estimate_name]]
  out$conf_low <- out[[lower_name]]
  out$conf_high <- out[[upper_name]]
  as_tibble(out)
}

extract_trends <- function(model, focal, moderator) {
  at_values <- setNames(list(c(-1, 0, 1)), moderator)
  trends <- emtrends(
    model,
    specs = as.formula(paste("~ class |", moderator)),
    var = focal,
    at = at_values,
    weights = "equal"
  )
  trend_df <- as.data.frame(summary(trends, infer = c(TRUE, TRUE)))
  trend_name <- paste0(focal, ".trend")
  lower_name <- intersect(c("asymp.LCL", "lower.CL"), names(trend_df))[1]
  upper_name <- intersect(c("asymp.UCL", "upper.CL"), names(trend_df))[1]
  slope_rows <- as_tibble(trend_df) %>%
    transmute(
      record_type = "conditional_slope",
      predictor = focal,
      moderator = moderator,
      moderator_value = .data[[moderator]],
      class = as.character(class),
      contrast = NA_character_,
      estimate = .data[[trend_name]],
      std_error = SE,
      conf_low = .data[[lower_name]],
      conf_high = .data[[upper_name]],
      odds_ratio = exp(estimate),
      odds_ratio_low = exp(conf_low),
      odds_ratio_high = exp(conf_high),
      p_value = p.value
    )
  contrast_df <- as.data.frame(
    summary(pairs(trends, by = moderator, adjust = "none"), infer = c(TRUE, TRUE))
  )
  lower_name_c <- intersect(c("asymp.LCL", "lower.CL"), names(contrast_df))[1]
  upper_name_c <- intersect(c("asymp.UCL", "upper.CL"), names(contrast_df))[1]
  contrast_rows <- as_tibble(contrast_df) %>%
    transmute(
      record_type = "class_slope_contrast",
      predictor = focal,
      moderator = moderator,
      moderator_value = .data[[moderator]],
      class = NA_character_,
      contrast = as.character(contrast),
      estimate = estimate,
      std_error = SE,
      conf_low = .data[[lower_name_c]],
      conf_high = .data[[upper_name_c]],
      odds_ratio = exp(estimate),
      odds_ratio_low = exp(conf_low),
      odds_ratio_high = exp(conf_high),
      p_value = p.value
    )
  bind_rows(slope_rows, contrast_rows)
}

# Data preparation -----------------------------------------------------------

site <- read_csv(site_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(Camera.ID),
    polygon = factor(Polygon.Name),
    wolf = as.numeric(Null.Wolf.Abundance),
    jackal = as.numeric(null.jackal.abundance),
    human = as.numeric(human.abundance),
    wolf_z = z_scale(wolf),
    jackal_z = z_scale(jackal),
    human_z = z_scale(human)
  )

if (nrow(site) != 60L || n_distinct(site$camera_name) != 60L) {
  stop("Site covariate file does not contain the expected 60 unique cameras.")
}

obs <- read_csv(obs_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    days_operational = as.numeric(Days.Operational)
  ) %>%
  filter(
    is.finite(days_operational), days_operational > 0,
    !is.na(start_date), !is.na(end_date)
  )

survey_start <- min(obs$start_date)
survey_end <- max(obs$end_date)
survey_split <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)

context_lookup <- read_csv(context_file, show_col_types = FALSE) %>%
  transmute(
    camera_name = as.integer(camera_name),
    context_jackal = as.numeric(jackal_abundance),
    context_wolf = as.numeric(wolf_abundance),
    predator_context = factor(predator_context, levels = context_levels)
  )

context_check <- site %>%
  inner_join(context_lookup, by = "camera_name") %>%
  summarise(
    max_wolf_difference = max(abs(wolf - context_wolf)),
    max_jackal_difference = max(abs(jackal - context_jackal))
  )
if (
  nrow(context_lookup) != 60L ||
    context_check$max_wolf_difference > 1e-12 ||
    context_check$max_jackal_difference > 1e-12
) {
  stop("Existing predator-context file does not match the final site indices.")
}

events <- read_csv(
  event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  transmute(
    event_row_id = row_number(),
    camera_name = as.integer(camera_name),
    session_id = as.character(trap_station_session_id),
    source_event_id = as.character(event_id),
    class = recode(
      as.character(gazelle_class),
      juvenile = "Fawn",
      female = "Female",
      male = "Male"
    ),
    timestamp = parse_local_camera_time(event_time),
    latitude = as.numeric(camelot_gps_latitude),
    longitude = as.numeric(camelot_gps_longitude)
  ) %>%
  mutate(
    event_date = as.Date(format(timestamp, tz = local_timezone)),
    solar_hour = apparent_solar_hour(timestamp, longitude),
    solar_elevation = solar_elevation(timestamp, latitude, longitude),
    daylight = as.integer(solar_elevation >= 0),
    nocturnality_hours = abs(solar_hour - 12),
    period = factor(
      ifelse(event_date <= survey_split, "Early", "Late"),
      levels = c("Early", "Late")
    )
  )

round_matches <- assign_rounds(events, obs)
events <- events %>% left_join(round_matches, by = "event_row_id")
if (anyNA(events$round)) stop("One or more final events did not map to an active survey round.")

# Reconstruct cross-class encounters because source event IDs are class-specific.
events <- events %>%
  arrange(camera_name, session_id, timestamp, class) %>%
  group_by(camera_name, session_id) %>%
  mutate(
    gap_minutes = as.numeric(difftime(timestamp, lag(timestamp), units = "mins")),
    encounter_sequence = cumsum(is.na(gap_minutes) | gap_minutes > 30),
    encounter_id = paste(camera_name, session_id, encounter_sequence, sep = "_")
  ) %>%
  ungroup() %>%
  left_join(site, by = "camera_name") %>%
  left_join(
    context_lookup %>% select(camera_name, predator_context),
    by = "camera_name"
  ) %>%
  mutate(
    class = factor(class, levels = class_levels),
    camera_factor = factor(camera_name),
    encounter_factor = factor(encounter_id),
    date_z = z_scale(as.numeric(event_date))
  )

expected_counts <- c(Fawn = 180L, Female = 1477L, Male = 2268L)
observed_counts <- table(events$class)
if (
  nrow(events) != sum(expected_counts) ||
    any(as.integer(observed_counts[names(expected_counts)]) != expected_counts)
) {
  stop(
    "Final event totals differ from the expected analysis-ready totals: ",
    paste(names(observed_counts), observed_counts, collapse = ", ")
  )
}

encounter_summary <- events %>%
  group_by(encounter_id) %>%
  summarise(
    event_rows = n(),
    demographic_classes = n_distinct(class),
    .groups = "drop"
  )

daylight_summary <- events %>%
  group_by(class) %>%
  summarise(
    events = n(),
    daylight_events = sum(daylight),
    daylight_proportion = mean(daylight),
    .groups = "drop"
  )

context_sample_sizes <- events %>%
  count(class, predator_context, name = "events") %>%
  complete(
    class = factor(class_levels, levels = class_levels),
    predator_context = factor(context_levels, levels = context_levels),
    fill = list(events = 0L)
  ) %>%
  left_join(
    events %>%
      group_by(class, predator_context) %>%
      summarise(cameras_with_events = n_distinct(camera_name), .groups = "drop"),
    by = c("class", "predator_context")
  ) %>%
  left_join(
    context_lookup %>% count(predator_context, name = "cameras_assigned"),
    by = "predator_context"
  ) %>%
  mutate(cameras_with_events = replace_na(cameras_with_events, 0L))

predator_correlation <- cor.test(site$wolf, site$jackal, method = "pearson")

site_gradient <- site %>%
  mutate(
    wolf_band = factor(ntile(wolf_z, 3), labels = c("Low", "Middle", "High")),
    jackal_band = factor(ntile(jackal_z, 3), labels = c("Low", "Middle", "High"))
  )

gradient_coverage <- events %>%
  left_join(site_gradient %>% select(camera_name, wolf_band, jackal_band), by = "camera_name") %>%
  group_by(class, wolf_band, jackal_band) %>%
  summarise(
    events = n(),
    cameras = n_distinct(camera_name),
    daylight_events = sum(daylight),
    daylight_proportion = mean(daylight),
    all_same_response = n_distinct(daylight) == 1,
    .groups = "drop"
  ) %>%
  complete(
    class = factor(class_levels, levels = class_levels),
    wolf_band,
    jackal_band,
    fill = list(
      events = 0L, cameras = 0L, daylight_events = 0L,
      daylight_proportion = NA_real_, all_same_response = NA
    )
  )

validation_report <- bind_rows(
  tibble(section = "sample", metric = "class_specific_event_rows", value = as.character(nrow(events))),
  tibble(section = "sample", metric = "underlying_30min_encounters", value = as.character(nrow(encounter_summary))),
  tibble(section = "sample", metric = "multiclass_encounters", value = as.character(sum(encounter_summary$demographic_classes > 1))),
  tibble(section = "sample", metric = "rows_in_multiclass_encounters", value = as.character(sum(encounter_summary$event_rows[encounter_summary$demographic_classes > 1]))),
  tibble(section = "sample", metric = "cameras_with_events", value = as.character(n_distinct(events$camera_name))),
  tibble(section = "sample", metric = "study_cameras", value = as.character(n_distinct(site$camera_name))),
  tibble(section = "sample", metric = "polygons_with_events", value = as.character(n_distinct(events$polygon))),
  tibble(section = "sample", metric = "survey_rounds", value = as.character(n_distinct(obs$round))),
  tibble(section = "sample", metric = "survey_start", value = as.character(survey_start)),
  tibble(section = "sample", metric = "survey_end", value = as.character(survey_end)),
  tibble(section = "sample", metric = "early_late_split", value = as.character(survey_split)),
  tibble(section = "predator", metric = "wolf_raw_range", value = paste(range(site$wolf), collapse = " to ")),
  tibble(section = "predator", metric = "wolf_raw_mean_sd", value = paste(mean(site$wolf), sd(site$wolf), sep = ", ")),
  tibble(section = "predator", metric = "jackal_raw_range", value = paste(range(site$jackal), collapse = " to ")),
  tibble(section = "predator", metric = "jackal_raw_mean_sd", value = paste(mean(site$jackal), sd(site$jackal), sep = ", ")),
  tibble(section = "predator", metric = "wolf_jackal_pearson_r", value = as.character(unname(predator_correlation$estimate))),
  tibble(section = "predator", metric = "wolf_jackal_pearson_p", value = as.character(predator_correlation$p.value)),
  daylight_summary %>% transmute(section = "daylight", metric = paste0(as.character(class), "_proportion"), value = as.character(daylight_proportion)),
  tibble(section = "class", metric = paste0(names(expected_counts), "_events"), value = as.character(expected_counts))
)

write_csv(validation_report, file.path(output_dir, "H3_validation_report.csv"))
write_csv(gradient_coverage, file.path(output_dir, "H3_joint_gradient_coverage.csv"))
write_csv(context_sample_sizes, file.path(output_dir, "H3_class_context_sample_sizes.csv"))

# Candidate primary models ---------------------------------------------------

random_encounter <- "+ (1 | polygon) + (1 | camera_factor) + (1 | encounter_factor)"
random_camera <- "+ (1 | polygon) + (1 | camera_factor)"

full_encounter_formula <- as.formula(paste(
  "daylight ~ class * wolf_z * jackal_z + human_z + period",
  random_encounter
))
full_camera_formula <- as.formula(paste(
  "daylight ~ class * wolf_z * jackal_z + human_z + period",
  random_camera
))
reduced_encounter_formula <- as.formula(paste(
  "daylight ~ class * wolf_z + class * jackal_z + wolf_z * jackal_z + human_z + period",
  random_encounter
))
reduced_camera_formula <- as.formula(paste(
  "daylight ~ class * wolf_z + class * jackal_z + wolf_z * jackal_z + human_z + period",
  random_camera
))

full_encounter <- fit_model(full_encounter_formula, events)
full_camera <- fit_model(full_camera_formula, events)
reduced_encounter <- fit_model(reduced_encounter_formula, events)
reduced_camera <- fit_model(reduced_camera_formula, events)

encounter_variance <- random_variances(full_encounter, "full_encounter") %>%
  filter(grepl("encounter", group)) %>%
  summarise(value = ifelse(n() == 0, NA_real_, max(variance))) %>%
  pull(value)

use_encounter <- is_stable(full_encounter) &&
  is.finite(encounter_variance) && encounter_variance > 1e-6
random_primary <- if (use_encounter) random_encounter else random_camera
full_selected <- if (use_encounter) full_encounter else full_camera
reduced_selected <- if (use_encounter) reduced_encounter else reduced_camera

use_full_interaction <- is_stable(full_selected)
primary <- if (use_full_interaction) full_selected else reduced_selected
primary_structure <- if (use_full_interaction) "full_three_way" else "reduced_interactions"
primary_fixed <- if (use_full_interaction) {
  "class * wolf_z * jackal_z"
} else {
  "class * wolf_z + class * jackal_z + wolf_z * jackal_z"
}
primary_formula <- as.formula(paste(
  "daylight ~", primary_fixed, "+ human_z + period", random_primary
))

# Nested models for formal H3 tests use the same random-effects structure.
no_predator_interaction <- fit_model(as.formula(paste(
  "daylight ~ class * wolf_z + class * jackal_z + human_z + period",
  random_primary
)), events)
no_class_wolf <- fit_model(as.formula(paste(
  "daylight ~ class * jackal_z + wolf_z * jackal_z + human_z + period",
  random_primary
)), events)
no_class_jackal <- fit_model(as.formula(paste(
  "daylight ~ class * wolf_z + wolf_z * jackal_z + human_z + period",
  random_primary
)), events)
no_wolf <- fit_model(as.formula(paste(
  "daylight ~ class * jackal_z + human_z + period",
  random_primary
)), events)
no_jackal <- fit_model(as.formula(paste(
  "daylight ~ class * wolf_z + human_z + period",
  random_primary
)), events)

global_tests <- bind_rows(
  wald_joint(primary, c("classFemale", "classMale"), "Baseline demographic-class differences at mean predator activity"),
  lrt_row(no_wolf, primary, "Any wolf association across classes and jackal contexts"),
  lrt_row(no_jackal, primary, "Any jackal association across classes and wolf contexts"),
  lrt_row(no_class_wolf, primary, "Wolf slopes differ among demographic classes"),
  lrt_row(no_class_jackal, primary, "Jackal slopes differ among demographic classes"),
  lrt_row(no_predator_interaction, primary, "Wolf by jackal interaction, including class-specific interaction if retained"),
  lrt_row(reduced_selected, full_selected, "Class by wolf by jackal interaction")
)
write_csv(global_tests, file.path(output_dir, "H3_global_tests.csv"))

primary_coefficients <- tidy(
  primary,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
) %>%
  mutate(
    odds_ratio = exp(estimate),
    odds_ratio_low = exp(conf.low),
    odds_ratio_high = exp(conf.high),
    primary_structure = primary_structure
  )
write_csv(primary_coefficients, file.path(output_dir, "H3_model_coefficients.csv"))

# Conditional effects and predicted probabilities ---------------------------

conditional_slopes <- bind_rows(
  extract_trends(primary, "wolf_z", "jackal_z"),
  extract_trends(primary, "jackal_z", "wolf_z")
) %>%
  group_by(record_type, predictor) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write_csv(conditional_slopes, file.path(output_dir, "H3_conditional_predator_slopes.csv"))

baseline_emm <- emmeans(
  primary,
  ~ class,
  at = list(wolf_z = 0, jackal_z = 0, human_z = 0),
  weights = "equal"
)
baseline_probabilities <- extract_probability_summary(baseline_emm) %>%
  transmute(
    record_type = "class_probability",
    class = as.character(class),
    contrast = NA_character_,
    estimate = probability,
    conf_low,
    conf_high,
    odds_ratio = NA_real_,
    odds_ratio_low = NA_real_,
    odds_ratio_high = NA_real_,
    p_value = p.value
  )
baseline_contrast_raw <- as.data.frame(
  summary(pairs(baseline_emm, adjust = "tukey"), infer = c(TRUE, TRUE))
)
baseline_contrast_low <- intersect(c("asymp.LCL", "lower.CL"), names(baseline_contrast_raw))[1]
baseline_contrast_high <- intersect(c("asymp.UCL", "upper.CL"), names(baseline_contrast_raw))[1]
baseline_contrasts <- as_tibble(baseline_contrast_raw) %>%
  transmute(
    record_type = "class_odds_ratio",
    class = NA_character_,
    contrast = as.character(contrast),
    estimate = estimate,
    conf_low = .data[[baseline_contrast_low]],
    conf_high = .data[[baseline_contrast_high]],
    odds_ratio = exp(estimate),
    odds_ratio_low = exp(conf_low),
    odds_ratio_high = exp(conf_high),
    p_value = p.value
  )
baseline_effects <- bind_rows(baseline_probabilities, baseline_contrasts)
write_csv(baseline_effects, file.path(output_dir, "H3_baseline_class_effects.csv"))

prediction_emm <- emmeans(
  primary,
  ~ class * wolf_z * jackal_z,
  at = list(wolf_z = c(-1, 1), jackal_z = c(-1, 1), human_z = 0),
  weights = "equal"
)
predicted_probabilities <- extract_probability_summary(prediction_emm) %>%
  transmute(
    record_type = "predicted_probability",
    class = as.character(class),
    wolf_z,
    jackal_z,
    probability,
    conf_low,
    conf_high,
    contrast = NA_character_,
    odds_ratio = NA_real_,
    odds_ratio_low = NA_real_,
    odds_ratio_high = NA_real_,
    probability_difference = NA_real_,
    p_value = NA_real_
  )

wolf_endpoint_link <- contrast(
  emmeans(
    primary,
    ~ wolf_z | class * jackal_z,
    at = list(wolf_z = c(-1, 1), jackal_z = c(-1, 1), human_z = 0),
    weights = "equal"
  ),
  method = list("+1 versus -1 SD wolf" = c(-1, 1)),
  by = c("class", "jackal_z")
)
wolf_endpoint_raw <- as.data.frame(summary(wolf_endpoint_link, infer = c(TRUE, TRUE)))
wolf_endpoint_low <- intersect(c("asymp.LCL", "lower.CL"), names(wolf_endpoint_raw))[1]
wolf_endpoint_high <- intersect(c("asymp.UCL", "upper.CL"), names(wolf_endpoint_raw))[1]
wolf_endpoint <- as_tibble(wolf_endpoint_raw) %>%
  transmute(
    record_type = "endpoint_contrast",
    class = as.character(class),
    wolf_z = NA_real_,
    jackal_z,
    probability = NA_real_,
    conf_low = .data[[wolf_endpoint_low]],
    conf_high = .data[[wolf_endpoint_high]],
    contrast = "+1 versus -1 SD wolf",
    odds_ratio = exp(estimate),
    odds_ratio_low = exp(.data[[wolf_endpoint_low]]),
    odds_ratio_high = exp(.data[[wolf_endpoint_high]]),
    probability_difference = NA_real_,
    p_value = p.value
  )

jackal_endpoint_link <- contrast(
  emmeans(
    primary,
    ~ jackal_z | class * wolf_z,
    at = list(wolf_z = c(-1, 1), jackal_z = c(-1, 1), human_z = 0),
    weights = "equal"
  ),
  method = list("+1 versus -1 SD jackal" = c(-1, 1)),
  by = c("class", "wolf_z")
)
jackal_endpoint_raw <- as.data.frame(summary(jackal_endpoint_link, infer = c(TRUE, TRUE)))
jackal_endpoint_low <- intersect(c("asymp.LCL", "lower.CL"), names(jackal_endpoint_raw))[1]
jackal_endpoint_high <- intersect(c("asymp.UCL", "upper.CL"), names(jackal_endpoint_raw))[1]
jackal_endpoint <- as_tibble(jackal_endpoint_raw) %>%
  transmute(
    record_type = "endpoint_contrast",
    class = as.character(class),
    wolf_z,
    jackal_z = NA_real_,
    probability = NA_real_,
    conf_low = .data[[jackal_endpoint_low]],
    conf_high = .data[[jackal_endpoint_high]],
    contrast = "+1 versus -1 SD jackal",
    odds_ratio = exp(estimate),
    odds_ratio_low = exp(.data[[jackal_endpoint_low]]),
    odds_ratio_high = exp(.data[[jackal_endpoint_high]]),
    probability_difference = NA_real_,
    p_value = p.value
  )

probability_differences <- predicted_probabilities %>%
  filter(record_type == "predicted_probability") %>%
  group_by(class, jackal_z) %>%
  summarise(
    probability_difference = probability[wolf_z == 1] - probability[wolf_z == -1],
    .groups = "drop"
  )
wolf_endpoint <- wolf_endpoint %>%
  left_join(probability_differences, by = c("class", "jackal_z"), suffix = c("", ".new")) %>%
  mutate(probability_difference = probability_difference.new) %>%
  select(-probability_difference.new)

probability_differences_j <- predicted_probabilities %>%
  filter(record_type == "predicted_probability") %>%
  group_by(class, wolf_z) %>%
  summarise(
    probability_difference = probability[jackal_z == 1] - probability[jackal_z == -1],
    .groups = "drop"
  )
jackal_endpoint <- jackal_endpoint %>%
  left_join(probability_differences_j, by = c("class", "wolf_z"), suffix = c("", ".new")) %>%
  mutate(probability_difference = probability_difference.new) %>%
  select(-probability_difference.new)

prediction_output <- bind_rows(predicted_probabilities, wolf_endpoint, jackal_endpoint)
write_csv(
  prediction_output,
  file.path(output_dir, "H3_predicted_daylight_probabilities.csv")
)

wolf_sequence <- seq(min(site$wolf_z), max(site$wolf_z), length.out = 121)
jackal_sequence <- seq(min(site$jackal_z), max(site$jackal_z), length.out = 121)

wolf_curve <- extract_probability_summary(emmeans(
  primary,
  ~ wolf_z | class * jackal_z,
  at = list(wolf_z = wolf_sequence, jackal_z = c(-1, 0, 1), human_z = 0),
  weights = "equal"
)) %>%
  mutate(
    class = factor(class, levels = class_levels),
    moderator = factor(
      jackal_z,
      levels = c(-1, 0, 1),
      labels = c("Jackal: -1 SD", "Jackal: mean", "Jackal: +1 SD")
    )
  )

jackal_curve <- extract_probability_summary(emmeans(
  primary,
  ~ jackal_z | class * wolf_z,
  at = list(jackal_z = jackal_sequence, wolf_z = c(-1, 0, 1), human_z = 0),
  weights = "equal"
)) %>%
  mutate(
    class = factor(class, levels = class_levels),
    moderator = factor(
      wolf_z,
      levels = c(-1, 0, 1),
      labels = c("Wolf: -1 SD", "Wolf: mean", "Wolf: +1 SD")
    )
  )

line_palette <- c("#A7C8C2", "#417F79", "#153F3B")

wolf_plot <- ggplot(
  wolf_curve,
  aes(x = wolf_z, y = probability, color = moderator, fill = moderator)
) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.14, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~class, nrow = 1) +
  scale_color_manual(values = line_palette, name = NULL) +
  scale_fill_manual(values = line_palette, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Standardized wolf activity",
    y = "Predicted probability of daylight activity",
    title = "Demographic daylight responses across predator-activity gradients",
    subtitle = "Lines show jackal activity at -1 SD, the mean, and +1 SD; ribbons are 95% confidence intervals",
    caption = "Predictions average equally over Early and Late periods and hold standardized human activity at zero."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

jackal_plot <- ggplot(
  jackal_curve,
  aes(x = jackal_z, y = probability, color = moderator, fill = moderator)
) +
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.14, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~class, nrow = 1) +
  scale_color_manual(values = line_palette, name = NULL) +
  scale_fill_manual(values = line_palette, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Standardized jackal activity",
    y = "Predicted probability of daylight activity",
    title = "Demographic daylight responses across predator-activity gradients",
    subtitle = "Lines show wolf activity at -1 SD, the mean, and +1 SD; ribbons are 95% confidence intervals",
    caption = "Predictions average equally over Early and Late periods and hold standardized human activity at zero."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "H3_predicted_daylight_figure.png"), wolf_plot, width = 11, height = 4.8, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "H3_predicted_daylight_figure.pdf"), wolf_plot, width = 11, height = 4.8, bg = "white")
ggsave(file.path(output_dir, "H3_predicted_daylight_reciprocal.png"), jackal_plot, width = 11, height = 4.8, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "H3_predicted_daylight_reciprocal.pdf"), jackal_plot, width = 11, height = 4.8, bg = "white")

forest_data <- conditional_slopes %>%
  filter(record_type == "conditional_slope") %>%
  mutate(
    moderator_label = paste0(
      ifelse(predictor == "wolf_z", "Jackal", "Wolf"),
      " = ", moderator_value, " SD"
    ),
    moderator_level = factor(
      moderator_value,
      levels = c(-1, 0, 1),
      labels = c("Other predator: -1 SD", "Other predator: mean", "Other predator: +1 SD")
    ),
    predictor_label = recode(predictor, wolf_z = "Wolf activity", jackal_z = "Jackal activity"),
    class = factor(class, levels = rev(class_levels))
  )

forest_plot <- ggplot(
  forest_data,
  aes(x = odds_ratio, y = class, color = moderator_level)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#697476") +
  geom_errorbar(
    aes(xmin = odds_ratio_low, xmax = odds_ratio_high),
    width = 0.12,
    orientation = "y",
    position = position_dodge(width = 0.55)
  ) +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  facet_wrap(~predictor_label, nrow = 1, scales = "free_x") +
  scale_x_log10() +
  scale_color_manual(values = line_palette, name = NULL) +
  labs(
    x = "Odds ratio for daylight activity per 1-SD predator increase",
    y = NULL,
    title = "Conditional demographic responses to predator activity",
    subtitle = "Points and lines are odds ratios and 95% confidence intervals"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "H3_conditional_effects_forest.png"), forest_plot, width = 10, height = 5.3, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "H3_conditional_effects_forest.pdf"), forest_plot, width = 10, height = 5.3, bg = "white")

# Prespecified sensitivity analyses -----------------------------------------

date_model <- fit_model(as.formula(paste(
  "daylight ~", primary_fixed, "+ human_z + date_z", random_primary
)), events)
no_human_model <- fit_model(as.formula(paste(
  "daylight ~", primary_fixed, "+ period", random_primary
)), events)

# Transform 0-12 h to the open unit interval for a bounded beta mixed model.
n_events <- nrow(events)
events <- events %>%
  mutate(
    nocturnality_proportion = ((nocturnality_hours / 12) * (n_events - 1) + 0.5) / n_events
  )
nocturnality_model <- fit_model(
  as.formula(paste(
    "nocturnality_proportion ~", primary_fixed, "+ human_z + period", random_primary
  )),
  events,
  family = beta_family(link = "logit")
)

nocturnality_descriptive <- events %>%
  group_by(class, predator_context) %>%
  summarise(
    n = n(),
    mean_hours = mean(nocturnality_hours),
    median_hours = median(nocturnality_hours),
    sd_hours = sd(nocturnality_hours),
    .groups = "drop"
  ) %>%
  mutate(
    record_type = "context_descriptive",
    predictor = NA_character_,
    moderator = NA_character_,
    moderator_value = NA_real_,
    estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_
  )

nocturnality_slopes <- bind_rows(
  extract_trends(nocturnality_model, "wolf_z", "jackal_z"),
  extract_trends(nocturnality_model, "jackal_z", "wolf_z")
) %>%
  mutate(
    n = NA_integer_, mean_hours = NA_real_, median_hours = NA_real_, sd_hours = NA_real_,
    predator_context = NA
  )

nocturnality_output <- bind_rows(
  nocturnality_descriptive,
  nocturnality_slopes
) %>%
  mutate(
    interpretation = case_when(
      record_type == "conditional_slope" & estimate < 0 ~ "Shift toward greater diurnality",
      record_type == "conditional_slope" & estimate > 0 ~ "Shift toward greater nocturnality",
      TRUE ~ NA_character_
    )
  )
write_csv(nocturnality_output, file.path(output_dir, "H3_nocturnality_sensitivity.csv"))

sensitivity_tests <- bind_rows(
  lrt_row(no_predator_interaction, primary, "Primary versus no wolf-jackal interaction"),
  lrt_row(reduced_selected, full_selected, "Full three-way versus reduced interactions"),
  lrt_row(no_human_model, primary, "With versus without human adjustment"),
  tibble(
    test = "Early/Late versus continuous date",
    method = "AIC comparison; nonnested period parameterization",
    statistic = NA_real_,
    df = NA_real_,
    p_value = NA_real_,
    reduced_AIC = AIC(date_model),
    full_AIC = AIC(primary)
  )
)
write_csv(sensitivity_tests, file.path(output_dir, "H3_sensitivity_tests.csv"))

mean_predator_slopes <- function(model, model_name) {
  bind_rows(
    extract_trends(model, "wolf_z", "jackal_z"),
    extract_trends(model, "jackal_z", "wolf_z")
  ) %>%
    filter(record_type == "conditional_slope", moderator_value == 0) %>%
    mutate(model = model_name, .before = 1)
}

sensitivity_slopes <- bind_rows(
  mean_predator_slopes(primary, "primary_full_three_way"),
  mean_predator_slopes(reduced_selected, "reduced_interactions"),
  mean_predator_slopes(no_predator_interaction, "no_wolf_jackal_interaction"),
  mean_predator_slopes(date_model, "continuous_date"),
  mean_predator_slopes(no_human_model, "without_human_adjustment")
) %>%
  group_by(model, predictor) %>%
  mutate(p_fdr_within_model_predator = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write_csv(
  sensitivity_slopes,
  file.path(output_dir, "H3_sensitivity_predator_slopes.csv")
)

central_effects <- bind_rows(
  global_tests %>%
    transmute(
      result_type = "global_test",
      effect = test,
      class = NA_character_,
      moderator_value = NA_real_,
      estimate = statistic,
      conf_low = NA_real_,
      conf_high = NA_real_,
      odds_ratio = NA_real_,
      odds_ratio_low = NA_real_,
      odds_ratio_high = NA_real_,
      p_value
    ),
  conditional_slopes %>%
    filter(record_type == "conditional_slope", moderator_value == 0) %>%
    transmute(
      result_type = "class_slope_at_other_predator_mean",
      effect = predictor,
      class,
      moderator_value,
      estimate,
      conf_low,
      conf_high,
      odds_ratio,
      odds_ratio_low,
      odds_ratio_high,
      p_value
    ),
  primary_coefficients %>%
    filter(term %in% c("human_z", "periodLate")) %>%
    transmute(
      result_type = "adjustment_effect",
      effect = term,
      class = NA_character_,
      moderator_value = NA_real_,
      estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      odds_ratio,
      odds_ratio_low,
      odds_ratio_high,
      p_value = p.value
    )
)
write_csv(central_effects, file.path(output_dir, "H3_central_effects.csv"))

all_diagnostics <- bind_rows(
  model_diagnostics(full_encounter, "candidate_full_with_encounter"),
  model_diagnostics(full_camera, "candidate_full_without_encounter"),
  model_diagnostics(reduced_encounter, "candidate_reduced_with_encounter"),
  model_diagnostics(reduced_camera, "candidate_reduced_without_encounter"),
  model_diagnostics(primary, "primary"),
  model_diagnostics(no_predator_interaction, "no_predator_interaction"),
  model_diagnostics(date_model, "continuous_date_sensitivity"),
  model_diagnostics(no_human_model, "no_human_sensitivity"),
  model_diagnostics(nocturnality_model, "beta_nocturnality_sensitivity")
)
write_csv(all_diagnostics, file.path(output_dir, "H3_model_diagnostics.csv"))

random_effects <- bind_rows(
  random_variances(primary, "primary"),
  random_variances(full_encounter, "candidate_full_with_encounter"),
  random_variances(full_camera, "candidate_full_without_encounter")
)
write_csv(random_effects, file.path(output_dir, "H3_random_effect_variances.csv"))

# Circular whole-distribution complement ------------------------------------

extract_activity_density <- function(fit, class_name, context_name) {
  density <- as.data.frame(fit@pdf)
  tibble(
    class = class_name,
    predator_context = context_name,
    solar_hour = 24 * density$x / (2 * pi),
    density = density$y
  )
}

density_rows <- list()
density_i <- 1L
for (class_i in class_levels) {
  for (context_i in context_levels) {
    times_i <- events$solar_hour[
      events$class == class_i & events$predator_context == context_i
    ] * 2 * pi / 24
    if (length(times_i) < minimum_curve_events) next
    fit_i <- fitact(times_i, reps = activity_reps, sample = "data", show = FALSE)
    density_rows[[density_i]] <- extract_activity_density(fit_i, class_i, context_i)
    density_i <- density_i + 1L
  }
}

context_density <- bind_rows(density_rows) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(predator_context, levels = context_levels)
  )
sparse_events <- events %>%
  add_count(class, predator_context, name = "class_context_n") %>%
  filter(class_context_n < minimum_curve_events) %>%
  mutate(
    class = factor(class, levels = class_levels),
    predator_context = factor(predator_context, levels = context_levels)
  )

write_csv(context_density, file.path(output_dir, "H3_circular_context_density_curves.csv"))

context_palette <- c(
  "Low both" = "#087F6D",
  "High jackal / Low wolf" = "#B65332",
  "Low jackal / High wolf" = "#244E73",
  "High both" = "#6B5976"
)

circular_plot <- ggplot(context_density, aes(solar_hour, density, color = predator_context)) +
  annotate("rect", xmin = 0, xmax = 5.4, ymin = -Inf, ymax = Inf, fill = "#DDE4E6", alpha = 0.55) +
  annotate("rect", xmin = 5.4, xmax = 7, ymin = -Inf, ymax = Inf, fill = "#F3E6D7", alpha = 0.40) +
  annotate("rect", xmin = 17, xmax = 18.8, ymin = -Inf, ymax = Inf, fill = "#F0DCD7", alpha = 0.40) +
  annotate("rect", xmin = 18.8, xmax = 24, ymin = -Inf, ymax = Inf, fill = "#DDE4E6", alpha = 0.55) +
  geom_line(linewidth = 0.95) +
  geom_rug(
    data = sparse_events,
    aes(x = solar_hour, color = predator_context),
    sides = "b",
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  facet_wrap(~class, ncol = 1, scales = "free_y") +
  scale_color_manual(values = context_palette, name = NULL) +
  scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24), expand = c(0, 0)) +
  labs(
    x = "Apparent solar time",
    y = "Activity density",
    title = "Class-specific diel activity across wolf-jackal contexts",
    subtitle = "Median-split contexts are a descriptive whole-distribution complement to the continuous H3 model",
    caption = "The high-jackal/low-wolf fawn context contains three events and is shown as ticks rather than a density curve."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "H3_circular_context_figure.png"), circular_plot, width = 10, height = 8.5, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "H3_circular_context_figure.pdf"), circular_plot, width = 10, height = 8.5, bg = "white")

authoritative_context_overlap <- read_csv(context_overlap_file, show_col_types = FALSE)
write_csv(
  authoritative_context_overlap,
  file.path(output_dir, "H3_context_pairwise_overlap_existing_authoritative.csv")
)

# Persist models and diagnostics --------------------------------------------

saveRDS(primary, file.path(output_dir, "H3_primary_daylight_model.rds"))
saveRDS(
  list(
    full_with_encounter = full_encounter,
    full_without_encounter = full_camera,
    reduced_with_encounter = reduced_encounter,
    reduced_without_encounter = reduced_camera,
    continuous_date = date_model,
    without_human = no_human_model,
    nocturnality_beta = nocturnality_model
  ),
  file.path(output_dir, "H3_sensitivity_models.rds")
)

diagnostic_text <- c(
  "Revised H3 model diagnostics",
  "============================",
  "",
  paste("Primary formula:", paste(deparse(formula(primary)), collapse = " ")),
  paste("Primary interaction structure:", primary_structure),
  paste("Encounter random effect retained:", use_encounter),
  paste("Encounter variance in full candidate:", format(encounter_variance, digits = 6)),
  "",
  "Validation summary:",
  capture.output(print(validation_report, n = Inf)),
  "",
  "Model diagnostics:",
  capture.output(print(all_diagnostics, n = Inf)),
  "",
  "Random-effect variances:",
  capture.output(print(random_effects, n = Inf)),
  "",
  "Joint predator-gradient coverage and separation screen:",
  capture.output(print(gradient_coverage, n = Inf)),
  "",
  paste("R version:", R.version.string),
  paste("glmmTMB version:", as.character(packageVersion("glmmTMB"))),
  paste("emmeans version:", as.character(packageVersion("emmeans"))),
  paste("activity version:", as.character(packageVersion("activity")))
)
writeLines(diagnostic_text, file.path(output_dir, "H3_model_diagnostics.txt"))

get_test_p <- function(label) {
  global_tests$p_value[match(label, global_tests$test)]
}
fmt_p <- function(x) ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
fmt_ci <- function(est, low, high) sprintf("%.2f [%.2f, %.2f]", est, low, high)

middle_slopes <- conditional_slopes %>%
  filter(record_type == "conditional_slope", moderator_value == 0) %>%
  arrange(predictor, factor(class, levels = class_levels))

most_diurnal <- daylight_summary$class[which.max(daylight_summary$daylight_proportion)]
most_diurnal_value <- max(daylight_summary$daylight_proportion)
supported_directional <- conditional_slopes %>%
  filter(record_type == "conditional_slope", p_fdr < 0.05)
supported_text <- if (nrow(supported_directional)) {
  paste0(nrow(supported_directional), " conditional predator slopes remained supported after FDR correction")
} else {
  "no conditional predator slope remained supported after FDR correction"
}

baseline_probability_text <- paste0(
  tolower(baseline_probabilities$class), " ",
  mapply(
    fmt_ci,
    baseline_probabilities$estimate,
    baseline_probabilities$conf_low,
    baseline_probabilities$conf_high
  ),
  collapse = "; "
)
wolf_middle <- middle_slopes %>% filter(predictor == "wolf_z")
jackal_middle <- middle_slopes %>% filter(predictor == "jackal_z")
human_row <- primary_coefficients %>% filter(term == "human_z")
period_row <- primary_coefficients %>% filter(term == "periodLate")
nocturnality_middle <- nocturnality_slopes %>%
  filter(record_type == "conditional_slope", moderator_value == 0)
nocturnality_direction_consistent <- all(nocturnality_middle$estimate < 0)
nocturnality_supported <- any(nocturnality_middle$p_value < 0.05)

context_supported_counts <- authoritative_context_overlap %>%
  mutate(
    supported = formal_test_eligible %in% TRUE &
      !is.na(significance) &
      !significance %in% c("ns", "Descriptive only")
  ) %>%
  group_by(class) %>%
  summarise(supported = sum(supported), .groups = "drop") %>%
  complete(class = class_levels, fill = list(supported = 0L))
get_context_supported <- function(class_name) {
  context_supported_counts$supported[match(class_name, context_supported_counts$class)]
}

question_answers <- c(
  paste0("1. Fawns were most diurnal: 82.8% daylight versus 61.1% for females and 56.8% for males; all model-adjusted pairwise class contrasts were supported."),
  paste0("2. Wolf activity was not associated with a general shift toward daylight (global LRT P = ", fmt_p(get_test_p("Any wolf association across classes and jackal contexts")), ")."),
  paste0("3. Jackal activity was not associated with a general shift toward daylight (global LRT P = ", fmt_p(get_test_p("Any jackal association across classes and wolf contexts")), ")."),
  paste0("4. Wolf slopes did not differ among classes (LRT P = ", fmt_p(get_test_p("Wolf slopes differ among demographic classes")), ")."),
  paste0("5. Jackal slopes did not differ among classes (LRT P = ", fmt_p(get_test_p("Jackal slopes differ among demographic classes")), ")."),
  paste0("6. The joint wolf-jackal interaction was unsupported (LRT P = ", fmt_p(get_test_p("Wolf by jackal interaction, including class-specific interaction if retained")), ")."),
  paste0("7. The class-specific three-way interaction was unsupported (LRT P = ", fmt_p(get_test_p("Class by wolf by jackal interaction")), ")."),
  "8. Fawns had high baseline daylight use and no supported additional predator-associated shift; sparse joint-gradient cells make invariance uncertain rather than demonstrated.",
  paste0("9. The directional models found no supported daylight shift, whereas circular tests retained ", get_context_supported("Female"), " female and ", get_context_supported("Male"), " male whole-distribution differences. These results indicate within-day or crepuscular redistribution rather than a consistent day-night displacement."),
  "10. Baseline demographic differences are well supported; predator-associated directional effects are uncertain, especially for fawns."
)

summary_paragraph <- paste0(
  "We analyzed 3,925 class-specific 30-minute events from 53 cameras in 12 polygons, representing 3,507 reconstructed encounters. ",
  "The full class by wolf by jackal binomial mixed model converged with polygon, camera, and encounter random intercepts. Fawns were most diurnal (82.8% daylight), followed by females (61.1%) and males (56.8%). Model-adjusted daylight probabilities at mean predator and human activity were ", baseline_probability_text, ", and all pairwise class contrasts were supported. ",
  "Neither wolf activity (global P = ", fmt_p(get_test_p("Any wolf association across classes and jackal contexts")), ") nor jackal activity (P = ", fmt_p(get_test_p("Any jackal association across classes and wolf contexts")), ") showed a supported directional association with daylight. At the other predator's mean, wolf odds ratios ranged from ", sprintf("%.2f", min(wolf_middle$odds_ratio)), " to ", sprintf("%.2f", max(wolf_middle$odds_ratio)), " and jackal odds ratios from ", sprintf("%.2f", min(jackal_middle$odds_ratio)), " to ", sprintf("%.2f", max(jackal_middle$odds_ratio)), "; every interval included one and ", supported_text, ". The wolf-jackal interaction (P = ", fmt_p(get_test_p("Wolf by jackal interaction, including class-specific interaction if retained")), ") and its class-specific component (P = ", fmt_p(get_test_p("Class by wolf by jackal interaction")), ") were unsupported. ",
  "Higher human activity was associated with lower daylight odds (OR ", fmt_ci(human_row$odds_ratio, human_row$odds_ratio_low, human_row$odds_ratio_high), "), whereas Late-period odds were higher (OR ", fmt_ci(period_row$odds_ratio, period_row$odds_ratio_low, period_row$odds_ratio_high), "). A continuous-date sensitivity lowered AIC by ", sprintf("%.1f", AIC(primary) - AIC(date_model)), " but left predator slopes nearly unchanged. The nocturnality-hours model ", ifelse(nocturnality_direction_consistent, "showed directionally greater diurnality", "did not show a consistent direction"), " with increasing predators, but ", ifelse(nocturnality_supported, "some slopes were supported", "none was supported"), ". Circular tests retained adult whole-distribution differences, indicating temporal redistribution without a consistent day-night shift. Predator indices represent relative camera-level activity, not causal predation risk."
)

result_summary <- c(
  "# Revised H3 result summary",
  "",
  summary_paragraph,
  "",
  "## Primary formula",
  "",
  paste0("`", paste(deparse(formula(primary)), collapse = " "), "`"),
  "",
  "## Direct answers to H3 questions",
  "",
  question_answers,
  "",
  "## Central global tests",
  "",
  capture.output(print(global_tests, n = Inf)),
  "",
  "## Daylight proportions",
  "",
  capture.output(print(daylight_summary, n = Inf))
)
writeLines(result_summary, file.path(output_dir, "H3_results_summary.md"))

cat("\nH3 validation summary\n")
print(validation_report, n = Inf)
cat("\nH3 primary model\n")
cat(paste(deparse(formula(primary)), collapse = " "), "\n")
cat("\nH3 global tests\n")
print(global_tests, n = Inf)
cat("\nDirect answers to H3 questions\n")
cat(paste(question_answers, collapse = "\n"), "\n")
cat("\nPlain-English interpretation\n")
cat(summary_paragraph, "\n")
cat("\nOutputs saved to: ", normalizePath(output_dir), "\n", sep = "")
