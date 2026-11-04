#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(tidyr)
})

set.seed(20260810)

site_file <- "SiteCovs_060324.csv"
obs_file <- "ObsCovs_050524.csv"
event_file <- "Gazella gazella_independent_events_30min.csv"
count_files <- c(
  Fawn = "gazelle_juvenile_abundance_60x7.csv",
  Female = "gazelle_female_abundance_60x7.csv",
  Male = "gazelle_male_abundance_60x7.csv"
)
prefix <- "exploratory_wolf_jackal_ratio"
local_timezone <- "Asia/Jerusalem"
bootstrap_reps <- 2000L
minimum_activity_events <- 5L

required_files <- c(site_file, obs_file, event_file, unname(count_files))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

z_score <- function(x) as.numeric(scale(x))

parse_local_camera_time <- function(x) {
  as.POSIXct(
    sub("Z$", "", x),
    format = "%Y-%m-%dT%H:%M:%S",
    tz = local_timezone
  )
}

# NOAA-style approximation, matching the unified spatiotemporal workflow.
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
  decl <-
    0.006918 -
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

site_raw <- read_csv(site_file, show_col_types = FALSE)
if (any(site_raw$Null.Wolf.Abundance <= 0, na.rm = TRUE) ||
    any(site_raw$null.jackal.abundance <= 0, na.rm = TRUE)) {
  stop("Wolf and jackal indices must both be positive for an unadjusted ratio.")
}

site <- site_raw %>%
  mutate(
    camera_name = as.integer(Camera.ID),
    wolf_jackal_ratio = Null.Wolf.Abundance / null.jackal.abundance,
    log_wolf_jackal_ratio = log(wolf_jackal_ratio),
    wolf_share = Null.Wolf.Abundance /
      (Null.Wolf.Abundance + null.jackal.abundance),
    predator_total_index = Null.Wolf.Abundance + null.jackal.abundance,
    log_wolf_jackal_ratio_z = z_score(log_wolf_jackal_ratio),
    wolf_share_z = z_score(wolf_share),
    predator_total_index_z = z_score(predator_total_index),
    open_spaces_z = z_score(Open.Spaces)
  )

site_enriched <- site %>%
  select(
    -camera_name, -log_wolf_jackal_ratio_z, -wolf_share_z,
    -predator_total_index_z, -open_spaces_z
  )
write_csv(site_enriched, "SiteCovs_060324_with_wolf_jackal_ratio.csv")

obs <- read.csv(
  obs_file,
  na.strings = c("NA", ""),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(
    camera_name = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    days_operational = as.numeric(`Days.Operational`)
  ) %>%
  filter(!is.na(camera_name), !is.na(round))

read_count_history <- function(file, class_name) {
  read_csv(file, na = c("NA", ""), show_col_types = FALSE) %>%
    rename(camera_name = camera_name) %>%
    pivot_longer(
      cols = matches("^X[1-7]$"),
      names_to = "round",
      values_to = "count"
    ) %>%
    mutate(
      camera_name = as.integer(camera_name),
      round = as.integer(sub("^X", "", round)),
      class = class_name,
      count = as.numeric(count)
    )
}

count_long <- bind_rows(lapply(names(count_files), function(class_i) {
  read_count_history(count_files[[class_i]], class_i)
})) %>%
  left_join(
    obs %>% select(camera_name, round, days_operational),
    by = c("camera_name", "round")
  )

class_rates <- count_long %>%
  filter(!is.na(count), is.finite(days_operational), days_operational > 0) %>%
  group_by(camera_name, class) %>%
  summarise(
    observed_count = sum(count),
    trap_days = sum(days_operational),
    count_rate_100_days = 100 * observed_count / trap_days,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = class,
    values_from = c(observed_count, trap_days, count_rate_100_days),
    names_glue = "{tolower(class)}_{.value}"
  )

camera_round_counts <- count_long %>%
  select(camera_name, round, class, count, days_operational) %>%
  pivot_wider(names_from = class, values_from = count) %>%
  filter(if_any(c(Fawn, Female, Male), ~ !is.na(.x))) %>%
  mutate(
    identified_count = rowSums(across(c(Fawn, Female, Male)), na.rm = TRUE)
  )

camera_counts <- camera_round_counts %>%
  group_by(camera_name) %>%
  summarise(
    total_identified_count = sum(identified_count),
    total_trap_days = sum(days_operational[is.finite(days_operational) & days_operational > 0]),
    observed_gazelle_rate_100_days = 100 * total_identified_count / total_trap_days,
    fawn_count = sum(Fawn, na.rm = TRUE),
    female_count = sum(Female, na.rm = TRUE),
    male_count = sum(Male, na.rm = TRUE),
    fawn_female_total = fawn_count + female_count,
    fawn_composition = if_else(
      fawn_female_total > 0,
      fawn_count / fawn_female_total,
      NA_real_
    ),
    fawns_per_100_females = if_else(
      female_count > 0,
      100 * fawn_count / female_count,
      NA_real_
    ),
    fawns_per_100_females_smoothed =
      100 * (fawn_count + 0.5) / (female_count + 0.5),
    .groups = "drop"
  )

events_raw <- read_csv(
  event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
) %>%
  mutate(
    event_row_id = row_number(),
    camera_name = as.integer(camera_name),
    class = recode(
      gazelle_class,
      juvenile = "Fawn",
      female = "Female",
      male = "Male"
    ),
    timestamp = parse_local_camera_time(event_time),
    event_date = as.Date(format(timestamp, tz = local_timezone)),
    solar_elevation = solar_elevation(
      timestamp,
      camelot_gps_latitude,
      camelot_gps_longitude
    ),
    daylight = solar_elevation >= 0
  )

active_event_ids <- events_raw %>%
  select(event_row_id, camera_name, event_date) %>%
  inner_join(
    obs %>% select(camera_name, round, start_date, end_date),
    by = "camera_name",
    relationship = "many-to-many"
  ) %>%
  filter(
    !is.na(start_date), !is.na(end_date),
    event_date >= start_date, event_date <= end_date
  ) %>%
  group_by(event_row_id) %>%
  slice_max(start_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(event_row_id)

events <- events_raw %>%
  inner_join(active_event_ids, by = "event_row_id") %>%
  filter(class %in% c("Fawn", "Female", "Male"), !is.na(daylight))

activity_all <- events %>%
  group_by(camera_name) %>%
  summarise(
    gazelle_daylight_events = sum(daylight),
    gazelle_dark_events = sum(!daylight),
    gazelle_activity_events = n(),
    gazelle_daylight_fraction = gazelle_daylight_events / gazelle_activity_events,
    gazelle_day_night_odds_smoothed =
      (gazelle_daylight_events + 0.5) / (gazelle_dark_events + 0.5),
    .groups = "drop"
  )

activity_class <- events %>%
  group_by(camera_name, class) %>%
  summarise(
    daylight_events = sum(daylight),
    dark_events = sum(!daylight),
    activity_events = n(),
    daylight_fraction = daylight_events / activity_events,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = class,
    values_from = c(
      daylight_events, dark_events, activity_events, daylight_fraction
    ),
    names_glue = "{tolower(class)}_{.value}"
  )

camera_data <- site %>%
  left_join(camera_counts, by = "camera_name") %>%
  left_join(class_rates, by = "camera_name") %>%
  left_join(activity_all, by = "camera_name") %>%
  left_join(activity_class, by = "camera_name") %>%
  mutate(
    gazelle_abundance_log_z = z_score(log(gazelleabundance)),
    observed_gazelle_rate_log_z = z_score(log1p(observed_gazelle_rate_100_days))
  )

write_csv(camera_data, paste0(prefix, "_camera_data.csv"))

bootstrap_spearman <- function(x, y, reps = bootstrap_reps) {
  n <- length(x)
  values <- replicate(reps, {
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  values <- values[is.finite(values)]
  if (length(values) < ceiling(0.8 * reps)) return(c(NA_real_, NA_real_))
  unname(quantile(values, c(0.025, 0.975)))
}

cor_specs <- tribble(
  ~scope, ~domain, ~variable, ~label, ~minimum_events_variable,
  "Primary", "Habitat", "Open.Spaces", "Open-space cover", NA_character_,
  "Primary", "Abundance", "gazelleabundance", "Gazelle abundance index", NA_character_,
  "Primary", "Composition", "fawn_composition", "Fawn share of fawn + female counts", "fawn_female_total",
  "Primary", "Activity", "gazelle_daylight_fraction", "All-gazelle daylight fraction", "gazelle_activity_events",
  "Secondary", "Abundance", "observed_gazelle_rate_100_days", "Observed identified count rate", NA_character_,
  "Secondary", "Abundance", "fawn_count_rate_100_days", "Fawn count rate", NA_character_,
  "Secondary", "Abundance", "female_count_rate_100_days", "Female count rate", NA_character_,
  "Secondary", "Abundance", "male_count_rate_100_days", "Male count rate", NA_character_,
  "Secondary", "Activity", "fawn_daylight_fraction", "Fawn daylight fraction", "fawn_activity_events",
  "Secondary", "Activity", "female_daylight_fraction", "Female daylight fraction", "female_activity_events",
  "Secondary", "Activity", "male_daylight_fraction", "Male daylight fraction", "male_activity_events",
  "Secondary", "Other pressures", "human.abundance", "Human abundance", NA_character_,
  "Secondary", "Other pressures", "null.boar.abundance", "Boar abundance", NA_character_,
  "Secondary", "Habitat", "Built.and.Roads", "Built and roads", NA_character_,
  "Secondary", "Habitat", "Nature.Reserves", "Nature reserves", NA_character_,
  "Diagnostic", "Ratio construction", "Null.Wolf.Abundance", "Wolf abundance", NA_character_,
  "Diagnostic", "Ratio construction", "null.jackal.abundance", "Jackal abundance", NA_character_
)

run_correlation <- function(spec) {
  data_i <- camera_data %>%
    filter(
      is.finite(log_wolf_jackal_ratio),
      is.finite(.data[[spec$variable]])
    )
  if (!is.na(spec$minimum_events_variable)) {
    threshold <- if (spec$domain == "Composition") 1L else minimum_activity_events
    data_i <- data_i %>%
      filter(.data[[spec$minimum_events_variable]] >= threshold)
  }
  test <- suppressWarnings(cor.test(
    data_i$log_wolf_jackal_ratio,
    data_i[[spec$variable]],
    method = "spearman",
    exact = FALSE
  ))
  interval <- bootstrap_spearman(
    data_i$log_wolf_jackal_ratio,
    data_i[[spec$variable]]
  )
  tibble(
    scope = spec$scope,
    domain = spec$domain,
    variable = spec$variable,
    outcome = spec$label,
    n_cameras = nrow(data_i),
    spearman_rho = unname(test$estimate),
    conf_low = interval[1],
    conf_high = interval[2],
    p_value = test$p.value
  )
}

correlations <- bind_rows(lapply(seq_len(nrow(cor_specs)), function(i) {
  run_correlation(cor_specs[i, ])
})) %>%
  group_by(scope) %>%
  mutate(
    p_fdr = if_else(scope == "Diagnostic", NA_real_, p.adjust(p_value, "BH"))
  ) %>%
  ungroup()

write_csv(correlations, paste0(prefix, "_correlations.csv"))

partial_spearman <- function(data, outcome, control) {
  data_i <- data %>%
    filter(
      is.finite(log_wolf_jackal_ratio),
      is.finite(.data[[outcome]]),
      is.finite(.data[[control]])
    )
  ratio_rank <- rank(data_i$log_wolf_jackal_ratio)
  outcome_rank <- rank(data_i[[outcome]])
  control_rank <- rank(data_i[[control]])
  ratio_residual <- residuals(lm(ratio_rank ~ control_rank))
  outcome_residual <- residuals(lm(outcome_rank ~ control_rank))
  rho <- cor(ratio_residual, outcome_residual)
  degrees_freedom <- nrow(data_i) - 3L
  statistic <- rho * sqrt(degrees_freedom / (1 - rho^2))
  tibble(
    check = paste0("Partial Spearman controlling ", control),
    outcome,
    excluded_camera = NA_integer_,
    n_cameras = nrow(data_i),
    spearman_rho = rho,
    p_value = 2 * pt(-abs(statistic), df = degrees_freedom)
  )
}

extreme_ratio_camera <- camera_data$camera_name[
  which.min(camera_data$wolf_jackal_ratio)
]

leave_one_out_spearman <- function(data, outcome, excluded_camera) {
  data_i <- data %>%
    filter(
      camera_name != excluded_camera,
      is.finite(log_wolf_jackal_ratio),
      is.finite(.data[[outcome]])
    )
  test <- suppressWarnings(cor.test(
    data_i$log_wolf_jackal_ratio,
    data_i[[outcome]],
    method = "spearman",
    exact = FALSE
  ))
  tibble(
    check = "Spearman excluding minimum-ratio camera",
    outcome,
    excluded_camera,
    n_cameras = nrow(data_i),
    spearman_rho = unname(test$estimate),
    p_value = test$p.value
  )
}

robustness_checks <- bind_rows(
  partial_spearman(camera_data, "gazelleabundance", "Open.Spaces"),
  partial_spearman(
    camera_data, "observed_gazelle_rate_100_days", "Open.Spaces"
  ),
  leave_one_out_spearman(
    camera_data, "gazelleabundance", extreme_ratio_camera
  ),
  leave_one_out_spearman(
    camera_data, "observed_gazelle_rate_100_days", extreme_ratio_camera
  )
) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH"))

write_csv(robustness_checks, paste0(prefix, "_robustness_checks.csv"))

extract_model_term <- function(
    model, model_name, term, scale_description, exponentiate = FALSE) {
  coefs <- summary(model)$coefficients
  estimate <- coefs[term, "Estimate"]
  std_error <- coefs[term, "Std. Error"]
  p_col <- grep("Pr\\(", colnames(coefs), value = TRUE)[1]
  p_value <- coefs[term, p_col]
  n_model <- nrow(stats::model.frame(model))
  tibble(
    model = model_name,
    term = term,
    scale = scale_description,
    estimate,
    std_error,
    conf_low = estimate - 1.96 * std_error,
    conf_high = estimate + 1.96 * std_error,
    odds_or_multiplicative_effect = if (exponentiate) exp(estimate) else NA_real_,
    effect_conf_low = if (exponentiate) exp(estimate - 1.96 * std_error) else NA_real_,
    effect_conf_high = if (exponentiate) exp(estimate + 1.96 * std_error) else NA_real_,
    p_value,
    n_cameras = n_model
  )
}

composition_data <- camera_data %>% filter(fawn_female_total > 0)
composition_model <- glm(
  cbind(fawn_count, female_count) ~ log_wolf_jackal_ratio_z,
  family = quasibinomial,
  data = composition_data
)

activity_data <- camera_data %>%
  filter(gazelle_activity_events >= minimum_activity_events)
activity_model <- glm(
  cbind(gazelle_daylight_events, gazelle_dark_events) ~
    log_wolf_jackal_ratio_z,
  family = quasibinomial,
  data = activity_data
)

abundance_model <- lm(
  gazelle_abundance_log_z ~ log_wolf_jackal_ratio_z + open_spaces_z,
  data = camera_data
)
observed_rate_model <- lm(
  observed_gazelle_rate_log_z ~ log_wolf_jackal_ratio_z + open_spaces_z,
  data = camera_data
)
abundance_share_model <- lm(
  gazelle_abundance_log_z ~ wolf_share_z + open_spaces_z,
  data = camera_data
)
observed_rate_share_model <- lm(
  observed_gazelle_rate_log_z ~ wolf_share_z + open_spaces_z,
  data = camera_data
)

model_results <- bind_rows(
  extract_model_term(
    composition_model,
    "Quasibinomial fawn composition",
    "log_wolf_jackal_ratio_z",
    "Odds ratio per 1 SD increase in log wolf:jackal ratio",
    exponentiate = TRUE
  ),
  extract_model_term(
    activity_model,
    "Quasibinomial daylight use",
    "log_wolf_jackal_ratio_z",
    "Odds ratio per 1 SD increase in log wolf:jackal ratio",
    exponentiate = TRUE
  ),
  extract_model_term(
    abundance_model,
    "Gazelle abundance adjusted for open space",
    "log_wolf_jackal_ratio_z",
    "Standardized log-abundance change per 1 SD log-ratio"
  ),
  extract_model_term(
    observed_rate_model,
    "Observed count rate adjusted for open space",
    "log_wolf_jackal_ratio_z",
    "Standardized log-rate change per 1 SD log-ratio"
  ),
  extract_model_term(
    abundance_share_model,
    "Gazelle abundance using bounded wolf share",
    "wolf_share_z",
    "Standardized log-abundance change per 1 SD wolf share"
  ),
  extract_model_term(
    observed_rate_share_model,
    "Observed count rate using bounded wolf share",
    "wolf_share_z",
    "Standardized log-rate change per 1 SD wolf share"
  )
) %>%
  mutate(p_fdr = p.adjust(p_value, "BH"))

write_csv(model_results, paste0(prefix, "_sensitivity_models.csv"))

balance_composition_model <- glm(
  cbind(fawn_count, female_count) ~
    wolf_share_z + predator_total_index_z + open_spaces_z,
  family = quasibinomial,
  data = composition_data
)
balance_activity_model <- glm(
  cbind(gazelle_daylight_events, gazelle_dark_events) ~
    wolf_share_z + predator_total_index_z + open_spaces_z,
  family = quasibinomial,
  data = activity_data
)
balance_abundance_model <- lm(
  gazelle_abundance_log_z ~
    wolf_share_z + predator_total_index_z + open_spaces_z,
  data = camera_data
)
balance_observed_rate_model <- lm(
  observed_gazelle_rate_log_z ~
    wolf_share_z + predator_total_index_z + open_spaces_z,
  data = camera_data
)

balance_load_results <- bind_rows(
  extract_model_term(
    balance_composition_model, "Fawn composition", "wolf_share_z",
    "Predator balance", exponentiate = TRUE
  ),
  extract_model_term(
    balance_composition_model, "Fawn composition", "predator_total_index_z",
    "Total predator index", exponentiate = TRUE
  ),
  extract_model_term(
    balance_activity_model, "All-gazelle daylight use", "wolf_share_z",
    "Predator balance", exponentiate = TRUE
  ),
  extract_model_term(
    balance_activity_model, "All-gazelle daylight use", "predator_total_index_z",
    "Total predator index", exponentiate = TRUE
  ),
  extract_model_term(
    balance_abundance_model, "Gazelle abundance index", "wolf_share_z",
    "Predator balance"
  ),
  extract_model_term(
    balance_abundance_model, "Gazelle abundance index", "predator_total_index_z",
    "Total predator index"
  ),
  extract_model_term(
    balance_observed_rate_model, "Observed identified count rate", "wolf_share_z",
    "Predator balance"
  ),
  extract_model_term(
    balance_observed_rate_model, "Observed identified count rate", "predator_total_index_z",
    "Total predator index"
  )
) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH"))

write_csv(
  balance_load_results,
  paste0(prefix, "_predator_balance_load_models.csv")
)

theme_exploratory <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#F7F3EA", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "#4E5D61"),
      axis.title = element_text(face = "bold")
    )
}

ratio_breaks <- log(c(0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8))
ratio_labels <- c("0.001", "0.01", "0.1", "0.5", "1", "2", "4", "8")

cor_subtitle <- function(target_variable) {
  row <- correlations %>%
    filter(.data$variable == .env$target_variable) %>%
    slice(1)
  sprintf(
    "Spearman rho = %.2f (95%% bootstrap CI %.2f to %.2f), FDR p = %s, n = %d",
    row$spearman_rho,
    row$conf_low,
    row$conf_high,
    ifelse(is.na(row$p_fdr), "NA", format.pval(row$p_fdr, digits = 2)),
    row$n_cameras
  )
}

base_scatter <- function(data, y, title, subtitle, y_label) {
  ggplot(data, aes(x = log_wolf_jackal_ratio, y = .data[[y]])) +
    geom_point(size = 2.1, alpha = 0.76, color = "#244F5A") +
    geom_text(
      data = data %>% filter(camera_name == extreme_ratio_camera),
      aes(label = paste("Camera", camera_name)),
      hjust = -0.08,
      vjust = -0.7,
      size = 2.8,
      color = "#6E3D32"
    ) +
    geom_smooth(
      method = "lm", formula = y ~ x, se = TRUE,
      color = "#B64B32", fill = "#DFA58E", linewidth = 0.8
    ) +
    scale_x_continuous(
      breaks = ratio_breaks,
      labels = ratio_labels,
      name = "Wolf:jackal abundance ratio (log scale)"
    ) +
    labs(y = y_label, title = title, subtitle = subtitle) +
    theme_exploratory()
}

p_open <- base_scatter(
  camera_data,
  "Open.Spaces",
  "A. Open-space cover",
  cor_subtitle("Open.Spaces"),
  "Open-space proportion"
) + scale_y_continuous(labels = label_percent(accuracy = 1))

p_abundance <- base_scatter(
  camera_data,
  "gazelleabundance",
  "B. Gazelle abundance index",
  cor_subtitle("gazelleabundance"),
  "Gazelle abundance index"
)

p_composition <- base_scatter(
  camera_data %>% filter(fawn_female_total > 0),
  "fawns_per_100_females_smoothed",
  "C. Camera-based fawn:female composition",
  cor_subtitle("fawn_composition"),
  "Smoothed fawns per 100 adult females"
)

p_daylight <- base_scatter(
  activity_data,
  "gazelle_daylight_fraction",
  "D. Diurnal versus nocturnal activity",
  cor_subtitle("gazelle_daylight_fraction"),
  "Proportion of independent events in daylight"
) + scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1))

primary_figure <- (p_open + p_abundance) / (p_composition + p_daylight) +
  plot_annotation(
    title = "Exploratory camera-level relationships with the wolf:jackal abundance ratio",
    subtitle = paste0(
      "Raw ratio shown on a log scale; correlations use ranks. ",
      "Activity panels require at least ", minimum_activity_events,
      " independent events per camera."
    ),
    caption = paste0(
      "Points are cameras. Lines are simple linear summaries on the log-ratio scale and are not causal models. ",
      "Fawn composition uses a 0.5 continuity correction only for visualization."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "#4E5D61"),
      plot.caption = element_text(color = "#5F696B", hjust = 0)
    )
  )

ggsave(
  paste0(prefix, "_primary_relationships.png"),
  primary_figure,
  width = 12,
  height = 9,
  dpi = 320
)
ggsave(
  paste0(prefix, "_primary_relationships.pdf"),
  primary_figure,
  width = 12,
  height = 9
)

forest_data <- correlations %>%
  filter(scope != "Diagnostic") %>%
  mutate(
    supported = if_else(!is.na(p_fdr) & p_fdr < 0.05, "FDR p < 0.05", "Not FDR-supported"),
    outcome = reorder(outcome, spearman_rho)
  )

correlation_figure <- ggplot(
  forest_data,
  aes(x = spearman_rho, y = outcome, color = supported)
) +
  geom_vline(xintercept = 0, color = "#788285", linewidth = 0.5) +
  geom_segment(
    aes(x = conf_low, xend = conf_high, yend = outcome),
    linewidth = 0.8
  ) +
  geom_point(size = 2.7) +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(
    values = c("FDR p < 0.05" = "#A54231", "Not FDR-supported" = "#285F73"),
    name = NULL
  ) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(
    x = "Spearman correlation with wolf:jackal abundance ratio",
    y = NULL,
    title = "Camera-level correlation screen",
    subtitle = "Points are rank correlations; intervals are 2,000 camera-bootstrap percentiles",
    caption = "FDR correction was applied separately to the four primary outcomes and the secondary exploratory screen."
  ) +
  theme_exploratory() +
  theme(
    legend.position = "top",
    strip.text.y = element_text(face = "bold", angle = 0),
    strip.background = element_rect(fill = "#E7DFCF", color = NA),
    panel.grid.major.y = element_blank()
  )

ggsave(
  paste0(prefix, "_correlation_summary.png"),
  correlation_figure,
  width = 10.5,
  height = 9.5,
  dpi = 320
)
ggsave(
  paste0(prefix, "_correlation_summary.pdf"),
  correlation_figure,
  width = 10.5,
  height = 9.5
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt(x, 3)))
}

ratio_summary <- camera_data %>%
  summarise(
    minimum = min(wolf_jackal_ratio),
    q1 = quantile(wolf_jackal_ratio, 0.25),
    median = median(wolf_jackal_ratio),
    q3 = quantile(wolf_jackal_ratio, 0.75),
    maximum = max(wolf_jackal_ratio)
  )

report <- c(
  "# Exploratory wolf:jackal abundance-ratio analysis",
  "",
  "## Ratio definition",
  "",
  paste0(
    "For each of 60 cameras, the raw covariate was `Null.Wolf.Abundance / null.jackal.abundance`. ",
    "Both indices were strictly positive, so no pseudocount was used. The raw ratio ranged from `",
    fmt(ratio_summary$minimum, 4), "` to `", fmt(ratio_summary$maximum),
    "` (median `", fmt(ratio_summary$median), "`). Analyses used its natural logarithm because the raw ratio was right-skewed; Spearman correlations are unchanged by this monotonic transformation."
  ),
  "",
  "A ratio above 1 indicates a wolf-dominated abundance index and a value below 1 a jackal-dominated index. This is a relative index, not a direct measure of predation risk or predator density.",
  "",
  "## Primary correlations",
  "",
  "| Outcome | Cameras | Spearman rho | 95% camera-bootstrap CI | Raw p | Primary FDR p |",
  "|---|---:|---:|---:|---:|---:|"
)

for (i in which(correlations$scope == "Primary")) {
  report <- c(report, sprintf(
    "| %s | %d | `%s` | `[%s, %s]` | `%s` | `%s` |",
    correlations$outcome[i], correlations$n_cameras[i],
    fmt(correlations$spearman_rho[i]), fmt(correlations$conf_low[i]),
    fmt(correlations$conf_high[i]), fmt_p(correlations$p_value[i]),
    fmt_p(correlations$p_fdr[i])
  ))
}

primary_open <- correlations %>% filter(variable == "Open.Spaces")
primary_abundance <- correlations %>% filter(variable == "gazelleabundance")
primary_composition <- correlations %>% filter(variable == "fawn_composition")
primary_daylight <- correlations %>%
  filter(variable == "gazelle_daylight_fraction")
adult_female_rate <- correlations %>%
  filter(variable == "female_count_rate_100_days")
adult_male_rate <- correlations %>%
  filter(variable == "male_count_rate_100_days")
balance_observed <- balance_load_results %>%
  filter(model == "Observed identified count rate", term == "wolf_share_z")
load_observed <- balance_load_results %>%
  filter(
    model == "Observed identified count rate",
    term == "predator_total_index_z"
  )

report <- c(
  report,
  "",
  "## Key interpretation",
  "",
  paste0(
    "- The wolf:jackal ratio was not detectably associated with open-space cover (`rho = ",
    fmt(primary_open$spearman_rho), "`, FDR `p = ",
    fmt_p(primary_open$p_fdr),
    "`), so the ratio was not simply an open-space gradient in this camera set."
  ),
  paste0(
    "- Cameras with more wolf-dominated predator indices had higher gazelle abundance (`rho = ",
    fmt(primary_abundance$spearman_rho), "`, 95% CI `",
    fmt(primary_abundance$conf_low), "-", fmt(primary_abundance$conf_high),
    "`, FDR `p = ", fmt_p(primary_abundance$p_fdr),
    "`). The same pattern was stronger for effort-standardized observed counts."
  ),
  paste0(
    "- The positive count-rate association was clearest for adult females (`rho = ",
    fmt(adult_female_rate$spearman_rho), "`) and males (`rho = ",
    fmt(adult_male_rate$spearman_rho),
    "`), whereas the fawn count-rate trend was weaker."
  ),
  paste0(
    "- Fawn representation (`rho = ", fmt(primary_composition$spearman_rho),
    "`, FDR `p = ", fmt_p(primary_composition$p_fdr),
    "`) and all-gazelle daylight use (`rho = ",
    fmt(primary_daylight$spearman_rho), "`, FDR `p = ",
    fmt_p(primary_daylight$p_fdr),
    "`) showed essentially no ratio association."
  ),
  paste0(
    "- In the model separating predator balance from total predator activity, observed gazelle count rate increased with wolf share (`beta = ",
    fmt(balance_observed$estimate), "`, FDR `p = ",
    fmt_p(balance_observed$p_fdr),
    "`) but not with the total wolf-plus-jackal index (`beta = ",
    fmt(load_observed$estimate), "`, FDR `p = ",
    fmt_p(load_observed$p_fdr),
    "`). This is an association with relative predator composition, not evidence that wolves cause higher gazelle abundance."
  )
)

report <- c(
  report,
  "",
  "## Secondary correlation screen",
  "",
  "| Domain | Outcome | Cameras | Spearman rho | 95% camera-bootstrap CI | Secondary FDR p |",
  "|---|---|---:|---:|---:|---:|"
)

for (i in which(correlations$scope == "Secondary")) {
  report <- c(report, sprintf(
    "| %s | %s | %d | `%s` | `[%s, %s]` | `%s` |",
    correlations$domain[i], correlations$outcome[i],
    correlations$n_cameras[i], fmt(correlations$spearman_rho[i]),
    fmt(correlations$conf_low[i]), fmt(correlations$conf_high[i]),
    fmt_p(correlations$p_fdr[i])
  ))
}

report <- c(
  report,
  "",
  "## Denominator-aware and partial-model sensitivities",
  "",
  "| Model | Cameras | Ratio effect | 95% CI | Odds ratio (quasibinomial only) | FDR p |",
  "|---|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(model_results))) {
  exponentiated_text <- if (
    is.finite(model_results$odds_or_multiplicative_effect[i])
  ) {
    paste0(
      fmt(model_results$odds_or_multiplicative_effect[i]), " [",
      fmt(model_results$effect_conf_low[i]), ", ",
      fmt(model_results$effect_conf_high[i]), "]"
    )
  } else {
    "not applicable"
  }
  report <- c(report, sprintf(
    "| %s | %d | `%s` | `[%s, %s]` | `%s` | `%s` |",
    model_results$model[i], model_results$n_cameras[i],
    fmt(model_results$estimate[i]), fmt(model_results$conf_low[i]),
    fmt(model_results$conf_high[i]), exponentiated_text,
    fmt_p(model_results$p_fdr[i])
  ))
}

report <- c(
  report,
  "",
  "## Abundance-signal robustness checks",
  "",
  "| Check | Outcome | Excluded camera | Cameras | Partial or Spearman rho | FDR p |",
  "|---|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(robustness_checks))) {
  report <- c(report, sprintf(
    "| %s | %s | `%s` | %d | `%s` | `%s` |",
    robustness_checks$check[i], robustness_checks$outcome[i],
    ifelse(
      is.na(robustness_checks$excluded_camera[i]),
      "none",
      robustness_checks$excluded_camera[i]
    ),
    robustness_checks$n_cameras[i],
    fmt(robustness_checks$spearman_rho[i]),
    fmt_p(robustness_checks$p_fdr[i])
  ))
}

report <- c(
  report,
  "",
  "## Predator balance versus total predator index",
  "",
  "These models include bounded wolf share, the standardized sum of wolf and jackal indices, and open-space cover simultaneously.",
  "",
  "| Outcome | Component | Effect per 1 SD | 95% CI | FDR p |",
  "|---|---|---:|---:|---:|"
)

for (i in seq_len(nrow(balance_load_results))) {
  report <- c(report, sprintf(
    "| %s | %s | `%s` | `[%s, %s]` | `%s` |",
    balance_load_results$model[i], balance_load_results$scale[i],
    fmt(balance_load_results$estimate[i]),
    fmt(balance_load_results$conf_low[i]),
    fmt(balance_load_results$conf_high[i]),
    fmt_p(balance_load_results$p_fdr[i])
  ))
}

diagnostic <- correlations %>% filter(scope == "Diagnostic")
report <- c(
  report,
  "",
  "## Construction diagnostics and cautions",
  "",
  paste0(
    "- As expected mathematically, the ratio correlated positively with wolf abundance (`rho = ",
    fmt(diagnostic$spearman_rho[diagnostic$variable == "Null.Wolf.Abundance"]),
    "`) and negatively with jackal abundance (`rho = ",
    fmt(diagnostic$spearman_rho[diagnostic$variable == "null.jackal.abundance"]), "`)."
  ),
  "- Ratio associations can arise from changes in wolves, jackals, or both. Any promising result should therefore be checked against models containing the two predator indices separately.",
  "- The ratio cannot distinguish low-wolf/low-jackal from high-wolf/high-jackal cameras with the same balance. The balance-plus-load models retain total predator activity as a separate covariate.",
  paste0(
    "- Camera `", extreme_ratio_camera,
    "` had the minimum ratio because its wolf index was almost zero. Rank correlations were not driven by this camera, but ordinary linear log-ratio models were sensitive to its leverage; bounded wolf share and rank-based checks are more stable summaries."
  ),
  "- The fawn:female outcome is a camera-based composition index, not recruitment or survival. The quasibinomial sensitivity accounts for unequal fawn-plus-female totals across cameras.",
  "- Daylight use is the fraction of 30-min independent events occurring at solar elevation >= 0 degrees. The quasibinomial sensitivity accounts for unequal event totals across cameras.",
  "- These tests do not account for spatial autocorrelation, camera detectability, or other confounding beyond the stated open-space partial models. They are hypothesis-generating rather than confirmatory.",
  "",
  "## Outputs",
  "",
  paste0("- Enriched site covariates: `SiteCovs_060324_with_wolf_jackal_ratio.csv`"),
  paste0("- Analysis-ready camera table: `", prefix, "_camera_data.csv`"),
  paste0("- Correlations: `", prefix, "_correlations.csv`"),
  paste0("- Sensitivity models: `", prefix, "_sensitivity_models.csv`"),
  paste0("- Robustness checks: `", prefix, "_robustness_checks.csv`"),
  paste0("- Predator balance/load models: `", prefix, "_predator_balance_load_models.csv`"),
  paste0("- Primary figure: `", prefix, "_primary_relationships.png`"),
  paste0("- Correlation summary: `", prefix, "_correlation_summary.png`")
)

writeLines(report, paste0(prefix, "_results.md"), useBytes = TRUE)

cat("Exploratory wolf:jackal ratio analysis complete.\n")
