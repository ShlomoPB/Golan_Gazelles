#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(broom.mixed)
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(mgcv)
  library(overlap)
  library(patchwork)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
cfg <- list(
  event_file = "Gazella gazella_independent_events_30min_analysis_ready.csv",
  site_cov_file = "SiteCovs_060324.csv",
  obs_cov_file = "ObsCovs_050524.csv",
  fawn_file = "gazelle_juvenile_abundance_60x7.csv",
  female_file = "gazelle_female_abundance_60x7.csv",
  timestamp_mode = "local",
  local_timezone = "Asia/Jerusalem",
  bootstrap_reps = "500",
  seed = "20260623",
  prefix = "publication_unified_spatiotemporal"
)

for (arg in args) {
  pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(pieces) == 2) {
    cfg[[sub("^--", "", pieces[1])]] <- pieces[2]
  }
}

if (!cfg$timestamp_mode %in% c("utc", "local")) {
  stop("timestamp_mode must be 'utc' or 'local'.")
}

bootstrap_reps <- as.integer(cfg$bootstrap_reps)
seed <- as.integer(cfg$seed)
set.seed(seed)

class_levels <- c("Fawn", "Female", "Male")
period_levels <- c("Early", "Late")
class_colors <- c(Fawn = "#149c8c", Female = "#c81d25", Male = "#17365d")
period_colors <- c(Early = "#d97706", Late = "#276678")

theme_publication <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = rel(1.15)),
      legend.position = "bottom"
    )
}

z_scale <- function(x) {
  x <- as.numeric(x)
  sd_x <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sd_x) || sd_x == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / sd_x
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, "<0.001", formatC(x, digits = 3, format = "f"))
  )
}

safe_fit <- function(expr, label) {
  tryCatch(
    expr,
    error = function(e) {
      warning(label, " failed: ", conditionMessage(e))
      NULL
    }
  )
}

model_diagnostics <- function(model, model_name) {
  if (is.null(model)) {
    return(data.frame(
      model = model_name,
      convergence_code = NA_integer_,
      positive_definite_hessian = NA,
      max_abs_gradient = NA_real_,
      AIC = NA_real_,
      nobs = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  gradient <- tryCatch(
    max(abs(model$sdr$gradient.fixed), na.rm = TRUE),
    error = function(e) NA_real_
  )
  data.frame(
    model = model_name,
    convergence_code = model$fit$convergence,
    positive_definite_hessian = isTRUE(model$sdr$pdHess),
    max_abs_gradient = gradient,
    AIC = stats::AIC(model),
    nobs = stats::nobs(model),
    stringsAsFactors = FALSE
  )
}

parse_event_time <- function(x, mode, local_timezone) {
  if (mode == "utc") {
    return(as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  }
  as.POSIXct(
    sub("Z$", "", x),
    format = "%Y-%m-%dT%H:%M:%S",
    tz = local_timezone
  )
}

equation_of_time <- function(doy) {
  b <- 2 * pi * (doy - 81) / 364
  9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
}

apparent_solar_hour <- function(ts, longitude) {
  doy <- as.integer(format(ts, "%j", tz = "UTC"))
  utc_hour <- as.numeric(format(ts, "%H", tz = "UTC")) +
    as.numeric(format(ts, "%M", tz = "UTC")) / 60 +
    as.numeric(format(ts, "%S", tz = "UTC")) / 3600
  (utc_hour + longitude / 15 + equation_of_time(doy) / 60) %% 24
}

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

assign_rounds <- function(events, obs) {
  event_index <- events %>%
    transmute(
      event_row_id = row_number(),
      camera_name,
      event_date
    )
  matches <- event_index %>%
    inner_join(
      obs %>% select(camera_name, round, start_date, end_date),
      by = "camera_name",
      relationship = "many-to-many"
    ) %>%
    filter(
      !is.na(start_date),
      event_date >= start_date,
      event_date <= end_date
    ) %>%
    group_by(event_row_id) %>%
    # Boundary dates can belong to adjacent rounds; retain the later round.
    slice_max(start_date, n = 1, with_ties = FALSE) %>%
    ungroup()
  matches$round[match(seq_len(nrow(events)), matches$event_row_id)]
}

read_count_long <- function(path, value_name) {
  raw <- read.csv(path, na.strings = c("NA", ""), check.names = FALSE)
  raw %>%
    pivot_longer(
      cols = all_of(paste0("X", 1:7)),
      names_to = "round_name",
      values_to = value_name
    ) %>%
    transmute(
      camera_name = as.integer(camera_name),
      round = as.integer(sub("^X", "", round_name)),
      !!value_name := as.numeric(.data[[value_name]])
    )
}

fixed_model_matrix <- function(model, newdata) {
  fixed_formula <- lme4::nobars(stats::formula(model))
  trm <- stats::delete.response(stats::terms(fixed_formula))
  mm <- stats::model.matrix(trm, newdata)
  beta_names <- names(glmmTMB::fixef(model)$cond)
  missing_cols <- setdiff(beta_names, colnames(mm))
  for (nm in missing_cols) {
    mm <- cbind(mm, setNames(data.frame(rep(0, nrow(mm))), nm))
  }
  mm[, beta_names, drop = FALSE]
}

linear_contrast <- function(model, contrast, label, class = NA_character_, period = NA_character_) {
  beta <- glmmTMB::fixef(model)$cond
  vc <- as.matrix(stats::vcov(model)$cond)
  estimate <- as.numeric(contrast %*% beta)
  se <- sqrt(as.numeric(contrast %*% vc %*% t(contrast)))
  data.frame(
    class = class,
    period = period,
    effect = label,
    estimate = estimate,
    std_error = se,
    lwr = estimate - 1.96 * se,
    upr = estimate + 1.96 * se,
    exp_estimate = exp(estimate),
    exp_lwr = exp(estimate - 1.96 * se),
    exp_upr = exp(estimate + 1.96 * se),
    stringsAsFactors = FALSE
  )
}

class_effect_contrasts <- function(model, base_row, classes, periods = NULL) {
  out <- list()
  idx <- 1
  effects <- c("wolf_z", "jackal_z", "open_z", "built_z", "human_z")
  period_values <- if (is.null(periods)) NA_character_ else periods
  for (class_i in classes) {
    for (period_i in period_values) {
      base <- base_row
      base$class <- factor(class_i, levels = class_levels)
      if (!is.null(periods)) {
        base$period <- factor(period_i, levels = period_levels)
      }
      for (effect_i in effects) {
        low <- high <- base
        low[[effect_i]] <- -0.5
        high[[effect_i]] <- 0.5
        contrast <- fixed_model_matrix(model, high) - fixed_model_matrix(model, low)
        out[[idx]] <- linear_contrast(
          model, contrast, effect_i, class_i,
          if (is.null(periods)) NA_character_ else period_i
        )
        idx <- idx + 1
      }
      interaction_points <- lapply(
        list(c(-0.5, -0.5), c(-0.5, 0.5), c(0.5, -0.5), c(0.5, 0.5)),
        function(vals) {
          x <- base
          x$wolf_z <- vals[1]
          x$jackal_z <- vals[2]
          x
        }
      )
      mm <- lapply(interaction_points, function(x) fixed_model_matrix(model, x))
      contrast <- mm[[4]] - mm[[3]] - mm[[2]] + mm[[1]]
      out[[idx]] <- linear_contrast(
        model, contrast, "wolf_z:jackal_z", class_i,
        if (is.null(periods)) NA_character_ else period_i
      )
      idx <- idx + 1
    }
  }
  bind_rows(out)
}

effect_vector <- function(model, base_row, class_i, effect_i, period_i = NULL) {
  base <- base_row
  base$class <- factor(class_i, levels = class_levels)
  if (!is.null(period_i)) {
    base$period <- factor(period_i, levels = period_levels)
  }
  if (effect_i != "wolf_z:jackal_z") {
    low <- high <- base
    low[[effect_i]] <- -0.5
    high[[effect_i]] <- 0.5
    return(
      fixed_model_matrix(model, high) -
        fixed_model_matrix(model, low)
    )
  }
  points <- lapply(
    list(c(-0.5, -0.5), c(-0.5, 0.5), c(0.5, -0.5), c(0.5, 0.5)),
    function(vals) {
      x <- base
      x$wolf_z <- vals[1]
      x$jackal_z <- vals[2]
      fixed_model_matrix(model, x)
    }
  )
  points[[4]] - points[[3]] - points[[2]] + points[[1]]
}

class_effect_differences <- function(model, base_row, classes) {
  effects <- c(
    "wolf_z", "jackal_z", "open_z", "built_z", "human_z",
    "wolf_z:jackal_z"
  )
  pairs <- combn(classes, 2, simplify = FALSE)
  rows <- list()
  idx <- 1
  for (effect_i in effects) {
    for (pair_i in pairs) {
      contrast <- effect_vector(model, base_row, pair_i[1], effect_i) -
        effect_vector(model, base_row, pair_i[2], effect_i)
      row <- linear_contrast(model, contrast, effect_i)
      row$comparison <- paste(pair_i[1], "-", pair_i[2])
      rows[[idx]] <- row
      idx <- idx + 1
    }
  }
  bind_rows(rows) %>%
    select(comparison, effect, everything(), -class, -period)
}

simulate_count_checks <- function(model, data, nsim = 200, seed = 1) {
  if (is.null(model)) {
    return(data.frame())
  }
  set.seed(seed)
  sims <- as.data.frame(simulate(model, nsim = nsim))
  stats_for <- function(x) {
    c(
      zeros = sum(x == 0),
      mean = mean(x),
      variance = stats::var(x),
      maximum = max(x)
    )
  }
  rows <- list()
  idx <- 1
  for (class_i in class_levels) {
    keep <- data$class == class_i
    observed <- stats_for(data$event_n[keep])
    simulated <- vapply(sims, function(x) stats_for(x[keep]), numeric(4))
    for (stat_i in names(observed)) {
      sim_i <- simulated[stat_i, ]
      rows[[idx]] <- data.frame(
        class = class_i,
        statistic = stat_i,
        observed = observed[stat_i],
        simulated_mean = mean(sim_i),
        simulated_lwr = unname(stats::quantile(sim_i, 0.025)),
        simulated_upr = unname(stats::quantile(sim_i, 0.975)),
        upper_tail_probability = mean(sim_i >= observed[stat_i]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  bind_rows(rows)
}

predict_glmmtmb <- function(model, newdata) {
  pred <- predict(
    model,
    newdata = newdata,
    type = "link",
    se.fit = TRUE,
    re.form = NA,
    allow.new.levels = TRUE
  )
  newdata %>%
    mutate(
      estimate_link = as.numeric(pred$fit),
      std_error_link = as.numeric(pred$se.fit),
      estimate = plogis(estimate_link),
      lwr = plogis(estimate_link - 1.96 * std_error_link),
      upr = plogis(estimate_link + 1.96 * std_error_link)
    )
}

predict_count_glmmtmb <- function(model, newdata) {
  pred <- predict(
    model,
    newdata = newdata,
    type = "link",
    se.fit = TRUE,
    re.form = NA,
    allow.new.levels = TRUE
  )
  newdata %>%
    mutate(
      estimate_link = as.numeric(pred$fit),
      std_error_link = as.numeric(pred$se.fit),
      estimate = exp(estimate_link),
      lwr = exp(estimate_link - 1.96 * std_error_link),
      upr = exp(estimate_link + 1.96 * std_error_link)
    )
}

cluster_overlap <- function(dat, class_a, class_b, period_i, cameras, reps) {
  dat_i <- if (period_i == "All") dat else filter(dat, period == period_i)
  a <- dat_i$time_rad[dat_i$class == class_a]
  b <- dat_i$time_rad[dat_i$class == class_b]
  overlap_type <- if (min(length(a), length(b)) < 75) "Dhat1" else "Dhat4"
  point <- tryCatch(
    unname(overlapEst(a, b, type = overlap_type)),
    error = function(e) NA_real_
  )
  boot <- rep(NA_real_, reps)
  for (r in seq_len(reps)) {
    sampled <- sample(cameras, length(cameras), replace = TRUE)
    a_boot <- unlist(lapply(
      sampled,
      function(cam) dat_i$time_rad[dat_i$camera_name == cam & dat_i$class == class_a]
    ))
    b_boot <- unlist(lapply(
      sampled,
      function(cam) dat_i$time_rad[dat_i$camera_name == cam & dat_i$class == class_b]
    ))
    if (length(a_boot) >= 10 && length(b_boot) >= 10) {
      boot[r] <- tryCatch(
        unname(overlapEst(a_boot, b_boot, type = overlap_type)),
        error = function(e) NA_real_
      )
    }
  }
  data.frame(
    period = period_i,
    class_a = class_a,
    class_b = class_b,
    n_events_a = length(a),
    n_events_b = length(b),
    n_cameras_a = n_distinct(dat_i$camera_name[dat_i$class == class_a]),
    n_cameras_b = n_distinct(dat_i$camera_name[dat_i$class == class_b]),
    estimator = overlap_type,
    overlap = point,
    lwr = unname(stats::quantile(boot, 0.025, na.rm = TRUE)),
    upr = unname(stats::quantile(boot, 0.975, na.rm = TRUE)),
    bootstrap_success = sum(is.finite(boot)),
    stringsAsFactors = FALSE
  )
}

events_raw <- read_csv(
  cfg$event_file,
  col_types = cols(event_time = col_character()),
  show_col_types = FALSE
)
site_raw <- read_csv(cfg$site_cov_file, show_col_types = FALSE)
obs_raw <- read_csv(cfg$obs_cov_file, show_col_types = FALSE)

site <- site_raw %>%
  transmute(
    camera_name = as.integer(Camera.ID),
    Polygon.Name = factor(Polygon.Name),
    wolf = as.numeric(Null.Wolf.Abundance),
    jackal = as.numeric(null.jackal.abundance),
    open = as.numeric(Open.Spaces),
    built = as.numeric(Built.and.Roads),
    human = as.numeric(human.abundance),
    reserve = as.numeric(Nature.Reserves),
    wolf_z = z_scale(wolf),
    jackal_z = z_scale(jackal),
    open_z = z_scale(open),
    built_z = z_scale(built),
    human_z = z_scale(human)
    ,
    reserve_z = z_scale(reserve)
  )

obs <- obs_raw %>%
  transmute(
    camera_name = as.integer(`Camera Number`),
    round = as.integer(`Round Number`),
    start_date = as.Date(`Start date`, format = "%d/%m/%Y"),
    end_date = as.Date(`End Date`, format = "%d/%m/%Y"),
    days_operational = as.numeric(Days.Operational),
    path_type = factor(Path.Type),
    camera_model = factor(Camera.Model),
    vegetation_structure = factor(Vegetation.Structure),
    area_behind_camera = factor(Area.Behind.Camera)
  ) %>%
  left_join(site, by = "camera_name") %>%
  mutate(
    midpoint_date = start_date + floor(as.numeric(end_date - start_date) / 2),
    camera_factor = factor(camera_name),
    round_factor = factor(round),
    camera_round_id = factor(paste(camera_name, round, sep = "_"))
  )

active_obs <- obs %>%
  filter(
    is.finite(days_operational),
    days_operational > 0,
    !is.na(start_date),
    !is.na(end_date)
  )

survey_start <- min(active_obs$start_date)
survey_end <- max(active_obs$end_date)
survey_split <- survey_start + floor(as.numeric(survey_end - survey_start) / 2)

active_obs <- active_obs %>%
  mutate(
    period = factor(
      ifelse(midpoint_date <= survey_split, "Early", "Late"),
      levels = period_levels
    ),
    log_days = log(days_operational)
  )

events <- events_raw %>%
  mutate(
    camera_name = as.integer(camera_name),
    class = recode(
      as.character(gazelle_class),
      juvenile = "Fawn",
      female = "Female",
      male = "Male"
    ),
    class = factor(class, levels = class_levels),
    timestamp = parse_event_time(event_time, cfg$timestamp_mode, cfg$local_timezone),
    event_date = as.Date(format(timestamp, tz = cfg$local_timezone)),
    solar_hour = apparent_solar_hour(timestamp, camelot_gps_longitude),
    solar_elevation = solar_elevation(
      timestamp,
      camelot_gps_latitude,
      camelot_gps_longitude
    ),
    time_rad = 2 * pi * solar_hour / 24,
    light_phase = factor(
      case_when(
        solar_elevation >= 0 ~ "Daylight",
        solar_elevation >= -6 ~ "Civil twilight",
        TRUE ~ "Night"
      ),
      levels = c("Daylight", "Civil twilight", "Night")
    )
  )

events$round <- assign_rounds(events, obs)

event_audit <- data.frame(
  metric = c(
    "raw_event_rows",
    "mapped_to_any_round",
    "mapped_to_active_round",
    "excluded_outside_active_survey",
    "survey_start",
    "survey_end",
    "survey_midpoint_split"
  ),
  value = c(
    nrow(events),
    sum(!is.na(events$round)),
    sum(
      paste(events$camera_name, events$round, sep = "_") %in%
        active_obs$camera_round_id
    ),
    sum(
      !paste(events$camera_name, events$round, sep = "_") %in%
        active_obs$camera_round_id
    ),
    as.character(survey_start),
    as.character(survey_end),
    as.character(survey_split)
  ),
  stringsAsFactors = FALSE
)

events_mapped <- events %>%
  filter(
    paste(camera_name, round, sep = "_") %in% active_obs$camera_round_id
  ) %>%
  select(
    camera_name, round, class, timestamp, event_date, event_count,
    solar_hour, solar_elevation, time_rad, light_phase
  ) %>%
  left_join(
    active_obs %>%
      select(
        camera_name, round, period, days_operational, log_days,
        Polygon.Name, wolf_z, jackal_z, open_z, built_z, human_z,
        reserve_z,
        camera_factor, round_factor, camera_round_id
      ),
    by = c("camera_name", "round")
  )

event_counts <- events_mapped %>%
  group_by(camera_name, round, class) %>%
  summarise(
    event_n = n(),
    individual_n = sum(event_count, na.rm = TRUE),
    .groups = "drop"
  )

space_data <- active_obs %>%
  select(
    camera_name, round, period, days_operational, log_days,
    Polygon.Name, wolf_z, jackal_z, open_z, built_z, human_z,
    reserve_z,
    camera_factor, round_factor, camera_round_id
  ) %>%
  crossing(class = factor(class_levels, levels = class_levels)) %>%
  left_join(event_counts, by = c("camera_name", "round", "class")) %>%
  mutate(
    event_n = replace_na(event_n, 0L),
    individual_n = replace_na(individual_n, 0),
    present = as.integer(event_n > 0),
    class = factor(class, levels = class_levels)
  )

space_full_formula <- event_n ~
  class * (wolf_z * jackal_z + open_z + built_z + human_z) +
  reserve_z +
  period +
  offset(log_days) +
  (1 | Polygon.Name / camera_factor)

space_common_formula <- event_n ~
  class +
  wolf_z * jackal_z +
  open_z + built_z + human_z +
  reserve_z +
  period +
  offset(log_days) +
  (1 | Polygon.Name / camera_factor)

space_full <- safe_fit(
  glmmTMB(
    space_full_formula,
    family = nbinom2,
    dispformula = ~ class,
    data = space_data
  ),
  "Full space-use NB model"
)

space_common <- safe_fit(
  glmmTMB(
    space_common_formula,
    family = nbinom2,
    dispformula = ~ class,
    data = space_data
  ),
  "Common-slope space-use NB model"
)

space_zinb <- safe_fit(
  glmmTMB(
    space_full_formula,
    family = nbinom2,
    ziformula = ~ class,
    dispformula = ~ class,
    data = space_data
  ),
  "Zero-inflated space-use NB model"
)

space_full_nb1 <- safe_fit(
  glmmTMB(
    space_full_formula,
    family = nbinom1,
    dispformula = ~ class,
    data = space_data
  ),
  "Full space-use NB1 model"
)

space_common_nb1 <- safe_fit(
  glmmTMB(
    space_common_formula,
    family = nbinom1,
    dispformula = ~ class,
    data = space_data
  ),
  "Common-slope space-use NB1 model"
)

encounter_full <- safe_fit(
  glmmTMB(
    update(space_full_formula, present ~ . - offset(log_days)),
    family = binomial,
    data = space_data
  ),
  "Encounter model"
)

positive_data <- filter(space_data, event_n > 0)
positive_intensity <- safe_fit(
  glmmTMB(
    space_full_formula,
    family = truncated_nbinom2,
    dispformula = ~ class,
    data = positive_data
  ),
  "Positive-event intensity model"
)

space_lrt <- if (!is.null(space_full) && !is.null(space_common)) {
  as.data.frame(anova(space_common, space_full)) %>%
    mutate(comparison = "Common versus class-specific ecological slopes")
} else {
  data.frame()
}

space_diagnostics <- bind_rows(
  model_diagnostics(space_common, "space_common_nb"),
  model_diagnostics(space_full, "space_class_specific_nb"),
  model_diagnostics(space_common_nb1, "space_common_nb1"),
  model_diagnostics(space_full_nb1, "space_class_specific_nb1"),
  model_diagnostics(space_zinb, "space_class_specific_zinb"),
  model_diagnostics(encounter_full, "space_encounter_binomial"),
  model_diagnostics(positive_intensity, "space_positive_truncated_nb")
)

base_space <- space_data[1, , drop = FALSE]
base_space$round <- 4L
base_space$round_factor <- factor(4, levels = levels(space_data$round_factor))
base_space$period <- factor("Early", levels = period_levels)
base_space$wolf_z <- 0
base_space$jackal_z <- 0
base_space$open_z <- 0
base_space$built_z <- 0
base_space$human_z <- 0
base_space$reserve_z <- 0
base_space$log_days <- 0

space_contrasts <- if (!is.null(space_full)) {
  class_effect_contrasts(space_full, base_space, class_levels)
} else {
  data.frame()
}

space_class_differences <- if (!is.null(space_full)) {
  class_effect_differences(space_full, base_space, class_levels)
} else {
  data.frame()
}

encounter_contrasts <- if (!is.null(encounter_full)) {
  class_effect_contrasts(encounter_full, base_space, class_levels)
} else {
  data.frame()
}

positive_intensity_contrasts <- if (!is.null(positive_intensity)) {
  class_effect_contrasts(positive_intensity, base_space, class_levels)
} else {
  data.frame()
}

space_simulation_checks <- simulate_count_checks(
  space_full,
  space_data,
  nsim = 200,
  seed = seed
)

space_nb1_contrasts <- if (!is.null(space_full_nb1)) {
  class_effect_contrasts(space_full_nb1, base_space, class_levels)
} else {
  data.frame()
}

space_nb1_simulation_checks <- simulate_count_checks(
  space_full_nb1,
  space_data,
  nsim = 200,
  seed = seed + 1
)

space_prediction_grid <- if (!is.null(space_full)) {
  expand_grid(
    class = factor(class_levels, levels = class_levels),
    wolf_z = seq(-1.5, 1.5, length.out = 61),
    jackal_z = seq(-1.5, 1.5, length.out = 61)
  ) %>%
    mutate(
      open_z = 0,
      built_z = 0,
      human_z = 0,
      reserve_z = 0,
      round = 4L,
      period = factor("Early", levels = period_levels),
      log_days = 0,
      Polygon.Name = base_space$Polygon.Name,
      camera_factor = base_space$camera_factor,
      round_factor = factor(4, levels = levels(space_data$round_factor)),
      camera_round_id = base_space$camera_round_id
    ) %>%
    predict_count_glmmtmb(space_full, .)
} else {
  data.frame()
}

fawn_long <- read_count_long(cfg$fawn_file, "fawn_n")
female_long <- read_count_long(cfg$female_file, "female_n")

composition_data <- active_obs %>%
  select(
    camera_name, round, period, Polygon.Name, wolf_z, jackal_z,
    open_z, built_z, human_z, camera_factor, round_factor
    , reserve_z
  ) %>%
  left_join(fawn_long, by = c("camera_name", "round")) %>%
  left_join(female_long, by = c("camera_name", "round")) %>%
  filter(is.finite(fawn_n), is.finite(female_n)) %>%
  mutate(total_classified = fawn_n + female_n) %>%
  filter(total_classified > 0)

composition_full <- safe_fit(
  glmmTMB(
    cbind(fawn_n, female_n) ~
      wolf_z * jackal_z +
      open_z + built_z + human_z +
      reserve_z +
      round_factor +
      (1 | Polygon.Name) +
      (1 | camera_factor),
    family = betabinomial,
    data = composition_data
  ),
  "Fawn:female composition model"
)

composition_no_predator_interaction <- safe_fit(
  glmmTMB(
    cbind(fawn_n, female_n) ~
      wolf_z + jackal_z +
      open_z + built_z + human_z +
      reserve_z +
      round_factor +
      (1 | Polygon.Name) +
      (1 | camera_factor),
    family = betabinomial,
    data = composition_data
  ),
  "Composition model without wolf-jackal interaction"
)

composition_lrt <- if (
  !is.null(composition_full) &&
    !is.null(composition_no_predator_interaction)
) {
  as.data.frame(anova(composition_no_predator_interaction, composition_full)) %>%
    mutate(comparison = "Additive versus interacting wolf-jackal composition model")
} else {
  data.frame()
}

composition_diagnostics <- bind_rows(
  model_diagnostics(composition_no_predator_interaction, "composition_additive"),
  model_diagnostics(composition_full, "composition_wolf_jackal_interaction")
)

base_composition <- composition_data[1, , drop = FALSE]
base_composition$round <- 4L
base_composition$round_factor <- factor(4, levels = levels(composition_data$round_factor))
base_composition$wolf_z <- 0
base_composition$jackal_z <- 0
base_composition$open_z <- 0
base_composition$built_z <- 0
base_composition$human_z <- 0
base_composition$reserve_z <- 0

composition_contrasts <- if (!is.null(composition_full)) {
  effects <- c("wolf_z", "jackal_z", "open_z", "built_z", "human_z", "reserve_z")
  rows <- lapply(effects, function(effect_i) {
    low <- high <- base_composition
    low[[effect_i]] <- -0.5
    high[[effect_i]] <- 0.5
    linear_contrast(
      composition_full,
      fixed_model_matrix(composition_full, high) -
        fixed_model_matrix(composition_full, low),
      effect_i
    )
  })
  interaction_points <- lapply(
    list(c(-0.5, -0.5), c(-0.5, 0.5), c(0.5, -0.5), c(0.5, 0.5)),
    function(vals) {
      x <- base_composition
      x$wolf_z <- vals[1]
      x$jackal_z <- vals[2]
      x
    }
  )
  mm <- lapply(interaction_points, function(x) fixed_model_matrix(composition_full, x))
  rows[[length(rows) + 1]] <- linear_contrast(
    composition_full,
    mm[[4]] - mm[[3]] - mm[[2]] + mm[[1]],
    "wolf_z:jackal_z"
  )
  bind_rows(rows)
} else {
  data.frame()
}

composition_prediction_grid <- if (!is.null(composition_full)) {
  expand_grid(
    wolf_z = seq(-1.5, 1.5, length.out = 61),
    jackal_z = seq(-1.5, 1.5, length.out = 61)
  ) %>%
    mutate(
      open_z = 0,
      built_z = 0,
      human_z = 0,
      reserve_z = 0,
      round = 4L,
      period = factor("Early", levels = period_levels),
      Polygon.Name = base_composition$Polygon.Name,
      camera_factor = base_composition$camera_factor,
      round_factor = factor(4, levels = levels(composition_data$round_factor))
    ) %>%
    predict_glmmtmb(composition_full, .)
} else {
  data.frame()
}

light_counts <- events_mapped %>%
  mutate(daylight = light_phase == "Daylight") %>%
  group_by(camera_name, round, class) %>%
  summarise(
    daylight_n = sum(daylight),
    dark_n = sum(!daylight),
    total_events = n(),
    .groups = "drop"
  ) %>%
  left_join(
    active_obs %>%
      select(
        camera_name, round, period, Polygon.Name, wolf_z, jackal_z,
        open_z, built_z, human_z, camera_factor, round_factor
        , reserve_z
      ),
    by = c("camera_name", "round")
  ) %>%
  mutate(class = factor(class, levels = class_levels))

diel_full <- safe_fit(
  glmmTMB(
    cbind(daylight_n, dark_n) ~
      class * (period + wolf_z * jackal_z + open_z + built_z + human_z) +
      reserve_z +
      (1 | Polygon.Name) +
      (1 | camera_factor),
    family = betabinomial,
    data = light_counts
  ),
  "Class-specific daylight-use model"
)

diel_common <- safe_fit(
  glmmTMB(
    cbind(daylight_n, dark_n) ~
      class +
      period +
      wolf_z * jackal_z +
      open_z + built_z + human_z +
      reserve_z +
      (1 | Polygon.Name) +
      (1 | camera_factor),
    family = betabinomial,
    data = light_counts
  ),
  "Common-slope daylight-use model"
)

diel_lrt <- if (!is.null(diel_full) && !is.null(diel_common)) {
  as.data.frame(anova(diel_common, diel_full)) %>%
    mutate(comparison = "Common versus class-specific daylight-use slopes")
} else {
  data.frame()
}

diel_diagnostics <- bind_rows(
  model_diagnostics(diel_common, "diel_daylight_common"),
  model_diagnostics(diel_full, "diel_daylight_class_specific")
)

base_diel <- light_counts[1, , drop = FALSE]
base_diel$wolf_z <- 0
base_diel$jackal_z <- 0
base_diel$open_z <- 0
base_diel$built_z <- 0
base_diel$human_z <- 0
base_diel$reserve_z <- 0

diel_contrasts <- if (!is.null(diel_full)) {
  class_effect_contrasts(diel_full, base_diel, class_levels)
} else {
  data.frame()
}

diel_period_contrasts <- if (!is.null(diel_full)) {
  bind_rows(lapply(class_levels, function(class_i) {
    early <- late <- base_diel
    early$class <- factor(class_i, levels = class_levels)
    late$class <- factor(class_i, levels = class_levels)
    early$period <- factor("Early", levels = period_levels)
    late$period <- factor("Late", levels = period_levels)
    linear_contrast(
      diel_full,
      fixed_model_matrix(diel_full, late) -
        fixed_model_matrix(diel_full, early),
      "Late versus early daylight use",
      class = class_i
    )
  }))
} else {
  data.frame()
}

diel_prediction_grid <- if (!is.null(diel_full)) {
  expand_grid(
    class = factor(class_levels, levels = class_levels),
    period = factor(period_levels, levels = period_levels),
    wolf_z = seq(-1.5, 1.5, length.out = 61),
    jackal_z = c(-1, 0, 1)
  ) %>%
    mutate(
      open_z = 0,
      built_z = 0,
      human_z = 0,
      reserve_z = 0,
      Polygon.Name = base_diel$Polygon.Name,
      camera_factor = base_diel$camera_factor,
      round_factor = base_diel$round_factor
    ) %>%
    predict_glmmtmb(diel_full, .)
} else {
  data.frame()
}

bin_width <- 2
hour_bins <- seq(0, 22, by = bin_width)
event_bins <- events_mapped %>%
  mutate(
    hour_bin = floor(solar_hour / bin_width) * bin_width,
    hour_mid = hour_bin + bin_width / 2
  ) %>%
  count(camera_name, round, class, hour_mid, name = "event_n")

diel_bin_data <- active_obs %>%
  select(
    camera_name, round, period, days_operational,
    camera_factor, round_factor
  ) %>%
  crossing(
    class = factor(class_levels, levels = class_levels),
    hour_mid = hour_bins + bin_width / 2
  ) %>%
  left_join(event_bins, by = c("camera_name", "round", "class", "hour_mid")) %>%
  mutate(
    event_n = replace_na(event_n, 0L),
    class_period = interaction(class, period, sep = " | ", drop = TRUE),
    log_exposure = log(days_operational * bin_width / 24)
  )

diel_gam <- safe_fit(
  gam(
    event_n ~
      class_period +
      s(hour_mid, by = class_period, bs = "cc", k = 8) +
      s(camera_factor, bs = "re") +
      s(round_factor, bs = "re") +
      offset(log_exposure),
    family = nb(),
    method = "REML",
    data = diel_bin_data,
    knots = list(hour_mid = c(0, 24))
  ),
  "Cyclic diel GAM"
)

diel_curve_data <- if (!is.null(diel_gam)) {
  expand_grid(
    class = factor(class_levels, levels = class_levels),
    period = factor(period_levels, levels = period_levels),
    hour_mid = seq(0, 24, length.out = 193)
  ) %>%
    mutate(
      class_period = interaction(class, period, sep = " | ", drop = TRUE),
      camera_factor = diel_bin_data$camera_factor[1],
      round_factor = diel_bin_data$round_factor[1],
      log_exposure = 0
    ) %>%
    bind_cols(
      as.data.frame(predict(diel_gam, newdata = ., type = "link", se.fit = TRUE)) %>%
        rename(estimate_link = fit, std_error_link = se.fit)
    ) %>%
    mutate(
      rate = exp(estimate_link),
      rate_lwr = exp(estimate_link - 1.96 * std_error_link),
      rate_upr = exp(estimate_link + 1.96 * std_error_link)
    ) %>%
    group_by(class, period) %>%
    mutate(
      density = rate / sum(rate),
      density_lwr = rate_lwr / sum(rate),
      density_upr = rate_upr / sum(rate)
    ) %>%
    ungroup()
} else {
  data.frame()
}

diel_gam_smooths <- if (!is.null(diel_gam)) {
  smooth_table <- as.data.frame(summary(diel_gam)$s.table)
  smooth_table$term <- rownames(smooth_table)
  rownames(smooth_table) <- NULL
  smooth_table %>% select(term, everything())
} else {
  data.frame()
}

camera_universe <- sort(unique(active_obs$camera_name))
overlap_rows <- list()
idx <- 1
for (period_i in c("All", period_levels)) {
  pairs <- combn(class_levels, 2, simplify = FALSE)
  for (pair_i in pairs) {
    overlap_rows[[idx]] <- cluster_overlap(
      events_mapped,
      pair_i[1],
      pair_i[2],
      period_i,
      camera_universe,
      bootstrap_reps
    )
    idx <- idx + 1
  }
}
overlap_results <- bind_rows(overlap_rows)

all_diagnostics <- bind_rows(
  space_diagnostics,
  composition_diagnostics,
  diel_diagnostics
)

space_coefficients <- if (!is.null(space_full)) {
  tidy(space_full, effects = "fixed", component = "cond", conf.int = TRUE)
} else {
  data.frame()
}

composition_coefficients <- if (!is.null(composition_full)) {
  tidy(composition_full, effects = "fixed", component = "cond", conf.int = TRUE)
} else {
  data.frame()
}

diel_coefficients <- if (!is.null(diel_full)) {
  tidy(diel_full, effects = "fixed", component = "cond", conf.int = TRUE)
} else {
  data.frame()
}

write_csv(event_audit, paste0(cfg$prefix, "_event_audit.csv"))
write_csv(space_data, paste0(cfg$prefix, "_camera_round_event_rates.csv"))
write_csv(all_diagnostics, paste0(cfg$prefix, "_model_diagnostics.csv"))
write_csv(space_lrt, paste0(cfg$prefix, "_space_class_interaction_lrt.csv"))
write_csv(space_coefficients, paste0(cfg$prefix, "_space_coefficients.csv"))
write_csv(space_contrasts, paste0(cfg$prefix, "_space_class_effects.csv"))
write_csv(space_class_differences, paste0(cfg$prefix, "_space_class_effect_differences.csv"))
write_csv(encounter_contrasts, paste0(cfg$prefix, "_space_encounter_effects.csv"))
write_csv(
  positive_intensity_contrasts,
  paste0(cfg$prefix, "_space_positive_intensity_effects.csv")
)
write_csv(space_simulation_checks, paste0(cfg$prefix, "_space_simulation_checks.csv"))
write_csv(space_nb1_contrasts, paste0(cfg$prefix, "_space_nb1_class_effects.csv"))
write_csv(
  space_nb1_simulation_checks,
  paste0(cfg$prefix, "_space_nb1_simulation_checks.csv")
)
write_csv(space_prediction_grid, paste0(cfg$prefix, "_space_prediction_grid.csv"))
write_csv(composition_lrt, paste0(cfg$prefix, "_composition_interaction_lrt.csv"))
write_csv(composition_coefficients, paste0(cfg$prefix, "_composition_coefficients.csv"))
write_csv(composition_contrasts, paste0(cfg$prefix, "_composition_effects.csv"))
write_csv(composition_prediction_grid, paste0(cfg$prefix, "_composition_prediction_grid.csv"))
write_csv(diel_lrt, paste0(cfg$prefix, "_diel_class_interaction_lrt.csv"))
write_csv(diel_coefficients, paste0(cfg$prefix, "_diel_coefficients.csv"))
write_csv(diel_contrasts, paste0(cfg$prefix, "_diel_class_effects.csv"))
write_csv(diel_period_contrasts, paste0(cfg$prefix, "_diel_period_effects.csv"))
write_csv(diel_prediction_grid, paste0(cfg$prefix, "_diel_prediction_grid.csv"))
write_csv(diel_curve_data, paste0(cfg$prefix, "_diel_curve_data.csv"))
write_csv(diel_gam_smooths, paste0(cfg$prefix, "_diel_gam_smooths.csv"))
write_csv(overlap_results, paste0(cfg$prefix, "_camera_clustered_overlap.csv"))

saveRDS(
  list(
    space_full = space_full,
    space_common = space_common,
    space_zinb = space_zinb,
    space_full_nb1 = space_full_nb1,
    space_common_nb1 = space_common_nb1,
    encounter_full = encounter_full,
    positive_intensity = positive_intensity,
    composition_full = composition_full,
    composition_no_predator_interaction = composition_no_predator_interaction,
    diel_full = diel_full,
    diel_common = diel_common,
    diel_gam = diel_gam,
    config = cfg
  ),
  paste0(cfg$prefix, "_model_fits.rds")
)
capture.output(
  sessionInfo(),
  file = paste0(cfg$prefix, "_session_info.txt")
)

if (nrow(space_prediction_grid) > 0) {
  space_heatmap <- space_prediction_grid %>%
    mutate(rate_plot = pmin(estimate, stats::quantile(estimate, 0.98, na.rm = TRUE))) %>%
    ggplot(aes(x = jackal_z, y = wolf_z, fill = rate_plot)) +
    geom_raster() +
    facet_wrap(~ class, nrow = 1) +
    scale_fill_viridis_c(option = "C", name = "Events per\ncamera-day") +
    labs(
      x = "Jackal abundance (z)",
      y = "Wolf abundance (z)",
      title = "Predicted gazelle independent-event rate"
    ) +
    theme_publication()
  ggsave(paste0(cfg$prefix, "_space_heatmap.png"), space_heatmap, width = 10, height = 4.2, dpi = 320)
  ggsave(paste0(cfg$prefix, "_space_heatmap.pdf"), space_heatmap, width = 10, height = 4.2)
}

if (nrow(space_contrasts) > 0) {
  space_forest <- space_contrasts %>%
    mutate(
      effect = recode(
        effect,
        wolf_z = "Wolf",
        jackal_z = "Jackal",
        open_z = "Open space",
        built_z = "Built and roads",
        human_z = "Human activity",
        `wolf_z:jackal_z` = "Wolf x jackal"
      )
    ) %>%
    ggplot(aes(x = estimate, y = effect, color = class)) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
    geom_errorbar(
      aes(xmin = lwr, xmax = upr),
      width = 0,
      orientation = "y",
      position = position_dodge(width = 0.55)
    ) +
    geom_point(position = position_dodge(width = 0.55), size = 2) +
    scale_color_manual(values = class_colors) +
    labs(
      x = "Effect on log independent-event rate per 1 SD",
      y = NULL,
      title = "Class-specific spatial-use effects"
    ) +
    theme_publication()
  ggsave(paste0(cfg$prefix, "_space_effects.png"), space_forest, width = 8.5, height = 5.3, dpi = 320)
  ggsave(paste0(cfg$prefix, "_space_effects.pdf"), space_forest, width = 8.5, height = 5.3)
}

if (nrow(composition_prediction_grid) > 0) {
  composition_heatmap <- composition_prediction_grid %>%
    ggplot(aes(x = jackal_z, y = wolf_z, fill = estimate)) +
    geom_raster() +
    scale_fill_viridis_c(
      option = "B",
      limits = c(0, max(composition_prediction_grid$estimate, na.rm = TRUE)),
      name = "Predicted fawn\nproportion"
    ) +
      labs(
        x = "Jackal abundance (z)",
        y = "Wolf abundance (z)",
        title = "Predicted fawn representation relative to adult females"
      ) +
      theme_publication() +
      theme(legend.position = "right")
  ggsave(paste0(cfg$prefix, "_composition_heatmap.png"), composition_heatmap, width = 6.8, height = 5.5, dpi = 320)
  ggsave(paste0(cfg$prefix, "_composition_heatmap.pdf"), composition_heatmap, width = 6.8, height = 5.5)
}

if (nrow(diel_curve_data) > 0) {
  diel_curves <- ggplot(
    diel_curve_data,
    aes(x = hour_mid, y = density, color = class, fill = class)
  ) +
    geom_ribbon(aes(ymin = density_lwr, ymax = density_upr), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    facet_grid(period ~ class, scales = "free_y") +
    scale_color_manual(values = class_colors) +
    scale_fill_manual(values = class_colors) +
    scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
    labs(
      x = "Apparent solar time",
      y = "Normalized event density",
      title = "Camera-aware diel activity curves by class and survey period"
    ) +
    theme_publication() +
    theme(legend.position = "none")
  ggsave(paste0(cfg$prefix, "_diel_curves.png"), diel_curves, width = 10, height = 6.5, dpi = 320)
  ggsave(paste0(cfg$prefix, "_diel_curves.pdf"), diel_curves, width = 10, height = 6.5)
}

if (nrow(diel_prediction_grid) > 0) {
  diel_lines <- diel_prediction_grid %>%
    mutate(
      jackal_level = factor(
        jackal_z,
        levels = c(-1, 0, 1),
        labels = c("Low jackal", "Mean jackal", "High jackal")
      )
    ) %>%
    ggplot(aes(x = wolf_z, y = estimate, color = jackal_level, fill = jackal_level)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_grid(period ~ class) +
    scale_color_manual(values = c("#2166ac", "#777777", "#b2182b")) +
    scale_fill_manual(values = c("#2166ac", "#777777", "#b2182b")) +
    labs(
      x = "Wolf abundance (z)",
      y = "Predicted probability an event occurred in daylight",
      color = NULL,
      fill = NULL,
      title = "Class-specific daylight use across predator gradients"
    ) +
    theme_publication()
  ggsave(paste0(cfg$prefix, "_daylight_predictions.png"), diel_lines, width = 10, height = 6.5, dpi = 320)
  ggsave(paste0(cfg$prefix, "_daylight_predictions.pdf"), diel_lines, width = 10, height = 6.5)
}

space_class_p <- if (nrow(space_lrt) > 0) tail(space_lrt$`Pr(>Chisq)`, 1) else NA_real_
composition_interaction_p <- if (nrow(composition_lrt) > 0) {
  tail(composition_lrt$`Pr(>Chisq)`, 1)
} else {
  NA_real_
}
diel_class_p <- if (nrow(diel_lrt) > 0) tail(diel_lrt$`Pr(>Chisq)`, 1) else NA_real_
nb1_delta_aic <- if (!is.null(space_full) && !is.null(space_full_nb1)) {
  stats::AIC(space_full_nb1) - stats::AIC(space_full)
} else {
  NA_real_
}

results_lines <- c(
  "# Unified spatiotemporal gazelle analysis",
  "",
  "## Scope and data treatment",
  "",
  paste0("- Timestamp interpretation: `", cfg$timestamp_mode, "`."),
  paste0("- Formal deployment window: `", survey_start, "` to `", survey_end, "`."),
  paste0("- Independent early/late split: `", survey_split, "`, the midpoint of the deployment window."),
  paste0("- Raw independent-event rows: `", nrow(events), "`."),
  paste0("- Event rows retained within active camera-rounds: `", nrow(events_mapped), "`."),
  paste0("- Event rows excluded outside active camera-rounds: `", nrow(events) - nrow(events_mapped), "`."),
  "- Spatial use was modeled as independent-event rate per operational camera-day, not latent animal abundance.",
  "- Juvenile:adult-female results are demographic composition indices, not direct recruitment or survival estimates.",
  "- Temporal overlap intervals used camera-level cluster bootstrap resampling.",
  "",
  "## Model comparisons",
  "",
  paste0(
    "- Evidence that ecological slopes differ among classes in the spatial event-rate model: likelihood-ratio `p = ",
    fmt_p(space_class_p), "`."
  ),
  paste0(
    "- Evidence for a wolf x jackal effect on fawn representation relative to females: likelihood-ratio `p = ",
    fmt_p(composition_interaction_p), "`."
  ),
  paste0(
    "- Evidence that daylight-use slopes differ among classes: likelihood-ratio `p = ",
    fmt_p(diel_class_p), "`."
  ),
  paste0(
    "- NB1 sensitivity model minus NB2 primary-model AIC: `",
    fmt_num(nb1_delta_aic, 1),
    "`; positive values favor NB2."
  ),
  "",
  "## Model diagnostics",
  "",
  "| Model | AIC | Convergence | Positive-definite Hessian | Maximum gradient |",
  "|---|---:|---:|---|---:|"
)

for (i in seq_len(nrow(all_diagnostics))) {
  results_lines <- c(
    results_lines,
    paste0(
      "| `", all_diagnostics$model[i], "` | `", fmt_num(all_diagnostics$AIC[i], 1),
      "` | `", all_diagnostics$convergence_code[i], "` | `",
      all_diagnostics$positive_definite_hessian[i], "` | `",
      fmt_num(all_diagnostics$max_abs_gradient[i], 4), "` |"
    )
  )
}

results_lines <- c(
  results_lines,
  "",
  "## Spatial model simulation checks",
  "",
  "| Class | Statistic | Observed | Simulated mean | Simulated 95% interval | Upper-tail probability |",
  "|---|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(space_simulation_checks))) {
  results_lines <- c(
    results_lines,
    paste0(
      "| `", space_simulation_checks$class[i], "` | `",
      space_simulation_checks$statistic[i], "` | `",
      fmt_num(space_simulation_checks$observed[i], 2), "` | `",
      fmt_num(space_simulation_checks$simulated_mean[i], 2), "` | `",
      fmt_num(space_simulation_checks$simulated_lwr[i], 2), " to ",
      fmt_num(space_simulation_checks$simulated_upr[i], 2), "` | `",
      fmt_num(space_simulation_checks$upper_tail_probability[i], 3), "` |"
    )
  )
}

results_lines <- c(
  results_lines,
  "",
  "## Class-specific spatial-use effects",
  "",
  "Effects are changes in log independent-event rate per 1 SD increase, evaluated at mean values of the interacting predator covariate.",
  "",
  "| Class | Effect | Estimate | 95% CI | Rate ratio |",
  "|---|---|---:|---:|---:|"
)

if (nrow(space_contrasts) > 0) {
  for (i in seq_len(nrow(space_contrasts))) {
    results_lines <- c(
      results_lines,
      paste0(
        "| `", space_contrasts$class[i], "` | `", space_contrasts$effect[i],
        "` | `", fmt_num(space_contrasts$estimate[i]), "` | `",
        fmt_num(space_contrasts$lwr[i]), " to ", fmt_num(space_contrasts$upr[i]),
        "` | `", fmt_num(space_contrasts$exp_estimate[i]), "` |"
      )
    )
  }
}

supported_class_differences <- space_class_differences %>%
  filter(lwr > 0 | upr < 0)

results_lines <- c(
  results_lines,
  "",
  "## Direct tests of class differences in spatial slopes",
  ""
)

if (nrow(supported_class_differences) == 0) {
  results_lines <- c(
    results_lines,
    "No individual pairwise class difference had a 95% confidence interval excluding zero, despite the omnibus class-by-environment test."
  )
} else {
  results_lines <- c(
    results_lines,
    "| Comparison | Effect | Difference | 95% CI |",
    "|---|---|---:|---:|"
  )
  for (i in seq_len(nrow(supported_class_differences))) {
    results_lines <- c(
      results_lines,
      paste0(
        "| `", supported_class_differences$comparison[i], "` | `",
        supported_class_differences$effect[i], "` | `",
        fmt_num(supported_class_differences$estimate[i]), "` | `",
        fmt_num(supported_class_differences$lwr[i]), " to ",
        fmt_num(supported_class_differences$upr[i]), "` |"
      )
    )
  }
}

results_lines <- c(
  results_lines,
  "",
  "## Fawn:adult-female composition effects",
  "",
  "| Effect | Log-odds estimate | 95% CI | Odds ratio |",
  "|---|---:|---:|---:|"
)

if (nrow(composition_contrasts) > 0) {
  for (i in seq_len(nrow(composition_contrasts))) {
    results_lines <- c(
      results_lines,
      paste0(
        "| `", composition_contrasts$effect[i], "` | `",
        fmt_num(composition_contrasts$estimate[i]), "` | `",
        fmt_num(composition_contrasts$lwr[i]), " to ",
        fmt_num(composition_contrasts$upr[i]), "` | `",
        fmt_num(composition_contrasts$exp_estimate[i]), "` |"
      )
    )
  }
}

results_lines <- c(
  results_lines,
  "",
  "## Early-to-late changes in daylight use",
  "",
  "| Class | Log-odds change | 95% CI | Odds ratio |",
  "|---|---:|---:|---:|"
)

if (nrow(diel_period_contrasts) > 0) {
  for (i in seq_len(nrow(diel_period_contrasts))) {
    results_lines <- c(
      results_lines,
      paste0(
        "| `", diel_period_contrasts$class[i], "` | `",
        fmt_num(diel_period_contrasts$estimate[i]), "` | `",
        fmt_num(diel_period_contrasts$lwr[i]), " to ",
        fmt_num(diel_period_contrasts$upr[i]), "` | `",
        fmt_num(diel_period_contrasts$exp_estimate[i]), "` |"
      )
    )
  }
}

results_lines <- c(
  results_lines,
  "",
  "## Camera-clustered temporal overlap",
  "",
  "| Period | Comparison | Events | Cameras | Overlap [95% cluster-bootstrap CI] |",
  "|---|---|---:|---:|---:|"
)

for (i in seq_len(nrow(overlap_results))) {
  results_lines <- c(
    results_lines,
    paste0(
      "| `", overlap_results$period[i], "` | `",
      overlap_results$class_a[i], " vs ", overlap_results$class_b[i],
      "` | `", overlap_results$n_events_a[i], " / ", overlap_results$n_events_b[i],
      "` | `", overlap_results$n_cameras_a[i], " / ", overlap_results$n_cameras_b[i],
      "` | `", fmt_num(overlap_results$overlap[i], 3), " [",
      fmt_num(overlap_results$lwr[i], 3), ", ",
      fmt_num(overlap_results$upr[i], 3), "]` |"
    )
  )
}

results_lines <- c(
  results_lines,
  "",
  "## Interpretation guardrails",
  "",
  "- Predator associations remain observational and cannot establish wolf-mediated suppression of jackals without temporally resolved predator data or a causal design.",
  "- The zero-inflated spatial model was retained only as a sensitivity fit and should be rejected if its Hessian is not positive definite.",
  "- NB1 reproduced the main coefficient directions and better matched some extreme-count statistics, but NB2 was retained when its information-criterion support was substantially stronger.",
  "- Rebuild `glmmTMB` against the installed `TMB` version before the archived publication run if R reports a binary-version mismatch.",
  "- The daylight model is conditional on successful demographic classification. Unclassified gazelle detections are needed to test whether fawns are harder to identify in infrared imagery.",
  "- Camera timestamps were confirmed to be Israel local time at installation, so solar-time and daylight analyses use `--timestamp_mode=local` as the primary interpretation.",
  "- The early/late comparison is a survey-period contrast. It should not be described as a temperature effect until external or field temperature data are joined."
)

writeLines(results_lines, paste0(cfg$prefix, "_results.md"))

cat("Unified spatiotemporal analysis complete.\n")
cat("Prefix:", cfg$prefix, "\n")
cat("Timestamp mode:", cfg$timestamp_mode, "\n")
cat("Survey split:", as.character(survey_split), "\n")
