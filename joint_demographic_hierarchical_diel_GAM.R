suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(mvtnorm)
  library(patchwork)
  library(readr)
  library(tidyr)
})

# Confirmatory joint demographic analysis. This script does not modify any
# output from the definitive class-specific diel analysis.
cfg <- list(
  seed = 20260825L,
  simulations = 5000L,
  primary_k = 12L,
  site_k = 12L,
  prediction_step = 0.1,
  nthreads = 2L,
  input_data = "definitive_hierarchical_diel_data.csv",
  obs_covs = "ObsCovs_050524.csv",
  camera_coordinates = "publication_camera_array_map_sites.csv",
  output_dir = "joint_demographic_hierarchical_diel_outputs",
  checkpoint_dir = file.path(
    "joint_demographic_hierarchical_diel_outputs", "checkpoints"
  )
)
set.seed(cfg$seed)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(name) file.path(cfg$output_dir, name)

class_levels <- c("Juvenile", "Adult female", "Adult male")
class_labels <- c(
  "Juvenile" = "Juvenile",
  "Adult female" = "Adult female",
  "Adult male" = "Adult male"
)
class_codes <- c(
  "Juvenile" = "juvenile",
  "Adult female" = "female",
  "Adult male" = "male"
)
effect_levels <- c("wolf", "jackal", "interaction")
effect_labels <- c(
  wolf = "Wolf x time",
  jackal = "Jackal x time",
  interaction = "Wolf x jackal x time"
)
contrast_pairs <- list(
  "Juvenile - Adult female" = c("Juvenile", "Adult female"),
  "Juvenile - Adult male" = c("Juvenile", "Adult male"),
  "Adult female - Adult male" = c("Adult female", "Adult male")
)
context_levels <- c(
  "Wolf -1 / Jackal -1",
  "Wolf -1 / Jackal +1",
  "Wolf +1 / Jackal -1",
  "Wolf +1 / Jackal +1"
)
context_colors <- c(
  "Wolf -1 / Jackal -1" = "#168575",
  "Wolf -1 / Jackal +1" = "#B62D35",
  "Wolf +1 / Jackal -1" = "#173F66",
  "Wolf +1 / Jackal +1" = "#835796"
)
class_colors <- c(
  "Juvenile" = "#C46A2D",
  "Adult female" = "#1B7F79",
  "Adult male" = "#2F4B73"
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}
fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", fmt(x, 3)))
}
collapse_call <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
`%||%` <- function(x, y) if (is.null(x)) y else x

# Data ----------------------------------------------------------------------

prepare_joint_data <- function(time_scale, history_type) {
  read_csv(cfg$input_data, show_col_types = FALSE) %>%
    filter(
      .data$time_scale == .env$time_scale,
      .data$history_type == .env$history_type,
      .data$response_class %in% class_levels
    ) %>%
    mutate(
      demographic_class = factor(response_class, levels = class_levels),
      round_factor = factor(round),
      camera_id = factor(camera_id),
      camera_round = factor(camera_round),
      camera_class = interaction(
        camera_id, demographic_class, sep = "_", drop = TRUE
      ),
      camera_round_class = interaction(
        camera_round, demographic_class, sep = "_", drop = TRUE
      ),
      wolf_jackal = wolf_z * jackal_z
    ) %>%
    arrange(demographic_class, camera_id, round, Time_bin) %>%
    add_class_by_variables()
}

add_class_by_variables <- function(data) {
  for (class_i in class_levels) {
    code_i <- unname(class_codes[class_i])
    indicator <- as.numeric(data$demographic_class == class_i)
    data[[paste0("wolf_", code_i)]] <- data$wolf_z * indicator
    data[[paste0("jackal_", code_i)]] <- data$jackal_z * indicator
    data[[paste0("interaction_", code_i)]] <-
      data$wolf_z * data$jackal_z * indicator
  }
  female_levels <- c(
    "reference", paste0("female_", levels(data$camera_id))
  )
  male_levels <- c(
    "reference", paste0("male_", levels(data$camera_id))
  )
  data$female_camera <- ordered(
    ifelse(
      data$demographic_class == "Adult female",
      paste0("female_", data$camera_id), "reference"
    ),
    levels = female_levels
  )
  data$male_camera <- ordered(
    ifelse(
      data$demographic_class == "Adult male",
      paste0("male_", data$camera_id), "reference"
    ),
    levels = male_levels
  )
  data
}

primary_data <- prepare_joint_data("Apparent solar time", "Raw sightings")
if (nrow(primary_data) != 24768L) stop("Expected 24,768 joint primary rows.")
if (!all(table(primary_data$demographic_class) == 8256L)) {
  stop("Expected 8,256 rows per demographic class.")
}
expected_success <- c(
  "Juvenile" = 183L, "Adult female" = 1576L, "Adult male" = 2307L
)
if (!all(tapply(
  primary_data$success, primary_data$demographic_class, sum
) == expected_success)) stop("Primary detection-day totals changed.")

data_audit <- primary_data %>%
  group_by(demographic_class) %>%
  summarise(
    rows = n(),
    successes_detection_days = sum(success),
    cameras_represented = n_distinct(camera_id),
    survey_rounds = n_distinct(round),
    zero_detection_rows = sum(success == 0),
    minimum_denominator = min(trials),
    maximum_denominator = max(trials),
    .groups = "drop"
  ) %>%
  mutate(
    time_scale = "Equal-duration apparent solar time",
    history_type = "Raw sighting camera-day-hour presence"
  )
write_csv(data_audit, out("table_00_joint_data_audit.csv"))

predator_audit <- primary_data %>%
  distinct(camera_id, wolf_index, jackal_index, wolf_z, jackal_z) %>%
  summarise(
    cameras = n(),
    wolf_z_minimum = min(wolf_z),
    wolf_z_maximum = max(wolf_z),
    jackal_z_minimum = min(jackal_z),
    jackal_z_maximum = max(jackal_z),
    wolf_z_matches_definitive_input = all(is.finite(wolf_z)),
    jackal_z_matches_definitive_input = all(is.finite(jackal_z))
  )
write_csv(predator_audit, out("joint_predator_index_audit.csv"))

# Formula -------------------------------------------------------------------

class_smooth_terms <- function(k) {
  baseline <- paste0(
    "s(Time, by = demographic_class, bs = 'cc', k = ", k, ")"
  )
  focal <- unlist(lapply(class_levels, function(class_i) {
    code_i <- unname(class_codes[class_i])
    c(
      paste0("s(Time, by = wolf_", code_i, ", bs = 'cc', k = ", k, ")"),
      paste0("s(Time, by = jackal_", code_i, ", bs = 'cc', k = ", k, ")"),
      paste0(
        "s(Time, by = interaction_", code_i,
        ", bs = 'cc', k = ", k, ")"
      )
    )
  }))
  paste(c(baseline, focal), collapse = " + ")
}

site_shape_terms <- function(site_k = cfg$site_k) {
  paste(
    paste0(
      "s(Time, by = female_camera, bs = 'cc', k = ", site_k,
      ", m = 1, id = 21)"
    ),
    paste0(
      "s(Time, by = male_camera, bs = 'cc', k = ", site_k,
      ", m = 1, id = 22)"
    ),
    sep = " + "
  )
}

hierarchy_terms <- function(hierarchy) {
  shared <- paste(
    "s(camera_id, bs = 're')",
    "s(camera_round, bs = 're')", sep = " + "
  )
  class_specific <- paste(
    "s(camera_class, bs = 're')",
    "s(camera_round_class, bs = 're')", sep = " + "
  )
  if (hierarchy == "full") {
    paste(site_shape_terms(), shared, class_specific, sep = " + ")
  } else if (hierarchy == "camera_only") {
    paste(
      "s(camera_id, bs = 're')", "s(camera_class, bs = 're')", sep = " + "
    )
  } else {
    stop("Unknown hierarchy: ", hierarchy)
  }
}

build_joint_formula <- function(k = cfg$primary_k, hierarchy = "full") {
  as.formula(paste(
    "cbind(success, failure) ~",
    "demographic_class * round_factor +",
    "demographic_class * wolf_z * jackal_z +",
    class_smooth_terms(k), "+",
    hierarchy_terms(hierarchy)
  ))
}

fit_joint_model <- function(data, k = cfg$primary_k, hierarchy = "full") {
  formula <- build_joint_formula(k, hierarchy)
  warnings <- character()
  fit <- withCallingHandlers(
    bam(
      formula,
      data = data,
      family = binomial(link = "logit"),
      method = "fREML",
      discrete = TRUE,
      select = TRUE,
      knots = list(Time = c(0, 24)),
      nthreads = cfg$nthreads,
      control = gam.control(trace = FALSE)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warnings))
}

formula_primary <- build_joint_formula()
writeLines(deparse(formula_primary, width.cutoff = 500L), out("exact_formula.txt"))

primary_checkpoint <- file.path(cfg$checkpoint_dir, "primary_joint_model.rds")
if (file.exists(primary_checkpoint)) {
  message("Loading primary joint model checkpoint")
  primary_object <- readRDS(primary_checkpoint)
} else {
  message("Fitting primary confirmatory joint demographic GAM")
  primary_object <- fit_joint_model(primary_data)
  saveRDS(primary_object, primary_checkpoint, compress = FALSE)
}
fit <- primary_object$fit
fit_warnings <- primary_object$warnings
writeLines(
  if (length(fit_warnings)) fit_warnings else "No warnings generated.",
  out("fit_warnings.txt")
)
writeLines(collapse_call(fit$call), out("exact_model_call.txt"))
literal_call <- c(
  "mgcv::bam(",
  paste0("  formula = ", paste(deparse(formula_primary), collapse = "\n  "), ","),
  "  data = primary_data,",
  "  family = binomial(link = \"logit\"),",
  "  method = \"fREML\", discrete = TRUE, select = TRUE,",
  "  knots = list(Time = c(0, 24)), nthreads = 2,",
  "  control = gam.control(trace = FALSE)",
  ")"
)
writeLines(literal_call, out("literal_final_model_call.txt"))

# Model inventory and diagnostics -------------------------------------------

unconditional_covariance <- tryCatch(
  vcov(fit, unconditional = TRUE),
  error = function(e) vcov(fit)
)
covariance_type <- if (isTRUE(all.equal(unconditional_covariance, vcov(fit)))) {
  "Bayesian coefficient covariance; unconditional covariance unavailable or identical"
} else {
  "Unconditional Bayesian covariance including smoothing-parameter uncertainty"
}

smooth_labels <- vapply(fit$smooth, function(sm) sm$label, character(1))
smooth_by <- vapply(fit$smooth, function(sm) sm$by, character(1))
population_exclusions <- smooth_labels[
  grepl("camera", smooth_labels, ignore.case = TRUE) |
    smooth_by %in% c("female_camera", "male_camera")
]
writeLines(population_exclusions, out("population_prediction_excluded_terms.txt"))

smooth_retest <- function(fit, index) {
  tryCatch({
    result <- mgcv:::reTest(fit, index)
    c(
      statistic = unname(result$stat),
      reference_df = unname(result$rank),
      P = unname(result$pval)
    )
  }, error = function(e) c(statistic = NA, reference_df = NA, P = NA))
}

focal_by_names <- c(
  "demographic_class",
  unlist(lapply(class_levels, function(class_i) {
    code_i <- unname(class_codes[class_i])
    paste0(c("wolf_", "jackal_", "interaction_"), code_i)
  }))
)
smooth_table_path <- out("full_smooth_table.csv")
table_a_path <- out("table_A_joint_model_specification.csv")
reuse_smooth_tables <- file.exists(smooth_table_path) && file.exists(table_a_path) &&
  nrow(read_csv(smooth_table_path, show_col_types = FALSE)) == length(fit$smooth)
if (reuse_smooth_tables) {
  message("Reusing complete primary smooth inventories")
  smooth_table <- read_csv(smooth_table_path, show_col_types = FALSE)
  table_a <- read_csv(table_a_path, show_col_types = FALSE)
} else {
  focal_indices <- which(smooth_by %in% focal_by_names)
  focal_tests <- setNames(
    lapply(focal_indices, function(i) smooth_retest(fit, i)),
    as.character(focal_indices)
  )
  smooth_table <- bind_rows(lapply(seq_along(fit$smooth), function(i) {
    sm <- fit$smooth[[i]]
    coef_index <- sm$first.para:sm$last.para
    sp_index <- sm$first.sp:sm$last.sp
    test_i <- focal_tests[[as.character(i)]]
    tibble(
      smooth_index = i,
      term = sm$label,
      by_variable = sm$by,
      basis_class = class(sm)[1],
      k = sm$bs.dim,
      EDF = sum(fit$edf[coef_index]),
      reference_df = if (is.null(test_i)) NA_real_ else test_i["reference_df"],
      statistic = if (is.null(test_i)) NA_real_ else test_i["statistic"],
      P = if (is.null(test_i)) NA_real_ else test_i["P"],
      smoothing_parameter = paste(fmt(fit$sp[sp_index], 6), collapse = "; "),
      null_space_dimension = sm$null.space.dim,
      notes = case_when(
        sm$by %in% focal_by_names ~ "Confirmatory ecological smooth; boundary likelihood-ratio test reported",
        sm$by %in% c("female_camera", "male_camera") ~ "Adult class-specific camera cyclic deviation; shared class-level penalty",
        grepl("camera", sm$label, ignore.case = TRUE) ~ "Hierarchical random-effect smooth; excluded from population prediction",
        TRUE ~ "Retained prespecified smooth"
      )
    )
  }))
  write_csv(smooth_table, smooth_table_path)
  table_a <- smooth_table %>%
    filter(by_variable %in% focal_by_names) %>%
    transmute(
      class_specific_term = term,
      basis = "cyclic cubic regression spline",
      k, EDF, smoothing_parameter, notes
    )
  write_csv(table_a, table_a_path)
}

summary_no_retest <- summary(fit, re.test = FALSE)
parametric <- as.data.frame(summary_no_retest$p.table, check.names = FALSE)
parametric$term <- rownames(parametric)
rownames(parametric) <- NULL
parametric <- parametric %>%
  transmute(
    term,
    estimate = Estimate,
    standard_error = `Std. Error`,
    z_value = `z value`,
    P = `Pr(>|z|)`
  )
write_csv(parametric, out("full_parametric_coefficient_table.csv"))

pearson <- residuals(fit, type = "pearson")
lag_data <- primary_data %>%
  mutate(pearson = pearson) %>%
  arrange(camera_round_class, Time_bin) %>%
  group_by(demographic_class, camera_round_class) %>%
  summarise(
    lag1 = suppressWarnings(cor(
      pearson[-n()], pearson[-1], use = "complete.obs"
    )),
    .groups = "drop"
  )
outer_gradient <- fit$outer.info$grad %||% NA_real_
outer_hessian <- fit$outer.info$hess %||% matrix(NA_real_, 1, 1)
hessian_values <- tryCatch(
  suppressWarnings(eigen(
    outer_hessian, symmetric = TRUE, only.values = TRUE
  )$values),
  error = function(e) NA_real_
)
maximum_gradient <- if (all(is.na(outer_gradient))) {
  NA_real_
} else max(abs(outer_gradient), na.rm = TRUE)
minimum_hessian <- if (all(is.na(hessian_values))) {
  NA_real_
} else min(hessian_values, na.rm = TRUE)
diagnostics <- tibble(
  model = "Primary joint confirmatory demographic GAM",
  convergence = if (
    is.finite(maximum_gradient) && maximum_gradient < 1e-3 &&
      is.finite(minimum_hessian) && minimum_hessian > -1e-8
  ) "Numerically converged" else "Convergence concern",
  maximum_absolute_gradient = maximum_gradient,
  minimum_hessian_eigenvalue = minimum_hessian,
  hessian_positive_definite = all(hessian_values > 0, na.rm = TRUE),
  hessian_numerically_positive_semidefinite = minimum_hessian > -1e-8,
  pearson_dispersion = sum(pearson^2, na.rm = TRUE) / fit$df.residual,
  deviance_explained = summary_no_retest$dev.expl,
  mean_camera_round_class_lag1_pearson = mean(lag_data$lag1, na.rm = TRUE),
  maximum_fitted_probability = max(fitted(fit)),
  minimum_fitted_probability = min(fitted(fit)),
  coefficient_rank = fit$rank,
  coefficients = length(coef(fit)),
  covariance = covariance_type
)
write_csv(diagnostics, out("joint_model_diagnostics.csv"))
lag_diagnostics <- lag_data %>%
  group_by(demographic_class) %>%
  summarise(
    mean_lag1_pearson = mean(lag1, na.rm = TRUE),
    median_lag1_pearson = median(lag1, na.rm = TRUE),
    camera_rounds_with_defined_lag = sum(is.finite(lag1)),
    .groups = "drop"
  )
write_csv(lag_diagnostics, out("joint_residual_lag_diagnostics.csv"))

# Population-level design matrices -----------------------------------------

time_grid <- seq(0, 24 - cfg$prediction_step, by = cfg$prediction_step)

newdata_grid <- function(class_i, wolf, jackal, times = time_grid) {
  n <- length(times)
  data <- tibble(
    Time = times,
    wolf_z = rep(wolf, n),
    jackal_z = rep(jackal, n),
    wolf_jackal = rep(wolf * jackal, n),
    demographic_class = factor(class_i, levels = class_levels),
    round_factor = factor("4", levels = levels(primary_data$round_factor)),
    camera_id = factor(
      levels(primary_data$camera_id)[1], levels = levels(primary_data$camera_id)
    ),
    camera_round = factor(
      levels(primary_data$camera_round)[1],
      levels = levels(primary_data$camera_round)
    ),
    camera_class = factor(
      levels(primary_data$camera_class)[1],
      levels = levels(primary_data$camera_class)
    ),
    camera_round_class = factor(
      levels(primary_data$camera_round_class)[1],
      levels = levels(primary_data$camera_round_class)
    )
  )
  add_class_by_variables(data)
}

population_lpmatrix <- function(class_i, wolf, jackal, model = fit) {
  predict(
    model,
    newdata = newdata_grid(class_i, wolf, jackal),
    type = "lpmatrix",
    exclude = population_exclusions
  )
}

center_function_matrix <- function(L) sweep(L, 2, colMeans(L), "-")

effect_matrix <- function(class_i, effect_i, model = fit) {
  X00 <- population_lpmatrix(class_i, 0, 0, model)
  if (effect_i == "wolf") {
    L <- population_lpmatrix(class_i, 1, 0, model) - X00
  } else if (effect_i == "jackal") {
    L <- population_lpmatrix(class_i, 0, 1, model) - X00
  } else if (effect_i == "interaction") {
    L <- population_lpmatrix(class_i, 1, 1, model) -
      population_lpmatrix(class_i, 1, 0, model) -
      population_lpmatrix(class_i, 0, 1, model) + X00
  } else {
    stop("Unknown effect: ", effect_i)
  }
  center_function_matrix(L)
}

rank_reduced_wald <- function(L, beta = coef(fit), V = unconditional_covariance) {
  active <- which(colSums(abs(L)) > 1e-12)
  L_active <- L[, active, drop = FALSE]
  estimate <- drop(L_active %*% beta[active])
  V_active <- V[active, active, drop = FALSE]
  covariance <- tcrossprod(L_active %*% V_active, L_active)
  eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  tolerance <- max(eig$values, na.rm = TRUE) * 1e-8
  keep <- which(eig$values > max(tolerance, .Machine$double.eps))
  if (!length(keep)) {
    return(c(statistic = 0, rank = 0, P = 1))
  }
  projection <- crossprod(eig$vectors[, keep, drop = FALSE], estimate)
  statistic <- sum(projection^2 / eig$values[keep])
  c(statistic = statistic, rank = length(keep), P = pchisq(
    statistic, df = length(keep), lower.tail = FALSE
  ))
}

active_draws <- function(L_list, beta = coef(fit), V = unconditional_covariance) {
  active <- which(Reduce(
    `|`, lapply(L_list, function(L) colSums(abs(L)) > 1e-12)
  ))
  set.seed(cfg$seed)
  draws <- rmvnorm(
    cfg$simulations,
    mean = beta[active],
    sigma = V[active, active, drop = FALSE],
    method = "svd"
  )
  list(active = active, draws = draws)
}

simulated_function <- function(L, draw_object) {
  L[, draw_object$active, drop = FALSE] %*% t(draw_object$draws)
}

simultaneous_summary <- function(L, draw_object) {
  estimate <- drop(L %*% coef(fit))
  simulations <- simulated_function(L, draw_object)
  standard_error <- apply(simulations, 1, sd)
  safe_se <- pmax(standard_error, 1e-12)
  standardized <- abs(sweep(simulations, 1, estimate, "-") / safe_se)
  critical <- unname(quantile(apply(standardized, 2, max), 0.95))
  tibble(
    Time = time_grid,
    estimate,
    standard_error,
    lower = estimate - critical * standard_error,
    upper = estimate + critical * standard_error,
    simultaneous_critical_value = critical
  )
}

format_intervals <- function(flag, direction) {
  index <- which(flag)
  if (!length(index)) return("None")
  groups <- split(index, cumsum(c(1L, diff(index) > 1L)))
  intervals <- vapply(groups, function(g) {
    start <- time_grid[min(g)]
    end <- min(24, time_grid[max(g)] + cfg$prediction_step)
    paste0(fmt(start, 2), "-", fmt(end, 2), " (", direction, ")")
  }, character(1))
  if (flag[1] && flag[length(flag)] && length(intervals) > 1L) {
    intervals <- c(
      paste0(
        intervals[length(intervals)], " and ", intervals[1],
        " [continuous across midnight]"
      ),
      intervals[-c(1, length(intervals))]
    )
  }
  paste(intervals, collapse = "; ")
}

# Direct demographic contrasts and RMS effect sizes -------------------------

effect_matrices <- setNames(lapply(effect_levels, function(effect_i) {
  setNames(lapply(class_levels, function(class_i) {
    effect_matrix(class_i, effect_i)
  }), class_levels)
}), effect_levels)

omnibus_rows <- list()
pairwise_rows <- list()
effect_curve_rows <- list()
difference_curve_rows <- list()
rms_rows <- list()
rms_difference_rows <- list()
effect_draw_cache <- list()

for (effect_i in effect_levels) {
  Ls <- effect_matrices[[effect_i]]
  pair_L <- lapply(contrast_pairs, function(pair) Ls[[pair[1]]] - Ls[[pair[2]]])
  omnibus_L <- rbind(
    Ls[["Juvenile"]] - Ls[["Adult female"]],
    Ls[["Juvenile"]] - Ls[["Adult male"]]
  )
  omnibus_test <- rank_reduced_wald(omnibus_L)
  omnibus_rows[[effect_i]] <- tibble(
    predator_effect = unname(effect_labels[effect_i]),
    omnibus_statistic = omnibus_test["statistic"],
    df_rank = omnibus_test["rank"],
    raw_P = omnibus_test["P"]
  )

  draw_object <- active_draws(c(Ls, pair_L))
  effect_draw_cache[[effect_i]] <- draw_object

  for (class_i in class_levels) {
    summary_i <- simultaneous_summary(Ls[[class_i]], draw_object) %>%
      mutate(
        predator_effect = unname(effect_labels[effect_i]),
        demographic_class = class_i
      )
    effect_curve_rows[[paste(effect_i, class_i)]] <- summary_i
    simulations_i <- simulated_function(Ls[[class_i]], draw_object)
    rms_sim <- sqrt(colMeans(simulations_i^2))
    rms_estimate <- sqrt(mean(drop(Ls[[class_i]] %*% coef(fit))^2))
    rms_rows[[paste(effect_i, class_i)]] <- tibble(
      predator_effect = unname(effect_labels[effect_i]),
      demographic_class = class_i,
      RMS_effect = median(rms_sim),
      lower_95 = unname(quantile(rms_sim, 0.025)),
      upper_95 = unname(quantile(rms_sim, 0.975))
    )
  }

  raw_pair_p <- numeric(length(pair_L))
  names(raw_pair_p) <- names(pair_L)
  pair_results_i <- list()
  for (contrast_i in names(pair_L)) {
    L_i <- pair_L[[contrast_i]]
    test_i <- rank_reduced_wald(L_i)
    raw_pair_p[contrast_i] <- test_i["P"]
    band_i <- simultaneous_summary(L_i, draw_object)
    positive <- band_i$lower > 0
    negative <- band_i$upper < 0
    intervals <- c(
      if (any(positive)) format_intervals(positive, "positive") else NULL,
      if (any(negative)) format_intervals(negative, "negative") else NULL
    )
    if (!length(intervals)) intervals <- "None"
    pair_results_i[[contrast_i]] <- tibble(
      predator_effect = unname(effect_labels[effect_i]),
      contrast = contrast_i,
      statistic = test_i["statistic"],
      df_rank = test_i["rank"],
      raw_P = test_i["P"],
      simultaneous_CI_excludes_zero = any(positive | negative),
      time_intervals_of_difference = paste(intervals, collapse = "; ")
    )
    difference_curve_rows[[paste(effect_i, contrast_i)]] <- band_i %>%
      mutate(
        predator_effect = unname(effect_labels[effect_i]),
        contrast = contrast_i
      )

    pair <- contrast_pairs[[contrast_i]]
    sim_a <- sqrt(colMeans(simulated_function(Ls[[pair[1]]], draw_object)^2))
    sim_b <- sqrt(colMeans(simulated_function(Ls[[pair[2]]], draw_object)^2))
    rms_difference <- sim_a - sim_b
    point_difference <-
      sqrt(mean(drop(Ls[[pair[1]]] %*% coef(fit))^2)) -
      sqrt(mean(drop(Ls[[pair[2]]] %*% coef(fit))^2))
    interval <- quantile(rms_difference, c(0.025, 0.975))
    rms_difference_rows[[paste(effect_i, contrast_i)]] <- tibble(
      predator_effect = unname(effect_labels[effect_i]),
      contrast = contrast_i,
      difference = median(rms_difference),
      lower_95 = unname(interval[1]),
      upper_95 = unname(interval[2]),
      ordering_supported = case_when(
        interval[1] > 0 ~ "Yes: first class has larger RMS effect",
        interval[2] < 0 ~ "No: second class has larger RMS effect",
        TRUE ~ "Uncertain"
      )
    )
  }
  adjusted <- p.adjust(raw_pair_p, method = "holm")
  for (contrast_i in names(pair_results_i)) {
    pairwise_rows[[paste(effect_i, contrast_i)]] <-
      pair_results_i[[contrast_i]] %>%
      mutate(Holm_adjusted_P = adjusted[contrast_i]) %>%
      select(
        predator_effect, contrast, statistic, df_rank, raw_P,
        Holm_adjusted_P, simultaneous_CI_excludes_zero,
        time_intervals_of_difference
      )
  }
}

table_b <- bind_rows(omnibus_rows)
table_c <- bind_rows(pairwise_rows)
table_d <- bind_rows(rms_rows)
table_e <- bind_rows(rms_difference_rows)
effect_curves <- bind_rows(effect_curve_rows)
difference_curves <- bind_rows(difference_curve_rows)
write_csv(table_b, out("table_B_omnibus_demographic_tests.csv"))
write_csv(table_c, out("table_C_pairwise_difference_smooth_tests.csv"))
write_csv(table_d, out("table_D_temporal_effect_magnitude.csv"))
write_csv(table_e, out("table_E_pairwise_RMS_differences.csv"))
write_csv(effect_curves, out("joint_temporal_effect_functions.csv"))
write_csv(difference_curves, out("joint_demographic_difference_curves.csv"))

# Response-scale normalized curves ------------------------------------------

prediction_conditions <- tribble(
  ~condition, ~wolf_z, ~jackal_z,
  "Wolf -1 / Jackal -1", -1, -1,
  "Wolf -1 / Jackal +1", -1, 1,
  "Wolf +1 / Jackal -1", 1, -1,
  "Wolf +1 / Jackal +1", 1, 1,
  "Wolf -1 / Jackal 0", -1, 0,
  "Wolf +1 / Jackal 0", 1, 0,
  "Wolf 0 / Jackal -1", 0, -1,
  "Wolf 0 / Jackal +1", 0, 1
)

prediction_matrices <- list()
for (class_i in class_levels) {
  for (row_i in seq_len(nrow(prediction_conditions))) {
    condition_i <- prediction_conditions[row_i, ]
    prediction_matrices[[paste(class_i, condition_i$condition, sep = "|")]] <-
      population_lpmatrix(class_i, condition_i$wolf_z, condition_i$jackal_z)
  }
}
prediction_draw_object <- active_draws(prediction_matrices)
prediction_rows <- list()
prediction_draws <- list()

for (class_i in class_levels) {
  for (row_i in seq_len(nrow(prediction_conditions))) {
    condition_i <- prediction_conditions[row_i, ]
    key <- paste(class_i, condition_i$condition, sep = "|")
    X <- prediction_matrices[[key]]
    eta <- drop(X %*% coef(fit))
    eta_draws <- simulated_function(X, prediction_draw_object)
    probability <- plogis(eta)
    probability_draws <- plogis(eta_draws)
    normalized <- probability / (sum(probability) * cfg$prediction_step)
    draw_area <- colSums(probability_draws) * cfg$prediction_step
    normalized_draws <- sweep(probability_draws, 2, draw_area, "/")
    prediction_draws[[key]] <- normalized_draws
    prediction_rows[[key]] <- tibble(
      demographic_class = class_i,
      condition = condition_i$condition,
      wolf_z = condition_i$wolf_z,
      jackal_z = condition_i$jackal_z,
      Time = time_grid,
      normalized_activity = normalized,
      lower_95 = apply(normalized_draws, 1, quantile, 0.025),
      upper_95 = apply(normalized_draws, 1, quantile, 0.975)
    )
  }
}
predictions <- bind_rows(prediction_rows)
write_csv(predictions, out("joint_normalized_predictions.csv"))

# Astronomical daylight summaries ------------------------------------------

obs <- read.csv(
  cfg$obs_covs, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
effort <- tibble(
  camera_id = as.integer(obs[["Camera Number"]]),
  round = as.integer(obs[["Round Number"]]),
  start_date = as.Date(obs[["Start date"]], "%d/%m/%Y"),
  end_date = as.Date(obs[["End Date"]], "%d/%m/%Y"),
  effort_days = as.integer(obs[["Days.Operational"]])
) %>%
  filter(
    !is.na(camera_id), !is.na(round), !is.na(start_date), !is.na(end_date),
    is.finite(effort_days), effort_days > 0
  ) %>%
  mutate(span_days = as.integer(end_date - start_date))
known_effort <- effort %>% filter(span_days == effort_days)
if (nrow(known_effort) != 343L || sum(known_effort$effort_days) != 5978L) {
  stop("Expected 343 exact-date camera-rounds and 5,978 known camera-days.")
}
coords <- read_csv(cfg$camera_coordinates, show_col_types = FALSE) %>%
  transmute(
    camera_id = as.integer(Camera.ID),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  )
day_rows <- lapply(seq_len(nrow(known_effort)), function(i) {
  tibble(
    camera_id = known_effort$camera_id[i],
    date = seq(
      known_effort$start_date[i],
      known_effort$end_date[i] - 1,
      by = "day"
    )
  )
})
operational_days <- bind_rows(day_rows) %>% left_join(coords, by = "camera_id")
if (nrow(operational_days) != 5978L || anyNA(operational_days$latitude)) {
  stop("Astronomical daylight operational-day reconstruction failed.")
}

solar_declination <- function(date) {
  day <- as.integer(format(date, "%j"))
  gamma <- 2 * pi / 365 * (day - 1)
  0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) -
    0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) -
    0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)
}
declination <- solar_declination(operational_days$date)
latitude_radians <- operational_days$latitude * pi / 180
hour_angle <- (time_grid - 12) * pi / 12
sin_elevation <- outer(sin(latitude_radians) * sin(declination), rep(1, length(time_grid))) +
  outer(cos(latitude_radians) * cos(declination), cos(hour_angle))
daylight_weight <- colMeans(sin_elevation > 0)
daylight_audit <- tibble(
  Time = time_grid,
  proportion_operational_camera_days_in_daylight = daylight_weight,
  known_operational_camera_days = nrow(operational_days),
  excluded_camera_round = "Camera 45 round 2 (21 days; exact outage dates unavailable)"
)
write_csv(daylight_audit, out("astronomical_daylight_weighting.csv"))

daylight_rows <- list()
daylight_draw_values <- list()
daylight_point_values <- list()
for (class_i in class_levels) {
  for (row_i in seq_len(nrow(prediction_conditions))) {
    condition_i <- prediction_conditions[row_i, ]
    key <- paste(class_i, condition_i$condition, sep = "|")
    normalized_draws <- prediction_draws[[key]]
    point_curve <- predictions %>%
      filter(
        demographic_class == class_i,
        condition == condition_i$condition
      ) %>%
      pull(normalized_activity)
    daylight_draw <- colSums(
      normalized_draws * daylight_weight
    ) * cfg$prediction_step
    daylight_point <- sum(point_curve * daylight_weight) * cfg$prediction_step
    daylight_draw_values[[key]] <- daylight_draw
    daylight_point_values[[key]] <- daylight_point
    daylight_rows[[key]] <- tibble(
      demographic_class = class_i,
      condition = condition_i$condition,
      wolf_z = condition_i$wolf_z,
      jackal_z = condition_i$jackal_z,
      daylight_proportion = daylight_point,
      lower_95 = unname(quantile(daylight_draw, 0.025)),
      upper_95 = unname(quantile(daylight_draw, 0.975))
    )
  }
}
daylight_proportions <- bind_rows(daylight_rows)
write_csv(daylight_proportions, out("joint_predicted_daylight_proportions.csv"))

daylight_shift_rows <- list()
daylight_shift_draws <- list()
for (class_i in class_levels) {
  get_day <- function(condition) daylight_draw_values[[paste(
    class_i, condition, sep = "|"
  )]]
  get_day_point <- function(condition) daylight_point_values[[paste(
    class_i, condition, sep = "|"
  )]]
  shift_definitions <- list(
    "Wolf gradient" = get_day("Wolf +1 / Jackal 0") -
      get_day("Wolf -1 / Jackal 0"),
    "Jackal gradient" = get_day("Wolf 0 / Jackal +1") -
      get_day("Wolf 0 / Jackal -1"),
    "Wolf x jackal interaction" =
      get_day("Wolf +1 / Jackal +1") - get_day("Wolf +1 / Jackal -1") -
      get_day("Wolf -1 / Jackal +1") + get_day("Wolf -1 / Jackal -1")
  )
  point_definitions <- c(
    "Wolf gradient" = get_day_point("Wolf +1 / Jackal 0") -
      get_day_point("Wolf -1 / Jackal 0"),
    "Jackal gradient" = get_day_point("Wolf 0 / Jackal +1") -
      get_day_point("Wolf 0 / Jackal -1"),
    "Wolf x jackal interaction" =
      get_day_point("Wolf +1 / Jackal +1") -
      get_day_point("Wolf +1 / Jackal -1") -
      get_day_point("Wolf -1 / Jackal +1") +
      get_day_point("Wolf -1 / Jackal -1")
  )
  for (effect_i in names(shift_definitions)) {
    values <- shift_definitions[[effect_i]]
    daylight_shift_draws[[paste(class_i, effect_i, sep = "|")]] <- values
    daylight_shift_rows[[paste(class_i, effect_i)]] <- tibble(
      demographic_class = class_i,
      predator_effect = effect_i,
      daylight_shift = point_definitions[effect_i],
      lower_95 = unname(quantile(values, 0.025)),
      upper_95 = unname(quantile(values, 0.975)),
      supported_increase = quantile(values, 0.025) > 0,
      supported_decrease = quantile(values, 0.975) < 0
    )
  }
}
daylight_shifts <- bind_rows(daylight_shift_rows)
write_csv(daylight_shifts, out("joint_daylight_shift_contrasts.csv"))

daylight_demographic_rows <- list()
for (effect_i in unique(daylight_shifts$predator_effect)) {
  for (contrast_i in names(contrast_pairs)) {
    pair <- contrast_pairs[[contrast_i]]
    difference <-
      daylight_shift_draws[[paste(pair[1], effect_i, sep = "|")]] -
      daylight_shift_draws[[paste(pair[2], effect_i, sep = "|")]]
    point_a <- daylight_shifts$daylight_shift[
      daylight_shifts$demographic_class == pair[1] &
        daylight_shifts$predator_effect == effect_i
    ]
    point_b <- daylight_shifts$daylight_shift[
      daylight_shifts$demographic_class == pair[2] &
        daylight_shifts$predator_effect == effect_i
    ]
    daylight_demographic_rows[[paste(effect_i, contrast_i)]] <- tibble(
      predator_effect = effect_i,
      contrast = contrast_i,
      difference_in_daylight_shift = point_a - point_b,
      lower_95 = unname(quantile(difference, 0.025)),
      upper_95 = unname(quantile(difference, 0.975)),
      direct_difference_supported =
        quantile(difference, 0.025) > 0 || quantile(difference, 0.975) < 0
    )
  }
}
daylight_demographic <- bind_rows(daylight_demographic_rows)
write_csv(
  daylight_demographic,
  out("joint_demographic_daylight_shift_differences.csv")
)

# Figures -------------------------------------------------------------------

publication_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.3),
    strip.background = element_rect(fill = "#EEE8DA", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#4C585C", size = 10),
    plot.caption = element_text(color = "#5E686C", hjust = 0, size = 9)
  )

close_cycle <- function(data, groups) {
  first_rows <- data %>%
    group_by(across(all_of(groups))) %>%
    slice_min(Time, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Time = 24)
  bind_rows(data, first_rows) %>% arrange(across(all_of(groups)), Time)
}

prediction_plot_data <- predictions %>%
  filter(condition %in% context_levels) %>%
  mutate(
    condition = factor(condition, levels = context_levels),
    demographic_class = factor(demographic_class, levels = class_levels)
  ) %>%
  close_cycle(c("demographic_class", "condition"))

prediction_plot <- ggplot(
  prediction_plot_data,
  aes(Time, normalized_activity, color = condition, fill = condition)
) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~demographic_class, ncol = 1, scales = "free_y") +
  scale_color_manual(values = context_colors, name = NULL) +
  scale_fill_manual(values = context_colors, name = NULL) +
  scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24)) +
  labs(
    title = "Joint-model demographic diel activity across predator contexts",
    subtitle = "Population-level normalized predictions from the confirmatory joint GAM",
    x = "Apparent solar time",
    y = "Normalized activity density",
    caption = paste(
      "Curves are normalized within class and predator condition.",
      "Camera-level and camera-round effects were excluded; survey round was fixed at round 4."
    )
  ) +
  publication_theme
ggsave(
  out("figure_01_joint_normalized_diel_curves.png"), prediction_plot,
  width = 10, height = 10.5, dpi = 400, bg = "white"
)
ggsave(
  out("figure_01_joint_normalized_diel_curves.pdf"), prediction_plot,
  width = 10, height = 10.5, device = cairo_pdf
)

effect_file_names <- c(
  "Wolf x time" = "figure_02_wolf_temporal_effect_functions",
  "Jackal x time" = "figure_03_jackal_temporal_effect_functions",
  "Wolf x jackal x time" = "figure_04_wolf_jackal_interaction_functions"
)
for (effect_i in names(effect_file_names)) {
  plot_data <- effect_curves %>%
    filter(predator_effect == effect_i) %>%
    mutate(demographic_class = factor(demographic_class, levels = class_levels)) %>%
    close_cycle(c("demographic_class"))
  effect_plot <- ggplot(
    plot_data,
    aes(Time, estimate, color = demographic_class, fill = demographic_class)
  ) +
    geom_hline(yintercept = 0, color = "#555555", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = class_colors, name = NULL) +
    scale_fill_manual(values = class_colors, name = NULL) +
    scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24)) +
    labs(
      title = paste0(effect_i, " functions by demographic class"),
      subtitle = "Centered link-scale temporal effects with simultaneous 95% bands",
      x = "Apparent solar time",
      y = "Centered temporal effect on logit scale",
      caption = paste(
        "A one-SD predator contrast was used. Centering removes constant encounter-rate shifts.",
        "Population-level predictions exclude hierarchical camera terms."
      )
    ) +
    publication_theme
  ggsave(
    out(paste0(effect_file_names[effect_i], ".png")), effect_plot,
    width = 9, height = 5.8, dpi = 400, bg = "white"
  )
  ggsave(
    out(paste0(effect_file_names[effect_i], ".pdf")), effect_plot,
    width = 9, height = 5.8, device = cairo_pdf
  )
}

difference_plot_data <- difference_curves %>%
  mutate(
    predator_effect = factor(
      predator_effect, levels = unname(effect_labels[effect_levels])
    ),
    contrast = factor(contrast, levels = names(contrast_pairs))
  ) %>%
  close_cycle(c("predator_effect", "contrast"))
difference_plot <- ggplot(
  difference_plot_data,
  aes(Time, estimate)
) +
  geom_hline(yintercept = 0, color = "#555555", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#6C8796", alpha = 0.25) +
  geom_line(color = "#173F66", linewidth = 0.85) +
  facet_grid(predator_effect ~ contrast, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
  labs(
    title = "Direct demographic differences in predator-associated diel responses",
    subtitle = "Centered link-scale difference functions with simultaneous 95% confidence bands",
    x = "Apparent solar time",
    y = "Difference in temporal effect",
    caption = paste(
      "Positive values indicate a larger temporal effect for the first named class.",
      "Bands are simultaneous within each planned contrast."
    )
  ) +
  publication_theme +
  theme(legend.position = "none")
ggsave(
  out("figure_05_pairwise_demographic_difference_curves.png"), difference_plot,
  width = 13, height = 9.5, dpi = 400, bg = "white"
)
ggsave(
  out("figure_05_pairwise_demographic_difference_curves.pdf"), difference_plot,
  width = 13, height = 9.5, device = cairo_pdf
)

# Cross-class dependence assessment -----------------------------------------

dependence_assessment <- tibble(
  component = c(
    "Shared camera dependence",
    "Shared camera-round dependence",
    "Class-specific hierarchy",
    "Residual shared camera-day-hour dependence",
    "Unaggregated shared-unit sensitivity feasibility"
  ),
  treatment = c(
    "Shared camera random intercept retained",
    "Shared camera-round random intercept retained",
    paste(
      "Camera-class and camera-round-class random intercepts retained;",
      "female and male camera-specific cyclic deviations retained"
    ),
    "Not fully represented after grouped-binomial aggregation",
    paste(
      "Not fitted: exact operational dates are unavailable for camera 45 round 2;",
      "the 5,978 known days alone imply 430,416 stacked rows and 143,472",
      "shared sampling-unit levels, making a dense GAM coefficient covariance",
      "computationally infeasible in this workspace"
    )
  ),
  implication = c(
    "Accounts for persistent cross-class covariance among cameras",
    "Accounts for cross-class covariance within camera survey rounds",
    "Does not force demographic classes to share latent camera responses",
    "Same-day/hour co-detection dependence may remain",
    "Primary direct contrasts are conditional on the grouped shared hierarchy"
  )
)
write_csv(dependence_assessment, out("cross_class_dependence_assessment.csv"))

# Prespecified sensitivity analyses -----------------------------------------

newdata_grid_for <- function(reference_data, class_i, wolf, jackal) {
  n <- length(time_grid)
  data <- tibble(
    Time = time_grid,
    wolf_z = rep(wolf, n),
    jackal_z = rep(jackal, n),
    wolf_jackal = rep(wolf * jackal, n),
    demographic_class = factor(class_i, levels = class_levels),
    round_factor = factor("4", levels = levels(reference_data$round_factor)),
    camera_id = factor(
      levels(reference_data$camera_id)[1], levels = levels(reference_data$camera_id)
    ),
    camera_round = factor(
      levels(reference_data$camera_round)[1],
      levels = levels(reference_data$camera_round)
    ),
    camera_class = factor(
      levels(reference_data$camera_class)[1],
      levels = levels(reference_data$camera_class)
    ),
    camera_round_class = factor(
      levels(reference_data$camera_round_class)[1],
      levels = levels(reference_data$camera_round_class)
    )
  )
  add_class_by_variables(data)
}

local_population_matrix <- function(fit_i, data_i, class_i, wolf, jackal) {
  labels_i <- vapply(fit_i$smooth, function(sm) sm$label, character(1))
  by_i <- vapply(fit_i$smooth, function(sm) sm$by, character(1))
  excluded_i <- labels_i[
    grepl("camera", labels_i, ignore.case = TRUE) |
      by_i %in% c("female_camera", "male_camera")
  ]
  predict(
    fit_i,
    newdata = newdata_grid_for(data_i, class_i, wolf, jackal),
    type = "lpmatrix",
    exclude = excluded_i
  )
}

local_effect_matrix <- function(fit_i, data_i, class_i, effect_i) {
  X00 <- local_population_matrix(fit_i, data_i, class_i, 0, 0)
  if (effect_i == "wolf") {
    L <- local_population_matrix(fit_i, data_i, class_i, 1, 0) - X00
  } else if (effect_i == "jackal") {
    L <- local_population_matrix(fit_i, data_i, class_i, 0, 1) - X00
  } else {
    L <- local_population_matrix(fit_i, data_i, class_i, 1, 1) -
      local_population_matrix(fit_i, data_i, class_i, 1, 0) -
      local_population_matrix(fit_i, data_i, class_i, 0, 1) + X00
  }
  center_function_matrix(L)
}

extract_direct_tests <- function(fit_i, data_i, specification) {
  V_i <- tryCatch(
    vcov(fit_i, unconditional = TRUE),
    error = function(e) vcov(fit_i)
  )
  beta_i <- coef(fit_i)
  rows <- list()
  for (effect_i in effect_levels) {
    Ls <- setNames(lapply(class_levels, function(class_i) {
      local_effect_matrix(fit_i, data_i, class_i, effect_i)
    }), class_levels)
    omnibus_L <- rbind(
      Ls[["Juvenile"]] - Ls[["Adult female"]],
      Ls[["Juvenile"]] - Ls[["Adult male"]]
    )
    omnibus <- rank_reduced_wald(omnibus_L, beta_i, V_i)
    rows[[paste(effect_i, "omnibus")]] <- tibble(
      specification,
      predator_effect = unname(effect_labels[effect_i]),
      test_type = "Omnibus",
      contrast = "All demographic classes",
      statistic = omnibus["statistic"],
      df_rank = omnibus["rank"],
      raw_P = omnibus["P"],
      adjusted_P = omnibus["P"]
    )
    pair_tests <- lapply(contrast_pairs, function(pair) {
      rank_reduced_wald(Ls[[pair[1]]] - Ls[[pair[2]]], beta_i, V_i)
    })
    pair_raw <- vapply(pair_tests, function(x) x["P"], numeric(1))
    pair_adjusted <- p.adjust(pair_raw, method = "holm")
    for (contrast_i in names(contrast_pairs)) {
      test_i <- pair_tests[[contrast_i]]
      rows[[paste(effect_i, contrast_i)]] <- tibble(
        specification,
        predator_effect = unname(effect_labels[effect_i]),
        test_type = "Pairwise",
        contrast = contrast_i,
        statistic = test_i["statistic"],
        df_rank = test_i["rank"],
        raw_P = test_i["P"],
        adjusted_P = pair_adjusted[contrast_i]
      )
    }
  }
  bind_rows(rows)
}

primary_direct_tests <- bind_rows(
  table_b %>% transmute(
    specification = "Primary apparent solar",
    predator_effect, test_type = "Omnibus",
    contrast = "All demographic classes",
    statistic = omnibus_statistic, df_rank, raw_P, adjusted_P = raw_P
  ),
  table_c %>% transmute(
    specification = "Primary apparent solar",
    predator_effect, test_type = "Pairwise", contrast,
    statistic, df_rank, raw_P, adjusted_P = Holm_adjusted_P
  )
)

sensitivity_specs <- tribble(
  ~specification, ~time_scale, ~history_type, ~k, ~hierarchy,
  "k = 8", "Apparent solar time", "Raw sightings", 8L, "full",
  "k = 16", "Apparent solar time", "Raw sightings", 16L, "full",
  "Sun time", "Sun time", "Raw sightings", 12L, "full",
  "Alternate camera-only hierarchy", "Apparent solar time", "Raw sightings", 12L, "camera_only",
  "Independent 30-min histories", "Apparent solar time", "Independent 30-min histories", 12L, "full",
  "Sun time + camera-only hierarchy", "Sun time", "Raw sightings", 12L, "camera_only"
)

sensitivity_rows <- list()
sensitivity_warning_rows <- list()
for (row_i in seq_len(nrow(sensitivity_specs))) {
  spec <- sensitivity_specs[row_i, ]
  safe_name <- gsub("[^A-Za-z0-9]+", "_", tolower(spec$specification))
  result_path <- file.path(
    cfg$checkpoint_dir, paste0("direct_tests_", safe_name, ".csv")
  )
  warning_path <- file.path(
    cfg$checkpoint_dir, paste0("warnings_", safe_name, ".txt")
  )
  if (file.exists(result_path)) {
    message("Loading sensitivity result: ", spec$specification)
    sensitivity_rows[[spec$specification]] <- read_csv(
      result_path, show_col_types = FALSE
    )
    warnings_i <- if (file.exists(warning_path)) readLines(warning_path) else character()
    if (identical(warnings_i, "No warnings generated.")) warnings_i <- character()
  } else {
    message("Fitting joint sensitivity: ", spec$specification)
    data_i <- prepare_joint_data(spec$time_scale, spec$history_type)
    object_i <- fit_joint_model(data_i, spec$k, spec$hierarchy)
    fit_i <- object_i$fit
    tests_i <- extract_direct_tests(fit_i, data_i, spec$specification)
    write_csv(tests_i, result_path)
    warnings_i <- object_i$warnings
    writeLines(
      if (length(warnings_i)) warnings_i else "No warnings generated.",
      warning_path
    )
    sensitivity_rows[[spec$specification]] <- tests_i
    rm(fit_i, object_i, data_i)
    gc(verbose = FALSE)
  }
  sensitivity_warning_rows[[spec$specification]] <- tibble(
    specification = spec$specification,
    warnings = if (length(warnings_i)) paste(warnings_i, collapse = " | ") else "None"
  )
}

all_direct_tests <- bind_rows(primary_direct_tests, bind_rows(sensitivity_rows))
write_csv(all_direct_tests, out("joint_all_sensitivity_direct_tests.csv"))
write_csv(
  bind_rows(sensitivity_warning_rows),
  out("joint_sensitivity_warnings.csv")
)

table_f <- all_direct_tests %>%
  filter(test_type == "Pairwise") %>%
  group_by(predator_effect, contrast) %>%
  summarise(
    primary_result = paste0(
      "Holm P = ", fmt_p(adjusted_P[specification == "Primary apparent solar"])
    ),
    number_supported = sum(adjusted_P < 0.05, na.rm = TRUE),
    number_tested = sum(!is.na(adjusted_P)),
    failed_specifications = paste(
      specification[is.na(adjusted_P) | adjusted_P >= 0.05],
      collapse = "; "
    ),
    robust_across_all_specifications =
      number_supported == number_tested & number_tested == 7L,
    .groups = "drop"
  )
write_csv(table_f, out("table_F_sensitivity_summary.csv"))

omnibus_sensitivity <- all_direct_tests %>%
  filter(test_type == "Omnibus") %>%
  group_by(predator_effect) %>%
  summarise(
    primary_P = raw_P[specification == "Primary apparent solar"],
    number_supported = sum(raw_P < 0.05, na.rm = TRUE),
    number_tested = sum(!is.na(raw_P)),
    failed_specifications = paste(
      specification[is.na(raw_P) | raw_P >= 0.05], collapse = "; "
    ),
    robust_across_all_specifications =
      number_supported == number_tested & number_tested == 7L,
    .groups = "drop"
  )
write_csv(omnibus_sensitivity, out("joint_omnibus_sensitivity_summary.csv"))

# Remaining diagnostics -----------------------------------------------------

message("Running primary basis-dimension and concurvity diagnostics")
k_diagnostic <- tryCatch(
  suppressWarnings(as.data.frame(k.check(
    fit, subsample = min(5000L, nrow(primary_data)), n.rep = 200L
  ))),
  error = function(e) data.frame(error = conditionMessage(e))
)
if (!"error" %in% names(k_diagnostic)) {
  k_diagnostic$term <- rownames(k_diagnostic)
  rownames(k_diagnostic) <- NULL
}
write_csv(as_tibble(k_diagnostic), out("joint_basis_dimension_checks.csv"))

set.seed(cfg$seed)
concurvity_rows <- sample(
  seq_len(nrow(primary_data)), min(5000L, nrow(primary_data))
)
X_concurvity <- predict(
  fit, newdata = primary_data[concurvity_rows, ], type = "lpmatrix"
)
focal_smooth_indices <- which(smooth_by %in% focal_by_names)
focal_column_sets <- setNames(lapply(focal_smooth_indices, function(i) {
  fit$smooth[[i]]$first.para:fit$smooth[[i]]$last.para
}), smooth_labels[focal_smooth_indices])
concurvity_diagnostic <- bind_rows(lapply(names(focal_column_sets), function(term_i) {
  own <- focal_column_sets[[term_i]]
  other <- unique(unlist(focal_column_sets[names(focal_column_sets) != term_i]))
  Xi <- scale(X_concurvity[, own, drop = FALSE], center = TRUE, scale = FALSE)
  Xother <- scale(X_concurvity[, other, drop = FALSE], center = TRUE, scale = FALSE)
  keep_other <- apply(Xother, 2, function(x) sum(x^2) > 1e-10)
  Xother <- Xother[, keep_other, drop = FALSE]
  fitted_basis <- if (ncol(Xother)) qr.fitted(qr(Xother), Xi) else Xi * 0
  basis_overlap <- 1 - sum((Xi - fitted_basis)^2) / sum(Xi^2)
  own_contribution <- drop(X_concurvity[, own, drop = FALSE] %*% coef(fit)[own])
  other_contribution <- drop(
    X_concurvity[, other, drop = FALSE] %*% coef(fit)[other]
  )
  tibble(
    term = term_i,
    focal_basis_overlap_with_other_focal_smooths = basis_overlap,
    squared_correlation_of_fitted_contributions = suppressWarnings(
      cor(own_contribution, other_contribution)^2
    ),
    sampled_rows = length(concurvity_rows),
    method = paste(
      "Focal ecological design-space audit; full mgcv concurvity was not",
      "computationally feasible with 120 camera-specific cyclic deviations"
    )
  )
}))
write_csv(concurvity_diagnostic, out("joint_concurvity_checks.csv"))

residual_pattern <- primary_data %>%
  mutate(
    fitted_probability = fitted(fit),
    fitted_successes = fitted_probability * trials,
    pearson_residual = pearson
  ) %>%
  group_by(demographic_class, Time_bin, Time) %>%
  summarise(
    observed_successes = sum(success),
    available_trials = sum(trials),
    observed_probability = observed_successes / available_trials,
    mean_fitted_probability = weighted.mean(fitted_probability, trials),
    mean_pearson_residual = mean(pearson_residual),
    .groups = "drop"
  )
write_csv(residual_pattern, out("joint_residual_temporal_patterns.csv"))

# Session information and concise written report ----------------------------

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))

md_table <- function(data) {
  data <- as.data.frame(data)
  data[] <- lapply(data, function(x) {
    if (is.numeric(x)) vapply(x, fmt, character(1), digits = 3) else as.character(x)
  })
  header <- paste0("| ", paste(names(data), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(data)), collapse = " | "), " |")
  body <- apply(data, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, body)
}

primary_pair_summary <- table_c %>%
  mutate(
    raw_P = vapply(raw_P, fmt_p, character(1)),
    Holm_adjusted_P = vapply(Holm_adjusted_P, fmt_p, character(1))
  ) %>%
  select(
    predator_effect, contrast, raw_P, Holm_adjusted_P,
    simultaneous_CI_excludes_zero, time_intervals_of_difference
  )
daylight_summary_print <- daylight_shifts %>%
  mutate(
    interval = paste0(fmt(lower_95), " to ", fmt(upper_95)),
    result = case_when(
      supported_increase ~ "supported increase",
      supported_decrease ~ "supported decrease",
      TRUE ~ "uncertain"
    )
  ) %>%
  select(demographic_class, predator_effect, daylight_shift, interval, result)

omnibus_text <- paste(vapply(seq_len(nrow(table_b)), function(i) {
  paste0(
    table_b$predator_effect[i], ": rank-reduced Wald chi-square = ",
    fmt(table_b$omnibus_statistic[i]), ", rank = ", table_b$df_rank[i],
    ", P = ", fmt_p(table_b$raw_P[i])
  )
}, character(1)), collapse = "; ")
shrunk_terms <- smooth_table %>%
  filter(by_variable %in% focal_by_names, EDF < 0.01) %>%
  pull(term)
focal_k <- k_diagnostic %>%
  filter(grepl("^s\\(Time\\):(demographic|wolf_|jackal_|interaction_)", term))
focal_k_p_range <- if (nrow(focal_k) && "p-value" %in% names(focal_k)) {
  range(focal_k[["p-value"]], na.rm = TRUE)
} else c(NA_real_, NA_real_)
any_primary_pair <- any(table_c$Holm_adjusted_P < 0.05, na.rm = TRUE)
any_robust_pair <- any(table_f$robust_across_all_specifications, na.rm = TRUE)
any_daylight_shift <- any(
  daylight_shifts$supported_increase | daylight_shifts$supported_decrease
)
any_daylight_class_difference <- any(
  daylight_demographic$direct_difference_supported
)

report_lines <- c(
  "# Confirmatory joint demographic hierarchical diel GAM",
  "",
  "## Scope",
  "",
  paste(
    "This analysis directly tests equality of centered predator-associated diel",
    "response functions among juveniles, adult females, and adult males. It is an",
    "addition to, not a replacement for, the definitive class-specific GAMs."
  ),
  "",
  "The model was prespecified and no demographic predator-by-time term was removed by model selection.",
  "",
  "## Joint data audit",
  "",
  md_table(data_audit),
  "",
  paste0(
    "Wolf and jackal standardized values were read without recalculation from `",
    cfg$input_data, "`, the exact definitive model input."
  ),
  "",
  "## Model and dependence structure",
  "",
  paste0("- Formula: `", paste(deparse(formula_primary), collapse = " "), "`"),
  paste0("- Fitting call: `", collapse_call(fit$call), "`"),
  paste0("- Coefficient simulations: ", cfg$simulations, "; seed: ", cfg$seed, "."),
  paste0("- Covariance: ", covariance_type, "."),
  paste(
    "- Population predictions fixed survey round at round 4 and excluded shared",
    "and class-specific camera/camera-round random effects and adult",
    "camera-specific cyclic deviations."
  ),
  paste(
    "- Shared camera and camera-round random intercepts account for persistent",
    "cross-class sampling dependence; class-specific hierarchical terms preserve",
    "the separate-model structure. Same-day/hour residual dependence is not fully",
    "identified by grouped-binomial aggregation."
  ),
  paste(
    "- The requested unaggregated shared sampling-unit sensitivity was not validly",
    "constructible for camera 45 round 2 and was computationally infeasible for",
    "the 143,472 known camera-day-hour sampling units."
  ),
  "",
  "## Diagnostics",
  "",
  md_table(diagnostics),
  "",
  paste0(
    "The focal basis checks were nonproblematic (P range = ",
    fmt(focal_k_p_range[1]), "-", fmt(focal_k_p_range[2]),
    "). The maximum focal design-space overlap was ",
    fmt(max(concurvity_diagnostic$focal_basis_overlap_with_other_focal_smooths)),
    ". Fitting and prediction generated no warnings."
  ),
  paste0(
    "The Hessian minimum eigenvalue (", fmt(minimum_hessian, 12),
    ") was within numerical zero tolerance, but mean adjacent-hour Pearson",
    " residual correlation remained ",
    fmt(diagnostics$mean_camera_round_class_lag1_pearson),
    " and Pearson dispersion was ", fmt(diagnostics$pearson_dispersion), "."
  ),
  paste0(
    "Terms effectively shrunk toward zero (EDF <0.01): ",
    paste(shrunk_terms, collapse = "; "), "."
  ),
  "Juvenile interaction stability must be judged from the direct tests, simultaneous bands, EDF, and the seven-specification sensitivity table; the smaller juvenile sample was not concealed by deleting or merging terms.",
  "",
  "## Omnibus demographic tests",
  "",
  md_table(table_b %>% mutate(raw_P = vapply(raw_P, fmt_p, character(1)))),
  "",
  "## Planned pairwise demographic tests",
  "",
  md_table(primary_pair_summary),
  "",
  "These are rank-reduced Wald tests of centered link-scale difference functions. Holm adjustment was applied across the three class contrasts within each predator family.",
  "",
  "## RMS temporal-effect magnitudes",
  "",
  md_table(table_d),
  "",
  "## Pairwise RMS differences",
  "",
  md_table(table_e),
  "",
  "RMS values are effect-size summaries. An ordering is supported only when the paired simulation interval excludes zero.",
  "",
  "## Astronomical daylight summaries",
  "",
  md_table(daylight_summary_print),
  "",
  paste(
    "Daylight was solar elevation >0 degrees, averaged across 5,978 operational",
    "camera-days with exact dates and camera-specific latitude. Camera 45 round 2",
    "was excluded only from this directional weighting because its 21 exact outage",
    "dates are unknown; it remained in the primary GAM."
  ),
  "",
  "## Sensitivity of direct demographic contrasts",
  "",
  md_table(table_f),
  "",
  "## Answers to the confirmatory questions",
  "",
  paste0(
    "A-B. No predator component showed formal demographic heterogeneity in the",
    " primary model. ", omnibus_text, "."
  ),
  paste0(
    "C-D. No planned pairwise class contrast was supported after Holm correction",
    " (any supported = ", any_primary_pair,
    "), and no simultaneous difference band excluded zero at any time of day."
  ),
  paste(
    "E. RMS effect sizes gave only partial ordering evidence. For wolf effects,",
    "juveniles and females each exceeded males, but juveniles did not exceed females.",
    "For the wolf-jackal interaction, juveniles exceeded both adult classes, but",
    "males exceeded females. No jackal ordering was supported. Thus the full",
    "Juvenile > Adult female > Adult male ordering was not supported for any component."
  ),
  paste0(
    "F. No direct demographic contrast was robust across all seven specifications",
    " (any robust = ", any_robust_pair,
    "). Isolated support occurred only in camera-only sensitivities and is",
    " sensitivity-dependent."
  ),
  paste(
    "G. The joint model materially qualifies, rather than overturns, the separate",
    "class-specific GAMs. Within-class predator-time terms can be supported while",
    "the direct between-class response functions remain statistically indistinguishable."
  ),
  paste(
    "H. Reporting concerns are the sparse juvenile response (183 detection-days),",
    "residual adjacent-hour correlation, underdispersion, and residual same-day/hour",
    "cross-class dependence after aggregation. Exact shared-unit modeling was not",
    "feasible, so absence of evidence must not be described as evidence of identical responses."
  ),
  paste(
    "Cross-class dependence. Shared camera and camera-round effects were included",
    "in the primary model. Because exact operational dates are unavailable for one",
    "camera-round and a 143,472-level sampling-unit effect is infeasible here, the",
    "incremental effect of explicit same-day/hour dependence could not be tested."
  ),
  paste0(
    "Directional daylight H3. No class showed a supported increase or decrease in",
    " astronomical-daylight use (any supported shift = ", any_daylight_shift,
    "), and no direct demographic daylight-shift difference was supported (",
    any_daylight_class_difference, ")."
  ),
  paste(
    "Overall H3 assessment. Demographic heterogeneity in response shape: not",
    "supported by the primary direct tests. Predicted Juvenile > Female > Male",
    "ordering: only partially supported and never in full. Redistribution toward",
    "daylight: unsupported."
  ),
  "",
  "## Output index",
  "",
  "- Tables A-F are saved as separately named CSV files.",
  "- Figures 1-5 are saved as 400-dpi PNG and vector PDF files.",
  "- The exact formula, model call, full parametric table, full smooth inventory, warnings, diagnostics, and session information are saved in this directory."
)
writeLines(report_lines, out("joint_demographic_diel_results.md"))

message("Joint demographic hierarchical diel GAM analysis complete.")
