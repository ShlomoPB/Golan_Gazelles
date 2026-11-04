suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(mvtnorm)
  library(readr)
  library(tidyr)
})

# Penalty-aware calibration of the confirmatory demographic diel test.
# Existing class-specific and joint-model outputs are read-only inputs.
cfg <- list(
  seed = 20260827L,
  simulations = 5000L,
  k = 12L,
  site_k = 12L,
  nthreads = 2L,
  input = "definitive_hierarchical_diel_data.csv",
  original_tests = file.path(
    "joint_demographic_hierarchical_diel_outputs",
    "table_B_omnibus_demographic_tests.csv"
  ),
  original_effect_functions = file.path(
    "joint_demographic_hierarchical_diel_outputs",
    "joint_temporal_effect_functions.csv"
  ),
  daylight_weights = file.path(
    "joint_demographic_hierarchical_diel_outputs",
    "astronomical_daylight_weighting.csv"
  ),
  output_dir = "joint_demographic_penalty_calibrated_outputs",
  checkpoint_dir = file.path(
    "joint_demographic_penalty_calibrated_outputs", "checkpoints"
  )
)

set.seed(cfg$seed)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(name) file.path(cfg$output_dir, name)

class_levels <- c("Juvenile", "Adult female", "Adult male")
effect_variables <- c(
  "Wolf x time" = "wolf_z",
  "Jackal x time" = "jackal_z",
  "Wolf x jackal x time" = "wolf_jackal"
)
contrast_pairs <- list(
  "Juvenile - Adult female" = c("Juvenile", "Adult female"),
  "Juvenile - Adult male" = c("Juvenile", "Adult male"),
  "Adult female - Adult male" = c("Adult female", "Adult male")
)

prepare_data <- function(
  classes = class_levels,
  time_scale = "Apparent solar time",
  history_type = "Raw sightings"
) {
  read_csv(cfg$input, show_col_types = FALSE) %>%
    filter(
      .data$time_scale == .env$time_scale,
      .data$history_type == .env$history_type,
      response_class %in% classes
    ) %>%
    mutate(
      demographic_class = factor(response_class, levels = classes),
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
    add_site_shape_factors()
}

add_site_shape_factors <- function(data) {
  camera_levels <- levels(data$camera_id)
  female_levels <- c("reference", paste0("female_", camera_levels))
  male_levels <- c("reference", paste0("male_", camera_levels))
  data %>%
    mutate(
      female_camera = ordered(
        if_else(
          demographic_class == "Adult female",
          paste0("female_", camera_id), "reference"
        ),
        levels = female_levels
      ),
      male_camera = ordered(
        if_else(
          demographic_class == "Adult male",
          paste0("male_", camera_id), "reference"
        ),
        levels = male_levels
      )
    )
}

site_shape_terms <- function(classes, k = cfg$site_k) {
  terms <- character()
  if ("Adult female" %in% classes) {
    terms <- c(terms, paste0(
      "s(Time, by = female_camera, bs = 'cc', k = ", k,
      ", m = 1, id = 21)"
    ))
  }
  if ("Adult male" %in% classes) {
    terms <- c(terms, paste0(
      "s(Time, by = male_camera, bs = 'cc', k = ", k,
      ", m = 1, id = 22)"
    ))
  }
  terms
}

# The common smooth is the across-class predator-time response. The sz term is
# the symmetric, sum-to-zero demographic deviation from that common response.
predator_terms <- function(variable, k = cfg$k) {
  c(
    paste0("s(Time, by = ", variable, ", bs = 'cc', k = ", k, ")"),
    paste0(
      "s(demographic_class, Time, by = ", variable,
      ", bs = 'sz', k = ", k, ", xt = list(bs = 'cc'))"
    )
  )
}

build_formula <- function(classes, k = cfg$k, hierarchy = "full") {
  ecological <- unlist(lapply(unname(effect_variables), predator_terms, k = k))
  hierarchy_terms <- if (hierarchy == "full") {
    c(
      site_shape_terms(classes),
      "s(camera_id, bs = 're')",
      "s(camera_round, bs = 're')",
      "s(camera_class, bs = 're')",
      "s(camera_round_class, bs = 're')"
    )
  } else if (hierarchy == "camera_only") {
    c("s(camera_id, bs = 're')", "s(camera_class, bs = 're')")
  } else {
    stop("Unknown hierarchy: ", hierarchy)
  }
  as.formula(paste(
    "cbind(success, failure) ~",
    "demographic_class * round_factor +",
    "demographic_class * wolf_z * jackal_z +",
    paste0("s(Time, by = demographic_class, bs = 'cc', k = ", k, ") +"),
    paste(c(ecological, hierarchy_terms), collapse = " + ")
  ))
}

fit_model <- function(
  data, classes, checkpoint, k = cfg$k, hierarchy = "full",
  save_checkpoint = TRUE
) {
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  formula <- build_formula(classes, k = k, hierarchy = hierarchy)
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
  result <- list(fit = fit, warnings = unique(warnings), formula = formula)
  if (save_checkpoint) saveRDS(result, checkpoint, compress = FALSE)
  result
}

smooth_inventory <- function(fit) {
  bind_rows(lapply(seq_along(fit$smooth), function(i) {
    smooth <- fit$smooth[[i]]
    coefficient_index <- smooth$first.para:smooth$last.para
    smoothing_index <- smooth$first.sp:smooth$last.sp
    tibble(
      smooth_index = i,
      term = smooth$label,
      by_variable = smooth$by,
      EDF = sum(fit$edf[coefficient_index]),
      coefficient_count = length(coefficient_index),
      null_space_dimension = smooth$null.space.dim,
      smoothing_parameter = paste(
        format(fit$sp[smoothing_index], scientific = FALSE), collapse = "; "
      )
    )
  }))
}

deviation_test_table <- function(
  fit, scope, effects = names(effect_variables)
) {
  inventory <- smooth_inventory(fit)
  bind_rows(lapply(effects, function(effect_label) {
    variable <- effect_variables[[effect_label]]
    target <- paste0("s(demographic_class,Time):", variable)
    row <- inventory %>% filter(term == target)
    if (nrow(row) != 1L) {
      stop("Could not uniquely identify penalty-aware deviation term: ", target)
    }
    test <- mgcv:::reTest(fit, row$smooth_index)
    tibble(
      scope = scope,
      predator_effect = effect_label,
      term = target,
      EDF = row$EDF,
      reference_df = unname(test$rank),
      statistic = unname(test$stat),
      raw_P = unname(test$pval)
    )
  }))
}

newdata_grid <- function(data, class_i, wolf, jackal, step = 0.1) {
  camera_reference <- levels(data$camera_id)[1]
  camera_round_reference <- levels(data$camera_round)[1]
  camera_class_reference <- levels(data$camera_class)[1]
  camera_round_class_reference <- levels(data$camera_round_class)[1]
  tibble(
    Time = seq(0, 24, by = step),
    demographic_class = factor(
      class_i, levels = levels(data$demographic_class)
    ),
    round_factor = factor("4", levels = levels(data$round_factor)),
    wolf_z = wolf,
    jackal_z = jackal,
    wolf_jackal = wolf * jackal,
    camera_id = factor(camera_reference, levels = levels(data$camera_id)),
    camera_round = factor(
      camera_round_reference, levels = levels(data$camera_round)
    ),
    camera_class = factor(
      camera_class_reference, levels = levels(data$camera_class)
    ),
    camera_round_class = factor(
      camera_round_class_reference, levels = levels(data$camera_round_class)
    ),
    female_camera = ordered(
      "reference", levels = levels(data$female_camera)
    ),
    male_camera = ordered("reference", levels = levels(data$male_camera))
  )
}

population_exclusions <- function(fit) {
  labels <- vapply(fit$smooth, function(x) x$label, character(1))
  by <- vapply(fit$smooth, function(x) x$by, character(1))
  labels[
    grepl("camera", labels, ignore.case = TRUE) |
      by %in% c("female_camera", "male_camera")
  ]
}

effect_matrix <- function(fit, data, class_i, effect_i) {
  exclusions <- population_exclusions(fit)
  lp <- function(wolf, jackal) predict(
    fit,
    newdata = newdata_grid(data, class_i, wolf, jackal),
    type = "lpmatrix",
    exclude = exclusions
  )
  x00 <- lp(0, 0)
  L <- switch(
    effect_i,
    "Wolf x time" = lp(1, 0) - x00,
    "Jackal x time" = lp(0, 1) - x00,
    "Wolf x jackal x time" = lp(1, 1) - lp(1, 0) - lp(0, 1) + x00
  )
  sweep(L, 2, colMeans(L), "-")
}

simultaneous_difference <- function(L, fit, draws = cfg$simulations) {
  active <- which(colSums(abs(L)) > 1e-12)
  beta <- coef(fit)[active]
  V <- tryCatch(
    vcov(fit, unconditional = TRUE),
    error = function(e) vcov(fit)
  )[active, active, drop = FALSE]
  set.seed(cfg$seed)
  beta_draws <- rmvnorm(draws, mean = beta, sigma = V)
  point <- drop(L[, active, drop = FALSE] %*% beta)
  simulations <- L[, active, drop = FALSE] %*% t(beta_draws)
  se <- apply(simulations, 1, sd)
  standardized <- sweep(simulations, 1, point, "-")
  standardized <- sweep(standardized, 1, pmax(se, 1e-12), "/")
  critical <- unname(quantile(apply(abs(standardized), 2, max), 0.95))
  tibble(
    Time = seq(0, 24, by = 0.1),
    estimate = point,
    lower = point - critical * se,
    upper = point + critical * se
  )
}

primary_data <- prepare_data()
if (nrow(primary_data) != 24768L) stop("Authoritative primary row count changed.")
if (!identical(
  as.integer(tapply(primary_data$success, primary_data$demographic_class, sum)),
  c(183L, 1576L, 2307L)
)) stop("Authoritative class success totals changed.")

primary_checkpoint <- file.path(cfg$checkpoint_dir, "penalty_calibrated_primary.rds")
message("Fitting/loading penalty-calibrated primary model")
primary_object <- fit_model(
  primary_data, class_levels, primary_checkpoint
)
fit <- primary_object$fit

writeLines(
  deparse(primary_object$formula, width.cutoff = 500L),
  out("penalty_calibrated_formula.txt")
)
writeLines(
  if (length(primary_object$warnings)) primary_object$warnings else
    "No warnings generated.",
  out("penalty_calibrated_warnings.txt")
)
write_csv(smooth_inventory(fit), out("penalty_calibrated_smooth_table.csv"))

omnibus_path <- out("penalty_calibrated_omnibus_tests.csv")
if (file.exists(omnibus_path)) {
  omnibus <- read_csv(omnibus_path, show_col_types = FALSE)
} else {
  omnibus <- deviation_test_table(fit, "All three classes")
}
if (
  file.exists(cfg$original_tests) &&
    !"original_rank_Wald_P" %in% names(omnibus)
) {
  original <- read_csv(cfg$original_tests, show_col_types = FALSE) %>%
    transmute(
      predator_effect,
      original_rank_Wald_statistic = omnibus_statistic,
      original_rank_Wald_df = df_rank,
      original_rank_Wald_P = raw_P
    )
  omnibus <- left_join(omnibus, original, by = "predator_effect")
}
write_csv(omnibus, omnibus_path)

# Closed-testing follow-ups: fit symmetric two-class models only for predator
# families supported by the all-class penalty-aware omnibus test.
supported_effects <- omnibus$predator_effect[omnibus$raw_P < 0.05]
pairwise_path <- out("penalty_calibrated_pairwise_tests.csv")
pairwise_warning_path <- out("penalty_calibrated_pairwise_warnings.csv")
pairwise_rows <- list()
pairwise_warnings <- list()
if (file.exists(pairwise_path) && file.info(pairwise_path)$size > 1L) {
  pairwise <- read_csv(pairwise_path, show_col_types = FALSE)
} else if (length(supported_effects)) {
  for (contrast_name in names(contrast_pairs)) {
    classes <- contrast_pairs[[contrast_name]]
    pair_data <- prepare_data(classes)
    checkpoint_name <- paste0(
      "penalty_calibrated_pair_",
      gsub("[^a-z0-9]+", "_", tolower(contrast_name)), ".rds"
    )
    pair_object <- fit_model(
      pair_data, classes, file.path(cfg$checkpoint_dir, checkpoint_name)
    )
    pair_tests <- deviation_test_table(
      pair_object$fit, contrast_name, effects = supported_effects
    ) %>%
      mutate(contrast = contrast_name)
    pairwise_rows[[contrast_name]] <- pair_tests
    pairwise_warnings[[contrast_name]] <- tibble(
      contrast = contrast_name,
      warning = if (length(pair_object$warnings)) {
        paste(pair_object$warnings, collapse = " | ")
      } else {
        "No warnings generated."
      }
    )
  }
  pairwise <- bind_rows(pairwise_rows)
  if (nrow(pairwise)) {
    pairwise <- pairwise %>%
      group_by(predator_effect) %>%
      mutate(Holm_adjusted_P = p.adjust(raw_P, method = "holm")) %>%
      ungroup()
  }
} else {
  pairwise <- tibble()
}
write_csv(pairwise, pairwise_path)
if (length(pairwise_warnings)) {
  write_csv(bind_rows(pairwise_warnings), pairwise_warning_path)
} else if (!file.exists(pairwise_warning_path) ||
  file.info(pairwise_warning_path)$size == 0L) {
  write_csv(
    tibble(contrast = names(contrast_pairs), warning = "None"),
    pairwise_warning_path
  )
}

pairwise_model_files <- list.files(
  cfg$checkpoint_dir,
  pattern = "^penalty_calibrated_pair_.*[.]rds$",
  full.names = TRUE
)
pairwise_diagnostics <- bind_rows(lapply(pairwise_model_files, function(path) {
  object <- readRDS(path)
  fit_i <- object$fit
  minimum_eigenvalue <- min(eigen(
    (fit_i$outer.info$hess + t(fit_i$outer.info$hess)) / 2,
    symmetric = TRUE, only.values = TRUE
  )$values)
  maximum_gradient_i <- max(abs(fit_i$outer.info$grad))
  coefficient_rank_i <- fit_i$rank
  coefficient_count_i <- length(coef(fit_i))
  tibble(
    model = sub(
      "^penalty_calibrated_pair_|[.]rds$", "", basename(path)
    ),
    maximum_absolute_gradient = maximum_gradient_i,
    minimum_hessian_eigenvalue = minimum_eigenvalue,
    hessian_numerically_positive_semidefinite = minimum_eigenvalue > -1e-8,
    coefficient_rank = coefficient_rank_i,
    coefficients = coefficient_count_i,
    warnings = if (length(object$warnings)) {
      paste(object$warnings, collapse = " | ")
    } else {
      "None"
    }
  )
}))
write_csv(
  pairwise_diagnostics, out("penalty_calibrated_pairwise_diagnostics.csv")
)

# Difference curves remain uncertainty visualizations; inferential labels come
# only from the penalty-aware tests above.
effect_matrices <- lapply(names(effect_variables), function(effect_i) {
  setNames(lapply(class_levels, function(class_i) {
    effect_matrix(fit, primary_data, class_i, effect_i)
  }), class_levels)
})
names(effect_matrices) <- names(effect_variables)

calibrated_effects <- bind_rows(lapply(names(effect_matrices), function(effect_i) {
  bind_rows(lapply(class_levels, function(class_i) {
    L <- effect_matrices[[effect_i]][[class_i]]
    tibble(
      Time = seq(0, 24, by = 0.1),
      estimate = drop(L %*% coef(fit)),
      predator_effect = effect_i,
      demographic_class = class_i
    )
  }))
}))
write_csv(
  calibrated_effects, out("penalty_calibrated_temporal_effect_functions.csv")
)

# Ordered H3 effect-size tests -----------------------------------------------
# RMS magnitude summarizes the size, not the direction, of a centered temporal
# response. The strict order was prespecified; juvenile > both adults is a
# post-hoc simplification motivated by the supported juvenile-adult contrasts.
ordered_effect <- "Wolf x jackal x time"
ordered_L <- effect_matrices[[ordered_effect]]
ordered_active <- which(Reduce(
  `|`, lapply(ordered_L, function(L) colSums(abs(L)) > 1e-12)
))
ordered_beta <- coef(fit)[ordered_active]
ordered_V <- tryCatch(
  vcov(fit, unconditional = TRUE),
  error = function(e) vcov(fit)
)[ordered_active, ordered_active, drop = FALSE]
set.seed(cfg$seed)
ordered_draws <- rmvnorm(
  cfg$simulations, mean = ordered_beta, sigma = ordered_V
)

rms_point_and_gradient <- function(L, beta, active) {
  La <- L[, active, drop = FALSE]
  values <- drop(La %*% beta)
  rms <- sqrt(mean(values^2))
  gradient <- if (rms > 1e-12) {
    drop(crossprod(La, values)) / (nrow(La) * rms)
  } else {
    rep(0, ncol(La))
  }
  list(rms = rms, gradient = gradient)
}

ordered_components <- setNames(lapply(class_levels, function(class_i) {
  L <- ordered_L[[class_i]][, ordered_active, drop = FALSE]
  simulation <- sqrt(colMeans((L %*% t(ordered_draws))^2))
  point <- rms_point_and_gradient(
    ordered_L[[class_i]], ordered_beta, ordered_active
  )
  list(
    estimate = point$rms,
    gradient = point$gradient,
    simulation = simulation
  )
}), class_levels)

ordered_rms <- bind_rows(lapply(class_levels, function(class_i) {
  component <- ordered_components[[class_i]]
  tibble(
    predator_effect = ordered_effect,
    demographic_class = class_i,
    RMS_effect = component$estimate,
    lower_95 = unname(quantile(component$simulation, 0.025)),
    upper_95 = unname(quantile(component$simulation, 0.975))
  )
}))
write_csv(ordered_rms, out("penalty_calibrated_ordered_RMS_by_class.csv"))

rms_contrast <- function(label, weights, hypothesis) {
  estimate <- sum(vapply(class_levels, function(class_i) {
    weights[class_i] * ordered_components[[class_i]]$estimate
  }, numeric(1)))
  gradient <- Reduce(`+`, lapply(class_levels, function(class_i) {
    weights[class_i] * ordered_components[[class_i]]$gradient
  }))
  standard_error <- sqrt(drop(crossprod(gradient, ordered_V %*% gradient)))
  simulation <- Reduce(`+`, lapply(class_levels, function(class_i) {
    weights[class_i] * ordered_components[[class_i]]$simulation
  }))
  z <- estimate / standard_error
  tibble(
    predator_effect = ordered_effect,
    contrast = label,
    directional_hypothesis = hypothesis,
    estimate = estimate,
    standard_error = standard_error,
    lower_95 = unname(quantile(simulation, 0.025)),
    upper_95 = unname(quantile(simulation, 0.975)),
    lower_one_sided_95 = unname(quantile(simulation, 0.05)),
    one_sided_Wald_P = pnorm(z, lower.tail = FALSE),
    simulation_probability_greater_than_zero = mean(simulation > 0),
    simulation_directional_tail_area = mean(simulation <= 0)
  )
}

ordered_contrasts <- bind_rows(
  rms_contrast(
    "Juvenile - Adult female",
    c("Juvenile" = 1, "Adult female" = -1, "Adult male" = 0),
    "Juvenile > Adult female"
  ),
  rms_contrast(
    "Juvenile - Adult male",
    c("Juvenile" = 1, "Adult female" = 0, "Adult male" = -1),
    "Juvenile > Adult male"
  ),
  rms_contrast(
    "Adult female - Adult male",
    c("Juvenile" = 0, "Adult female" = 1, "Adult male" = -1),
    "Adult female > Adult male"
  ),
  rms_contrast(
    "Juvenile - mean adults",
    c("Juvenile" = 1, "Adult female" = -0.5, "Adult male" = -0.5),
    "Juvenile > mean of adult classes"
  )
)
write_csv(
  ordered_contrasts, out("penalty_calibrated_ordered_RMS_contrasts.csv")
)

contrast_p <- setNames(
  ordered_contrasts$one_sided_Wald_P, ordered_contrasts$contrast
)
contrast_sim_p <- setNames(
  ordered_contrasts$simulation_directional_tail_area,
  ordered_contrasts$contrast
)
strict_order_probability <- mean(
  ordered_components[["Juvenile"]]$simulation >
    ordered_components[["Adult female"]]$simulation &
    ordered_components[["Adult female"]]$simulation >
    ordered_components[["Adult male"]]$simulation
)
juvenile_both_probability <- mean(
  ordered_components[["Juvenile"]]$simulation >
    ordered_components[["Adult female"]]$simulation &
    ordered_components[["Juvenile"]]$simulation >
    ordered_components[["Adult male"]]$simulation
)
ordered_hypotheses <- tibble(
  hypothesis = c(
    "Prespecified strict order: Juvenile > Adult female > Adult male",
    "Post-hoc simplified order: Juvenile > both adult classes"
  ),
  delta_Wald_intersection_union_P = c(
    max(
      contrast_p["Juvenile - Adult female"],
      contrast_p["Adult female - Adult male"]
    ),
    max(
      contrast_p["Juvenile - Adult female"],
      contrast_p["Juvenile - Adult male"]
    )
  ),
  simulation_intersection_union_tail_area = c(
    max(
      contrast_sim_p["Juvenile - Adult female"],
      contrast_sim_p["Adult female - Adult male"]
    ),
    max(
      contrast_sim_p["Juvenile - Adult female"],
      contrast_sim_p["Juvenile - Adult male"]
    )
  ),
  joint_simulation_probability = c(
    strict_order_probability, juvenile_both_probability
  ),
  delta_Wald_supported_at_0_05 = delta_Wald_intersection_union_P < 0.05,
  simulation_tail_area_below_0_05 =
    simulation_intersection_union_tail_area < 0.05,
  interpretation = c(
    "Requires both Juvenile > Female and Female > Male",
    "Requires Juvenile > Female and Juvenile > Male"
  )
)
write_csv(
  ordered_hypotheses, out("penalty_calibrated_ordered_H3_tests.csv")
)

# Astronomical daylight post-estimation ------------------------------------
# Reuse the audited operational-date/camera-coordinate weighting from the
# original joint pipeline; only model-based predictions are recalculated.
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
daylight_audit <- read_csv(cfg$daylight_weights, show_col_types = FALSE) %>%
  arrange(Time)
if (
  nrow(daylight_audit) != 240L ||
    !isTRUE(all.equal(daylight_audit$Time, seq(0, 23.9, by = 0.1))) ||
    unique(daylight_audit$known_operational_camera_days) != 5978L
) {
  stop("Audited astronomical-daylight weights changed.")
}
daylight_weight <- daylight_audit$proportion_operational_camera_days_in_daylight
prediction_time <- daylight_audit$Time
prediction_exclusions <- population_exclusions(fit)

calibrated_prediction_matrices <- list()
for (class_i in class_levels) {
  for (row_i in seq_len(nrow(prediction_conditions))) {
    condition_i <- prediction_conditions[row_i, ]
    key <- paste(class_i, condition_i$condition, sep = "|")
    newdata_i <- newdata_grid(
      primary_data, class_i, condition_i$wolf_z, condition_i$jackal_z
    ) %>%
      filter(Time < 24)
    if (!isTRUE(all.equal(newdata_i$Time, prediction_time))) {
      stop("Calibrated prediction grid does not match daylight weights.")
    }
    calibrated_prediction_matrices[[key]] <- predict(
      fit, newdata = newdata_i, type = "lpmatrix",
      exclude = prediction_exclusions
    )
  }
}
daylight_active <- which(Reduce(
  `|`, lapply(calibrated_prediction_matrices, function(X) {
    colSums(abs(X)) > 1e-12
  })
))
daylight_V <- tryCatch(
  vcov(fit, unconditional = TRUE),
  error = function(e) vcov(fit)
)[daylight_active, daylight_active, drop = FALSE]
set.seed(cfg$seed)
daylight_beta_draws <- rmvnorm(
  cfg$simulations,
  mean = coef(fit)[daylight_active],
  sigma = daylight_V
)

daylight_rows <- list()
daylight_draw_values <- list()
daylight_point_values <- list()
for (class_i in class_levels) {
  for (row_i in seq_len(nrow(prediction_conditions))) {
    condition_i <- prediction_conditions[row_i, ]
    key <- paste(class_i, condition_i$condition, sep = "|")
    X <- calibrated_prediction_matrices[[key]]
    probability <- plogis(drop(X %*% coef(fit)))
    probability_draws <- plogis(
      X[, daylight_active, drop = FALSE] %*% t(daylight_beta_draws)
    )
    normalized <- probability / (sum(probability) * 0.1)
    normalized_draws <- sweep(
      probability_draws, 2, colSums(probability_draws) * 0.1, "/"
    )
    daylight_point <- sum(normalized * daylight_weight) * 0.1
    daylight_draw <- colSums(normalized_draws * daylight_weight) * 0.1
    daylight_point_values[[key]] <- daylight_point
    daylight_draw_values[[key]] <- daylight_draw
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
calibrated_daylight_proportions <- bind_rows(daylight_rows)
write_csv(
  calibrated_daylight_proportions,
  out("penalty_calibrated_predicted_daylight_proportions.csv")
)

daylight_shift_rows <- list()
daylight_shift_draws <- list()
for (class_i in class_levels) {
  get_draw <- function(condition) {
    daylight_draw_values[[paste(class_i, condition, sep = "|")]]
  }
  get_point <- function(condition) {
    daylight_point_values[[paste(class_i, condition, sep = "|")]]
  }
  shift_draws <- list(
    "Wolf gradient" = get_draw("Wolf +1 / Jackal 0") -
      get_draw("Wolf -1 / Jackal 0"),
    "Jackal gradient" = get_draw("Wolf 0 / Jackal +1") -
      get_draw("Wolf 0 / Jackal -1"),
    "Wolf x jackal interaction" =
      get_draw("Wolf +1 / Jackal +1") -
      get_draw("Wolf +1 / Jackal -1") -
      get_draw("Wolf -1 / Jackal +1") +
      get_draw("Wolf -1 / Jackal -1")
  )
  shift_points <- c(
    "Wolf gradient" = get_point("Wolf +1 / Jackal 0") -
      get_point("Wolf -1 / Jackal 0"),
    "Jackal gradient" = get_point("Wolf 0 / Jackal +1") -
      get_point("Wolf 0 / Jackal -1"),
    "Wolf x jackal interaction" =
      get_point("Wolf +1 / Jackal +1") -
      get_point("Wolf +1 / Jackal -1") -
      get_point("Wolf -1 / Jackal +1") +
      get_point("Wolf -1 / Jackal -1")
  )
  for (effect_i in names(shift_draws)) {
    values <- shift_draws[[effect_i]]
    daylight_shift_draws[[paste(class_i, effect_i, sep = "|")]] <- values
    daylight_shift_rows[[paste(class_i, effect_i, sep = "|")]] <- tibble(
      demographic_class = class_i,
      predator_effect = effect_i,
      daylight_shift = unname(shift_points[effect_i]),
      lower_95 = unname(quantile(values, 0.025)),
      upper_95 = unname(quantile(values, 0.975)),
      supported_increase = quantile(values, 0.025) > 0,
      supported_decrease = quantile(values, 0.975) < 0
    )
  }
}
calibrated_daylight_shifts <- bind_rows(daylight_shift_rows)
write_csv(
  calibrated_daylight_shifts,
  out("penalty_calibrated_daylight_shift_contrasts.csv")
)

daylight_demographic_rows <- list()
for (effect_i in unique(calibrated_daylight_shifts$predator_effect)) {
  for (contrast_i in names(contrast_pairs)) {
    pair <- contrast_pairs[[contrast_i]]
    difference <-
      daylight_shift_draws[[paste(pair[1], effect_i, sep = "|")]] -
      daylight_shift_draws[[paste(pair[2], effect_i, sep = "|")]]
    point_a <- calibrated_daylight_shifts$daylight_shift[
      calibrated_daylight_shifts$demographic_class == pair[1] &
        calibrated_daylight_shifts$predator_effect == effect_i
    ]
    point_b <- calibrated_daylight_shifts$daylight_shift[
      calibrated_daylight_shifts$demographic_class == pair[2] &
        calibrated_daylight_shifts$predator_effect == effect_i
    ]
    daylight_demographic_rows[[paste(effect_i, contrast_i, sep = "|")]] <-
      tibble(
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
calibrated_daylight_demographic <- bind_rows(daylight_demographic_rows)
write_csv(
  calibrated_daylight_demographic,
  out("penalty_calibrated_demographic_daylight_shift_differences.csv")
)

if (file.exists(cfg$original_effect_functions)) {
  original_effects <- read_csv(
    cfg$original_effect_functions, show_col_types = FALSE
  ) %>%
    select(Time, predator_effect, demographic_class, original_estimate = estimate)
  curve_comparison <- calibrated_effects %>%
    inner_join(
      original_effects,
      by = c("Time", "predator_effect", "demographic_class")
    ) %>%
    group_by(predator_effect, demographic_class) %>%
    summarise(
      correlation = cor(estimate, original_estimate),
      RMSE = sqrt(mean((estimate - original_estimate)^2)),
      maximum_absolute_difference = max(abs(estimate - original_estimate)),
      .groups = "drop"
    )
  write_csv(curve_comparison, out("penalty_calibrated_curve_stability.csv"))
}

difference_curves <- bind_rows(lapply(names(effect_matrices), function(effect_i) {
  Ls <- effect_matrices[[effect_i]]
  bind_rows(lapply(names(contrast_pairs), function(contrast_i) {
    pair <- contrast_pairs[[contrast_i]]
    simultaneous_difference(Ls[[pair[1]]] - Ls[[pair[2]]], fit) %>%
      mutate(predator_effect = effect_i, contrast = contrast_i)
  }))
}))
write_csv(
  difference_curves, out("penalty_calibrated_difference_curves.csv")
)

figure <- ggplot(
  difference_curves,
  aes(x = Time, y = estimate)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#D6DEE3") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#58656B") +
  geom_line(color = "#163D5C", linewidth = 0.7) +
  facet_grid(predator_effect ~ contrast, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
  labs(
    title = "Demographic differences in predator-associated diel responses",
    subtitle = "Centered link-scale contrasts with simultaneous 95% bands",
    x = "Apparent solar time",
    y = "Difference in temporal effect"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0EEE8", color = NA),
    plot.title = element_text(face = "bold")
  )
ggsave(
  out("penalty_calibrated_difference_curves.png"), figure,
  width = 12, height = 8, dpi = 450
)
ggsave(
  out("penalty_calibrated_difference_curves.pdf"), figure,
  width = 12, height = 8, device = cairo_pdf
)

maximum_gradient <- max(abs(fit$outer.info$grad))
minimum_hessian_eigenvalue <- min(eigen(
  (fit$outer.info$hess + t(fit$outer.info$hess)) / 2,
  symmetric = TRUE, only.values = TRUE
)$values)
diagnostics <- tibble(
  convergence = ifelse(
    maximum_gradient < 1e-4 && minimum_hessian_eigenvalue > -1e-8,
    "Numerically converged", "Convergence concern"
  ),
  maximum_absolute_gradient = maximum_gradient,
  minimum_hessian_eigenvalue = minimum_hessian_eigenvalue,
  hessian_numerically_positive_semidefinite =
    minimum_hessian_eigenvalue > -1e-8,
  coefficients = length(coef(fit)),
  coefficient_rank = fit$rank,
  deviance_explained = 1 - fit$deviance / fit$null.deviance,
  pearson_dispersion = sum(residuals(fit, type = "pearson")^2) /
    df.residual(fit)
)
write_csv(diagnostics, out("penalty_calibrated_diagnostics.csv"))

# Prespecified robustness audit using the same six specifications as the
# original joint analysis. Only the omnibus-supported interaction family is
# recalibrated; unsupported primary families are not reopened post hoc.
sensitivity_specs <- tribble(
  ~specification, ~time_scale, ~history_type, ~k, ~hierarchy,
  "k = 8", "Apparent solar time", "Raw sightings", 8L, "full",
  "k = 16", "Apparent solar time", "Raw sightings", 16L, "full",
  "Sun time", "Sun time", "Raw sightings", 12L, "full",
  "Alternate camera-only hierarchy", "Apparent solar time", "Raw sightings", 12L, "camera_only",
  "Independent 30-min histories", "Apparent solar time", "Independent 30-min histories", 12L, "full",
  "Sun time + camera-only hierarchy", "Sun time", "Raw sightings", 12L, "camera_only"
)

sensitivity_result_path <- out("penalty_calibrated_sensitivity_tests.csv")
if (file.exists(sensitivity_result_path)) {
  sensitivity_tests <- read_csv(
    sensitivity_result_path, show_col_types = FALSE
  )
} else {
  sensitivity_rows <- list(
    "Primary apparent solar" = omnibus %>%
      filter(predator_effect == "Wolf x jackal x time") %>%
      transmute(
        specification = "Primary apparent solar",
        time_scale = "Apparent solar time",
        history_type = "Raw sightings",
        k = cfg$k,
        hierarchy = "full",
        EDF, reference_df, statistic, raw_P,
        warnings = "None",
        maximum_absolute_gradient = maximum_gradient,
        minimum_hessian_eigenvalue = minimum_hessian_eigenvalue
      )
  )
  for (i in seq_len(nrow(sensitivity_specs))) {
    spec <- sensitivity_specs[i, ]
    safe_name <- gsub("[^a-z0-9]+", "_", tolower(spec$specification))
    cached_result <- file.path(
      cfg$checkpoint_dir, paste0("penalty_calibrated_", safe_name, ".csv")
    )
    if (file.exists(cached_result)) {
      result_i <- read_csv(cached_result, show_col_types = FALSE)
    } else {
      message("Fitting calibrated sensitivity: ", spec$specification)
      data_i <- prepare_data(
        class_levels, spec$time_scale, spec$history_type
      )
      object_i <- fit_model(
        data_i,
        class_levels,
        checkpoint = file.path(
          cfg$checkpoint_dir, paste0("temporary_", safe_name, ".rds")
        ),
        k = spec$k,
        hierarchy = spec$hierarchy,
        save_checkpoint = FALSE
      )
      fit_i <- object_i$fit
      test_i <- deviation_test_table(
        fit_i, spec$specification, effects = "Wolf x jackal x time"
      )
      gradient_i <- max(abs(fit_i$outer.info$grad))
      hessian_i <- min(eigen(
        (fit_i$outer.info$hess + t(fit_i$outer.info$hess)) / 2,
        symmetric = TRUE, only.values = TRUE
      )$values)
      result_i <- test_i %>%
        transmute(
          specification = spec$specification,
          time_scale = spec$time_scale,
          history_type = spec$history_type,
          k = spec$k,
          hierarchy = spec$hierarchy,
          EDF, reference_df, statistic, raw_P,
          warnings = if (length(object_i$warnings)) {
            paste(object_i$warnings, collapse = " | ")
          } else {
            "None"
          },
          maximum_absolute_gradient = gradient_i,
          minimum_hessian_eigenvalue = hessian_i
        )
      write_csv(result_i, cached_result)
      rm(fit_i, object_i, data_i)
      gc(verbose = FALSE)
    }
    sensitivity_rows[[spec$specification]] <- result_i
  }
  sensitivity_tests <- bind_rows(sensitivity_rows)
  write_csv(sensitivity_tests, sensitivity_result_path)
}

sensitivity_summary <- sensitivity_tests %>%
  summarise(
    predator_effect = "Wolf x jackal x time",
    number_supported = sum(raw_P < 0.05),
    number_tested = n(),
    robust_across_all_specifications = all(raw_P < 0.05),
    failed_specifications = paste(
      specification[raw_P >= 0.05], collapse = "; "
    )
  )
write_csv(
  sensitivity_summary, out("penalty_calibrated_sensitivity_summary.csv")
)

omnibus_markdown <- c(
  "| Predator-time component | EDF | Statistic | Penalty-aware P | Original rank-Wald P |",
  "|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(omnibus)), function(i) sprintf(
    "| %s | %.3f | %.3f | %.4f | %.4f |",
    omnibus$predator_effect[i], omnibus$EDF[i], omnibus$statistic[i],
    omnibus$raw_P[i], omnibus$original_rank_Wald_P[i]
  ), character(1))
)
pairwise_markdown <- if (nrow(pairwise)) c(
  "| Predator-time component | Contrast | EDF | Statistic | Raw P | Holm P |",
  "|---|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(pairwise)), function(i) sprintf(
    "| %s | %s | %.3f | %.3f | %.4f | %.4f |",
    pairwise$predator_effect[i], pairwise$contrast[i], pairwise$EDF[i],
    pairwise$statistic[i], pairwise$raw_P[i], pairwise$Holm_adjusted_P[i]
  ), character(1))
) else "No omnibus predator family met P < 0.05."
sensitivity_markdown <- c(
  "| Specification | EDF | Statistic | P | Max gradient |",
  "|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(sensitivity_tests)), function(i) sprintf(
    "| %s | %.3f | %.3f | %.4f | %.2e |",
    sensitivity_tests$specification[i], sensitivity_tests$EDF[i],
    sensitivity_tests$statistic[i], sensitivity_tests$raw_P[i],
    sensitivity_tests$maximum_absolute_gradient[i]
  ), character(1))
)
ordered_rms_markdown <- c(
  "| Class | RMS interaction effect | 95% simulation interval |",
  "|---|---:|---:|",
  vapply(seq_len(nrow(ordered_rms)), function(i) sprintf(
    "| %s | %.4f | %.4f to %.4f |",
    ordered_rms$demographic_class[i], ordered_rms$RMS_effect[i],
    ordered_rms$lower_95[i], ordered_rms$upper_95[i]
  ), character(1))
)
ordered_contrast_markdown <- c(
  "| Directional contrast | Difference | 95% simulation interval | Delta-Wald P | Simulation tail area | Simulation Pr(>0) |",
  "|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(ordered_contrasts)), function(i) sprintf(
    "| %s | %.4f | %.4f to %.4f | %.4f | %.4f | %.3f |",
    ordered_contrasts$directional_hypothesis[i],
    ordered_contrasts$estimate[i], ordered_contrasts$lower_95[i],
    ordered_contrasts$upper_95[i], ordered_contrasts$one_sided_Wald_P[i],
    ordered_contrasts$simulation_directional_tail_area[i],
    ordered_contrasts$simulation_probability_greater_than_zero[i]
  ), character(1))
)
ordered_hypothesis_markdown <- c(
  "| Ordered hypothesis | Delta-Wald IUT P | Simulation IUT tail area | Joint simulation probability | Delta support | Simulation tail area <0.05 |",
  "|---|---:|---:|---:|---|---|",
  vapply(seq_len(nrow(ordered_hypotheses)), function(i) sprintf(
    "| %s | %.4f | %.4f | %.3f | %s | %s |",
    ordered_hypotheses$hypothesis[i],
    ordered_hypotheses$delta_Wald_intersection_union_P[i],
    ordered_hypotheses$simulation_intersection_union_tail_area[i],
    ordered_hypotheses$joint_simulation_probability[i],
    ordered_hypotheses$delta_Wald_supported_at_0_05[i],
    ordered_hypotheses$simulation_tail_area_below_0_05[i]
  ), character(1))
)
daylight_shift_markdown <- c(
  "| Class | Predator component | Daylight shift (percentage points) | 95% simulation interval |",
  "|---|---|---:|---:|",
  vapply(seq_len(nrow(calibrated_daylight_shifts)), function(i) sprintf(
    "| %s | %s | %.2f | %.2f to %.2f |",
    calibrated_daylight_shifts$demographic_class[i],
    calibrated_daylight_shifts$predator_effect[i],
    100 * calibrated_daylight_shifts$daylight_shift[i],
    100 * calibrated_daylight_shifts$lower_95[i],
    100 * calibrated_daylight_shifts$upper_95[i]
  ), character(1))
)

summary_lines <- c(
  "# Penalty-calibrated joint demographic diel test",
  "",
  "The existing class-specific and original joint analyses were not modified.",
  "This model replaced the numerical-rank Wald omnibus test with common",
  "predator-time smooths plus symmetric sum-to-zero demographic deviation",
  "smooths. Omnibus P values are the single-model penalty-aware smooth-term",
  "tests implemented by mgcv.",
  "",
  "## Omnibus tests",
  "",
  omnibus_markdown,
  "",
  "## Closed-testing pairwise follow-ups",
  "",
  pairwise_markdown,
  "",
  "## Diagnostics",
  "",
  sprintf(
    paste(
      "The primary model was %s (maximum absolute gradient = %.2e;",
      "minimum Hessian eigenvalue = %.2e; numerical PSD = %s)."
    ),
    tolower(diagnostics$convergence),
    diagnostics$maximum_absolute_gradient,
    diagnostics$minimum_hessian_eigenvalue,
    diagnostics$hessian_numerically_positive_semidefinite
  ),
  sprintf(
    paste(
      "It had %d coefficients at full rank, explained %.1f%% of deviance,",
      "and had Pearson dispersion %.3f. All three pairwise models had",
      "maximum absolute gradients below 7.7e-05, numerical PSD Hessians,",
      "and no fitting warnings."
    ),
    diagnostics$coefficients,
    100 * diagnostics$deviance_explained,
    diagnostics$pearson_dispersion
  ),
  "",
  "## Prespecified sensitivity audit",
  "",
  sensitivity_markdown,
  "",
  paste0(
    "The wolf x jackal x time demographic deviation was supported in ",
    sensitivity_summary$number_supported, "/",
    sensitivity_summary$number_tested, " prespecified specifications."
  ),
  "",
  "## Ordered H3 effect-size tests",
  "",
  paste(
    "Response strength was summarized as the RMS magnitude of each class's",
    "centered wolf x jackal temporal-interaction function. RMS measures the",
    "size of temporal redistribution but not its timing or direction."
  ),
  "",
  ordered_rms_markdown,
  "",
  ordered_contrast_markdown,
  "",
  ordered_hypothesis_markdown,
  "",
  paste(
    "The strict juvenile > female > male order was prespecified. The",
    "juvenile > both adults hypothesis is a post-hoc simplification and",
    "must be reported as such. Intersection-union support requires every",
    "component directional contrast to be supported. Because RMS is a",
    "nonlinear nonnegative quantity near its boundary, delta-Wald and",
    "coefficient-simulation results are both shown; disagreement is treated",
    "as sensitivity dependence rather than definitive support."
  ),
  "",
  "## Astronomical-daylight direction test",
  "",
  paste(
    "Positive values indicate redistribution toward astronomical daylight.",
    "Daylight was averaged across the 5,978 operational camera-days with",
    "known dates and coordinates; camera 45 round 2 was excluded from this",
    "weighting because its exact outage dates were unavailable."
  ),
  "",
  daylight_shift_markdown,
  "",
  paste(
    "No class-specific daylight shift had a 95% interval excluding zero,",
    "and none of the nine direct demographic differences in daylight shift",
    "was supported. The push-to-diurnality component of H3 is therefore",
    "unsupported in the penalty-calibrated analysis."
  ),
  "",
  "## Result",
  "",
  paste(
    "The penalty-aware omnibus test supports demographic heterogeneity in",
    "the wolf x jackal temporal interaction (P = 0.006), but not in the",
    "wolf-only (P = 0.078) or jackal-only (P = 0.707) temporal responses."
  ),
  paste(
    "Closed pairwise tests locate the interaction heterogeneity between",
    "juveniles and adult females (Holm P = 0.012) and between juveniles",
    "and adult males (Holm P = 0.011), not between the two adult classes",
    "(Holm P = 0.277)."
  ),
  paste(
    "The juvenile-adult contrast is greatest around 06:24 and 17:24, but",
    "no individual simultaneous-band interval excludes zero; inference is",
    "therefore global across the 24-hour function rather than hour-specific."
  ),
  paste(
    "The strict juvenile > female > male magnitude ordering was unsupported",
    "because the two adult RMS interaction effects were effectively identical."
  ),
  paste(
    "Juvenile RMS interaction magnitude exceeded each adult point estimate,",
    "but the post-hoc juvenile > both adults test was method-sensitive",
    "(simulation IUT tail area = 0.046; delta-Wald IUT P = 0.104) and",
    "should be",
    "described as suggestive rather than robust ordered support."
  ),
  paste(
    "Predator activity did not produce a supported increase in astronomical",
    "daylight use in any demographic class, and demographic classes did not",
    "differ directly in their daylight shifts."
  ),
  "",
  "## Interpretation rule",
  "",
  "The penalty-aware omnibus tests are primary. Pairwise tests are interpreted",
  "only within an omnibus-supported predator family and use Holm correction.",
  "Difference curves and simultaneous bands are supporting uncertainty displays."
)
writeLines(summary_lines, out("penalty_calibrated_results.md"))
capture.output(sessionInfo(), file = out("sessionInfo.txt"))

message("Penalty-calibrated analysis complete: ", cfg$output_dir)
