#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mgcv)
  library(readr)
})

checkpoint_dir <- "definitive_hierarchical_diel_checkpoints_apparent_primary"
selection_file <- "definitive_hierarchical_diel_selection_summary.csv"
model_summary_file <- "definitive_hierarchical_diel_model_summary.csv"
legacy_summary_file <- "definitive_hierarchical_diel_suntime_primary_model_summary.csv"

candidate_output <- "definitive_hierarchical_diel_audit_candidate_models.csv"
smooth_output <- "definitive_hierarchical_diel_audit_smooth_statistics.csv"
parametric_output <- "definitive_hierarchical_diel_audit_parametric_statistics.csv"
implementation_output <- "definitive_hierarchical_diel_audit_model_implementation.csv"
robustness_output <- "definitive_hierarchical_diel_audit_robustness.csv"

classes <- c("Juvenile", "Adult female", "Adult male", "Pooled species")
model_labels <- c(
  A = "baseline",
  B = "wolf-only",
  C = "jackal-only",
  D = "additive wolf + jackal",
  E = "wolf x jackal"
)
term_effects <- c(
  "s(Time):wolf_z" = "wolf x time",
  "s(Time):jackal_z" = "jackal x time",
  "s(Time):wolf_jackal" = "wolf x jackal x time"
)
term_effect_keys <- c(
  "s(Time):wolf_z" = "wolf",
  "s(Time):jackal_z" = "jackal",
  "s(Time):wolf_jackal" = "interaction"
)

collapse_code <- function(x) {
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}

checkpoint_path <- function(response_class, model_id, selected = FALSE) {
  class_key <- gsub(" ", "_", response_class, fixed = TRUE)
  stem <- if (selected) "final_fREML" else "ecology_fREML"
  file.path(
    checkpoint_dir,
    paste0("apparent_primary__", stem, "__", class_key, "__", model_id, ".rds")
  )
}

selection <- read_csv(selection_file, show_col_types = FALSE)
stopifnot(
  identical(selection$response_class, classes),
  all(selection$selected_ecology %in% names(model_labels))
)

# Candidate ecological models ------------------------------------------------

candidate_rows <- vector("list", length(classes) * length(model_labels))
row_index <- 0L
for (response_i in classes) {
  selected_id <- selection$selected_ecology[selection$response_class == response_i]
  hierarchy_i <- selection$selected_hierarchy[selection$response_class == response_i]
  for (model_i in names(model_labels)) {
    row_index <- row_index + 1L
    path_i <- checkpoint_path(response_i, model_i)
    if (!file.exists(path_i)) stop("Missing candidate checkpoint: ", path_i)
    fit_i <- readRDS(path_i)
    likelihood_i <- logLik(fit_i)
    candidate_rows[[row_index]] <- data.frame(
      response_class = response_i,
      candidate_model = model_i,
      candidate_structure = unname(model_labels[model_i]),
      hierarchy = hierarchy_i,
      formula = collapse_code(formula(fit_i)),
      fitting_function = "mgcv::bam",
      smoothing_method = fit_i$method,
      family = fit_i$family$family,
      link = fit_i$family$link,
      logLik = unname(likelihood_i),
      df = attr(likelihood_i, "df"),
      AIC = AIC(fit_i),
      selected = model_i == selected_id,
      stringsAsFactors = FALSE
    )
    rm(fit_i)
    gc(verbose = FALSE)
  }
}
candidate_table <- do.call(rbind, candidate_rows)
candidate_table$delta <- ave(
  candidate_table$AIC, candidate_table$response_class,
  FUN = function(x) x - min(x)
)
candidate_table <- candidate_table[, c(
  "response_class", "candidate_model", "candidate_structure", "hierarchy",
  "formula", "fitting_function", "smoothing_method", "family", "link",
  "logLik", "df", "AIC", "delta", "selected"
)]
write_csv(candidate_table, candidate_output)

# Exact selected-model summaries --------------------------------------------

smooth_rows <- list()
parametric_rows <- list()
implementation_rows <- list()

completed_summaries_exist <- all(file.exists(c(
  smooth_output, parametric_output, implementation_output
))) && setequal(
  unique(read_csv(smooth_output, show_col_types = FALSE)$response_class), classes
)

if (completed_summaries_exist) {
  message("Reusing completed selected-model summary exports.")
} else for (response_i in classes) {
  selected_id <- selection$selected_ecology[selection$response_class == response_i]
  hierarchy_i <- selection$selected_hierarchy[selection$response_class == response_i]
  path_i <- checkpoint_path(response_i, selected_id, selected = TRUE)
  if (!file.exists(path_i)) stop("Missing selected-model checkpoint: ", path_i)
  message("Summarizing selected model: ", response_i)
  fit_i <- readRDS(path_i)
  summary_i <- summary(fit_i)

  smooth_i <- as.data.frame(summary_i$s.table, check.names = FALSE)
  smooth_i$term <- rownames(smooth_i)
  rownames(smooth_i) <- NULL
  smooth_objects <- fit_i$smooth[match(smooth_i$term, vapply(
    fit_i$smooth, function(x) x$label, character(1)
  ))]
  smooth_table_i <- data.frame(
    response_class = response_i,
    selected_model = selected_id,
    selected_structure = unname(model_labels[selected_id]),
    hierarchy = hierarchy_i,
    term = smooth_i$term,
    EDF = smooth_i$edf,
    reference_df = smooth_i$Ref.df,
    statistic = smooth_i$Chi.sq,
    statistic_label = "Chi.sq",
    P = smooth_i$`p-value`,
    null_space_dimension = vapply(
      smooth_objects, function(x) x$null.space.dim, numeric(1)
    ),
    test_method = vapply(smooth_objects, function(x) {
      if (!x$fixed && x$null.space.dim == 0L && !is.null(fit_i$R)) {
        "mgcv reTest likelihood-ratio test with boundary reference distribution"
      } else {
        "mgcv approximate smooth significance test"
      }
    }, character(1)),
    stringsAsFactors = FALSE
  )
  smooth_rows[[response_i]] <- smooth_table_i

  parametric_i <- as.data.frame(summary_i$p.table, check.names = FALSE)
  parametric_i$term <- rownames(parametric_i)
  rownames(parametric_i) <- NULL
  parametric_rows[[response_i]] <- data.frame(
    response_class = response_i,
    selected_model = selected_id,
    selected_structure = unname(model_labels[selected_id]),
    hierarchy = hierarchy_i,
    term = parametric_i$term,
    estimate = parametric_i$Estimate,
    std_error = parametric_i$`Std. Error`,
    z_value = parametric_i$`z value`,
    P = parametric_i$`Pr(>|z|)`,
    test_method = "Approximate Wald z test using mgcv Bayesian covariance",
    stringsAsFactors = FALSE
  )

  implementation_rows[[response_i]] <- data.frame(
    response_class = response_i,
    selected_model = selected_id,
    selected_structure = unname(model_labels[selected_id]),
    hierarchy = hierarchy_i,
    formula = collapse_code(formula(fit_i)),
    fitting_call = collapse_code(fit_i$call),
    fitting_function = "mgcv::bam",
    smoothing_method = fit_i$method,
    family = fit_i$family$family,
    link = fit_i$family$link,
    select = TRUE,
    discrete = TRUE,
    knots = "list(Time = c(0, 24))",
    nthreads = 2L,
    stringsAsFactors = FALSE
  )

  write_csv(do.call(rbind, smooth_rows), smooth_output)
  write_csv(do.call(rbind, parametric_rows), parametric_output)
  write_csv(do.call(rbind, implementation_rows), implementation_output)
  rm(fit_i, summary_i, smooth_i, parametric_i)
  gc(verbose = FALSE)
}

# Existing seven-specification robustness rule ------------------------------

model_summary <- read_csv(model_summary_file, show_col_types = FALSE)
current_stages <- c(
  "Full E fREML", "Sensitivity k8", "Sensitivity k16",
  "Sensitivity sun_time", "Sensitivity alternate_hierarchy",
  "Sensitivity independent"
)
current_effects <- model_summary[
  model_summary$stage %in% current_stages &
    model_summary$term %in% names(term_effects),
]
# The current file also appends the legacy cross-sensitivity with the same
# stage label. Keep the first current-pipeline row for each stage and term.
current_effects <- current_effects[
  !duplicated(current_effects[, c("response_class", "stage", "term")]),
]
current_effects$specification <- ifelse(
  current_effects$stage == "Full E fREML",
  "Primary apparent-solar full model",
  sub("Sensitivity ", "", current_effects$stage, fixed = TRUE)
)

legacy_summary <- read_csv(legacy_summary_file, show_col_types = FALSE)
legacy_effects <- legacy_summary[
  legacy_summary$stage == "Sensitivity alternate_hierarchy" &
    legacy_summary$term %in% names(term_effects),
]
legacy_effects$specification <- "sun_time + alternate_hierarchy"
all_effects <- rbind(
  current_effects[, c("response_class", "term", "P", "specification")],
  legacy_effects[, c("response_class", "term", "P", "specification")]
)

stability <- read_csv(
  "definitive_hierarchical_diel_effect_stability.csv", show_col_types = FALSE
)
robustness_rows <- list()
for (response_i in classes) {
  for (term_i in names(term_effects)) {
    subset_i <- all_effects[
      all_effects$response_class == response_i & all_effects$term == term_i,
    ]
    stability_i <- stability[
      stability$response_class == response_i &
        stability$effect == unname(term_effect_keys[term_i]),
    ]
    failed_i <- subset_i$specification[is.na(subset_i$P) | subset_i$P >= 0.05]
    robustness_rows[[paste(response_i, term_i)]] <- data.frame(
      response_class = response_i,
      effect = unname(term_effects[term_i]),
      primary_P = subset_i$P[
        subset_i$specification == "Primary apparent-solar full model"
      ][1],
      primary_supported = subset_i$P[
        subset_i$specification == "Primary apparent-solar full model"
      ][1] < 0.05,
      supported_specifications = sum(subset_i$P < 0.05, na.rm = TRUE),
      tested_specifications = sum(!is.na(subset_i$P)),
      failed_specifications = if (length(failed_i)) {
        paste(failed_i, collapse = "; ")
      } else {
        "None"
      },
      term_retained = stability_i$term_retained,
      robustness_status = stability_i$robustness_status,
      stringsAsFactors = FALSE
    )
  }
}
robustness_table <- do.call(rbind, robustness_rows)
rownames(robustness_table) <- NULL
write_csv(robustness_table, robustness_output)

message("Definitive hierarchical diel model audit exports completed.")
