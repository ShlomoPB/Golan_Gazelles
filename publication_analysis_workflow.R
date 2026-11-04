#!/usr/bin/env Rscript

# Staged, retained-analysis-only publication workflow.
#
# Usage:
#   Rscript publication_analysis_workflow.R --stage=list
#   Rscript publication_analysis_workflow.R --stage=audit
#   Rscript publication_analysis_workflow.R --stage=figures
#   Rscript publication_analysis_workflow.R --stage=unmarked
#   Rscript publication_analysis_workflow.R --stage=models --confirm-expensive=true
#   Rscript publication_analysis_workflow.R --stage=all --confirm-expensive=true

args <- commandArgs(trailingOnly = TRUE)
stage <- "audit"
confirm_expensive <- FALSE
for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  key <- if (length(parts) == 2) sub("^--", "", parts[1]) else ""
  if (key == "stage") {
    stage <- tolower(parts[2])
  } else if (key == "confirm-expensive") {
    confirm_expensive <- identical(tolower(parts[2]), "true")
  } else {
    stop(
      "Expected --stage=<list|audit|figures|unmarked|models|all> and optional ",
      "--confirm-expensive=true."
    )
  }
}

valid_stages <- c("list", "audit", "figures", "unmarked", "models", "all")
if (!stage %in% valid_stages) {
  stop("Unknown stage: ", stage, ". Valid stages: ", paste(valid_stages, collapse = ", "))
}

run_expensive <- identical(
  tolower(Sys.getenv("RUN_EXPENSIVE_MODELS", unset = "false")),
  "true"
) || confirm_expensive

required_anchor <- "publication_analysis_manifest.md"
if (!file.exists(required_anchor)) {
  stop("Run this workflow from the project root; missing ", required_anchor, ".")
}

run_step <- function(script, args = character()) {
  if (!file.exists(script)) stop("Missing workflow script: ", script)
  command_args <- c(script, args)
  message("Running: Rscript ", paste(command_args, collapse = " "))
  status <- system2("Rscript", args = command_args)
  if (!identical(status, 0L)) stop("Workflow step failed: ", script)
}

run_data_corrections <- function() {
  run_step("apply_camera_timestamp_corrections.R")
  run_step("apply_event_classification_corrections.R")
}

run_audit <- function() {
  run_data_corrections()
  run_step("audit_publication_analysis_phase1.R")
  run_step("audit_publication_model_objects.R")
  run_step("audit_expanded_nmix_models.R")
  run_step("consolidate_publication_authoritative_tables.R")
  run_step("assemble_publication_manuscript.R")
  run_step("freeze_publication_reproducibility.R")
}

run_unmarked <- function() {
  run_step("fit_unmarked_wolf_jackal_sensitivity.R")
  run_step("consolidate_publication_authoritative_tables.R")
}

run_figures <- function() {
  run_data_corrections()
  run_step("make_spabundance_relative_wolf_jackal_figure.R")
  run_step("summarize_fawn_female_beta_binomial_coefficients.R")
  run_step("summarize_fawn_female_beta_binomial_interaction_sensitivity.R")
  run_step("make_fawn_female_beta_binomial_coefficient_plot.R")
  run_step("make_overall_gazelle_diel_activity.R")
  run_step(
    "make_predator_context_class_diel_2x2.R",
    c(
      "--lower_quantile=0.5",
      "--upper_quantile=0.5",
      "--output_prefix=publication_predator_context_class_diel_2x2_median",
      "--seed=20260810"
    )
  )
  run_step("make_predator_context_diel_by_class.R")
  run_step("analyze_predator_context_diel_pairwise_overlap.R")
  run_step("consolidate_publication_authoritative_tables.R")
}

run_models <- function() {
  if (!run_expensive) {
    stop(
      "The models stage refits long Bayesian chains. Set ",
      "--confirm-expensive=true to confirm that this is intentional."
    )
  }

  run_data_corrections()

  class_data <- list(
    fawn = c("gazelle_juvenile_abundance_60x7.csv", "gazelle_fawn_spNMix_data.rds"),
    male = c("gazelle_male_abundance_60x7.csv", "gazelle_male_spNMix_data.rds"),
    female = c("gazelle_female_abundance_60x7.csv", "gazelle_female_spNMix_data.rds")
  )
  for (class_i in names(class_data)) {
    run_step(
      "prepare_single_species_spabundance_data.R",
      c(
        paste0("--species=", class_i),
        paste0("--count_file=", class_data[[class_i]][1]),
        paste0("--output_file=", class_data[[class_i]][2])
      )
    )
  }

  run_step("analyze_wolf_jackal_interaction_nmix.R")
  run_step(
    "prepare_multispecies_nmix_data.R",
    c(
      "--fawn_data=gazelle_fawn_spNMix_data.rds",
      "--male_data=gazelle_male_spNMix_data.rds",
      "--female_data=gazelle_female_spNMix_data.rds",
      "--output_file=gazelle_multispecies_NMix_data.rds"
    )
  )
  run_step(
    "fit_focused_multiclass_nmix.R",
    c(
      "--data=gazelle_multispecies_NMix_data.rds",
      "--model=lf2",
      "--class_order=fawn,male,female",
      "--n_batch=2000",
      "--batch_length=25",
      "--n_burn=5000",
      "--n_thin=10",
      "--n_chains=3",
      "--seed=20260819"
    )
  )
  run_step(
    "evaluate_focused_multiclass_nmix.R",
    c(
      "--structures=lf2",
      "--run_tag=2000b_25bl_10th_3ch",
      "--class_order=fawn,male,female",
      "--output_prefix=publication_focused_lfmsnmix",
      "--run_waic=false"
    )
  )
  run_step("analyze_fawn_female_beta_binomial_wolf_jackal.R")
  run_figures()
}

if (stage == "list") {
  cat(
    paste(
      "audit     Validate data, model objects, rejected models, and inventories.",
      "figures   Regenerate retained tables and figures from saved fits/data.",
      "unmarked  Refit and archive the matched frequentist sensitivity models.",
      "models    Refit retained Bayesian and beta-binomial models; guarded.",
      "all       Run models, figures, unmarked archive, and final audit; guarded.",
      sep = "\n"
    ),
    "\n"
  )
} else if (stage == "audit") {
  run_audit()
} else if (stage == "figures") {
  run_figures()
} else if (stage == "unmarked") {
  run_unmarked()
} else if (stage == "models") {
  run_models()
} else if (stage == "all") {
  run_models()
  run_unmarked()
  run_audit()
}

message("Completed publication workflow stage: ", stage)
