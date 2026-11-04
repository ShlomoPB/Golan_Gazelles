#!/usr/bin/env Rscript

# Consolidate existing tidy outputs. No model is fitted by this script.

options(stringsAsFactors = FALSE)

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

sp_all <- read.csv(
  "publication_wolf_jackal_interaction_nmix_coefficients.csv",
  check.names = FALSE
)

sp_terms <- c(
  "Null.Wolf.Abundance.z",
  "null.jackal.abundance.z",
  "Null.Wolf.Abundance.z:null.jackal.abundance.z"
)
sp_term_labels <- c(
  "Null.Wolf.Abundance.z" = "Wolf index (1 SD)",
  "null.jackal.abundance.z" = "Jackal index (1 SD)",
  "Null.Wolf.Abundance.z:null.jackal.abundance.z" = "Wolf x jackal interaction"
)

primary_sp <- sp_all[
  sp_all$framework == "spAbundance" & sp_all$term %in% sp_terms,
  c("species", "term", "estimate", "std_error", "lwr", "upr", "rhat", "ess")
]
primary_sp$predictor <- unname(sp_term_labels[primary_sp$term])
primary_sp$summary_statistic <- "Posterior mean; 95% equal-tail credible interval"
primary_sp$credible_interval_excludes_zero <-
  primary_sp$lwr > 0 | primary_sp$upr < 0
primary_sp <- primary_sp[c(
  "species", "predictor", "term", "estimate", "std_error", "lwr", "upr",
  "rhat", "ess", "credible_interval_excludes_zero", "summary_statistic"
)]
write_csv(primary_sp, "publication_primary_spabundance_coefficients.csv")

unmarked <- read.csv(
  "publication_unmarked_wolf_jackal_sensitivity_ecological_coefficients.csv",
  check.names = FALSE
)
unmarked$predictor <- c(
  wolf_z = "Wolf index (1 SD)",
  jackal_z = "Jackal index (1 SD)",
  `wolf_z:jackal_z` = "Wolf x jackal interaction"
)[unmarked$term]
unmarked$direction <- ifelse(unmarked$estimate > 0, "Positive", "Negative")
unmarked <- unmarked[c(
  "species", "predictor", "term", "estimate", "std_error", "z_value",
  "p_value", "conf_low", "conf_high", "direction"
)]
write_csv(unmarked, "publication_unmarked_sensitivity_coefficients.csv")

diel <- read.csv(
  "publication_predator_context_diel_pairwise_overlap_median_results.csv",
  check.names = FALSE
)
diel_out <- diel[c(
  "class", "context_a", "context_b", "n_a", "n_b", "cameras_with_events_a",
  "cameras_with_events_b", "estimator", "overlap", "conf_low", "conf_high",
  "camera_bootstrap_successes", "formal_test_eligible", "compareCkern_p_raw",
  "compareCkern_p", "compareCkern_p_fdr_within_class", "significance"
)]
names(diel_out)[names(diel_out) == "context_a"] <- "context_1"
names(diel_out)[names(diel_out) == "context_b"] <- "context_2"
names(diel_out)[names(diel_out) == "n_a"] <- "n_context_1"
names(diel_out)[names(diel_out) == "n_b"] <- "n_context_2"
write_csv(diel_out, "publication_authoritative_diel_context_comparisons.csv")

analysis_status <- data.frame(
  analysis = c(
    "Separate-class spAbundance NMix",
    "Expanded human/landscape single-class NMix",
    "Matched unmarked pcount",
    "Focused two-factor lfMsNMix residual correlations",
    "Broad shared-environment lfMsNMix",
    "Observed camera-based fawn:female ratio",
    "Season-adjusted additive beta-binomial composition model",
    "Beta-binomial wolf x jackal and season interaction models",
    "Overall gazelle class diel activity",
    "Gazelle diel activity across median wolf-jackal contexts",
    "Early/Late gazelle diel models",
    "Gazelle-predator temporal overlap",
    "Domestic-dog recurrent-event PAMM",
    "spOccupancy class models",
    "Unified independent-event spatial rate models",
    "Season-specific spAbundance wolf-jackal NMix",
    "Categorical predator-context spatial count models",
    "Wolf:jackal ratio exploration",
    "N-mixture-derived fawn:female ratio surfaces",
    "lfMsAbund event-rate latent-factor models",
    "Older broad msNMix/lfMsNMix models"
  ),
  scientific_purpose = c(
    "Primary class-specific spatial predator responses",
    "Assess additional human and landscape adjustment",
    "Framework sensitivity for primary spatial slopes",
    "Residual spatial correlation among classes",
    "Sensitivity of residual correlations to broader covariates",
    "Descriptive demographic-composition index",
    "Inference on fawn representation relative to adult females",
    "Sensitivity to predator and survey-period interactions",
    "Describe differences among gazelle demographic classes",
    "Test within-class timing differences among predator contexts",
    "Assess survey-period changes in gazelle timing",
    "Test direct temporal overlap with predators",
    "Directional timing around dog detections",
    "Earlier occupancy/site-use exploration",
    "Independent-event camera-use-rate sensitivity analysis",
    "Earlier seasonal predator-response exploration",
    "Earlier categorical spatial exploration",
    "Explore predator-balance indices",
    "Earlier model-derived ratio exploration",
    "Alternative event-rate joint model",
    "Earlier broad multispecies model development"
  ),
  status = c(
    "PRIMARY", "REJECTED/DIAGNOSTIC", "SENSITIVITY", "SECONDARY",
    "SUPERSEDED/SENSITIVITY", "DESCRIPTIVE", "PRIMARY", "SENSITIVITY",
    "PRIMARY DESCRIPTIVE", "PRIMARY", "EXCLUDED", "LEGACY/EXCLUDED",
    "EXPLORATORY/EXCLUDED", "EXPLORATORY/SUPERSEDED",
    "SENSITIVITY/EXPLORATORY", "EXPLORATORY/SUPERSEDED",
    "EXPLORATORY/SUPERSEDED", "EXPLORATORY/EXCLUDED",
    "EXPLORATORY/REJECTED", "EXPLORATORY/EXCLUDED",
    "SUPERSEDED/REJECTED"
  ),
  main_text = c(
    "Yes", "No", "Brief robustness statement", "Yes, concise", "No", "Yes, concise",
    "Yes", "No", "Yes", "Yes", "No", "No", "No", "No", "No", "No",
    "No", "No", "No", "No", "No"
  ),
  supporting_information = c(
    "Yes", "Yes, failure diagnostics", "Yes", "Yes", "Optional sensitivity",
    "Yes", "Yes", "Yes", "Yes", "Yes", "No", "No", "No", "No",
    "Only if explicitly retained", "No", "No", "No", "No", "No", "No"
  ),
  authoritative_script = c(
    "analyze_wolf_jackal_interaction_nmix.R",
    "fit_single_species_spnmix.R; analyze_multiclass_lfmsnmix_class_threats.R",
    "fit_unmarked_wolf_jackal_sensitivity.R; analyze_wolf_jackal_interaction_nmix.R",
    "fit_focused_multiclass_nmix.R; evaluate_focused_multiclass_nmix.R",
    "audit_multispecies_shared_env_lfmsnmix.R",
    "analyze_fawn_female_ratio.R (descriptive outputs require reconciliation)",
    "analyze_fawn_female_beta_binomial_wolf_jackal.R; summarize_fawn_female_beta_binomial_coefficients.R",
    "summarize_fawn_female_beta_binomial_interaction_sensitivity.R; test_fawn_female_beta_binomial_fixed_season.R",
    "make_overall_gazelle_diel_activity.R",
    "make_predator_context_class_diel_2x2.R; analyze_predator_context_diel_pairwise_overlap.R",
    "run_unified_spatiotemporal_analysis.R",
    "analyze_gazelle_predator_overlap_context_contrasts.R",
    "analyze_gazelle_pressure_recurrent_events.R",
    "fit_fawn_spoccupancy_model.R and related outputs",
    "run_unified_spatiotemporal_analysis.R",
    "analyze_wolf_jackal_interaction_spabundance_by_season.R",
    "analyze_predator_context_spatial_models.R",
    "analyze_wolf_jackal_ratio_exploratory.R",
    "summarize_lfmsnmix_wolf_jackal_ratio.R",
    "fit_focused_multiclass_event_rate.R",
    "run_msNMix_model.R; model_comparison.R and related scripts"
  ),
  saved_model_object = c(
    "publication_wolf_jackal_interaction_nmix_spabundance_fits.rds",
    "Multiple class-specific *_fit.rds files",
    "publication_unmarked_wolf_jackal_sensitivity_fits.rds",
    "gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds",
    "gazelle_lfMsNMix_2f_NB_shared_env_1200b_25bl_5th_3ch_fit.rds",
    "None required",
    "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds",
    "Contained in primary beta-binomial fit object and sensitivity RDS files",
    "None; activity densities and overlap are regenerated from event data",
    "None; activity densities/results regenerated from event data",
    "publication_unified_spatiotemporal_model_fits.rds",
    "Existing legacy output objects/tables; excluded",
    "publication_gazelle_dog_recurrent_events_30min_model_fits.rds",
    "Multiple *_spPGOcc_*_fit.rds files",
    "publication_unified_spatiotemporal_model_fits.rds",
    "publication_wolf_jackal_interaction_spabundance_by_season_fits.rds",
    "None or glmmTMB objects not retained consistently",
    "None required",
    "publication_fawn_female_posterior_ratio_wolf_jackal_posterior_draws.rds",
    "Smoke-test objects only",
    "Multiple gazelle_msNMix* and gazelle_lfMsNMix* objects"
  ),
  authoritative_table = c(
    "publication_primary_spabundance_coefficients.csv; publication_spabundance_wolf_jackal_relative_conditional_effect_ratios.csv",
    "publication_expanded_nmix_rejection_diagnostics.csv",
    "publication_unmarked_sensitivity_coefficients.csv; publication_unmarked_wolf_jackal_sensitivity_diagnostics.csv",
    "publication_focused_lfmsnmix_class_correlations.csv; publication_focused_lfmsnmix_diagnostics.csv; publication_focused_lfmsnmix_ppc.csv",
    "publication_multispecies_shared_env_audit_correlations.csv",
    "publication_descriptive_fawn_female_ratio.csv",
    "publication_fawn_female_beta_binomial_coefficient_details.csv",
    "publication_fawn_female_beta_binomial_interaction_sensitivity_model_comparison.csv; publication_fawn_female_beta_binomial_interaction_sensitivity_coefficients.csv",
    "publication_overall_gazelle_diel_activity_event_counts.csv; publication_overall_gazelle_diel_activity_class_overlap.csv",
    "publication_authoritative_diel_context_comparisons.csv",
    "EXCLUDED: no manuscript table",
    "EXCLUDED: no manuscript table",
    "EXCLUDED: no manuscript table",
    "EXCLUDED: no manuscript table",
    "publication_unified_spatiotemporal_space_coefficients.csv",
    "EXCLUDED: no manuscript table",
    "EXCLUDED: no manuscript table",
    "EXCLUDED: no manuscript table",
    "EXCLUDED: no manuscript table",
    "EXCLUDED: no manuscript table",
    "EXCLUDED: no manuscript table"
  ),
  authoritative_figure = c(
    "publication_spabundance_wolf_jackal_relative_interaction_curves.png/.pdf",
    "SI diagnostic figure not yet designated",
    "SI table only; no final figure designated",
    "Main text table/prose; no figure required",
    "No figure retained",
    "No figure required",
    "publication_fawn_female_beta_binomial_coefficients.png/.pdf",
    "SI table only; no primary figure",
    "publication_overall_gazelle_diel_activity.png/.pdf",
    "publication_predator_context_diel_by_class_median_combined.png/.pdf",
    "EXCLUDED: no manuscript figure",
    "EXCLUDED: no manuscript figure",
    "EXCLUDED: no manuscript figure",
    "EXCLUDED: no manuscript figure",
    "No figure retained unless sensitivity is selected",
    "EXCLUDED: no manuscript figure",
    "EXCLUDED: no manuscript figure",
    "EXCLUDED: no manuscript figure",
    "EXCLUDED: no manuscript figure",
    "EXCLUDED: no manuscript figure",
    "EXCLUDED: no manuscript figure"
  ),
  known_limitation = c(
    "Abundance and detection intercepts failed convergence; closure and individual-identification assumptions limit absolute-abundance interpretation",
    "All 18 audited saved variants failed through Rhat > 1.1, ESS < 100, or both",
    "Detection model differs from spAbundance and is treated only as framework sensitivity",
    "Intercepts and replicate-level fit are poor; use correlations only",
    "Superseded covariate structure",
    "Repeated unknown individuals; not recruitment or survival",
    "Sparse factor levels and generated predator covariates require explicit caveats",
    "Interaction terms unsupported",
    "Overall overlap is descriptive",
    "Median split is observational; compareCkern is event-level while intervals are camera-clustered",
    "Excluded by the final analysis decision; do not report seasonal/daylight diel models",
    "Explicitly excluded from manuscript and SI",
    "Outside current manuscript question",
    "Superseded by abundance and diel hierarchy",
    "Different estimand from latent relative abundance",
    "Season splitting weakens information and is no longer primary",
    "Categorization discards continuous information",
    "Ratio is denominator-sensitive and added little beyond continuous models",
    "Inadequate convergence for demographic-ratio inference",
    "Not selected for residual class correlations",
    "Obsolete formulas and/or convergence problems"
  ),
  stringsAsFactors = FALSE
)
write_csv(analysis_status, "publication_analysis_status_table.csv")

figures <- data.frame(
  figure = c(
    "Class-specific wolf x jackal relative-abundance curves",
    "Fawn-representation beta-binomial coefficient plot",
    "Class-specific gazelle diel activity across wolf-jackal contexts",
    "Overall gazelle class diel activity",
    "Relative-abundance 3 x 3 heatmap",
    "Within-class diel context overlap forest plot"
  ),
  placement = c("Main", "Main", "Main", "SI", "SI", "SI"),
  status = c(
    "PRESENT", "PRESENT", "PRESENT", "PRESENT", "PRESENT", "PRESENT"
  ),
  script = c(
    "make_spabundance_relative_wolf_jackal_figure.R",
    "make_fawn_female_beta_binomial_coefficient_plot.R",
    "make_predator_context_diel_by_class.R",
    "make_overall_gazelle_diel_activity.R",
    "make_spabundance_relative_wolf_jackal_figure.R",
    "analyze_predator_context_diel_pairwise_overlap.R"
  ),
  png = c(
    "publication_spabundance_wolf_jackal_relative_interaction_curves.png",
    "publication_fawn_female_beta_binomial_coefficients.png",
    "publication_predator_context_diel_by_class_median_combined.png",
    "publication_overall_gazelle_diel_activity.png",
    "publication_spabundance_wolf_jackal_relative_heatmap.png",
    "publication_predator_context_diel_pairwise_overlap_median.png"
  ),
  pdf = c(
    "publication_spabundance_wolf_jackal_relative_interaction_curves.pdf",
    "publication_fawn_female_beta_binomial_coefficients.pdf",
    "publication_predator_context_diel_by_class_median_combined.pdf",
    "publication_overall_gazelle_diel_activity.pdf",
    "publication_spabundance_wolf_jackal_relative_heatmap.pdf",
    "publication_predator_context_diel_pairwise_overlap_median.pdf"
  ),
  stringsAsFactors = FALSE
)
figures$file_check <- vapply(seq_len(nrow(figures)), function(i) {
  paths <- na.omit(c(figures$png[i], figures$pdf[i]))
  if (!length(paths)) return("Missing")
  if (all(file.exists(paths))) "Present" else "Incomplete"
}, character(1))
write_csv(figures, "publication_authoritative_figure_inventory.csv")

repro_paths <- c(
  "publication_spabundance_wolf_jackal_relative_interaction_curves.png",
  "publication_spabundance_wolf_jackal_relative_interaction_curves.pdf",
  "publication_fawn_female_beta_binomial_coefficients.png",
  "publication_fawn_female_beta_binomial_coefficients.pdf",
  "publication_predator_context_diel_by_class_median_combined.png",
  "publication_predator_context_diel_by_class_median_combined.pdf",
  "publication_overall_gazelle_diel_activity.png",
  "publication_overall_gazelle_diel_activity.pdf",
  "publication_spabundance_wolf_jackal_relative_heatmap.png",
  "publication_spabundance_wolf_jackal_relative_heatmap.pdf",
  "publication_predator_context_diel_pairwise_overlap_median.png",
  "publication_predator_context_diel_pairwise_overlap_median.pdf",
  "publication_primary_spabundance_coefficients.csv",
  "publication_spabundance_wolf_jackal_relative_conditional_effect_ratios.csv",
  "publication_unmarked_sensitivity_coefficients.csv",
  "publication_unmarked_wolf_jackal_sensitivity_fits.rds",
  "publication_unmarked_wolf_jackal_sensitivity_diagnostics.csv",
  "publication_expanded_nmix_rejection_diagnostics.csv",
  "publication_focused_lfmsnmix_class_correlations.csv",
  "publication_focused_lfmsnmix_diagnostics.csv",
  "publication_focused_lfmsnmix_ppc.csv",
  "publication_descriptive_fawn_female_ratio.csv",
  "publication_fawn_female_beta_binomial_coefficient_details.csv",
  "publication_fawn_female_beta_binomial_coefficients.csv",
  "publication_authoritative_diel_context_comparisons.csv",
  "publication_overall_gazelle_diel_activity_class_overlap.csv"
)
repro_generators <- c(
  rep("make_spabundance_relative_wolf_jackal_figure.R", 2),
  rep("make_fawn_female_beta_binomial_coefficient_plot.R", 2),
  rep("make_predator_context_diel_by_class.R", 2),
  rep("make_overall_gazelle_diel_activity.R", 2),
  rep("make_spabundance_relative_wolf_jackal_figure.R", 2),
  rep("analyze_predator_context_diel_pairwise_overlap.R", 2),
  "consolidate_publication_authoritative_tables.R",
  "make_spabundance_relative_wolf_jackal_figure.R",
  "consolidate_publication_authoritative_tables.R",
  "fit_unmarked_wolf_jackal_sensitivity.R",
  "fit_unmarked_wolf_jackal_sensitivity.R",
  "audit_expanded_nmix_models.R",
  "evaluate_focused_multiclass_nmix.R",
  "evaluate_focused_multiclass_nmix.R",
  "evaluate_focused_multiclass_nmix.R",
  "audit_publication_analysis_phase1.R",
  "summarize_fawn_female_beta_binomial_coefficients.R",
  "make_fawn_female_beta_binomial_coefficient_plot.R",
  "consolidate_publication_authoritative_tables.R",
  "make_overall_gazelle_diel_activity.R"
)
repro <- data.frame(
  artifact = repro_paths,
  generator = repro_generators,
  exists = file.exists(repro_paths),
  audited_on = as.character(Sys.Date()),
  known_issue = "",
  stringsAsFactors = FALSE
)
write_csv(repro, "publication_reproducibility_checklist.csv")

cat("Authoritative publication tables consolidated without model refitting.\n")
