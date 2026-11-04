#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

input_files <- c(
  "ObsCovs_050524.csv",
  "SiteCovs_060324.csv",
  "gazelle_juvenile_abundance_60x7.csv",
  "gazelle_female_abundance_60x7.csv",
  "gazelle_male_abundance_60x7.csv",
  "Gazella gazella_independent_events_30min.csv",
  "camera_timestamp_corrections.csv",
  "event_classification_exclusions.csv"
)

script_files <- c(
  "publication_analysis_workflow.R",
  "apply_camera_timestamp_corrections.R",
  "apply_event_classification_corrections.R",
  "prepare_single_species_spabundance_data.R",
  "analyze_wolf_jackal_interaction_nmix.R",
  "fit_unmarked_wolf_jackal_sensitivity.R",
  "prepare_multispecies_nmix_data.R",
  "fit_focused_multiclass_nmix.R",
  "evaluate_focused_multiclass_nmix.R",
  "analyze_fawn_female_beta_binomial_wolf_jackal.R",
  "make_spabundance_relative_wolf_jackal_figure.R",
  "make_fawn_female_beta_binomial_coefficient_plot.R",
  "make_overall_gazelle_diel_activity.R",
  "make_predator_context_diel_by_class.R",
  "analyze_predator_context_diel_pairwise_overlap.R",
  "consolidate_publication_authoritative_tables.R",
  "assemble_publication_manuscript.R",
  "freeze_publication_reproducibility.R",
  "make_conbio_word_doc.ps1"
)

model_files <- c(
  "publication_wolf_jackal_interaction_nmix_spabundance_fits.rds",
  "publication_unmarked_wolf_jackal_sensitivity_fits.rds",
  "gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds",
  "gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_metadata.rds",
  "publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds"
)

result_files <- c(
  "publication_primary_spabundance_coefficients.csv",
  "publication_spabundance_wolf_jackal_relative_conditional_effect_ratios.csv",
  "publication_unmarked_sensitivity_coefficients.csv",
  "publication_unmarked_wolf_jackal_sensitivity_diagnostics.csv",
  "publication_focused_lfmsnmix_class_correlations.csv",
  "publication_focused_lfmsnmix_diagnostics.csv",
  "publication_focused_lfmsnmix_ppc.csv",
  "publication_fawn_female_beta_binomial_coefficients.csv",
  "publication_fawn_female_beta_binomial_wolf_jackal_fit_diagnostics.csv",
  "publication_overall_gazelle_diel_activity_class_overlap.csv",
  "publication_authoritative_diel_context_comparisons.csv",
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
  "manuscript_spabundance_methods_results.md",
  "manuscript_fawn_female_beta_binomial_methods_results.md",
  "manuscript_diel_activity_methods_results.md",
  "publication_analysis_manifest.md",
  "publication_code_map.md",
  "publication_SI_inventory.md",
  "publication_analysis_status_table.csv",
  "publication_reproducibility_checklist.csv",
  "publication_package_versions.csv",
  "conbio_methods_results_raw_package.md",
  "Golan_Gazelles_ConBio_Methods_Results.docx"
)

groups <- list(
  input = input_files,
  script = script_files,
  model_object = model_files,
  result = result_files
)

manifest <- do.call(rbind, lapply(names(groups), function(group_i) {
  paths <- groups[[group_i]]
  info <- file.info(paths)
  data.frame(
    group = group_i,
    file = paths,
    exists = file.exists(paths),
    bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}))

write_csv(manifest, "publication_final_artifact_checksums.csv")

model_metadata <- data.frame(
  analysis = c(
    "Separate-class spAbundance NMix",
    "Matched unmarked pcount sensitivity",
    "Focused two-factor lfMsNMix correlations",
    "Fawn representation beta-binomial",
    "Overall and predator-context diel activity"
  ),
  family = c(
    "Negative binomial",
    "Negative binomial",
    "Negative binomial",
    "Beta-binomial",
    "Circular kernel density and overlap"
  ),
  abundance_or_fixed_formula = c(
    "~ Null.Wolf.Abundance.z * null.jackal.abundance.z",
    "~ wolf_z * jackal_z",
    "~ Open.Spaces.z + Null.Wolf.Abundance.z * null.jackal.abundance.z",
    paste(
      "~ wolf_z + jackal_z + season + open_z + built_z + reserve_z +",
      "human_z + path_type + camera_model + vegetation_structure +",
      "area_behind_camera + (1 | Polygon.Name / camera_id)"
    ),
    "Class-specific activity distributions; median wolf-jackal contexts"
  ),
  detection_formula = c(
    "~ log_days_operational.z + cattle_sightings.z",
    "~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera",
    "~ log_days_operational.z + cattle_sightings.z",
    NA_character_,
    NA_character_
  ),
  chains = c(3L, NA_integer_, 3L, NA_integer_, NA_integer_),
  samples_per_chain = c(50000L, NA_integer_, 50000L, NA_integer_, NA_integer_),
  burn_in = c(5000L, NA_integer_, 5000L, NA_integer_, NA_integer_),
  thinning = c(10L, NA_integer_, 10L, NA_integer_, NA_integer_),
  seed = c("20260623", "20260623", "20260819", "20260804", "20260810; 20260811"),
  inferential_role = c("Primary", "Sensitivity", "Secondary", "Primary", "Primary"),
  stringsAsFactors = FALSE
)
write_csv(model_metadata, "publication_final_model_specifications.csv")

writeLines(capture.output(sessionInfo()), "publication_final_session_info.txt")

if (any(!manifest$exists)) {
  stop(
    "Reproducibility freeze failed; missing files: ",
    paste(manifest$file[!manifest$exists], collapse = ", ")
  )
}

cat("Publication reproducibility freeze complete for", nrow(manifest), "artifacts.\n")
