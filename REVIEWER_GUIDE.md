# Reviewer guide

This repository is a complete computational archive. It intentionally retains
exploratory analyses, sensitivity checks, saved model objects, diagnostics, and
derived outputs, but those files are not all part of the primary manuscript
workflow. This guide identifies the shortest path through the analysis.

## Core inputs

- `ObsCovs_050524.csv`: camera-round observation covariates and effort.
- `SiteCovs_060324.csv`: camera-site covariates.
- `gazelle_male_abundance_60x7.csv`, `gazelle_female_abundance_60x7.csv`, and
  `gazelle_juvenile_abundance_60x7.csv`: class-specific count histories.
- `Gazella gazella_independent_events_30min_analysis_ready.csv`: corrected
  class-specific independent-event data.
- `camera_timestamp_corrections.csv`: documented camera-clock corrections.

## Manuscript-critical workflows

1. `apply_camera_timestamp_corrections.R` and
   `apply_event_classification_corrections.R` document event-level corrections.
2. `prepare_spabundance_data.R` prepares the spatial count-model inputs.
3. `analyze_multiclass_lfmsnmix_class_threats.R` and the
   `publication_*spatial*` outputs contain the focal spatial demographic models
   and their diagnostics.
4. `analyze_fawn_female_beta_binomial_wolf_jackal.R` contains the published
   juvenile:adult-female ratio model.
5. `definitive_hierarchical_diel_GAM.R` contains the primary class-specific and
   pooled diel GAM workflow.
6. `calibrate_joint_demographic_diel_GAM.R` contains the confirmatory joint
   demographic GAM workflow.
7. `R/generate_supplementary_tables.R` and
   `supplementary_tables/supplementary_materials_final.qmd` assemble the
   supplementary tables from saved results.

## Authoritative review outputs

- `supplementary_tables/`: final CSV/XLSX table sources and assembly source.
- `definitive_hierarchical_diel_*`: primary diel results, diagnostics, model
  summaries, predictions, and session information.
- `joint_demographic_penalty_calibrated_outputs/`: confirmatory joint-GAM
  estimates and diagnostics.
- `publication_*`: publication-focused spatial outputs and sensitivity checks.
- `juvenile_female_seasonal_rate_analysis/` and
  `juvenile_female_seasonal_audit/`: explicitly exploratory follow-up analyses.

## Files that are not primary evidence

Files or folders labelled `smoke`, `exploratory`, `audit`, `checkpoint`, or
`in_progress` preserve development, diagnostics, or recovery state. They should
not be treated as manuscript results unless a primary script or final table
explicitly references them. Saved `.rds` objects are retained to verify fitted
models without requiring costly refits and are managed through Git LFS.

No model should be rerun merely to inspect the reported estimates: the relevant
coefficient tables, diagnostics, predictions, and session-information files are
already included alongside the scripts.
