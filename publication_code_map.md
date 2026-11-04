# Publication Code Map

This file lists the main scripts that reproduce the class-specific and joint-model results used for publication.

## Data Preparation

- [prepare_single_species_spabundance_data.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/prepare_single_species_spabundance_data.R:1)  
  Builds class-specific `spAbundance` data objects for fawns, males, and females, including site covariates, detection covariates, and projected camera coordinates.

- [prepare_multispecies_nmix_data.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/prepare_multispecies_nmix_data.R:1)  
  Combines the three class-specific data objects into a single multi-class `lfMsNMix` data object.

## Primary Class-specific Models

- [fit_single_species_spnmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_single_species_spnmix.R:1)  
  Fits class-specific `NMix` or `spNMix` models. This is the main fitting script for the recommended non-spatial negative-binomial abundance analysis.

## Sensitivity Occupancy Models

- [prepare_fawn_spoccupancy_data.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/prepare_fawn_spoccupancy_data.R:1)  
  Builds the fawn occupancy dataset.

- [fit_fawn_spoccupancy_model.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_fawn_spoccupancy_model.R:1)  
  Fits the spatial occupancy sensitivity models.

## Exploratory Joint Models

- [fit_multiclass_lfmsnmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_multiclass_lfmsnmix.R:1)  
  Fits the shared-environment latent-factor negative-binomial N-mixture model.

- [fit_multiclass_lfmsnmix_class_threats.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_multiclass_lfmsnmix_class_threats.R:1)  
  Fits the class-threat latent-factor negative-binomial N-mixture models with design-based detection and optional effort.

- [analyze_multiclass_lfmsnmix_class_threats.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_multiclass_lfmsnmix_class_threats.R:1)  
  Summarizes the class-threat joint models, including class-specific coefficients, `kappa`, manual posterior predictive checks, and cross-class latent mean-abundance correlations.

## Workflow Entry Point

- [publication_analysis_workflow.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_analysis_workflow.R:1)  
  Staged retained-analysis-only driver for audit, figure, `unmarked`, and guarded long-model refits. Excluded and superseded analyses are not invoked.

- [fit_unmarked_wolf_jackal_sensitivity.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_unmarked_wolf_jackal_sensitivity.R:1)  
  Refits and archives the matched negative-binomial `unmarked::pcount` sensitivity models, complete coefficient tables, optimizer diagnostics, and session metadata.

- [make_overall_gazelle_diel_activity.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_overall_gazelle_diel_activity.R:1)  
  Generates the authoritative overall gazelle-class activity curves and camera-clustered class-overlap intervals from local camera time converted to apparent solar time.

- [make_publication_camera_array_map.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_publication_camera_array_map.R:1)  
  Generates the provisional 60-camera, 12-stratum sampling-layout map in WGS 84 / UTM Zone 36N without inferring an unavailable boundary layer.

- [make_conbio_word_doc.ps1](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_conbio_word_doc.ps1:1)  
  Builds [Golan_Gazelles_ConBio_Methods_Results.docx](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/Golan_Gazelles_ConBio_Methods_Results.docx) from the raw methods/results package and embeds the main figures and summary tables.

## Directional Recurrent-Event Models

- [analyze_gazelle_pressure_recurrent_events.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_gazelle_pressure_recurrent_events.R:1)  
  Fits configurable directional recurrent-event PAMMs from the complete Camelot export, using camera sessions as the sampling frame. The current publication workflow selects domestic dogs and 30-min independent events, fits six dog-gazelle directions, and performs 999 camera-session circular-shift randomization tests with false-discovery-rate correction.

- [publication_gazelle_dog_recurrent_events_30min_results.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_gazelle_dog_recurrent_events_30min_results.md:1)  
  Reports the independent-event audit, dog-gazelle PAMM contrasts, randomization effects, multiplicity-adjusted inference, interpretation guardrails, and output manifest.

## Predator-context Diel Activity

- [make_predator_context_class_diel_2x2.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_predator_context_class_diel_2x2.R:1)  
  Produces four-corner wolf-jackal diel figures from 30-min independent events, confirmed local-clock timestamps converted to apparent solar time, and camera-clustered overlap intervals. The default [33rd/67th-percentile figure](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_class_diel_2x2.png) emphasizes extreme contexts; the [median-split figure](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_class_diel_2x2_median.png) assigns all 60 cameras. Sparse class-context combinations are explicitly flagged rather than assigned unsupported curves.

- [make_predator_context_diel_by_class.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_predator_context_diel_by_class.R:1)  
  Reorients the all-camera median-split densities so each demographic class is held constant and the four predator contexts are overlaid. It produces a [combined class-specific figure](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_diel_by_class_median_combined.png) and separate [fawn](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_diel_by_class_median_fawn.png), [female](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_diel_by_class_median_female.png), and [male](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_diel_by_class_median_male.png) figures.

- [make_gazelle_predator_diel_context_grid.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_gazelle_predator_diel_context_grid.R:1)  
  Exploratory only; not retained in the manuscript workflow. It overlays context-specific wolf and jackal activity on each gazelle panel.

- [make_gazelle_predator_diel_main_figure.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_gazelle_predator_diel_main_figure.R:1)  
  Exploratory only; not retained in the manuscript workflow. It produces a gazelle-predator activity overlay with one predator panel above each context.

- [analyze_gazelle_predator_overlap_context_contrasts.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_gazelle_predator_overlap_context_contrasts.R:1)  
  Exploratory only; not retained in the manuscript workflow. It tests whether each class's overlap with wolves or jackals differs among median-split contexts.

- [analyze_predator_context_diel_pairwise_overlap.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_predator_context_diel_pairwise_overlap.R:1)  
  Quantifies all six pairwise temporal-overlap contrasts among the four median-split predator contexts within each demographic class. It selects `Dhat1` when the smaller sample has fewer than 75 events and `Dhat4` otherwise, computes camera-cluster bootstrap intervals, and runs `activity::compareCkern` with 999 null replicates. The retained outputs report unadjusted P-values and include the [ranked overlap figure](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_diel_pairwise_overlap_median.png), [complete result table](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_diel_pairwise_overlap_median_results.csv), and [analysis report](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_diel_pairwise_overlap_median_results.md).

- [analyze_predator_context_nocturnality.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_predator_context_nocturnality.R:1)  
  Calculates the Procko et al. (2023) event-level nocturnality metric as absolute hours from apparent solar noon and fits class-specific Gaussian mixed models with predator context as a fixed effect and camera as a random intercept. Outputs include the [nocturnality figure](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_nocturnality_median.png), [omnibus tests](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_nocturnality_median_omnibus_tests.csv), [raw pairwise contrasts](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_nocturnality_median_pairwise_contrasts.csv), and [analysis report](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_nocturnality_median_results.md).

## Exploratory Predator Balance

- [analyze_wolf_jackal_ratio_exploratory.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_wolf_jackal_ratio_exploratory.R:1)  
  Creates [SiteCovs_060324_with_wolf_jackal_ratio.csv](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/SiteCovs_060324_with_wolf_jackal_ratio.csv) without modifying the original site file. It calculates the raw wolf:jackal ratio, log-ratio, bounded wolf share, and total predator index; tests camera-level associations with habitat, gazelle abundance, class-specific count rates, fawn composition, and daylight use; and adds partial-rank, extreme-camera, denominator-aware, and predator-balance-versus-load sensitivity checks. Main outputs are the [results report](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/exploratory_wolf_jackal_ratio_results.md), [primary figure](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/exploratory_wolf_jackal_ratio_primary_relationships.png), and [correlation summary](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/exploratory_wolf_jackal_ratio_correlation_summary.png).

## Manuscript Methods and Results

- [manuscript_spabundance_methods_results.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/manuscript_spabundance_methods_results.md:1)  
  Provides manuscript-ready Methods, Results, table and figure captions, convergence limits, and references. It leads with the finalized separate-class wolf-by-jackal N-mixture models and uses the focused two-factor `lfMsNMix` only for secondary residual class correlations.

- [focused_multiclass_spabundance_model_review.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/focused_multiclass_spabundance_model_review.md:1)  
  Documents the literature-based model review, focused `Open.Spaces + wolf * jackal` N-mixture screening, convergence and predictive-fit limitations, and the decision to retain the two-factor model only for residual class correlations.

- [fit_focused_multiclass_nmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_focused_multiclass_nmix.R:1) and [evaluate_focused_multiclass_nmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/evaluate_focused_multiclass_nmix.R:1)  
  Fit and compare no-factor, one-factor, and two-factor negative-binomial multi-class N-mixture models under the focused habitat and predator-interaction formula. Outputs include full diagnostics, multiple posterior predictive checks, coefficient summaries, latent correlations, and common latent-abundance integration ceilings.

- [audit_multispecies_shared_env_lfmsnmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/audit_multispecies_shared_env_lfmsnmix.R:1)  
  Audits the superseded broad shared-environment `lfMsNMix` fit. It remains available for provenance and sensitivity comparison but is no longer designated as the primary joint model.

- [manuscript_diel_activity_methods_results.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/manuscript_diel_activity_methods_results.md:1)  
  Provides manuscript-ready Methods, Results, figure captions, and references for overall class activity, all-camera predator contexts, within-class overlap tests, and Procko-style nocturnality models. Direct gazelle-predator temporal overlap is excluded.

- [make_spabundance_relative_wolf_jackal_figure.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_spabundance_relative_wolf_jackal_figure.R:1)  
  Produces posterior relative-abundance curves and a 3 by 3 heatmap from the converged predator slopes. The predictions exclude the nonconverged abundance intercept and are normalized to mean wolf and jackal activity.

- [analyze_wolf_jackal_interaction_spabundance_by_season.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_wolf_jackal_interaction_spabundance_by_season.R:1)  
  Fits separate Early and Late negative-binomial `spAbundance` N-mixture models for each class using a calendar-midpoint split defined from camera deployment dates. It exports complete MCMC diagnostics, relative response surfaces, conditional wolf effects at three jackal levels, and approximate posterior Late-minus-Early contrasts. Primary outputs are the [seasonal results report](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_wolf_jackal_interaction_spabundance_by_season_results.md), [relative heatmap](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_wolf_jackal_interaction_spabundance_by_season_relative_heatmap.png), and [conditional wolf-effect figure](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_wolf_jackal_interaction_spabundance_by_season_conditional_wolf_effects.png).

## Demographic Ratio Visualization

- [manuscript_fawn_female_beta_binomial_methods_results.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/manuscript_fawn_female_beta_binomial_methods_results.md:1)  
  Provides manuscript-ready, first-person Methods and Results for the primary season-adjusted additive beta-binomial composition model and the otherwise identical wolf-by-jackal interaction sensitivity analysis. It includes the main figure caption, methodological rationale, exact verified estimates, and literature references.

- [analyze_fawn_female_beta_binomial_wolf_jackal.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_fawn_female_beta_binomial_wolf_jackal.R:1)  
  Fits the primary adjusted camera-round beta-binomial model set with continuous wolf, jackal, season, landscape and camera-design controls, and nested camera random effects. The retained specification uses additive predator slopes and a common Early-Late shift; wolf-by-jackal and season-specific predator interactions were unsupported.

- [summarize_fawn_female_beta_binomial_coefficients.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/summarize_fawn_female_beta_binomial_coefficients.R:1)  
  Produces the preferred-model coefficient table, interaction-model simple slopes, factor-level support audit, random-effect summary, and [coefficient interpretation report](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_coefficient_interpretation.md).

- [summarize_fawn_female_beta_binomial_interaction_sensitivity.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/summarize_fawn_female_beta_binomial_interaction_sensitivity.R:1)  
  Compares the primary season-adjusted additive model with an otherwise identical model containing a wolf-by-jackal interaction. It exports the interaction coefficient, conditional predator slopes at `-1`, `0`, and `+1` SD of the other predator, convergence checks, and the nested model comparison in a dedicated [sensitivity report](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_interaction_sensitivity_results.md).

- [fit_fawn_female_beta_binomial_focused.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_fawn_female_beta_binomial_focused.R:1)  
  Fits a reduced-model sensitivity analysis with only standardized wolf activity, standardized jackal activity, their interaction, and standardized open-space cover as fixed effects. It preserves the original beta-binomial response and nested polygon/camera random intercepts. Outputs include the [model report](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_focused_results.md), [coefficient table](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_focused_coefficients.csv), interaction likelihood-ratio test, conditional predator slopes, diagnostics, and saved fits.

- [test_fawn_female_beta_binomial_fixed_season.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/test_fawn_female_beta_binomial_fixed_season.R:1)  
  Evaluates Early versus Late survey period as a fixed-effect sensitivity analysis; two season levels are not treated as a random-effect distribution. It compares the focused interaction model with and without season and retests the wolf-by-jackal interaction after season adjustment. Outputs include the [sensitivity report](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_fixed_season_sensitivity_results.md), coefficients, nested model tests, diagnostics, and saved fits.

- [make_fawn_female_beta_binomial_coefficient_plot.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_fawn_female_beta_binomial_coefficient_plot.R:1)  
  Verifies the retained season-adjusted additive fit against the saved coefficient output and produces the primary demographic-composition [forest plot](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_coefficients.png) and [underlying coefficient table](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_coefficients.csv). The figure displays season, open-space cover, wolf abundance, and jackal abundance as odds ratios with 95% confidence intervals.

- [summarize_pooled_fawn_female_beta_binomial.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/summarize_pooled_fawn_female_beta_binomial.R:1)  
  Summarizes the full-survey model without a season term, compares additive and wolf-by-jackal structures, tests the loss of the common seasonal shift, and produces the [pooled additive heatmap](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_pooled_heatmap.png), [interaction curves](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_pooled_interaction_curves.png), and [pooled results report](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_beta_binomial_pooled_results.md).

- [fit_lfmsnmix_wolf_jackal_ratio.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_lfmsnmix_wolf_jackal_ratio.R:1) and [summarize_lfmsnmix_wolf_jackal_ratio.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/summarize_lfmsnmix_wolf_jackal_ratio.R:1)  
  Fit and summarize paired-draw, two-factor negative-binomial `lfMsNMix` ratios for Overall, Early, and Late periods. The method correctly preserves fawn-female posterior covariance, but grid-cell log-ratio convergence was inadequate (`Rhat` up to `1.87`), so [publication_fawn_female_posterior_ratio_wolf_jackal_heatmap.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_posterior_ratio_wolf_jackal_heatmap.png) is exploratory only.

- [make_fawn_female_count_ratio_heatmaps.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_fawn_female_count_ratio_heatmaps.R:1)  
  Produces pooled-count and equal-camera/session descriptive heatmaps. Use [publication_fawn_female_count_ratio_wolf_jackal_combined.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_count_ratio_wolf_jackal_combined.png) only as a supplementary empirical-support diagnostic.

- [make_fawn_female_wolf_jackal_ratio_heatmaps.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_fawn_female_wolf_jackal_ratio_heatmaps.R:1)  
  Converts saved separate-class posterior mean prediction grids into ratios. Retained for provenance, but superseded because dividing posterior means does not propagate ratio uncertainty or preserve cross-class posterior covariance.

- [publication_fawn_female_ratio_wolf_jackal_results.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_ratio_wolf_jackal_results.md:1)  
  Reports all 27 predicted abundance ratios and the interpretation limits of ratios derived from separate posterior mean predictions.

## Preliminary Results Export Manifest

This section maps the figures and tables in `Golan Gazelles Prelim Results.docx` to the current reproducible outputs in this workspace. The six embedded PNG figures in the Word file were exact hash matches to workspace PNGs at the time of audit.

### Figures

| Word media | Workspace output | Producing script | Status |
|---|---|---|---|
| `image1.png` | [publication_class_coefficient_figure.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_class_coefficient_figure.png) | [make_publication_coefficient_figure.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_publication_coefficient_figure.R:1) | Exact match |
| `image3.png` | [publication_wolf_jackal_interaction_spabundance_by_season_interaction_curves.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_wolf_jackal_interaction_spabundance_by_season_interaction_curves.png) | Legacy version of [analyze_wolf_jackal_interaction_spabundance_by_season.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_wolf_jackal_interaction_spabundance_by_season.R:1) | Exact match to preliminary document, but superseded by the audited relative curves and heatmap because abundance intercepts did not converge |
| `image4.png` | [publication_wolf_jackal_interaction_nmix_heatmap_spabundance_only.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_wolf_jackal_interaction_nmix_heatmap_spabundance_only.png) | [analyze_wolf_jackal_interaction_nmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_wolf_jackal_interaction_nmix.R:1) | Exact match |
| `image5.png` | [publication_temporal_overall_activity.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_temporal_overall_activity.png) | [make_temporal_publication_plots.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_temporal_publication_plots.R:1) | Exact match, older local-time plot |
| `image6.png` | [gazelle_temporal_activity_solar_season_overall_class_activity_curves.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/gazelle_temporal_activity_solar_season_overall_class_activity_curves.png) | [make_temporal_solar_season_plots.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_temporal_solar_season_plots.R:1) | Exact match, older solar-season plot |
| `image7.png` | [publication_predator_context_class_overlap_q33_67.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_class_overlap_q33_67.png) | [make_predator_context_class_overlap_plots.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_predator_context_class_overlap_plots.R:1) | Exact match |
| `image2.emf` | No matching workspace `.emf` found | Unknown/Word vector object | Not currently mapped; likely nonessential unless it is a manuscript figure |

The final retained temporal figures are [publication_predator_context_diel_by_class_median_combined.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_diel_by_class_median_combined.png), [publication_overall_gazelle_diel_activity.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_overall_gazelle_diel_activity.png), [publication_predator_context_nocturnality_median.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_nocturnality_median.png), and [publication_predator_context_diel_pairwise_overlap_median.png](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_predator_context_diel_pairwise_overlap_median.png). Unified Early/Late and daylight outputs are excluded.

### Tables And Text Results

| Preliminary result | Workspace output | Producing script | Status |
|---|---|---|---|
| Wolf x jackal `unmarked` and `spAbundance` coefficient table | [publication_wolf_jackal_interaction_nmix_results.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_wolf_jackal_interaction_nmix_results.md:1) and [publication_wolf_jackal_interaction_nmix_coefficients.csv](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_wolf_jackal_interaction_nmix_coefficients.csv) | [analyze_wolf_jackal_interaction_nmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_wolf_jackal_interaction_nmix.R:1) | Reproducible |
| Seasonal spAbundance wolf x jackal coefficient table | Historical seasonal outputs | [analyze_wolf_jackal_interaction_spabundance_by_season.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_wolf_jackal_interaction_spabundance_by_season.R:1) | Excluded from manuscript and SI |
| Class-specific coefficient comparison table | [publication_class_comparison_table.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_class_comparison_table.md:1) and [publication_class_comparison_table.csv](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_class_comparison_table.csv) | [fit_single_species_spnmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_single_species_spnmix.R:1) and [make_publication_coefficient_figure.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/make_publication_coefficient_figure.R:1) | Reproducible |
| Latent class correlations | [publication_focused_lfmsnmix_class_correlations.csv](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_focused_lfmsnmix_class_correlations.csv) and [publication_multispecies_shared_env_audit_results.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_multispecies_shared_env_audit_results.md:1) | [fit_focused_multiclass_nmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/fit_focused_multiclass_nmix.R:1) and [evaluate_focused_multiclass_nmix.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/evaluate_focused_multiclass_nmix.R:1) | Reproducible secondary joint result; inference is restricted to residual correlations because intercepts and absolute abundance are weakly identified |
| Juvenile:adult-female ratio and open-space effect | [publication_fawn_female_ratio_results.md](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_ratio_results.md:1), [publication_fawn_female_ratio_by_round.csv](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_ratio_by_round.csv), and [publication_fawn_female_ratio_model_coefficients.csv](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/publication_fawn_female_ratio_model_coefficients.csv) | [analyze_fawn_female_ratio.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/analyze_fawn_female_ratio.R:1) | Reproducible |
| Unified Early/Late and daylight synthesis | Historical unified outputs | [run_unified_spatiotemporal_analysis.R](C:/Users/preissbloom/Documents/Github%20Repositories/Golan_Gazelles/Gazelle%20VSCode/run_unified_spatiotemporal_analysis.R:1) | Excluded from manuscript and SI |
