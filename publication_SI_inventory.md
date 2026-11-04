# Proposed Supporting Information Inventory

**Audit date:** 2026-08-18  
**Rule:** Only retained analyses and diagnostics are included. Gazelle-predator temporal-overlap analyses are excluded from all appendices.

## Appendix S1. Gazelle Classification and Data Processing

### Methods

- Demographic definitions for fawn, adult female, and adult male
- Treatment of `Not Adult Male`, `Sub-adult`, missing, too-far, and unidentifiable records
- Thirty-minute independence rule and event-level group-count rule
- Camera ID and survey-round conventions
- Preservation of missing camera-rounds as `NA`
- Event-to-round allocation rule, including boundary dates

### Tables

| Proposed item | Current artifact | Status |
|---|---|---|
| Data-audit summary | `publication_data_audit_summary.csv` | Ready, pending camera 45 correction |
| Count and population provenance | `publication_count_provenance.csv` | Ready |
| Descriptive fawn:female ratio populations | `publication_descriptive_fawn_female_ratio.csv` | Ready |
| Raw classification audit | `publication_raw_gazelle_classification_audit.csv` | Ready; category definition needs approval |
| Subadult/not-adult provenance | `publication_subadult_male_provenance.csv` | Ready |
| Events outside active rounds | `publication_events_excluded_from_active_rounds.csv` | Ready |
| Count-history/event reconciliation | `publication_count_history_event_reconciliation.csv` | Ready; discrepancies unresolved |

### Figures

- Gazelle classification plate: **not yet designated**
- Additional diagnostic classification photographs: **not yet designated**

## Appendix S2. Covariate Derivation and Camera Sampling

### Methods

- Camera deployment and operational-day calculation
- Site-level wolf and jackal index derivation
- Human activity-index derivation
- Open-space, built-and-road, and nature-reserve covariate derivation
- Cattle-sighting detection covariate
- Standardization procedure and interpretation of 1-SD effects
- Early/Late split date and rationale
- Camera clocks as Israel local civil time and apparent-solar-time conversion

### Tables

| Proposed item | Current artifact | Status |
|---|---|---|
| Variable and terminology crosswalk | `publication_variable_crosswalk.csv` | Ready |
| Input hashes and timestamps | `publication_input_checksums.csv` | Ready |
| Package versions | `publication_package_versions.csv` | Ready; `glmmTMB`/`TMB` mismatch flagged |
| Camera effort details | `ObsCovs_050524.csv` plus data-audit summary | Camera 45 correction required |
| Predictor distributions and correlations | Existing exploratory outputs require selection/regeneration | Not designated |

### Figures

- Camera-array and sampling-strata map: `publication_camera_array_map.png` and `.pdf`; provisional until an authoritative regional boundary layer is supplied or deemed unnecessary
- Covariate distributions and predictor-correlation panel: optional, not designated

## Appendix S3. Spatial Abundance Diagnostics and Sensitivity

### Primary N-Mixture Methods

- Negative-binomial `spAbundance::NMix` specification
- Priors, initial values, MCMC settings, and seed
- Class-specific count histories
- Abundance and detection formulas
- Convergence criteria and restriction to relative effects
- Paired-posterior conditional-effect calculation

### Tables

| Proposed item | Current artifact | Status |
|---|---|---|
| Full primary focal coefficients | `publication_primary_spabundance_coefficients.csv` | Ready |
| All primary MCMC diagnostics | `publication_wolf_jackal_interaction_nmix_spabundance_diagnostics.csv` | Ready |
| Conditional wolf/jackal abundance ratios | `publication_spabundance_wolf_jackal_relative_conditional_effect_ratios.csv` | Ready |
| Expanded-model rejection diagnostics | `publication_expanded_nmix_rejection_diagnostics.csv` | Ready; priors/seeds unavailable from objects |
| `unmarked` sensitivity coefficients | `publication_unmarked_sensitivity_coefficients.csv` | Ready |
| `unmarked` fit diagnostics | `publication_unmarked_wolf_jackal_sensitivity_diagnostics.csv`; fit object `publication_unmarked_wolf_jackal_sensitivity_fits.rds` | Ready |

### Figures

- Relative-abundance 3 x 3 wolf-jackal heatmap: `publication_spabundance_wolf_jackal_relative_heatmap.png` and `.pdf`
- Expanded-model trace/diagnostic examples: optional, not designated
- `unmarked` sensitivity figure: not required if the coefficient table is sufficient

## Appendix S4. Demographic Composition and Latent Class Association

### Beta-Binomial Methods and Tables

| Proposed item | Current artifact | Status |
|---|---|---|
| Complete primary coefficients | `publication_fawn_female_beta_binomial_coefficient_details.csv` | Ready |
| Focal coefficient subset | `publication_fawn_female_beta_binomial_coefficients.csv` | Ready |
| Factor-level support | `publication_fawn_female_beta_binomial_coefficient_factor_support.csv` | Ready |
| Random-effect summary | `publication_fawn_female_beta_binomial_coefficient_random_effects.csv` | Ready |
| Fit diagnostics | `publication_fawn_female_beta_binomial_wolf_jackal_fit_diagnostics.csv` | Ready; software mismatch noted |
| Wolf x jackal model comparison | `publication_fawn_female_beta_binomial_interaction_sensitivity_model_comparison.csv` | Ready |
| Interaction coefficients and conditional slopes | `publication_fawn_female_beta_binomial_interaction_sensitivity_coefficients.csv`; simple-effect tables | Ready |
| Fixed-season sensitivity | `publication_fawn_female_beta_binomial_fixed_season_sensitivity_*.csv` | Ready; optional |

### Latent-Factor Methods and Tables

| Proposed item | Current artifact | Status |
|---|---|---|
| Residual class correlations | `publication_focused_lfmsnmix_class_correlations.csv` | Ready |
| Factor loadings and coefficients | `publication_focused_lfmsnmix_coefficients.csv` | Ready |
| Parameter diagnostics | `publication_focused_lfmsnmix_diagnostics.csv`; diagnostic summary | Ready |
| Posterior predictive checks | `publication_focused_lfmsnmix_ppc.csv` | Ready |
| One- versus two-factor comparison | `publication_focused_multiclass_nmix_comparison_*.csv` | Available; designation should be limited to technical support |

### Figures

- Full beta-binomial coefficient plot including adjustment terms: optional, not designated
- Raw descriptive fawn:female heatmaps: optional only if clearly labeled descriptive
- Latent loading/correlation figure: not necessary if a compact table is used

Do not include N-mixture-derived fawn:female ratio surfaces because their source fits are not adequate for demographic-ratio inference.

## Appendix S5. Gazelle Diel Activity

### Methods

- Active-round event filtering and event sample sizes
- Israel local-clock interpretation and apparent-solar-time conversion
- Daylight definition from solar elevation
- Class-specific activity-density estimation
- Coefficient of overlap and `Dhat1`/`Dhat4` rule
- Camera-cluster bootstrap with 500 replicates
- Median wolf-jackal context thresholds and camera counts
- Minimum five events for drawing a curve
- Minimum 20 events per group for formal tests
- `activity::compareCkern` with 999 null replicates and +1 correction
- Benjamini-Hochberg correction within demographic class
- Limitations of event-level tests versus camera-cluster uncertainty

### Tables

| Proposed item | Current artifact | Status |
|---|---|---|
| Overall class-class overlap | `publication_overall_gazelle_diel_activity_class_overlap.csv` | Ready |
| Complete context comparison table | `publication_authoritative_diel_context_comparisons.csv` | Ready |
| Context camera and event counts | Embedded in context outputs; consolidate if journal requires a separate table | Available |

### Figures

- Within-class context-overlap forest plot: `publication_predator_context_diel_pairwise_overlap_median.png` and `.pdf`
- Corrected overall gazelle-only activity curves: `publication_overall_gazelle_diel_activity.png` and `.pdf`
- Additional gazelle-only context panels: optional

## Explicit SI Exclusions

The following must not appear in the Supporting Information:

- Gazelle-wolf temporal-overlap coefficients or tests
- Gazelle-jackal temporal-overlap coefficients or tests
- Predator-overlay activity figures
- Predator-overlap forest plots
- Context contrasts in gazelle-predator overlap
- Domestic-dog PAMM analyses
- Wolf:jackal ratio exploration
- spOccupancy analyses
- Early/Late and daylight diel models
- Season-specific N-mixture models
- Smoke-test or failed model outputs except the concise expanded-N-mixture rejection diagnostics

## Suggested SI Tables and Figures

| Identifier | Proposed content |
|---|---|
| Table S1 | Data audit, class totals, event populations, camera-round effort, and missingness |
| Table S2 | Variable derivation and terminology crosswalk |
| Table S3 | Primary N-mixture coefficients and MCMC diagnostics |
| Table S4 | Conditional wolf and jackal effects at -1, 0, and +1 SD |
| Table S5 | Expanded N-mixture convergence failures |
| Table S6 | `unmarked` sensitivity coefficients |
| Table S7 | `lfMsNMix` loadings, correlations, diagnostics, and PPCs |
| Table S8 | Complete beta-binomial coefficients and factor support |
| Table S9 | Beta-binomial interaction and period sensitivities |
| Table S10 | Overall class overlap and complete gazelle context comparisons |
| Table S11 | Early/Late gazelle-only diel results, if retained |
| Figure S1 | Gazelle classification plate |
| Figure S2 | Predictor distributions/correlations, if needed |
| Figure S3 | Relative-abundance 3 x 3 interaction heatmap |
| Figure S4 | N-mixture diagnostics or representative traces |
| Figure S5 | Full beta-binomial coefficient plot, if useful |
| Figure S6 | Within-class context-overlap forest plot |
| Figure S7 | Corrected overall and/or Early/Late gazelle-only activity curves |
