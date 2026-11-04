# Publication Analysis Cleanup Report

**Audit date:** 2026-08-18  
**Scope:** Non-destructive audit and reproducibility consolidation. No source data were edited, no files were moved or deleted, no Bayesian model was rerun, and manuscript prose was not changed. The matched `unmarked` sensitivity models were refitted and archived.

## A. Confirmed Authoritative Analyses

1. The primary spatial analysis is the full-survey, separate-class negative-binomial `spAbundance::NMix` analysis with `wolf_z * jackal_z` in abundance and operational days plus cattle sightings in detection.
2. The matched `unmarked::pcount` fits are sensitivity analyses only. Their predator-gradient conclusions agree qualitatively with the primary models.
3. The focused two-factor `lfMsNMix` fit is secondary and is retained only for residual class correlations after open-space cover and wolf-jackal context are considered.
4. The pooled fawn:adult-female ratio is a descriptive count-history index: 199 / 1,951 = 0.102, or 10.2 fawns per 100 adult females.
5. The primary demographic-composition analysis is the season-adjusted additive beta-binomial mixed model. Late survey period and open-space cover are supported; wolf and jackal main effects are not.
6. The primary temporal analysis concerns gazelles only: overall class activity and within-class comparisons across median-split wolf-jackal camera contexts.
7. The authoritative context figure parses camera timestamps as Israel local civil time before converting them to apparent solar time.

The full specifications and artifacts are in `publication_analysis_manifest.md` and `publication_analysis_status_table.csv`.

## B. Superseded Analyses

The following should no longer be cited as primary manuscript evidence:

- Broad shared-environment `lfMsNMix` models containing nature-reserve, open-space, and built-and-road covariates together.
- Season-specific `spAbundance` wolf-jackal N-mixture models.
- Categorical predator-context spatial count models.
- Older pooled or round-factor fawn:female models in `analyze_fawn_female_ratio.R`.
- Independent-event spatial rate models where they duplicate the primary latent relative-abundance question.
- spOccupancy analyses of fawns or other classes.
- Old `msNMix`, broad `lfMsNMix`, and smoke-test model objects.
- Older overall diel figures that interpreted the stored timestamps as UTC before converting to local or solar time.

These files were retained in place because moving them now would risk breaking historical scripts.

## C. Rejected Analyses

Eighteen saved nonspatial class-specific N-mixture variants containing class-specific landscape, human, boar, wolf, or jackal combinations were audited from their saved objects. Every fit failed the audit criterion through `Rhat > 1.1`, ESS < 100, or both.

| Class | Saved fits | Worst Rhat | Lowest ESS |
|---|---:|---:|---:|
| Fawn | 6 | 2.28 | 6.78 |
| Adult female | 6 | 3.09 | 4.98 |
| Adult male | 6 | 2.06 | 6.64 |

The complete table is `publication_expanded_nmix_rejection_diagnostics.csv`. These models must not support ecological inference. Priors and seeds are not embedded in their fit objects, which is an additional provenance limitation.

N-mixture-derived fawn:female ratio surfaces are also rejected because their component abundance models did not provide adequate convergence for demographic-ratio inference.

## D. Legacy and Excluded Analyses

All gazelle-predator temporal-overlap work is **LEGACY / EXCLUDED FROM MANUSCRIPT AND SUPPORTING INFORMATION**. This includes:

- `analyze_gazelle_predator_overlap_context_contrasts.R`
- `publication_gazelle_predator_overlap_context_contrasts*`
- `publication_gazelle_predator_diel_by_context_median_gazelle_predator_overlap.csv`
- Predator-overlay diel figures and density tables
- Gazelle-wolf and gazelle-jackal overlap coefficients
- Predator-overlap forest plots, context tests, and FDR results

Domestic-dog recurrent-event PAMM outputs, wolf:jackal ratio exploration, and `lfMsAbund` event-rate smoke fits are also excluded from the current manuscript.

`publication_analysis_workflow.R` is now a staged retained-analysis-only driver. Long model refits require the explicit `RUN_EXPENSIVE_MODELS=true` guard; dog PAMM, predator-overlap, ratio exploration, occupancy, and superseded seasonal spatial models are absent.

## E. Inconsistencies Found

### Counts and analysis populations

| Quantity | Current value | Competing value | Audit conclusion |
|---|---:|---:|---|
| Fawns | 199 | 195 | The three retained count-history/beta analyses use 199. No reproducible current source for 195 was found. |
| Adult females | 1,951 | 2,135 | The final count history and beta model use 1,951. No reproducible current source for 2,135 was found. |
| Adult females in all prepared event data | 1,953 | 1,951 | Two-individual difference reflects nonidentical source-event and final count-history populations. |
| Adult females in date-resolved active diel rounds | 1,931 | 1,951 | Diel filtering removes events outside mapped active rounds; event counts are not interchangeable with count-history totals. |
| Adult males in final count history | 2,355 | 2,346 active-round event count | Different round allocation/filtering populations. |
| Not-adult/subadult males | 48 events / 62 event maxima | 51 | Current raw exports reproduce 43 `Not Adult Male` events / 55 event maxima, or 48 events / 62 event maxima when `Sub-adult` life stage is added. The origin and unit of 51 remain unresolved. |

The authoritative survey and N-mixture totals should be 199 fawns, 1,951 adult females, and 2,355 adult males, subject to final count-history regeneration. The authoritative diel sample sizes should be reported as event rows: 180, 1,466, and 2,260, respectively.

### Effort and camera 45

`ObsCovs_050524.csv` contains 344 positive-effort rounds totaling 5,999 days, but only 343 have valid dates and these total 5,970 days. File line 313, camera 45 round 4, has 29 operational days and blank dates. The adult-male count history assigns two males to round 4, while the two corresponding October-November events map to round 2 under the dates currently stored on file lines 310-316.

This is not merely a reporting difference: it affects event-to-round reconciliation and the choice between 5,999 recorded camera-days and 5,970 date-resolved camera-days. It must be corrected in the source data or explicitly documented before final submission.

### Count-history versus date-mapped events

Direct aggregation of active-round independent events did not reproduce every count-history cell:

| Class | Differing class-camera-round cells | Sum of absolute cell differences |
|---|---:|---:|
| Fawn | 10 | 12 |
| Adult female | 50 | 122 |
| Adult male | 76 | 123 |

Exact rows are in `publication_count_history_event_reconciliation.csv`. Likely causes include boundary-date allocation, manual curation, and the camera 45 discrepancy, but these causes have not been proven for every row. The count histories and event file must therefore remain separate authoritative inputs for their current analyses until a single generation script reproduces them.

### Model and output provenance

- Primary `spAbundance` formulas are stored symbolically in the fit calls, but the authoritative script and metadata confirm the resolved formulas and MCMC settings.
- The matched `unmarked` models were rebuilt with `fit_unmarked_wolf_jackal_sensitivity.R`; fit objects, full coefficients, ecological coefficients, optimizer diagnostics, and session metadata are now saved. All optimizers returned convergence code 0 and ecological estimates reproduced the prior table to rounding.
- The focused `lfMsNMix` manuscript draft contains slightly stale site-grouped PPC values (0.427, 0.402, 0.326). Current saved output gives 0.429, 0.397, 0.330. Use the current output or round consistently to 0.43, 0.40, 0.33.
- The primary N-mixture intercepts remain nonconverged despite well-mixed focal slopes. Absolute abundance predictions that retain those intercepts are not authoritative.
- The beta-binomial fit object is valid, but loading `glmmTMB 1.1.13` warns that it was built with `TMB 1.9.18` while `TMB 1.9.17` is installed.

### Figures

- The retained analytical figures are present: spatial interaction curves, beta-binomial coefficient plot, gazelle-only context diel figure, and corrected overall-class diel figure.
- `make_publication_camera_array_map.R` provides a provisional array-and-strata map. It does not infer a regional or political boundary because no authoritative boundary layer is stored in the project.
- `make_overall_gazelle_diel_activity.R` is the authoritative standalone overall-class figure and parses camera clocks as Israel local time before apparent-solar-time conversion.
- Legacy predator-overlay figures remain in the repository and are excluded.

### Portability and software environment

- The retained authoritative scripts use project-relative paths. The only detected hard-coded workstation path is in obsolete `unmarked attempt.R`.
- R and package versions are recorded in `publication_audit_session_info.txt` and `publication_package_versions.csv`.
- Seeds are recorded for the primary Bayesian fits and simulation/bootstrap procedures in their scripts or metadata.
- `publication_reproducibility_checklist.csv` verifies that all 29 designated retained artifacts currently exist and identifies their generator scripts.

## F. Remaining Decisions Requiring Scientific Judgment

1. Correct or document camera 45 round dates and decide whether its two adult-male detections belong to round 2 or round 4.
2. Choose the reported effort total after that correction: 5,999 recorded days or the corrected date-resolved total.
3. Decide whether the curated count histories remain authoritative as-is or must be regenerated from the independent-event data with a documented round-boundary rule.
4. Define the manuscript's unidentifiable-gazelle category. Raw sex, life-stage, and distance/visibility labels do not map to one unambiguous total.
5. Confirm that subadult/not-adult males remain excluded from every focal demographic class and remove the unsupported number 51 unless its source is supplied.
6. Decide whether overall class-class diel overlap remains in the main Results or moves entirely to SI.
7. Decide whether Early/Late diel models support a necessary claim; otherwise retain them only in SI or omit them.
8. Decide whether the provisional camera-array map is sufficient or provide an authoritative regional boundary layer.
9. Rerun the archived `unmarked` sensitivity fits only if resolving camera 45 changes a focal count history.
10. Approve an environment repair or lockfile strategy for the `glmmTMB`/`TMB` mismatch before final model regeneration.

## G. Manuscript Numbers to Update

| Manuscript quantity | Remove or verify | Replacement |
|---|---:|---:|
| Fawns in retained count analyses | 195 | 199 |
| Adult females in retained count analyses | 2,135 | 1,951 |
| Adult males in retained count analyses | Any event-derived total | 2,355 |
| Not-adult males | 51 | Do not report until unit/provenance is resolved; current alternatives are 43 events/55 event maxima or 48 events/62 event maxima |
| Fawn:adult-female ratio | Any ratio using old totals | 0.102, or 10.2 per 100 adult females |
| Diel fawn sample | Individual count total | 180 independent events |
| Diel female sample | Individual count total | 1,466 independent events |
| Diel male sample | Individual count total | 2,260 independent events |
| `lfMsNMix` grouped PPC | 0.427, 0.402, 0.326 | 0.429, 0.397, 0.330, or 0.43, 0.40, 0.33 when rounded |

Do not finalize a camera-night value until the camera 45 discrepancy is resolved.

## H. Supporting Information Inventory

The detailed organization is in `publication_SI_inventory.md`. Core retained SI artifacts are:

- Data audit, count provenance, variable crosswalk, event exclusions, and round reconciliation
- Primary N-mixture coefficients, conditional effects, and MCMC diagnostics
- Expanded-model rejection diagnostics
- Full `unmarked` sensitivity coefficients
- Focused `lfMsNMix` loadings, residual correlations, diagnostics, and PPCs
- Complete beta-binomial coefficients and model-comparison sensitivities
- Overall gazelle class overlap and complete within-class context comparison table
- Early/Late gazelle-only results if retained

Gazelle-predator temporal-overlap files are excluded from every SI appendix.

## Repository Organization Recommendation

Do not move files until a retained-only driver uses explicit relative paths. The proposed structure is:

```text
analysis/01_data_audit/
analysis/02_spatial_primary_nmix/
analysis/03_spatial_sensitivity_unmarked/
analysis/04_class_correlations_lfmsnmix/
analysis/05_fawn_female_beta_binomial/
analysis/06_diel_activity/
figures/main/
figures/supplementary/
tables/main/
tables/supplementary/
models/final/
models/exploratory/
models/rejected/
archive/superseded_scripts/
archive/superseded_outputs/
archive/legacy_predator_overlap/
reports/
```

`publication_proposed_file_move_manifest.csv` records the non-destructive proposed mapping.

## Main Text Versus SI Recommendation

Keep the main Methods and Results centered on four components: separate-class `spAbundance` predator responses, the three secondary residual class correlations, the additive beta-binomial fawn-representation model, and gazelle-only diel distributions across predator contexts. Retain one brief `unmarked` robustness statement and the pooled descriptive ratio.

Move full coefficients, conditional effects, convergence diagnostics, rejected expanded models, all `unmarked` details, latent-factor diagnostics, complete beta-binomial adjustment terms, all context contrasts, and optional Early/Late diel analyses to SI. Exclude all gazelle-predator temporal overlap, dog PAMM, ratio exploration, occupancy, and superseded model-derived demographic surfaces from both.
