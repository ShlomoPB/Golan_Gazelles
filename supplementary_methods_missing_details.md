# Supplementary methods missing-details audit

Audit date: 2026-09-17. No model was refitted. Values were read from final
scripts, saved model objects, final CSV outputs, and the archived predator-index
workflow. `VERIFIED FINAL` means that both implementation and output are
traceable. `LIKELY FINAL BUT NEEDS CONFIRMATION` identifies a provenance or
script conflict. `OLD/OBSOLETE` identifies superseded work. `NOT FOUND` means
that the requested information is absent from the inspected files.

## Manuscript and supplementary-table alignment

| Check | Finding | Source | Location | Status |
|---|---|---|---|---|
| Appendix/table structure | The simplified supplement is internally aligned: Appendix S1 = Tables S1-S4; Appendix S2 = Tables S5-S6; Appendix S3 = Tables S7-S10. | `supplementary_tables/supplementary_materials_final.qmd`; `supplementary_table_qc.md` | Appendix and table headings | VERIFIED FINAL |
| Numerical provenance | Saved final objects reproduce the reported S1-S10 values to manuscript rounding. | `supplementary_table_qc.md` | `Numerical validation` | VERIFIED FINAL |
| Pooled selected-GAM EDF | The selected pooled wolf-only model has wolf x time EDF = 3.869643 (reported as 3.870) and P = 0.000740. EDF = 3.041 belongs to the full interaction robustness fit, not the selected pooled model. | `definitive_hierarchical_diel_audit_model_summary.csv`; `supplementary_tables/table_S08_diel_effects_and_robustness.csv` | Pooled selected model / Table S8 | VERIFIED FINAL; current manuscript value needs correction if it still says 3.04 |
| Early/Late wording | Code assigns camera-round midpoints on or before 2020-10-17 to Early and later midpoints to Late. No midpoint occurred exactly on the boundary, so the fitted result is unaffected, but wording such as "before 17 October / 17 October onward" does not exactly describe the code. | `analyze_fawn_female_beta_binomial_wolf_jackal.R` | `build_season_lookup()`, lines 36-65 | VERIFIED FINAL; wording needs alignment |
| Joint-GAM robustness | The all-class wolf x jackal x time omnibus was audited in 7/7 specifications. Pairwise demographic contrasts were fitted only for the primary apparent-solar model. | `calibrate_joint_demographic_diel_GAM.R`; `joint_demographic_penalty_calibrated_outputs/penalty_calibrated_sensitivity_tests.csv` | lines 357-401 and 905-1008 | VERIFIED FINAL |
| Current DOCX access | The supplied DOCX was exclusively locked by another process during this audit. The alignment statement above uses the prior successful manuscript extraction recorded in `supplementary_table_qc.md`, not a new extraction of the open file. | Current DOCX; `supplementary_table_qc.md` | file-access audit | LIKELY FINAL BUT NEEDS CONFIRMATION |

# Appendix S1 - Spatial analyses

## A. Predator activity indices

### Final model and index construction

| Requested item | Final value / setting | Source file | Function / object and location | Status |
|---|---|---|---|---|
| Wolf function and distribution | `unmarked::pcount(..., mixture = "ZIP")`; zero-inflated Poisson. | External `Part 2 - N-mixture Models.R` under `Golan Wolf Study/Paper 2 - Trophic Thunder/Golan Codes and Dataframes for Github` | `wm1`, `models_detection`, and `Wolf_models_state`, lines 48-71 and 231-296 | LIKELY FINAL BUT NEEDS CONFIRMATION |
| Jackal function and distribution | `unmarked::pcount()` with no `mixture` argument; the `unmarked` default is Poisson. | Same external script | `jm1`, `jackal_models_detection`, and `Jackal_models_state`, lines 82-105 and 313-358 | LIKELY FINAL BUT NEEDS CONFIRMATION |
| Final state formula | `~ 1` for both species. | Same external script | `Wolf_models_state$Null`, line 234; `Jackal_models_state$jNull`, line 315 | LIKELY FINAL BUT NEEDS CONFIRMATION |
| Final detection formula | `~ Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera` for both species. | Same external script; saved detection-selection CSVs | `wm10`, line 71; `jm10`, line 104; state nulls at lines 234 and 315 | VERIFIED by saved output; archived script needs correction |
| Candidate-selection method | Explicit detection candidate sets were fitted with a constant state formula and ranked with `unmarked::fitList()` plus `AICcmodavg::modSel()`. This is an explicit candidate set, not sequential stepwise selection. | Same external script | lines 51-78 and 85-111 | VERIFIED by saved output |
| Selected wolf detection model | `wm10`, AIC 1803.803091, delta AIC 0, weight 0.449624. Runner-up `wm2`, AIC 1804.833953, delta 1.030862. All saved candidate rows have `Converge = 0`. | External `Wolf_Detection_Functions_All.csv` | rows `wm10` and `wm2` | VERIFIED FINAL |
| Selected jackal detection model | `jm10`, AIC 9304.758208, delta AIC 0, weight 1.000. Runner-up `jm9`, AIC 9567.482461, delta 262.724253. All saved candidate rows have `Converge = 0`. | External `Jackal_Detection_Functions_All.csv` | rows `jm10` and `jm9` | VERIFIED FINAL |
| Empirical-Bayes extraction | Intended code: `re <- ranef(model)`; `posteriorSamples(re, nsims = 500)`; site posterior mean, SD, and 2.5/97.5% quantiles. | Same external script | `estimate_site_activity()`, lines 184-217 | LIKELY FINAL BUT NEEDS CONFIRMATION |
| Pre-z transformation | Site posterior means were divided by the species-specific maximum: `activity_norm <- mean / max(mean)`. The archived source CSVs call these fields `wolf_abundance_stand` and `jackal_abundance_stand`; they are max-normalized 0-1 indices, not z scores. | External script; external `Null Wolf abundance with SE.csv` and `Null Jackal abundance with SE.csv`; `make_camera_wolf_jackal_abundance_coordinate_plot.R` | external lines 193-213; validation lines 32-82 and 110-152 | VERIFIED FINAL |
| Downstream z-standardization | `(x - mean(x)) / sd(x)` separately for wolf and jackal. The 60-site pre-z means/SDs are wolf 0.578505/0.201986 and jackal 0.500347/0.226118, yielding sample mean 0 and sample SD 1 after scaling. | `analyze_wolf_jackal_interaction_nmix.R`; `definitive_hierarchical_diel_GAM.R`; `SiteCovs_060324.csv` | `z_scale()`, lines 48-54 and site mutation lines 117-122; GAM lines 220-224 | VERIFIED FINAL |
| Wolf-jackal correlation | Pearson r = -0.1826892895, t = -1.415136, df = 58, P = 0.1623756, n = 60. | `SiteCovs_060324.csv`; independently reproduced in this audit | columns `Null.Wolf.Abundance`, `null.jackal.abundance` | VERIFIED FINAL |
| Sampling dimensions | 60 sites and seven repeated count occasions for each predator. | External `Wolf Detection History for Manuscript.csv`; `Canis_aureus_Abundance.csv` | 60 rows x 7 count columns | VERIFIED FINAL |
| Other settings | No `K`, optimizer, starting values, or goodness-of-fit settings are supplied in the archived `pcount()` calls; package defaults were therefore used. The `unmarked` package version used for the original fit is not recorded. | External script | all predator `pcount()` calls | NOT FOUND for exact package version/default realization |

### Exact detection candidate sets

The state side of every detection candidate was `~ 1`. All candidates and AIC
values are retained in the two external `*_Detection_Functions_All.csv` files.

| Species | Candidate names and detection formulas | Source / location | Status |
|---|---|---|---|
| Wolf | `Null: ~1`; `Agri: ~Agricultural.Presence`; `Cattle: ~Cattle.Sightings`; `Model: ~Camera.Model`; `Delay: ~Capture.Delay`; `Veg: ~Vegetation.Structure`; `Path: ~Path.Type`; `BehindCam: ~Area.Behind.Camera`; `Bottleneck: ~Bottleneck.`; `wm2: ~Path.Type + Camera.Model + Bottleneck.`; `wm3: ~Path.Type + Camera.Model`; `wm4/wm7: ~Path.Type + Camera.Model + Cattle.Sightings`; `wm5: ~Path.Type + Camera.Model + Area.Behind.Camera`; `wm6: ~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera`; `wm8: ~Cattle.Sightings + Camera.Model + Bottleneck.`; `wm85: ~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck.`; `wm87: ~Cattle.Sightings + Camera.Model + Agricultural.Presence + Bottleneck.`; `wm9: ~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck.`; `wm10: ~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera`. | External script lines 48-72; full AICs in `Wolf_Detection_Functions_All.csv` | VERIFIED by saved output |
| Jackal | `Null: ~1`; `Cattle: ~Cattle.Sightings`; `Model: ~Camera.Model`; `Sensitivity: ~Camera.Sensitivity`; `Delay: ~Capture.Delay`; `Veg: ~Vegetation.Structure`; `Path: ~Path.Type`; `BehindCam: ~Area.Behind.Camera`; `Bottleneck: ~Bottleneck.`; and `jm2`-`jm10`/`jm85`/`jm87` matching the corresponding multi-covariate wolf formulas above. | External script lines 82-105; full AICs in `Jackal_Detection_Functions_All.csv` | VERIFIED by saved output |

### Predator-index provenance conflicts

| Conflict | Evidence | Status |
|---|---|---|
| Wolf selection call is not executable as archived | Line 75 reads `fitList(c(Null = list(wm1) models_detection))`, missing a comma. The saved AIC CSV proves that the intended set was fitted in another/corrected execution. | LIKELY FINAL BUT NEEDS CONFIRMATION |
| Jackal empirical-Bayes call is not executable as archived | The state-list object is `Jackal_models_state` and its null member is `jNull`, but line 222 calls `jackal_models_state$Null`. The saved index CSV and exact downstream match prove that a corrected object was used somewhere, but that corrected call is not present. | LIKELY FINAL BUT NEEDS CONFIRMATION |

### Citations

| Citation | Exact project entry | Source / location | Status |
|---|---|---|---|
| Royle 2004 | Royle, J. A. 2004. N-mixture models for estimating population size from spatially replicated counts. *Biometrics* 60:108-115. DOI 10.1111/j.0006-341X.2004.00142.x. | `conbio_methods_results_raw_package.md:110` | VERIFIED FINAL |
| Barker et al. 2018 | Barker, R. J., M. R. Schofield, W. A. Link, and J. R. Sauer. 2018. On the reliability of N-mixture models for count data. *Biometrics* 74:369-377. DOI 10.1111/biom.12734. | `conbio_methods_results_raw_package.md:102` | VERIFIED FINAL |
| Kellner et al. 2023 | Kellner, K. F., et al. 2023. The unmarked R package: Twelve years of advances in occurrence and abundance modelling in ecology. *Methods in Ecology and Evolution* 14:1408-1415. DOI 10.1111/2041-210X.14123. | `conbio_methods_results_raw_package.md:106` | VERIFIED FINAL |

## B. Joint demographic `lfMsNMix` model

| Requested item | Final value / setting | Source file | Function / object and location | Status |
|---|---|---|---|---|
| Abundance formula | `~ Open.Spaces.z + Null.Wolf.Abundance.z * null.jackal.abundance.z`. | `fit_focused_multiclass_nmix.R` | `abund_formula`, lines 53-56 | VERIFIED FINAL |
| Detection formula | `~ log_days_operational.z + cattle_sightings.z`. | same | `det_formula`, line 56 | VERIFIED FINAL |
| Distribution | Negative binomial, `family = "NB"`, with class-specific `kappa`. | same | `common_args`, lines 127-142 | VERIFIED FINAL |
| Candidate structures | No-factor: `msNMix()` (`q = 0` conceptually); one-factor: `lfMsNMix(n.factors = 1)`; two-factor: `lfMsNMix(n.factors = 2)`. | same | lines 145-152 | VERIFIED FINAL |
| Table S1 model | Focused two-factor model, class order `fawn-male-female`; object `gazelle_lf2_NMix_NB_open_wolf_jackal_interaction_order-fawn-male-female_2000b_25bl_10th_3ch_fit.rds`. | model RDS; `publication_focused_lfmsnmix_class_correlations.csv` | saved final object / Table S1 generator | VERIFIED FINAL |
| Why one factor was rejected | With three classes and one factor, residual correlations have absolute value 1 by construction; the two free loadings switched signs among chains, R-hat 22.68-23.35. | `focused_multiclass_spabundance_model_review.md` | line 43 | VERIFIED FINAL |
| Priors | Community abundance coefficients Normal(0, variance 100); community detection coefficients Normal(0, variance 2.72); abundance and detection hypervariances inverse-gamma(a = 0.1, b = 0.1); class-specific `kappa` Uniform(0,10). | `fit_focused_multiclass_nmix.R` | `priors`, lines 68-83 | VERIFIED FINAL |
| Initialization | All coefficient vectors 0; `tau.sq.beta = tau.sq.alpha = 0.5`; `kappa = 1`; latent N initialized to each class-site observed maximum. | same | `inits`, lines 85-96 | VERIFIED FINAL |
| MCMC | 3 chains; 2,000 batches x 25 = 50,000 iterations/chain; burn-in 5,000; thinning 10; 13,500 retained pooled draws. Seeds: no-factor 20260817, one-factor 20260818, two-factor 20260819. | script plus metadata RDS files | lines 6-17, 41-49, 122-142, metadata | VERIFIED FINAL |
| Tuning / computing | `beta = alpha = kappa = 0.5`; latent-factor fits also `lambda = w = 0.5`; one OpenMP thread; `n.report = 200`. | same | lines 98-102 and 127-142 | VERIFIED FINAL |
| Diagnostics used | Saved `spAbundance` R-hat and ESS for community, class, kappa, and free loading parameters. | `evaluate_focused_multiclass_nmix.R` | lines 90-120 | VERIFIED FINAL |
| Posterior predictive check | `ppcAbund(..., type = "marginal")`; Freeman-Tukey at groups 0 and 1 and chi-squared at group 2; Bayesian P = `colMeans(fit.y.rep > fit.y)`. | same | lines 200-228 | VERIFIED FINAL |
| WAIC | Disabled by default (`run_waic = false`) and not used. The review records that weakly identified latent N required very high integration ceilings, making full common-Nmax WAIC impractical and scientifically difficult to defend. No final WAIC file exists. | `evaluate_focused_multiclass_nmix.R`; `focused_multiclass_spabundance_model_review.md` | lines 11, 39, 177-197; review line 47 | VERIFIED FINAL |
| Retained-model limitations | Overall R-hat range 1.0002-1.8151 and ESS 9.29-13,500. Class abundance intercept R-hat 1.345-1.581 (ESS 21.0-30.4); detection intercept R-hat 1.497-1.815 (ESS 9.29-16.68); hypervariance R-hat up to 1.279. Free loadings mixed adequately. | `publication_focused_lfmsnmix_diagnostics.csv` | final LF2 rows | VERIFIED FINAL |

### Final LF2 PPC values

All rows below are `VERIFIED FINAL` from
`publication_focused_lfmsnmix_ppc.csv`, created by
`evaluate_focused_multiclass_nmix.R:200-228`.

| Statistic / group | Juvenile | Adult female | Adult male |
|---|---:|---:|---:|
| Freeman-Tukey, group 0 | 0.046519 | 0.002000 | 0.001852 |
| Freeman-Tukey, group 1 | 0.407778 | 0.318296 | 0.402370 |
| Chi-squared, group 2 | 0.130741 | 0.147630 | 0.027259 |

### Free-loading diagnostics and Table S1 extraction

| Parameter / result | Value | Source / code | Status |
|---|---|---|---|
| `male-1` loading | R-hat 1.002074; ESS 709.846 | `publication_focused_lfmsnmix_diagnostics.csv` | VERIFIED FINAL |
| `female-1` loading | R-hat 1.000819; ESS 393.761 | same | VERIFIED FINAL |
| `female-2` loading | R-hat 1.014726; ESS 322.860 | same | VERIFIED FINAL |
| Residual-correlation calculation | For every posterior draw, reshape `lambda.samples`, compute `tcrossprod(lambda_i)`, then `cov2cor()`. | `evaluate_focused_multiclass_nmix.R:123-166` | VERIFIED FINAL |
| Juvenile-adult male | mean 0.657368; 95% CrI 0.435371-0.806227 | `publication_focused_lfmsnmix_class_correlations.csv` | VERIFIED FINAL |
| Juvenile-adult female | mean 0.863347; 95% CrI 0.523338-0.998449 | same | VERIFIED FINAL |
| Adult male-adult female | mean 0.906536; 95% CrI 0.712481-0.999566 | same | VERIFIED FINAL |

Doser citation: Doser, J. W., A. O. Finley, M. Kery, and E. F. Zipkin. 2024.
spAbundance: An R package for single-species and multi-species spatially
explicit abundance models. *Methods in Ecology and Evolution* 15:1024-1033.
DOI 10.1111/2041-210X.14332 (`conbio_methods_results_raw_package.md:104`).
Status: `VERIFIED FINAL`.

## C. Class-specific `spAbundance::NMix` models

### Shared implementation

| Requested item | Final value / setting | Source / location | Status |
|---|---|---|---|
| Function and formulas | One `spAbundance::NMix()` per class; abundance `~ Null.Wolf.Abundance.z * null.jackal.abundance.z`; detection `~ log_days_operational.z + cattle_sightings.z`. | `analyze_wolf_jackal_interaction_nmix.R:181-213`, loop over class definitions | VERIFIED FINAL |
| Distribution | `family = "NB"`; package-defined negative-binomial model with one class-specific `kappa`. No custom variance parameterization was coded. | same, lines 198-212 | VERIFIED FINAL |
| Priors | Abundance coefficients Normal(0, variance 100); detection coefficients Normal(0, variance 2.72); `kappa` Uniform(0,10). | same, lines 185-191 | VERIFIED FINAL |
| MCMC | 3 chains; 50,000 iterations/chain (2,000 x 25); burn-in 5,000; thinning 10; 13,500 retained draws/class; seed 20260623. | same, lines 23-35, 142, 197-212 | VERIFIED FINAL |
| Initialization / tuning | Coefficients 0, `kappa = 1`, latent N = site maximum; tuning `beta = 0.1`, `alpha = 0.1`, `kappa = 0.5`; one OpenMP thread; `verbose = FALSE`. | same, lines 185-212 | VERIFIED FINAL |
| Posterior-draw source | `publication_wolf_jackal_interaction_nmix_spabundance_fits.rds`, each class object's `beta.samples`. | same line 312; `extract_final_spabundance_predator_combinations.R:8-15,62-84` | VERIFIED FINAL |

### Exact R-hat and ESS

All rows are `VERIFIED FINAL` from
`publication_wolf_jackal_interaction_nmix_spabundance_diagnostics.csv`.

| Class | Component / term | R-hat | ESS |
|---|---|---:|---:|
| Juvenile | abundance intercept | 1.207430 | 29.16 |
| Juvenile | detection intercept | 1.230864 | 23.86 |
| Juvenile | wolf | 1.002421 | 5454.74 |
| Juvenile | jackal | 1.004166 | 3846.93 |
| Juvenile | wolf x jackal | 1.003836 | 3548.63 |
| Adult female | abundance intercept | 2.885785 | 10.72 |
| Adult female | detection intercept | 3.491553 | 7.47 |
| Adult female | wolf | 1.004164 | 9092.53 |
| Adult female | jackal | 1.000744 | 8124.88 |
| Adult female | wolf x jackal | 1.001827 | 6834.82 |
| Adult male | abundance intercept | 1.291654 | 14.35 |
| Adult male | detection intercept | 1.379977 | 9.90 |
| Adult male | wolf | 1.001038 | 9537.76 |
| Adult male | jackal | 1.002949 | 8105.73 |
| Adult male | wolf x jackal | 1.001287 | 7188.02 |

Focal-slope range: R-hat 1.000744-1.004166; ESS 3548.63-9537.76.
Intercepts did not mix adequately; interpretation was therefore restricted to
well-mixed slopes and intercept-free relative predictions.

### Conditional contrasts, normalization, and Figure 3

| Requested item | Exact implementation | Source / location | Status |
|---|---|---|---|
| Prediction combinations | Wolf and jackal each evaluated at -1 and +1 SD. | `extract_final_spabundance_predator_combinations.R:32-38` | VERIFIED FINAL |
| Relative abundance | Intercept omitted; for each draw `R(w,j) = exp(beta_w*w + beta_j*j + beta_wj*w*j)`. Because all design terms are zero at `(w,j)=(0,0)`, this is exactly relative to the mean-predator condition, which equals 1 draw by draw. | same, lines 62-96 | VERIFIED FINAL |
| Jackal -1 to +1 contrast at wolf level `w` | `R(w,-1)/R(w,+1) = exp[-2(beta_j + w*beta_wj)]`. Draws are divided directly before taking medians and 95% CrIs. | same, lines 106-123 and 132-135 | VERIFIED FINAL |
| Ratio-of-ratios | `[R(-1,-1)/R(-1,+1)] / [R(+1,-1)/R(+1,+1)] = exp(4*beta_wj)`. | same, lines 130-153 | VERIFIED FINAL |
| Figure 3 sources | Fit RDS above; saved grid `publication_spabundance_wolf_jackal_relative_prediction_grid.csv`; figure script `make_spabundance_relative_wolf_jackal_figure.R`. | figure script lines 9-39 and 109-216 | VERIFIED FINAL |

### Exploratory habitat/human fits

| Requested item | Finding | Source / location | Status |
|---|---|---|---|
| Audited fits | Eighteen saved superseded `NMix` fits (six/class) were audited without refitting. | `audit_expanded_nmix_models.R:26-88`; `publication_expanded_nmix_rejection_diagnostics.csv` | OLD/OBSOLETE |
| Juvenile formulas | Base `~ Open.Spaces.z`; base + human; base + boar; base + human + boar; and `~ Open.Spaces.z + null.jackal.abundance.z + null.boar.abundance.z + human.abundance.z + Null.Wolf.Abundance.z`; one duplicate base file. | rejection CSV, juvenile rows | OLD/OBSOLETE |
| Adult-female formulas | Same templates, replacing base with `Nature.Reserves.z`. | rejection CSV, female rows | OLD/OBSOLETE |
| Adult-male formulas | Same templates, replacing base with `Built.and.Roads.z`. | rejection CSV, male rows | OLD/OBSOLETE |
| Detection formula | `~ log_days_operational.z + cattle_sightings.z` in all audited fits. | rejection CSV | OLD/OBSOLETE |
| Rejection reason | Every saved fit had at least one R-hat >1.1 or ESS <100 under the audit rule, generally abundance and/or detection intercepts. Across all 18 fits, max R-hat reached 3.088 and minimum ESS 4.982. Even the female all-predictor fit had intercept ESS <100. PPC was not retained. | `audit_expanded_nmix_models.R:35-81`; rejection CSV | VERIFIED diagnostic rejection |
| Role in final selection | No AIC, WAIC, DIC, LOO, or PPC comparison selected the final predator-only model. The expanded models were exploratory and rejected/superseded on diagnostics. | audit script and rejection CSV | VERIFIED FINAL |

## D. Matched `unmarked` sensitivity

| Requested item | Final value / setting | Source / location | Status |
|---|---|---|---|
| Function/distribution | `unmarked::pcount(..., mixture = "NB")`, fitted separately for juvenile, adult female, and adult male. | `fit_unmarked_wolf_jackal_sensitivity.R:132-173` | VERIFIED FINAL |
| State formula | `~ wolf_z * jackal_z`. | same, lines 163-172 | VERIFIED FINAL |
| Detection formula | `~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera`. | same, lines 163-172 | VERIFIED FINAL |
| Differences from Bayesian fits | Maximum-likelihood `unmarked`; camera-design detection formula instead of effort + cattle; `K = max(150, max(y)+50)`; no priors/MCMC; Wald normal CIs/P values. | same, lines 91-109 and 153-200 | VERIFIED FINAL |
| Fit diagnostics | Juvenile K=150, AIC=667.6463; female K=246, AIC=2547.3209; male K=150, AIC=2618.0323. Each used 60 sites x 7 repeats and convergence code 0. | `publication_unmarked_wolf_jackal_sensitivity_diagnostics.csv` | VERIFIED FINAL |

All Table S4 values below are `VERIFIED FINAL` from
`publication_unmarked_wolf_jackal_sensitivity_ecological_coefficients.csv`.

| Class | Term | Estimate | SE | 95% CI | z | P |
|---|---|---:|---:|---:|---:|---:|
| Juvenile | wolf | 0.265344 | 0.224593 | -0.174849, 0.705538 | 1.181447 | 0.237425 |
| Juvenile | jackal | -0.742420 | 0.293881 | -1.318416, -0.166424 | -2.526263 | 0.011528 |
| Juvenile | wolf x jackal | 0.485236 | 0.288202 | -0.079629, 1.050100 | 1.683667 | 0.092246 |
| Adult female | wolf | 0.354697 | 0.160464 | 0.040194, 0.669200 | 2.210451 | 0.027074 |
| Adult female | jackal | -0.546524 | 0.171691 | -0.883032, -0.210016 | -3.183182 | 0.001457 |
| Adult female | wolf x jackal | 0.444205 | 0.186118 | 0.079421, 0.808989 | 2.386691 | 0.017001 |
| Adult male | wolf | 0.318499 | 0.137890 | 0.048240, 0.588757 | 2.309809 | 0.020899 |
| Adult male | jackal | -0.541790 | 0.144677 | -0.825353, -0.258228 | -3.744819 | 0.000181 |
| Adult male | wolf x jackal | 0.422511 | 0.155286 | 0.118157, 0.726866 | 2.720860 | 0.006511 |

# Appendix S2 - Juvenile:female ratio

## Data, family, and coding

| Requested item | Final value / setting | Source / location | Status |
|---|---|---|---|
| Response | `cbind(fawn, female)` at camera-round level; rows with `fawn + female == 0` excluded. It models juvenile representation among classified juvenile plus adult-female counts. | `analyze_fawn_female_beta_binomial_wolf_jackal.R:23-34,68-73,106-124` | VERIFIED FINAL |
| Sample | 221 informative camera-rounds at 50 cameras in 12 polygons. | model data; `mountain_gazelle_manuscript_final_verification.md:84` | VERIFIED FINAL |
| Family/link | `glmmTMB::betabinomial(link = "logit")`. | analysis script lines 147-188 | VERIFIED FINAL |
| Random effects | `(1 | Polygon.Name / camera_id)`, equivalent to polygon and camera-within-polygon random intercepts. | same | VERIFIED FINAL |
| Fixed adjustment set | `open_z + built_z + reserve_z + human_z + path_type + camera_model + vegetation_structure + area_behind_camera`; held in all six candidates. | same, lines 140-188 | VERIFIED FINAL |
| Standardized variables | `wolf_z`, `jackal_z`, `open_z`, `built_z`, `reserve_z`, `human_z`, each constructed by sample mean/SD scaling. | same, lines 75-85 | VERIFIED FINAL |
| Survey-period coding | Camera-round midpoint; study midpoint `split_date`; `mid_date <= 2020-10-17` = `Early (hotter)` reference, later dates = `Late (mild-cold)`. No midpoint fell exactly on 17 October. | `build_season_lookup()`, lines 36-65; final QMD Appendix S2 | VERIFIED FINAL |
| Requested verbal split | "Early before 17 Oct; Late 17 Oct onward" is not the literal code. Use "Early on or before 17 Oct; Late after 17 Oct" unless an a priori definition is intentionally being stated separately. | same | LIKELY FINAL BUT NEEDS CONFIRMATION |

## Six candidate models

All formulas use the same random effects and adjustment set shown above. Values
are from the saved `publication_fawn_female_beta_binomial_wolf_jackal_model_fits.rds`.
All six fits had convergence code 0 and positive-definite Hessians. Status for
all rows: `VERIFIED FINAL`.

| Object | Exact varying part of formula | df | logLik | AIC | delta AIC |
|---|---|---:|---:|---:|---:|
| `overall_fit` | `wolf_z * jackal_z` | 24 | -198.456261 | 444.912522 | 3.837036 |
| `overall_additive_fit` | `wolf_z + jackal_z` | 23 | -198.611202 | 443.222404 | 2.146918 |
| `seasonal_fit` | `wolf_z * jackal_z * season` | 28 | -194.097040 | 444.194079 | 3.118593 |
| `seasonal_no_predator_interaction_fit` | `(wolf_z + jackal_z) * season` | 26 | -196.052413 | 444.104827 | 3.029340 |
| `seasonal_stable_response_fit` | `wolf_z * jackal_z + season` | 25 | -196.360881 | 442.721763 | 1.646276 |
| `seasonal_additive_stable_fit` | `wolf_z + jackal_z + season` | 24 | -196.537743 | 441.075486 | 0.000000 |

Exact full calls are in
`analyze_fawn_female_beta_binomial_wolf_jackal.R:147-188`; the first additive
fit is created by `update(overall_fit, . ~ . - wolf_z:jackal_z)` at line 155.

## Nested tests and retained-model quantities

| Requested result | Comparison / value | Source / location | Status |
|---|---|---|---|
| Common Early/Late shift | Add `season` to additive model: LR = 4.146918, df = 1, P = 0.041710; larger-minus-smaller delta AIC = -2.146918. | `publication_fawn_female_beta_binomial_wolf_jackal_likelihood_ratio_tests.csv`; script lines 212-256 | VERIFIED FINAL |
| Wolf x jackal interaction | Add interaction to retained season-adjusted additive model: LR = 0.353724, df = 1, P = 0.552013; delta AIC = +1.646276. | `publication_fawn_female_beta_binomial_interaction_sensitivity_model_comparison.csv` | VERIFIED FINAL |
| Period-specific additive slopes | Compare common additive slopes with `(wolf_z + jackal_z) * season`: LR = 0.970660, df = 2, P = 0.615494; delta AIC = +3.029340. | likelihood-ratio CSV | VERIFIED FINAL |
| Dispersion | Beta-binomial dispersion = 41.813302. | fit diagnostics CSV, retained row | VERIFIED FINAL |
| Random effects | Camera within polygon variance 0.182006, SD 0.426621; polygon variance 1.008e-10, SD 1.004e-5. | `publication_fawn_female_beta_binomial_coefficient_random_effects.csv` | VERIFIED FINAL |
| Interaction OR sensitivity | `seasonal_stable_response_fit`; interaction beta = -0.151229, SE = 0.245838, OR = 0.859651, 95% CI 0.530959-1.391822, P = 0.538452. | `publication_fawn_female_beta_binomial_interaction_sensitivity_coefficients.csv` | VERIFIED FINAL |
| Other retained sensitivity | The six-model comparison and the two interaction/period-slope LRT extensions above are the retained final sensitivity workflow. Files labelled `fixed_season_sensitivity` and earlier pooled/focused/posterior-ratio analyses are not the final selected workflow. | final script and final QMD | VERIFIED FINAL / OLD files excluded |

## Full retained fixed-effect table

All values are `VERIFIED FINAL` from
`publication_fawn_female_beta_binomial_coefficient_details.csv`, model
`Preferred additive-season`.

| Term (reference where applicable) | Beta | SE | OR | 95% CI | P |
|---|---:|---:|---:|---:|---:|
| Intercept | -4.085386 | 0.629914 | 0.016817 | 0.004893-0.057801 | <0.000001 |
| Wolf activity | -0.073787 | 0.143856 | 0.928869 | 0.700651-1.231424 | 0.608005 |
| Jackal activity | -0.059050 | 0.204763 | 0.942660 | 0.631041-1.408162 | 0.773055 |
| Late vs Early | 0.450592 | 0.224285 | 1.569240 | 1.011053-2.435595 | 0.044536 |
| Open-space cover | 0.440523 | 0.204757 | 1.553520 | 1.039977-2.320649 | 0.031441 |
| Built-and-road cover | 0.005407 | 0.360555 | 1.005422 | 0.495950-2.038258 | 0.988035 |
| Nature-reserve cover | 0.009663 | 0.233680 | 1.009710 | 0.638681-1.596279 | 0.967016 |
| Human activity | 0.002340 | 0.178575 | 1.002343 | 0.706335-1.422402 | 0.989543 |
| Animal path without vehicles vs 4x4 road | 0.411292 | 0.832066 | 1.508765 | 0.295361-7.707091 | 0.621093 |
| Dirt path <1 m vs 4x4 road | 1.129674 | 0.802136 | 3.094647 | 0.642420-14.907439 | 0.159033 |
| Dirt path >1 m vs 4x4 road | 0.172504 | 0.365184 | 1.188276 | 0.580853-2.430908 | 0.636660 |
| Paved road vs 4x4 road | -0.126905 | 0.523877 | 0.880818 | 0.315466-2.459344 | 0.808593 |
| Browning PRO X long vs Browning MAX | 0.984662 | 0.628609 | 2.676908 | 0.780820-9.177320 | 0.117251 |
| Browning PRO X normal vs Browning MAX | 0.911431 | 0.562178 | 2.487879 | 0.826597-7.487980 | 0.104964 |
| Bushnell Trophy Cam vs Browning MAX | 0.101836 | 0.556670 | 1.107202 | 0.371860-3.296659 | 0.854847 |
| Shrubs <2 m vs grassland | 2.459745 | 0.712133 | 11.701828 | 2.897835-47.253478 | 0.000552 |
| Shrubs >2 m vs grassland | 0.666588 | 0.635578 | 1.947580 | 0.560377-6.768775 | 0.294275 |
| Overhead tree cover vs grassland | 0.331179 | 0.489890 | 1.392609 | 0.533122-3.637741 | 0.499023 |
| Fence vs clear behind camera | 0.766953 | 0.621766 | 2.153195 | 0.636539-7.283526 | 0.217387 |
| Impassable topography vs clear | 0.805444 | 0.571847 | 2.237689 | 0.729515-6.863814 | 0.158984 |
| Vegetation obstruction vs clear | 0.370296 | 0.447448 | 1.448163 | 0.602479-3.480911 | 0.407913 |

Brooks citation: Brooks, M. E., et al. 2017. glmmTMB balances speed and
flexibility among packages for zero-inflated generalized linear mixed
modeling. *The R Journal* 9:378-400. DOI 10.32614/RJ-2017-066
(`conbio_methods_results_raw_package.md:165`). Status: `VERIFIED FINAL`.

# Appendix S3 - Diel activity analyses

## A. Solar-time conversion

| Requested item | Exact implementation | Source / location | Status |
|---|---|---|---|
| Input timezone | Camera-clock strings parsed as `Asia/Jerusalem`; date/clock corrections were applied before time-scale conversion. | `definitive_hierarchical_diel_GAM.R:126-138`; config lines 17-35 | VERIFIED FINAL |
| Equation of time | `B = 2*pi*(DOY-81)/364`; `E = 9.87*sin(2B) - 7.53*cos(B) - 1.5*sin(B)` minutes. | same, `equation_of_time()`, lines 141-144 | VERIFIED FINAL |
| Apparent-solar formula | `hour = (UTC decimal hour + longitude/15 + E/60) %% 24`. Longitude is degrees east and contributes four minutes per degree. | same, `apparent_solar_hour()`, lines 146-152 | VERIFIED FINAL |
| Noon alignment | The longitude and equation-of-time corrections define apparent solar noon as 12:00; the transformation is additive and preserves equal clock-hour durations. | same | VERIFIED FINAL |
| Package | No package performs the primary conversion; it is explicit base-R arithmetic. | same | VERIFIED FINAL |
| Sunrise/sunset alternative | Local clock hour is converted to radians, passed with POSIX timestamp and WGS84 camera coordinates to `overlap::sunTime()`, then converted back by `*12/pi`. | same, `sun_time_hour()`, lines 154-163 | VERIFIED FINAL |

Nouvellet citation: Nouvellet, P., G. S. A. Rasmussen, D. W. Macdonald, and
F. Courchamp. 2012. Noisy clocks and silent sunrises: measurement methods of
daily activity pattern. *Journal of Zoology* 286:179-184. DOI
10.1111/j.1469-7998.2011.00864.x (`conbio_methods_results_raw_package.md:231`).
Status: `VERIFIED FINAL`.

## B. Day/night and kernel overlap

| Requested item | Exact implementation | Source / location | Status |
|---|---|---|---|
| Solar elevation | Custom NOAA-style calculation of fractional-year angle, equation of time, declination, true solar minutes, hour angle, and zenith; no external solar package. | `demographic_nocturnality_astronomical.R`, `solar_elevation()`, lines 51-75 | VERIFIED FINAL |
| Night definition | `solar_elevation_degrees < 0`; elevation 0 or greater is daylight. | same, lines 281-290 | VERIFIED FINAL |
| Camera bootstrap | Construct camera x class total/night matrices for all 60 positive-effort cameras, sample 60 cameras with replacement, retain all detections within sampled cameras, and calculate class proportions and pairwise differences. Percentile 2.5/97.5% intervals. | same, lines 414-527 | VERIFIED FINAL |
| Day/night replicates | 10,000; seed 20260831. | same, lines 11-21 and 453-466 | VERIFIED FINAL |
| Formal night model | Retained `night ~ demographic_class + round_factor + (1|camera_factor) + (1|camera_round)`, binomial logit in `glmmTMB`. It was retained over camera-only after convergence, variance, and AIC checks. | same, lines 561-628 | VERIFIED FINAL |
| Pairwise night contrasts | `emmeans(model, ~ demographic_class, weights="equal")`; three planned log-odds contrasts; raw P values adjusted with `p.adjust(method="holm")`. | same, `extract_pairwise()`, lines 185-212 | VERIFIED FINAL |
| Kernel density | `activity::fitact(class_times, reps=199, sample="data")` for class density curves. | `make_overall_gazelle_diel_activity.R:210-214` | VERIFIED FINAL |
| Overlap / Dhat4 | `overlap::overlapEst(time_a,time_b,type=estimator)`; `Dhat1` if the smaller class has <75 events, otherwise `Dhat4`. All final class pairs used `Dhat4`. | same, lines 97-105; final overlap CSV | VERIFIED FINAL |
| Dhat4 CI bootstrap | Sample the 60 study cameras with replacement 500 times, concatenate all class events from each sampled camera, recalculate the same overlap estimator, and use 2.5/97.5% percentiles. | same, lines 101-142; config lines 14-38 | VERIFIED FINAL |

Final Dhat4 values (`VERIFIED FINAL` from
`publication_overall_gazelle_diel_activity_class_overlap.csv`): juvenile-female
0.787259 (0.695279-0.855186), juvenile-male 0.742062
(0.650439-0.782342), female-male 0.947562 (0.858529-0.960255); 500 successful
bootstrap replicates each.

Relevant citations (`VERIFIED FINAL`): Ridout, M. S., and M. Linkie. 2009.
*Journal of Agricultural, Biological, and Environmental Statistics* 14:322-337,
DOI 10.1198/jabes.2009.08038 (`conbio_methods_results_raw_package.md:235`);
Rowcliffe, J. M., et al. 2014. *Methods in Ecology and Evolution* 5:1170-1179,
DOI 10.1111/2041-210X.12278 (`conbio_methods_results_raw_package.md:239`).

## C. Class-specific hierarchical GAMs

### Data construction and fitting

| Requested item | Final implementation | Source / location | Status |
|---|---|---|---|
| Camera-day-hour | Corrected raw sighting rows were assigned to camera, survey round, operational date, and apparent-solar hour. Within each class, `distinct(camera_id, round, event_date, hour)` defines positive binary camera-day-hours. Pooled uses any gazelle sighting. | `definitive_hierarchical_diel_GAM.R`, `binary_keys()` and `make_successes()`, lines 279-285 and 342-405 | VERIFIED FINAL |
| Aggregation / zeros | A complete crossing of each positive-effort camera-round and 24 hours is joined to successes; missing successes become 0. `success` is the number of dates positive in that hour; `failure = Days.Operational - success`; `trials = Days.Operational`; zeros retained. | same, lines 408-431 | VERIFIED FINAL |
| Sample sizes | Each final class/pooled model: 8,256 rows = 344 camera-rounds x 24 hours; 143,976 binomial trials. Positive camera-day-hours: juvenile 183, female 1,576, male 2,307, pooled 3,857. | `definitive_hierarchical_diel_data_audit.csv`; Table S10 | VERIFIED FINAL |
| Fitting call | `mgcv::bam(formula, family=binomial(link="logit"), method="fREML", discrete=TRUE, select=TRUE, knots=list(Time=c(0,24)), nthreads=2, control=gam.control(trace=FALSE))`. | script `fit_model()`, lines 507-530 | VERIFIED FINAL |
| Baseline cyclic smooth | `s(Time, bs="cc", k=12)` with cyclic knots 0 and 24. | `ecological_terms()`, lines 482-497 | VERIFIED FINAL |
| Survey round | `round_factor` categorical fixed effect in every candidate. | `build_formula()`, lines 500-504 | VERIFIED FINAL |
| Camera hierarchy | Juvenile: `s(camera_round,bs="re") + s(camera_id,bs="re")`. Adults/pooled additionally: `s(Time,bs="cc",k=12,by=camera_id,m=1,id=20)`, with shared smoothing parameter, plus both random intercepts. | `hierarchy_terms()`, lines 464-479; final implementation CSV | VERIFIED FINAL |
| Weights/offset | No separate weights or offset. The grouped `cbind(success,failure)` response carries operational-day trials. | final call/formulas | VERIFIED FINAL |
| AR structure | No `rho`, `AR.start`, or other AR1 term was fitted. Lag-1 Pearson residual correlation was diagnostic only. | final call; `extract_diagnostics()` | VERIFIED FINAL |

### Exact selected formulas

All formulas below are `VERIFIED FINAL` from
`definitive_hierarchical_diel_audit_model_implementation.csv`.

| Response | Formula beyond `cbind(success,failure) ~ round_factor +` |
|---|---|
| Juvenile | `wolf_z * jackal_z + s(Time,cc,k=12) + s(Time,by=wolf_z,cc,k=12) + s(Time,by=jackal_z,cc,k=12) + s(Time,by=wolf_jackal,cc,k=12) + s(camera_round,re) + s(camera_id,re)` |
| Adult female | Juvenile ecological terms plus `s(Time,by=camera_id,cc,k=12,m=1,id=20) + s(camera_round,re) + s(camera_id,re)` |
| Adult male | Same as adult female. |
| Pooled | `wolf_z + s(Time,cc,k=12) + s(Time,by=wolf_z,cc,k=12) + s(Time,by=camera_id,cc,k=12,m=1,id=20) + s(camera_round,re) + s(camera_id,re)` |

### Ecological candidate formulas

`VERIFIED FINAL` from `definitive_hierarchical_diel_GAM.R:482-504`:

| ID | Parametric and varying-coefficient terms, in addition to baseline `s(Time,cc,k)` |
|---|---|
| A baseline | no predator terms |
| B wolf-only | `wolf_z + s(Time,by=wolf_z,cc,k)` |
| C jackal-only | `jackal_z + s(Time,by=jackal_z,cc,k)` |
| D additive | `wolf_z + jackal_z` plus both corresponding time smooths |
| E interaction | `wolf_z * jackal_z` plus wolf, jackal, and `wolf_z*jackal_z` time smooths |

### Candidate AICs and selections

All values are `VERIFIED FINAL` from
`definitive_hierarchical_diel_ecological_model_comparison.csv` and Table S7.

| Response | A baseline | B wolf | C jackal | D additive | E interaction | Selected under least-complex-within-2 rule |
|---|---:|---:|---:|---:|---:|---|
| Juvenile | 1336.322 (14.966) | 1335.187 (13.831) | 1331.744 (10.388) | 1331.021 (9.666) | 1321.356 (0) | E |
| Adult female | 5995.249 (10.439) | 5990.347 (5.537) | 5991.937 (7.127) | 5987.305 (2.495) | 5984.810 (0) | E |
| Adult male | 8589.608 (6.000) | 8587.376 (3.768) | 8588.870 (5.263) | 8586.869 (3.262) | 8583.608 (0) | E |
| Pooled | 11276.834 (2.152) | 11274.682 (0) | 11276.994 (2.313) | 11274.878 (0.196) | 11276.252 (1.570) | B |

Numbers in parentheses are delta AIC. The rule yields the manuscript models.
For pooled gazelles, wolf-only is also the minimum-AIC model, not merely the
simplest model within two units.

### Citations

| Citation | Verified project source | Status |
|---|---|---|
| Cunningham et al. 2019, *Ecography* 42:2157-2168, DOI 10.1111/ecog.04485 | `manuscript_definitive_hierarchical_diel_methods.md:39` | VERIFIED FINAL |
| Iannarilli et al. 2025, *Journal of Animal Ecology* 94:182-194, DOI 10.1111/1365-2656.14213 | same, line 41 | VERIFIED FINAL |
| Burnham & Anderson 2004 | An in-text occurrence was recovered from the DOCX audit, but the exact bibliography entry is not present in accessible repository text and the current DOCX was locked. | LIKELY FINAL BUT NEEDS CONFIRMATION |

## D. Seven sensitivity specifications

Final implementation: `definitive_hierarchical_diel_GAM.R:811-1003` plus the
explicitly loaded cached combined sensitivity at lines 950-968. Status for all
rows: `VERIFIED FINAL`, with the provenance caveat shown for specification 7.

| # | Specification | Difference from primary |
|---:|---|---|
| 1 | Primary apparent solar | Full interaction model E, k=12, raw-sighting camera-day-hour histories, selected hierarchy. |
| 2 | k=8 | Primary model with predator-time and baseline cyclic basis dimension 8. |
| 3 | k=16 | Primary model with basis dimension 16. |
| 4 | Sun time | `overlap::sunTime()` response scale; otherwise k=12, raw histories, selected hierarchy. |
| 5 | Alternate hierarchy | Apparent solar and raw histories; if selected hierarchy is camera-only use camera-round, otherwise use camera-only. For the selected adult/pooled site-shape fits, this removes site-shape and camera-round terms and retains camera. |
| 6 | Independent 30-min histories | Apparent solar, k=12, selected hierarchy; binary histories reconstructed from the frozen independent-event file (pooled independently reconstructed to 3,688 encounters). |
| 7 | Sun time + camera-only hierarchy | Sun time and alternate camera-only hierarchy together. Class-specific values are imported from cached `definitive_hierarchical_diel_suntime_primary_model_summary.csv`, not refitted by the final apparent-solar script. |

### Support counts

All counts are `VERIFIED FINAL` from
`definitive_hierarchical_diel_audit_robustness.csv`.

| Response | Wolf x time | Jackal x time | Wolf x jackal x time |
|---|---:|---:|---:|
| Juvenile | 0/7 | 0/7 | 7/7 |
| Adult female | 7/7 | 2/7 | 6/7 |
| Adult male | 4/7 | 1/7 | 6/7 |
| Pooled full model | 7/7 | 2/7 | 7/7 |

For pooled gazelles, jackal and interaction are not retained in the selected
wolf-only model, irrespective of their full-model sensitivity counts.

## E. Prediction normalization

| Requested item | Exact implementation | Source / location | Status |
|---|---|---|---|
| Grid | Time centers 0.05-23.95 by 0.1 h; all four wolf/jackal combinations at -1/+1 SD; `round_factor` fixed at round 4; placeholder first camera levels supplied. | `definitive_hierarchical_diel_GAM.R:1021-1043` | VERIFIED FINAL |
| Random effects | The prediction lpmatrix excludes every smooth whose label contains `camera`; population predator/baseline smooths and the round-4 fixed effect remain. | same, lines 1045-1049 | VERIFIED FINAL |
| Simulation | 1,000 multivariate-normal coefficient draws from `vcov(fit, unconditional=TRUE)` with ordinary `vcov` fallback; `mvtnorm::rmvnorm(..., method="svd")`. | config lines 17-35; lines 1050-1053 | VERIFIED FINAL |
| Normalization | For each predator context and each draw, logistic probabilities are divided by rectangle-rule area `sum(p)*0.1`. Thus every curve integrates numerically to 1 over 24 h. | same, lines 1054-1063 | VERIFIED FINAL |
| Bands | Pointwise 2.5% and 97.5% quantiles of the normalized simulated curves. | same, lines 1065-1074 | VERIFIED FINAL |
| Peak time | Grid time at `which.max(normalized_activity)`. | same, lines 1087-1095 | VERIFIED FINAL |

## F. Confirmatory joint GAM

| Requested item | Exact implementation | Source / location | Status |
|---|---|---|---|
| Formula | Full saved formula is in `joint_demographic_penalty_calibrated_outputs/penalty_calibrated_formula.txt`; it includes class x round, class x wolf x jackal parametric terms, class baseline curves, common predator-time curves, demographic deviations, adult site shapes, and four random-effect smooths. | formula file; `calibrate_joint_demographic_diel_GAM.R:134-155` | VERIFIED FINAL |
| Sum-to-zero deviations | For each of `wolf_z`, `jackal_z`, and `wolf_jackal`, common `s(Time,by=variable,bs="cc",k=12)` plus `s(demographic_class,Time,by=variable,bs="sz",k=12,xt=list(bs="cc"))`. | script `predator_terms()`, lines 122-132 | VERIFIED FINAL |
| Adult site shapes | Ordered reference factors for female and male cameras; separate `cc` deviations with k=12, m=1, linked penalties `id=21` and `id=22`. | same, lines 82-119 | VERIFIED FINAL |
| Random effects | `s(camera_id,re) + s(camera_round,re) + s(camera_class,re) + s(camera_round_class,re)`. | same, lines 134-155 | VERIFIED FINAL |
| Fit settings | `bam`, binomial logit, fREML, `discrete=TRUE`, `select=TRUE`, cyclic knots 0/24, two threads, no AR1. | same, lines 158-184 | VERIFIED FINAL |
| Omnibus tests | Identify each `bs="sz"` demographic-deviation smooth and call `mgcv:::reTest()`; report summed EDF, weighted-reference rank, statistic, and boundary-test P value. Reference rank is not conventional chi-square df. | same, lines 187-227 | VERIFIED FINAL |
| Pairwise tests | Only for an omnibus-supported family: refit each two-class model with the same common + symmetric deviation construction, apply the same `reTest()`, then Holm-adjust three raw P values within the family. | same, lines 357-401 | VERIFIED FINAL |
| Pairwise time scale | Primary apparent-solar, raw-history model only; pairwise models were not repeated across all seven sensitivities. | same, `prepare_data()` defaults lines 54-58 and pairwise block | VERIFIED FINAL |
| Omnibus robustness | Refit interaction demographic deviation under primary, k8, k16, sun time, camera-only hierarchy, independent histories, and sun-time + camera-only; all seven P<0.01. | same, lines 905-1008; sensitivity CSVs | VERIFIED FINAL |

### Joint-model test values

All rows are `VERIFIED FINAL` from the penalty-calibrated output directory.

| Scope / term | EDF | Reference rank | Statistic | Raw P | Holm P / sensitivity |
|---|---:|---:|---:|---:|---:|
| Omnibus wolf x time | 1.006554 | 30 | 2.080503 | 0.077787 | primary unsupported; not reopened |
| Omnibus jackal x time | 0.000414 | 30 | 0.000160 | 0.706565 | primary unsupported; not reopened |
| Omnibus wolf x jackal x time | 1.575468 | 30 | 6.109909 | 0.005964 | 7/7 |
| Juvenile - adult female, interaction | 1.503706 | 20 | 5.127298 | 0.006037 | 0.012074 |
| Juvenile - adult male, interaction | 1.629777 | 20 | 6.724420 | 0.003704 | 0.011112 |
| Adult female - adult male, interaction | 0.000329 | 20 | 0.000136 | 0.276733 | 0.276733 |

## G. GAM diagnostics actually used

### Class-specific selected models

The final script computed convergence labels, deviance explained, an adjusted
R-squared diagnostic, Pearson dispersion, observed/predicted encounter
frequency, mean within-camera-round lag-1 Pearson correlation, `k.check()` with
300 simulations, and `concurvity()`. For site-shape models, k-index and
concurvity were intentionally returned as NA because the large factor-by-camera
smooths made those summaries impractical. It did not fit an AR model.

All values below are `VERIFIED FINAL` from
`definitive_hierarchical_diel_diagnostics.csv`, rows `Final fREML`.

| Model | Convergence | Deviance explained | Dispersion | Mean lag-1 Pearson | Min k-index (P) | Max concurvity | Observed / predicted frequency |
|---|---|---:|---:|---:|---:|---:|---:|
| Juvenile | no warning | 0.465714 | 0.413479 | 0.729110 | 0.852007 (0.020) | 0.183389 | 0.001271045 / 0.001271045 |
| Adult female | no warning | 0.573552 | 0.634084 | 0.312642 | NA | NA | 0.010946269 / 0.010946269 |
| Adult male | no warning | 0.466244 | 0.784112 | 0.172726 | NA | NA | 0.016023504 / 0.016023504 |
| Pooled | no warning | 0.547133 | 0.874321 | 0.158126 | NA | NA | 0.026789187 / 0.026789187 |

### Confirmatory joint model

| Diagnostic | Value | Source | Status |
|---|---:|---|---|
| Maximum absolute gradient | 6.902e-5 | `penalty_calibrated_diagnostics.csv` | VERIFIED FINAL |
| Minimum Hessian eigenvalue | -1.406e-13; treated as numerical zero | same | VERIFIED FINAL |
| Hessian numerically PSD | TRUE | same | VERIFIED FINAL |
| Coefficients / rank | 3008 / 3008 | same | VERIFIED FINAL |
| Deviance explained | 0.563573 | same | VERIFIED FINAL |
| Pearson dispersion | 0.682957 | same | VERIFIED FINAL |
| Pairwise fits | max gradients <7.61e-5; numerical PSD Hessians; full coefficient rank; no warnings | `penalty_calibrated_pairwise_diagnostics.csv` | VERIFIED FINAL |
| Diagnostics not run | No final joint-GAM `k.check`, concurvity, residual ACF, or AR1 fit is recorded. | calibration script/output set | NOT FOUND |

# Details that genuinely need to appear in supplementary methods

- Predator indices: final `pcount` distributions, constant state formula,
  selected camera-design detection formula, max-normalization before downstream
  z scaling, 60 x 7 sampling dimensions, and the Pearson correlation.
- Joint LF model: two-factor structure, priors/MCMC/tuning, loading-to-correlation
  calculation, PPC grouping/statistics, non-use of WAIC, and the intercept
  convergence limitation.
- Class-specific NMix: exact priors/MCMC/tuning, intercept diagnostic limitation,
  intercept-free relative normalization, and draw-level contrast formulas.
- Ratio model: literal midpoint coding, exact six-model candidate set,
  beta-binomial dispersion, random-effect variances, and LRT definitions.
- Diel methods: explicit apparent-solar formula, astronomical elevation
  calculation, camera-cluster bootstrap algorithms, exact GAM hierarchy and
  fitting settings, cached provenance of sensitivity 7, population-level
  exclusion of camera smooths, and numerical normalization rule.
- Joint GAM: exact `sz` construction, penalty-aware `reTest()` interpretation,
  and the fact that pairwise contrasts were primary-model only.

# Details already stated in the main manuscript

- The ecological purpose and cautious interpretation of N-mixture responses as
  relative local abundance/site use, the final class-specific predator formula,
  and the final joint LF abundance formula.
- The beta-binomial response concept, core adjustment covariates, and selected
  additive Early/Late model.
- Apparent solar time as the primary diel scale, binomial camera-day-hour
  histories, candidate ecological structures, least-complex-within-two-AIC
  rule, selected class hierarchies, seven-sensitivity concept, normalized curve
  interpretation, and headline joint-GAM results.
- These items should be referenced rather than repeated verbatim in a concise
  supplement; only implementation details needed to reproduce S1-S10 should be
  added.

# Unresolved items requiring Shlomo's decision

1. Locate or reconstruct a clean, executable final predator-index script. The
   archived script's wolf `fitList()` call is missing a comma and its jackal
   empirical-Bayes call uses the wrong object/member name. Saved AIC and index
   files verify the intended results, but the executable provenance is incomplete.
2. Confirm the original `unmarked` version and any implicit `pcount()` defaults
   (especially K) used for the predator-index fits; they are not recorded.
3. Align the manuscript's pooled selected-model EDF to 3.870 if it still reports
   3.04, and align the Early/Late verbal boundary to the implemented midpoint rule.
4. Decide whether to disclose the cached provenance of the seventh class-specific
   GAM sensitivity explicitly; it was imported from the legacy sun-time-primary
   output rather than refitted in the final apparent-solar script.
5. Confirm the exact Burnham & Anderson (2004) bibliography entry after the DOCX
   is closed or an unlocked copy is available.
