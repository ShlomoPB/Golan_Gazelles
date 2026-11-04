# Mountain Gazelle Manuscript: Final Verification Report

Verified against the frozen publication outputs and corrected source data on 19 August 2026.

## 1. Survey and sample counts

- The final demographic analyses used **3,925 independent 30-minute class-events**: **2,268 adult-male**, **1,477 adult-female**, and **180 fawn** events.
- These events contained summed maximum counts of **2,355 adult males**, **1,951 adult females**, and **199 fawns**. These are repeated classified-animal detections, not unique individuals.
- The fawn:female analysis used **1,951**, rather than 1,477, because it modeled summed maximum animals per camera-round; 1,477 is the number of independent female event rows.
- The raw-export audit identified **48 subadult/not-adult-male events** containing 62 recorded animals. These were excluded from the adult-male count history, fawn:female analysis, and diel analysis.
- It identified **391 unidentifiable/nonfocal 30-minute event groups** containing 463 animals; these were excluded from all retained class-specific analyses.
- The raw export contains **3,688 species-level gazelle encounters** before demographic splitting. This is not directly additive with class-event counts because one encounter can contain multiple demographic classes.
- Survey effort comprised **60 cameras**, **7 nominal rounds**, **344 positive-effort camera-rounds**, and **5,999 camera-nights**. Positive-effort cameras by round were 51, 56, 60, 57, 56, 55, and 9.

## 2. Camera survey details

First positive-effort assignments were **46 Browning DARK OPS HD PRO X** cameras (20 normal-sensitivity and 26 long-sensitivity settings), **12 Bushnell TROPHY CAM HD AGGRESSOR NO-GLOW**, and **2 Browning DARK OPS HD MAX** cameras. Replacements and configuration changes produced the following recorded exposures.

| Recorded camera model | Camera-rounds | Sites ever using model | Camera-nights |
|---|---:|---:|---:|
| Browning DARK OPS HD PRO X - Normal Sensitivity | 195 | 40 | 3,503 |
| Browning DARK OPS HD PRO X - Long Sensitivity | 68 | 28 | 1,130 |
| Browning DARK OPS HD MAX | 58 | 16 | 1,049 |
| Bushnell TROPHY CAM HD AGGRESSOR NO-GLOW | 23 | 12 | 317 |

Recorded operation spans **20 August-15 December 2020**. Cameras were deployed asynchronously and did not all operate over identical dates. Any Methods statement specifying **30 August-14 December, 5,997 nights, six occasions, all cameras concurrent, only BTC-6HDPX cameras, or a 20-minute independence rule** is inconsistent with the finalized gazelle data. The retained manuscript analyses use seven survey rounds and a 30-minute event rule.

## 3. Land-use and human covariates

Open-space, infrastructure, cropland, and orchard rasters came from Hamaarag, Israel's National Ecosystem Assessment Program; reserve boundaries came from the Israel Nature and Parks Authority; and firing-zone and minefield layers came from the Israel Defense Forces. Open space represents undeveloped natural cover, including grassland, scrubland, woodland, and planted forest. `Built.and.Roads` represents settlements, principal roads, industrial areas, quarries, solar facilities, and construction. `Nature.Reserves` represents protected reserve cover.

The retained beta-binomial model used `Open.Spaces`, `Built.and.Roads`, and `Nature.Reserves`; the focused `lfMsNMix` model used only `Open.Spaces`; and the primary separate-class N-mixture models contained no landscape covariate. Cropland, firing zones, minefields, orchards, and proximity to water remain in the site-covariate file but were not retained in the primary analyses.

Binary 25-m rasters were weighted inside a 7-km circular kernel using `w = exp[-0.01 x sqrt(dx^2 + dy^2)]`, where distance is measured in 25-m cells. Equivalently, `w = exp(-distance_m/2500)`, truncated at 7 km. Weighted raster values were summed, extracted at camera locations, rescaled to 0-1, and subsequently z-standardized for modeling. Sources and implementation are documented by Preiss-Bloom et al. (2025) and in the published source code:

- Preiss-Bloom, S., H. Shamon, Z. Amir, D. Ben-Ami, et al. 2025. Human disturbance thresholds determine the ecological role of an apex predator. *Scientific Reports* 15:37604. https://doi.org/10.1038/s41598-025-21402-x
- Source code: https://github.com/zacha46/golan_wolves_SEM

The correct citation for the decay-method reference is:

- Shamoon, H., T. Dayan, and D. Saltz. 2017. Cattle grazing effects on mountain gazelles in Mediterranean natural landscapes. *Journal of Wildlife Management* 81:1351-1362. https://doi.org/10.1002/jwmg.21323

The human variable is exactly `human.abundance` (`human_z` after standardization). Camera-round human counts were standardized to 15 operational days as `round(count / days_operational x 15)`, modeled with a Poisson `unmarked::pcount` model using state `~ 1` and detection `~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera`. Empirical-Bayes best unbiased predictors from `ranef()` and `bup()` were divided by their maximum to produce the 0-1 index.

## 4. Wolf and jackal indices

| Index | Model family | State formula | Detection formula |
|---|---|---|---|
| Wolf | Zero-inflated Poisson | `~ 1` | `~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera` |
| Jackal | Poisson | `~ 1` | `~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera` |

The retained `SiteCovs` values exactly match the empirical-Bayes BUP columns in `Null Wolf abundance with SE.csv` and `Null Jackal abundance with SE.csv`, divided by each species' maximum BUP. They are **relative site-level activity indices**, not density or population-abundance estimates. The wolf and jackal count histories were not normalized to 15 operational days.

Wolf and jackal indices were weakly negatively correlated: **Pearson r = -0.183, P = 0.162, n = 60**. Median thresholds for the diel contexts were **wolf = 0.556** and **jackal = 0.495**.

## 5. Final primary spAbundance results

| Class | Predictor | Beta [95% CrI] | Rhat | ESS |
|---|---|---:|---:|---:|
| Fawn | Wolf | 0.265 [-0.274, 0.855] | 1.002 | 5,455 |
| Fawn | Jackal | -0.803 [-1.533, -0.126] | 1.004 | 3,847 |
| Fawn | Wolf x jackal | 0.497 [-0.208, 1.304] | 1.004 | 3,549 |
| Adult female | Wolf | 0.540 [0.132, 0.972] | 1.004 | 9,093 |
| Adult female | Jackal | -0.721 [-1.186, -0.273] | 1.001 | 8,125 |
| Adult female | Wolf x jackal | 0.521 [0.039, 1.020] | 1.002 | 6,835 |
| Adult male | Wolf | 0.352 [0.007, 0.712] | 1.001 | 9,538 |
| Adult male | Jackal | -0.607 [-0.959, -0.260] | 1.003 | 8,106 |
| Adult male | Wolf x jackal | 0.395 [0.025, 0.769] | 1.001 | 7,188 |

These values match `manuscript_spabundance_methods_results.md`. Abundance and detection intercepts did not converge, so inference is restricted to relative predator slopes rather than absolute abundance.

## 6. Residual class correlations

| Class pair | Posterior correlation | 95% CrI |
|---|---:|---:|
| Fawn-adult male | 0.657 | 0.435-0.806 |
| Fawn-adult female | 0.863 | 0.523-0.998 |
| Adult male-adult female | 0.907 | 0.712-1.000 |

Manuscript rounding should be **0.66 [0.44-0.81]**, **0.86 [0.52-1.00]**, and **0.91 [0.71-1.00]**. These corrected-data focused-model values supersede **0.60/0.86/0.88**, which came from the older broad shared-environment model. The focused model controls for open-space cover, wolf activity, jackal activity, and wolf x jackal activity.

## 7. Fawn:female analysis

The pooled descriptive ratio is **199/1,951 = 0.102**, or **10.20 fawns per 100 adult females**. The primary model used **221 informative camera-rounds from 50 cameras** in 12 sampling polygons.

| Predictor | Odds ratio [95% CI] | P |
|---|---:|---:|
| Late vs Early | 1.57 [1.01, 2.44] | 0.045 |
| Open-space cover | 1.55 [1.04, 2.32] | 0.031 |
| Wolf activity | 0.93 [0.70, 1.23] | 0.608 |
| Jackal activity | 0.94 [0.63, 1.41] | 0.773 |

Adding wolf x jackal gave **OR = 0.86 [0.53, 1.39], P = 0.538**, did not improve fit (**likelihood-ratio P = 0.552**), and increased AIC by **1.65**. The season-adjusted additive model remains primary.

## 8. Diel analysis

Daylight percentages based on independent event rows were **fawns 82.8% (149/180)**, **adult females 61.1% (903/1,477)**, and **adult males 56.8% (1,288/2,268)**. Daylight was defined as solar elevation greater than or equal to 0 degrees.

Overall activity overlap was **fawn-female 0.787 [0.695-0.855]**, **fawn-male 0.742 [0.650-0.782]**, and **female-male 0.948 [0.859-0.960]**.

| Wolf-jackal context | Cameras | Fawn events | Female events | Male events |
|---|---:|---:|---:|---:|
| Low both | 14 | 67 | 531 | 791 |
| High jackal / Low wolf | 16 | 3 | 69 | 163 |
| Low jackal / High wolf | 16 | 72 | 601 | 770 |
| High both | 14 | 38 | 276 | 544 |

Supported context comparisons after within-class false-discovery-rate correction were **0 of 3 eligible comparisons for fawns**, **4 of 6 for females**, and **4 of 6 for males**.

- The lowest female overlap was low-jackal/high-wolf versus high-both: **Dhat4 = 0.788, 95% camera-bootstrap CI = 0.652-0.883, FDR P = 0.002**.
- The lowest supported male overlap was high-jackal/low-wolf versus high-both: **Dhat4 = 0.838, 95% camera-bootstrap CI = 0.652-0.894, FDR P = 0.012**.
- All eligible fawn comparisons were unsupported. The lowest eligible fawn overlap was low-jackal/high-wolf versus high-both: **Dhat1 = 0.796, 95% camera-bootstrap CI = 0.477-0.854, FDR P = 0.269**. Comparisons involving the three-event high-jackal/low-wolf fawn group were data deficient.

## 9. Sensitivity and supporting-information checks

- Matched negative-binomial `unmarked` models recovered the same qualitative spatial pattern: negative jackal effects for all classes and positive wolf and wolf x jackal effects for adult females and males.
- All 18 expanded human/landscape N-mixture variants failed retained convergence or effective-sample-size criteria. Only failure diagnostics should be included in Supporting Information; their ecological coefficients should not be interpreted.
- Gazelle-predator temporal-overlap analyses are excluded from both the manuscript and Supporting Information.

## 10. Numbers or statements in the current manuscript that must be changed

- Replace any `30 August-14 December`, `5,997 camera-nights`, six-round, all-concurrent, single-camera-model, or 20-minute-event statements with the finalized survey details above.
- Replace stale `3,932/3,906` event totals with **3,925 retained focal class-events**.
- Replace stale daylight values of `901/1,466` females and `1,286/2,260` males with **903/1,477** and **1,288/2,268**, respectively.
- Update the results-outline correlation intervals from `0.43-0.81`, `0.51-1.00`, and `0.73-1.00` to **0.44-0.81**, **0.52-1.00**, and **0.71-1.00**.
- Describe 199, 1,951, and 2,355 as **summed classified-animal counts**, not independent events.
- Describe wolf and jackal covariates as **relative site-level activity indices**, not density estimates.

## Authoritative local sources

- `publication_analysis_manifest.md`
- `publication_data_audit_summary.csv`
- `publication_count_provenance.csv`
- `publication_final_model_specifications.csv`
- `publication_primary_spabundance_coefficients.csv`
- `publication_focused_lfmsnmix_class_correlations.csv`
- `publication_fawn_female_beta_binomial_coefficients.csv`
- `publication_fawn_female_beta_binomial_interaction_sensitivity_model_comparison.csv`
- `publication_gazelle_daylight_proportions.csv`
- `publication_overall_gazelle_diel_activity_class_overlap.csv`
- `publication_authoritative_diel_context_comparisons.csv`
