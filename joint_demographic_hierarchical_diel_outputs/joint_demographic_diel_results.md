# Confirmatory joint demographic hierarchical diel GAM

## Scope

This analysis directly tests equality of centered predator-associated diel response functions among juveniles, adult females, and adult males. It is an addition to, not a replacement for, the definitive class-specific GAMs.

The model was prespecified and no demographic predator-by-time term was removed by model selection.

## Joint data audit

| demographic_class | rows | successes_detection_days | cameras_represented | survey_rounds | zero_detection_rows | minimum_denominator | maximum_denominator | time_scale | history_type |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Juvenile | 8256.000 | 183.000 | 60.000 | 7.000 | 8097.000 | 5.000 | 33.000 | Equal-duration apparent solar time | Raw sighting camera-day-hour presence |
| Adult female | 8256.000 | 1576.000 | 60.000 | 7.000 | 7281.000 | 5.000 | 33.000 | Equal-duration apparent solar time | Raw sighting camera-day-hour presence |
| Adult male | 8256.000 | 2307.000 | 60.000 | 7.000 | 6723.000 | 5.000 | 33.000 | Equal-duration apparent solar time | Raw sighting camera-day-hour presence |

Wolf and jackal standardized values were read without recalculation from `definitive_hierarchical_diel_data.csv`, the exact definitive model input.

## Model and dependence structure

- Formula: `cbind(success, failure) ~ demographic_class * round_factor +      demographic_class * wolf_z * jackal_z + s(Time, by = demographic_class,      bs = "cc", k = 12) + s(Time, by = wolf_juvenile, bs = "cc",      k = 12) + s(Time, by = jackal_juvenile, bs = "cc", k = 12) +      s(Time, by = interaction_juvenile, bs = "cc", k = 12) + s(Time,      by = wolf_female, bs = "cc", k = 12) + s(Time, by = jackal_female,      bs = "cc", k = 12) + s(Time, by = interaction_female, bs = "cc",      k = 12) + s(Time, by = wolf_male, bs = "cc", k = 12) + s(Time,      by = jackal_male, bs = "cc", k = 12) + s(Time, by = interaction_male,      bs = "cc", k = 12) + s(Time, by = female_camera, bs = "cc",      k = 12, m = 1, id = 21) + s(Time, by = male_camera, bs = "cc",      k = 12, m = 1, id = 22) + s(camera_id, bs = "re") + s(camera_round,      bs = "re") + s(camera_class, bs = "re") + s(camera_round_class,      bs = "re")`
- Fitting call: `bam(formula = formula, family = binomial(link = "logit"), data = data, method = "fREML", control = gam.control(trace = FALSE), select = TRUE, knots = list(Time = c(0, 24)), discrete = TRUE, nthreads = cfg$nthreads)`
- Coefficient simulations: 5000; seed: 20260825.
- Covariance: Unconditional Bayesian covariance including smoothing-parameter uncertainty.
- Population predictions fixed survey round at round 4 and excluded shared and class-specific camera/camera-round random effects and adult camera-specific cyclic deviations.
- Shared camera and camera-round random intercepts account for persistent cross-class sampling dependence; class-specific hierarchical terms preserve the separate-model structure. Same-day/hour residual dependence is not fully identified by grouped-binomial aggregation.
- The requested unaggregated shared sampling-unit sensitivity was not validly constructible for camera 45 round 2 and was computationally infeasible for the 143,472 known camera-day-hour sampling units.

## Diagnostics

| model | convergence | maximum_absolute_gradient | minimum_hessian_eigenvalue | hessian_positive_definite | hessian_numerically_positive_semidefinite | pearson_dispersion | deviance_explained | mean_camera_round_class_lag1_pearson | maximum_fitted_probability | minimum_fitted_probability | coefficient_rank | coefficients | covariance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Primary joint confirmatory demographic GAM | Numerically converged | 0.000 | -0.000 | FALSE | TRUE | 0.684 | 0.564 | 0.407 | 0.501 | 0.000 | 2975.000 | 2975.000 | Unconditional Bayesian covariance including smoothing-parameter uncertainty |

The focal basis checks were nonproblematic (P range = 0.825-0.895). The maximum focal design-space overlap was 0.094. Fitting and prediction generated no warnings.
The Hessian minimum eigenvalue (-0.000000000000) was within numerical zero tolerance, but mean adjacent-hour Pearson residual correlation remained 0.407 and Pearson dispersion was 0.684.
Terms effectively shrunk toward zero (EDF <0.01): s(Time):jackal_juvenile; s(Time):jackal_female; s(Time):interaction_female; s(Time):wolf_male; s(Time):jackal_male.
Juvenile interaction stability must be judged from the direct tests, simultaneous bands, EDF, and the seven-specification sensitivity table; the smaller juvenile sample was not concealed by deleting or merging terms.

## Omnibus demographic tests

| predator_effect | omnibus_statistic | df_rank | raw_P |
| --- | --- | --- | --- |
| Wolf x time | 3.795 | 20.000 | 1.000 |
| Jackal x time | 0.000 | 20.000 | 1.000 |
| Wolf x jackal x time | 8.147 | 20.000 | 0.991 |

## Planned pairwise demographic tests

| predator_effect | contrast | raw_P | Holm_adjusted_P | simultaneous_CI_excludes_zero | time_intervals_of_difference |
| --- | --- | --- | --- | --- | --- |
| Wolf x time | Juvenile - Adult female | 0.987 | 1.000 | FALSE | None |
| Wolf x time | Juvenile - Adult male | 1.000 | 1.000 | FALSE | None |
| Wolf x time | Adult female - Adult male | 0.989 | 1.000 | FALSE | None |
| Jackal x time | Juvenile - Adult female | 1.000 | 1.000 | FALSE | None |
| Jackal x time | Juvenile - Adult male | 1.000 | 1.000 | FALSE | None |
| Jackal x time | Adult female - Adult male | 1.000 | 1.000 | FALSE | None |
| Wolf x jackal x time | Juvenile - Adult female | 0.771 | 1.000 | FALSE | None |
| Wolf x jackal x time | Juvenile - Adult male | 0.905 | 1.000 | FALSE | None |
| Wolf x jackal x time | Adult female - Adult male | 0.998 | 1.000 | FALSE | None |

These are rank-reduced Wald tests of centered link-scale difference functions. Holm adjustment was applied across the three class contrasts within each predator family.

## RMS temporal-effect magnitudes

| predator_effect | demographic_class | RMS_effect | lower_95 | upper_95 |
| --- | --- | --- | --- | --- |
| Wolf x time | Juvenile | 0.156 | 0.047 | 0.336 |
| Wolf x time | Adult female | 0.139 | 0.060 | 0.249 |
| Wolf x time | Adult male | 0.003 | 0.001 | 0.008 |
| Jackal x time | Juvenile | 0.001 | 0.000 | 0.003 |
| Jackal x time | Adult female | 0.001 | 0.000 | 0.003 |
| Jackal x time | Adult male | 0.001 | 0.000 | 0.001 |
| Wolf x jackal x time | Juvenile | 0.338 | 0.136 | 0.557 |
| Wolf x jackal x time | Adult female | 0.010 | 0.002 | 0.032 |
| Wolf x jackal x time | Adult male | 0.105 | 0.041 | 0.209 |

## Pairwise RMS differences

| predator_effect | contrast | difference | lower_95 | upper_95 | ordering_supported |
| --- | --- | --- | --- | --- | --- |
| Wolf x time | Juvenile - Adult female | 0.015 | -0.134 | 0.213 | Uncertain |
| Wolf x time | Juvenile - Adult male | 0.153 | 0.044 | 0.333 | Yes: first class has larger RMS effect |
| Wolf x time | Adult female - Adult male | 0.136 | 0.057 | 0.245 | Yes: first class has larger RMS effect |
| Jackal x time | Juvenile - Adult female | -0.000 | -0.002 | 0.002 | Uncertain |
| Jackal x time | Juvenile - Adult male | 0.001 | -0.000 | 0.003 | Uncertain |
| Jackal x time | Adult female - Adult male | 0.001 | -0.000 | 0.003 | Uncertain |
| Wolf x jackal x time | Juvenile - Adult female | 0.325 | 0.124 | 0.546 | Yes: first class has larger RMS effect |
| Wolf x jackal x time | Juvenile - Adult male | 0.229 | 0.009 | 0.467 | Yes: first class has larger RMS effect |
| Wolf x jackal x time | Adult female - Adult male | -0.094 | -0.198 | -0.026 | No: second class has larger RMS effect |

RMS values are effect-size summaries. An ordering is supported only when the paired simulation interval excludes zero.

## Astronomical daylight summaries

| demographic_class | predator_effect | daylight_shift | interval | result |
| --- | --- | --- | --- | --- |
| Juvenile | Wolf gradient | 0.002 | -0.066 to 0.077 | uncertain |
| Juvenile | Jackal gradient | 0.000 | -0.001 to 0.001 | uncertain |
| Juvenile | Wolf x jackal interaction | 0.098 | -0.116 to 0.305 | uncertain |
| Adult female | Wolf gradient | -0.016 | -0.120 to 0.095 | uncertain |
| Adult female | Jackal gradient | 0.001 | -0.001 to 0.003 | uncertain |
| Adult female | Wolf x jackal interaction | -0.001 | -0.028 to 0.026 | uncertain |
| Adult male | Wolf gradient | -0.001 | -0.006 to 0.005 | uncertain |
| Adult male | Jackal gradient | 0.001 | -0.000 to 0.002 | uncertain |
| Adult male | Wolf x jackal interaction | 0.111 | -0.083 to 0.302 | uncertain |

Daylight was solar elevation >0 degrees, averaged across 5,978 operational camera-days with exact dates and camera-specific latitude. Camera 45 round 2 was excluded only from this directional weighting because its 21 exact outage dates are unknown; it remained in the primary GAM.

## Sensitivity of direct demographic contrasts

| predator_effect | contrast | primary_result | number_supported | number_tested | failed_specifications | robust_across_all_specifications |
| --- | --- | --- | --- | --- | --- | --- |
| Jackal x time | Adult female - Adult male | Holm P = 1.000 | 0.000 | 7.000 | Primary apparent solar; k = 8; k = 16; Sun time; Alternate camera-only hierarchy; Independent 30-min histories; Sun time + camera-only hierarchy | FALSE |
| Jackal x time | Juvenile - Adult female | Holm P = 1.000 | 2.000 | 7.000 | Primary apparent solar; k = 8; k = 16; Sun time; Independent 30-min histories | FALSE |
| Jackal x time | Juvenile - Adult male | Holm P = 1.000 | 0.000 | 7.000 | Primary apparent solar; k = 8; k = 16; Sun time; Alternate camera-only hierarchy; Independent 30-min histories; Sun time + camera-only hierarchy | FALSE |
| Wolf x jackal x time | Adult female - Adult male | Holm P = 1.000 | 1.000 | 7.000 | Primary apparent solar; k = 8; k = 16; Sun time; Alternate camera-only hierarchy; Independent 30-min histories | FALSE |
| Wolf x jackal x time | Juvenile - Adult female | Holm P = 1.000 | 0.000 | 7.000 | Primary apparent solar; k = 8; k = 16; Sun time; Alternate camera-only hierarchy; Independent 30-min histories; Sun time + camera-only hierarchy | FALSE |
| Wolf x jackal x time | Juvenile - Adult male | Holm P = 1.000 | 0.000 | 7.000 | Primary apparent solar; k = 8; k = 16; Sun time; Alternate camera-only hierarchy; Independent 30-min histories; Sun time + camera-only hierarchy | FALSE |
| Wolf x time | Adult female - Adult male | Holm P = 1.000 | 2.000 | 7.000 | Primary apparent solar; k = 8; k = 16; Sun time; Independent 30-min histories | FALSE |
| Wolf x time | Juvenile - Adult female | Holm P = 1.000 | 0.000 | 7.000 | Primary apparent solar; k = 8; k = 16; Sun time; Alternate camera-only hierarchy; Independent 30-min histories; Sun time + camera-only hierarchy | FALSE |
| Wolf x time | Juvenile - Adult male | Holm P = 1.000 | 0.000 | 7.000 | Primary apparent solar; k = 8; k = 16; Sun time; Alternate camera-only hierarchy; Independent 30-min histories; Sun time + camera-only hierarchy | FALSE |

## Answers to the confirmatory questions

A-B. No predator component showed formal demographic heterogeneity in the primary model. Wolf x time: rank-reduced Wald chi-square = 3.795, rank = 20, P = 1.000; Jackal x time: rank-reduced Wald chi-square = 0.000, rank = 20, P = 1.000; Wolf x jackal x time: rank-reduced Wald chi-square = 8.147, rank = 20, P = 0.991.
C-D. No planned pairwise class contrast was supported after Holm correction (any supported = FALSE), and no simultaneous difference band excluded zero at any time of day.
E. RMS effect sizes gave only partial ordering evidence. For wolf effects, juveniles and females each exceeded males, but juveniles did not exceed females. For the wolf-jackal interaction, juveniles exceeded both adult classes, but males exceeded females. No jackal ordering was supported. Thus the full Juvenile > Adult female > Adult male ordering was not supported for any component.
F. No direct demographic contrast was robust across all seven specifications (any robust = FALSE). Isolated support occurred only in camera-only sensitivities and is sensitivity-dependent.
G. The joint model materially qualifies, rather than overturns, the separate class-specific GAMs. Within-class predator-time terms can be supported while the direct between-class response functions remain statistically indistinguishable.
H. Reporting concerns are the sparse juvenile response (183 detection-days), residual adjacent-hour correlation, underdispersion, and residual same-day/hour cross-class dependence after aggregation. Exact shared-unit modeling was not feasible, so absence of evidence must not be described as evidence of identical responses.
Cross-class dependence. Shared camera and camera-round effects were included in the primary model. Because exact operational dates are unavailable for one camera-round and a 143,472-level sampling-unit effect is infeasible here, the incremental effect of explicit same-day/hour dependence could not be tested.
Directional daylight H3. No class showed a supported increase or decrease in astronomical-daylight use (any supported shift = FALSE), and no direct demographic daylight-shift difference was supported (FALSE).
Overall H3 assessment. Demographic heterogeneity in response shape: not supported by the primary direct tests. Predicted Juvenile > Female > Male ordering: only partially supported and never in full. Redistribution toward daylight: unsupported.

## Output index

- Tables A-F are saved as separately named CSV files.
- Figures 1-5 are saved as 400-dpi PNG and vector PDF files.
- The exact formula, model call, full parametric table, full smooth inventory, warnings, diagnostics, and session information are saved in this directory.
