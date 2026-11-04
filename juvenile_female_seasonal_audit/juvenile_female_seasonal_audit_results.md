# Exploratory audit of the seasonal juvenile:female ratio increase

## Scope and provenance

This audit investigates the published Early/Late association without changing
the primary beta-binomial analysis. It reconstructs the response and period
definition from `analyze_fawn_female_beta_binomial_wolf_jackal.R`, the three
60 x 7 class-count histories, `ObsCovs_050524.csv`, and
`SiteCovs_060324.csv`. Event-level summaries use
`Gazella gazella_independent_events_30min_analysis_ready.csv`; the raw export
was consulted only to determine whether finer age or image metadata existed.

The source code assigns each camera-round from its midpoint date and defines
Early as midpoint date on or before 2020-10-17. Thus, Late technically begins
on October 18 rather than October 17 as stated in the shorthand description;
no camera-round midpoint falls exactly on October 17, so this wording
difference does not change the fitted data. Rounds 3 and 4 contain
camera-rounds in both periods. The primary model contains 221 informative
camera-rounds from 50 cameras and exactly reproduces the saved season effect
(Late vs Early odds ratio = 1.569, 95% CI 1.011-2.436, P = 0.0445).
Exploratory refits reported below preserve the primary fixed and random effects
unless explicitly identified as focused fallback models.

## Counts, effort, and observed composition

Across all positive-effort camera-rounds, the study contributed 344
camera-rounds and 5,999 camera-nights. The 221 model rows contain all 199
juvenile and 1,951 adult-female counts but exclude 123 operational
camera-rounds in which both counts were zero. Encounter rates below therefore
use all documented operational effort; the beta-binomial model instead
conditions on informative rows with at least one juvenile or adult female.

| Period | Cameras | Camera-nights | Valid camera-rounds | Informative/model rows | Juveniles | Adult females | Juveniles/100 nights | Adult females/100 nights | Juveniles/100 females |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Early | 60 | 2,701 | 175 | 111 | 74 | 934 | 2.74 | 34.58 | 7.92 |
| Late | 59 | 3,298 | 169 | 110 | 125 | 1,017 | 3.79 | 30.84 | 12.29 |

From Early to Late, the effort-standardized juvenile encounter rate increased
by 38.3%, the adult-female rate decreased by 10.8%, and the observed
juvenile:female ratio increased by 55.1%. Thus, both a higher numerator rate
and a lower denominator rate contributed, with the juvenile rate showing the
larger proportional change.

| Round | Period coverage | Cameras | Camera-nights | Informative rows | Juveniles | Adult females | Juveniles/100 nights | Adult females/100 nights | Juveniles/100 females |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Early | 51 | 658 | 32 | 14 | 198 | 2.13 | 30.09 | 7.07 |
| 2 | Early | 56 | 1,077 | 35 | 38 | 434 | 3.53 | 40.30 | 8.76 |
| 3 | Early/Late | 60 | 874 | 36 | 18 | 263 | 2.06 | 30.09 | 6.84 |
| 4 | Early/Late | 57 | 1,170 | 41 | 29 | 377 | 2.48 | 32.22 | 7.69 |
| 5 | Late | 56 | 1,127 | 36 | 48 | 402 | 4.26 | 35.67 | 11.94 |
| 6 | Late | 55 | 909 | 36 | 38 | 238 | 4.18 | 26.18 | 15.97 |
| 7 | Late | 9 | 184 | 5 | 14 | 39 | 7.61 | 21.20 | 35.90 |

The observed ratio remained between 6.84 and 8.76 during rounds 1-4 and rose
during rounds 5-7. Round 7 had only nine operational cameras and should not be
interpreted as a survey-wide estimate by itself.

## Sampling and influence audit

Camera turnover was limited: 59 of 60 cameras operated in both periods. In the
model subset, 47 cameras were represented in each period and 44 occurred in
both. Among the 59 operational paired cameras, the mean within-camera change
was +1.02 juvenile detections and -3.77 adult-female detections per 100
camera-nights. However, median changes were zero, and juvenile rates increased
at only 19 cameras, reflecting the sparse juvenile response rather than a
uniform site-level shift.

Juvenile detections were not more concentrated in Late. Early detections
occurred at 18 cameras and the five largest contributors supplied 62.2% of
juveniles; Late detections occurred at 28 cameras and the five largest supplied
47.2%. The largest single camera supplied 14.9% of Early juveniles and 13.6%
of Late juveniles. Omitting each camera from the pooled observed ratio left all
50 ratio-of-ratios above 1 (range 1.326-1.676).

The fitted-model checks support a stable direction but not uniformly stable
threshold significance:

| Exploratory fit | Late/Early OR | 95% CI | P | Interpretation |
|---|---:|---:|---:|---|
| Published primary model | 1.569 | 1.011-2.436 | 0.0445 | Primary result |
| Exclude largest juvenile camera-round | 1.542 | 0.988-2.406 | 0.0566 | Same direction; interval includes 1 |
| Exclude three largest juvenile camera-rounds | 1.611 | 1.014-2.559 | 0.0436 | Same direction and supported |
| Cameras represented in both periods | 1.565 | 1.005-2.436 | 0.0476 | Same direction and supported |
| Exclude round 7 | 1.550 | 0.991-2.425 | 0.0548 | Same direction; interval includes 1 |
| Exclude rounds 6 and 7 | 1.467 | 0.887-2.428 | 0.1356 | Same direction; less precise |

The full adjustment set became computationally singular when the highest-count
camera or three highest-count cameras were removed because sparse factor
levels and nested random effects no longer supported that parameterization.
A focused fallback retaining wolf activity, jackal activity, open-space cover,
period, and the nested spatial random effects gave OR = 1.558 (95% CI
1.047-2.319) after excluding the highest-count camera and OR = 1.611 (95% CI
1.016-2.555) after excluding the three highest-count cameras.

Of 50 full leave-one-camera-out fits, 49 returned finite season estimates and
47 had positive-definite Hessians without fitting warnings. All 47 clean odds
ratios exceeded 1 (range 1.436-1.731), but only 30 had 95% intervals excluding
1. The direction is therefore not attributable to one camera, while the
strength of conventional statistical support is sensitive to the limited
juvenile information.

## Timing of the change

Replacing the binary term with a continuous round effect produced an odds
ratio of 1.139 per round (95% CI 0.990-1.310, P = 0.068; likelihood-ratio
P = 0.066). Treating round as a seven-level factor was not supported overall
(likelihood-ratio P = 0.300). In contrast, the published binary term improved
fit over the otherwise identical no-time model (likelihood-ratio P = 0.0417;
AIC 441.08 vs 443.22).

Alternative cutoffs yielded Late/Early odds ratios of 0.873 for October 3,
1.446 for October 10, 1.569 for October 17, 1.804 for October 24, and 1.901 for
October 31. Their P values were 0.590, 0.098, 0.045, 0.006, and 0.003,
respectively. This pattern is more consistent with a change building late in
the survey than with an abrupt biological transition at the prespecified
October 17 boundary. Statistical support is sensitive to the period definition.

## Event-level grouping and detectability indicators

The authoritative independent-event file contains 180 juvenile event rows
(199 animals) and 1,477 adult-female event rows (1,951 animals), all mapped to
documented active camera dates. Cross-class episodes were reconstructed within
camera sessions by joining adjacent class-event starts separated by no more
than 30 minutes.

| Period | Juvenile episodes | With adult female | Without adult female | Proportion with adult female | Adult-female-only episodes |
|---|---:|---:|---:|---:|---:|
| Early | 64 | 42 | 22 | 65.6% | 669 |
| Late | 116 | 65 | 51 | 56.0% | 701 |

Juvenile episode rows increased from 64 to 116, whereas mean detected juvenile
group size changed little (1.14 to 1.09 animals; median 1 in both periods).
Adult-female mean group size was also similar (1.30 Early and 1.34 Late;
median 1). The proportion of juvenile episodes containing an adult-female
event declined by 9.6 percentage points. This is consistent with a change in
grouping, movement, or class-specific detectability, but a camera episode
without an adult-female detection does not establish that a juvenile was
unaccompanied or maternally independent.

## Competing biological explanations

| Explanation | Assessment from available data |
|---|---|
| Juvenile maturation and changing activity | Consistent with more juvenile event rows and no increase in juvenile group size, but not directly testable because the export contains no body-size notes or finer juvenile age classes. |
| Autumn social reorganization | Consistent with the lower fraction of juvenile episodes containing adult females and the lower adult-female encounter rate, but camera co-detection cannot establish social or maternal association. |
| Autumn birth pulse | Neither supported nor contradicted. The available labels cannot identify neonates or estimate juvenile age, and no reproducibly linked image filenames or paths are present for a targeted photograph review. |
| Increased juvenile survival or recruitment | Not estimable from these camera-based composition data. The analysis does not track individuals, births, deaths, or population size. |
| Classification change as juveniles mature | Possible but unquantifiable. Broad `Juvenile` and rare `Sub-adult` labels exist, but individuals cannot be followed. Reclassification of maturing juveniles out of the juvenile class would tend to reduce, not create, later juvenile counts, although seasonal classification or detectability changes cannot be excluded. |

The corrected class histories treat juvenile and adult-female responses as
mutually exclusive: records marked Female and Juvenile are retained as
juveniles and excluded from the adult-female count. Double-counting the same
classification in both response columns therefore does not explain the result.

## Recommended Discussion interpretation

Juvenile representation relative to adult females was higher later in the
survey, but this shift should be interpreted as a seasonal change in
camera-based demographic composition rather than as direct evidence of
recruitment, survival, or population growth. The increase reflected both a
higher effort-standardized juvenile encounter rate and a lower adult-female
encounter rate, developed mainly during rounds 5-7, persisted in direction
across camera influence analyses, and was not attributable to wholesale camera
turnover or concentration at a single station. Event-level patterns were
consistent with increased juvenile detectability or changing grouping as the
season progressed, but the available classifications and metadata cannot
distinguish maturation, autumn social reorganization, a late birth pulse, or
other detection processes.

## Reproducible outputs

The audit script is `investigate_juvenile_female_seasonal_increase.R`. All new
outputs are confined to `juvenile_female_seasonal_audit/`. The principal
tables are `juvenile_female_summary_by_round.csv`,
`juvenile_female_summary_by_period.csv`, `camera_period_composition.csv`,
`juvenile_camera_concentration.csv`, `season_effect_model_sensitivity.csv`,
`temporal_model_tests.csv`, `temporal_model_coefficients.csv`,
`period_cutoff_sensitivity.csv`, `cross_class_episodes_by_round.csv`,
`cross_class_episodes_by_period.csv`, `event_group_sizes_by_round.csv`, and
`event_group_sizes_by_period.csv`. Diagnostic figures are provided as both PNG
and PDF in `juvenile_female_round_diagnostics.*` and
`paired_camera_period_diagnostics.*`.
