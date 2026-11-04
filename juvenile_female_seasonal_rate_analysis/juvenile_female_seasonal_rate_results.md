# Exploratory seasonal decomposition of the juvenile:adult-female ratio

## Scope and existing inputs

This analysis is exploratory and does not modify the primary beta-binomial model or manuscript outputs. It reuses `gazelle_juvenile_abundance_60x7.csv`, `gazelle_female_abundance_60x7.csv`, `ObsCovs_050524.csv`, `SiteCovs_060324.csv`, and the clock- and classification-corrected `Gazella gazella_independent_events_30min_analysis_ready.csv`. The published ratio workflow is `analyze_fawn_female_beta_binomial_wolf_jackal.R`; the earlier influence and metadata audit is `investigate_juvenile_female_seasonal_increase.R`.

The source-code definition assigns a camera-round to Early when its midpoint is on or before 17 October 2020 and to Late otherwise. Thus, Late technically begins 18 October; no camera-round midpoint occurs exactly on 17 October, so the common shorthand does not change the fitted data.

## Reconstruction

The count histories contain 199 juvenile and 1951 adult-female individuals across 344 positive-effort camera-rounds from 60 cameras. The independent-event file contains 180 juvenile and 1477 adult-female 30-min events. Summing each event's maximum simultaneously observed individuals reproduces the overall totals exactly. However, 46 camera-round rows differ because the count histories retain session-based allocation whereas the event file is assigned by corrected event date. The project cleanup audit identifies boundary-date allocation and manual curation as likely causes but does not provide a complete historical session-to-round crosswalk. Counts are repeated detections, not known unique animals.

## Descriptive decomposition

Early contributed 2701 camera-nights, 74 juvenile individuals, and 934 adult-female individuals; Late contributed 3298 camera-nights, 125 juvenile individuals, and 1017 adult-female individuals. The pooled juvenile rate increased by 38.3%, the pooled adult-female rate changed by -10.8%, and the observed juvenile:female ratio increased by 55.1%.

The increase was expressed primarily as more juvenile encounters rather than larger juvenile groups: juvenile event rows increased from 64 to 116, while mean detected group size changed from 1.14 to 1.09.

## Count-model results

Negative-binomial mixed models included every positive-effort camera-round, retained zero detections, used log camera-nights as an offset, and accounted for repeated cameras. In models using the same stable additive adjustment set as the published ratio analysis, the Late/Early juvenile count-rate ratio was 1.76 (95% CI 1.04-2.99, P = 0.035), whereas the adult-female rate ratio was 1.20 (95% CI 0.92-1.57, P = 0.181).

The event-count models gave the same interpretation: juvenile independent events increased (rate ratio 1.64, 95% CI 1.03-2.62, P = 0.037), but adult-female events did not show a supported change (rate ratio 1.03, 95% CI 0.82-1.30, P = 0.773). The direct joint class-by-period contrast was close to, but did not cross, the conventional threshold (ratio of Late/Early rate ratios = 1.43, 95% CI 0.99-2.06, Wald P = 0.053).

Restricting the analysis to the 59 cameras represented in both periods left all effect directions essentially unchanged. Camera turnover therefore does not explain the pattern, although the limited juvenile information makes interval estimates imprecise.

## Biological interpretation

- **Additional autumn births:** neither supported nor excluded. The available labels contain no newborn category, body-size notes, or finer juvenile ages.
- **Increasing juvenile mobility or detectability:** consistent with the increase in juvenile event frequency, the absence of increased detected group size, and the spread of juvenile detections across more cameras in the prior influence audit. It is not directly proven.
- **Adult-female movement during the rut:** the pooled adult-female rate declined, but mixed models did not support a decline after accounting for camera and effort. This mechanism remains possible but receives weak direct support here.
- **Early juvenile mortality:** not estimable from these repeated-detection data and not, by itself, an explanation for an increase during the survey. It remains relevant background for interpreting predator associations.

These are camera-based detection rates and composition indices. They do not estimate survival, recruitment, births, or absolute abundance.

## Recommendation

Add one cautious sentence to Section 3.2 and place the three-panel figure in the supplement. A defensible sentence is: 'The later increase in juvenile representation was driven primarily by a higher frequency of juvenile camera detections rather than larger detected juvenile groups; adult-female rates showed no supported decline in mixed models, and the available classifications could not distinguish increased juvenile detectability from additional births.'

## Reproducible outputs

- `juvenile_female_counts_rates_by_round.csv`: requested seven-round descriptive table with camera-bootstrap intervals.
- `juvenile_female_counts_rates_by_period.csv`: Early/Late descriptive decomposition.
- `seasonal_rate_model_results.csv`: effect sizes, intervals, tests, and convergence diagnostics.
- `event_count_reconstruction_audit.csv`: count-history versus independent-event verification.
- `event_count_reconstruction_mismatches.csv`: camera-round allocation differences retained for audit.
- `juvenile_female_seasonal_rates_supplement.png` and `.pdf`: supplementary figure.
- `run_juvenile_female_seasonal_rate_analysis.R`: complete analysis script.
