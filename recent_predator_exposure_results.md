# Recent local predator exposure analysis

## 1. PAMM audit

The existing wolf-jackal PAMM was `analyze_gazelle_pressure_recurrent_events.R`, with results in `publication_gazelle_pressure_recurrent_events_results.md`. It modeled the recurrent-event rate of a demographic gazelle class after each focal predator timestamp at the same camera and Camelot camera session. Exact duplicate export rows at a timestamp were collapsed, but the retained wolf-jackal run used every distinct timestamp rather than 30-minute-independent events. Follow-up began at a wolf or jackal timestamp, ended at the next same-predator event, camera-session end, or 14 days, and could cross calendar dates. Thus, the 0-6 h interval necessarily represented the same or immediately following date, but the complete risk interval was not restricted to one date.

The response was the number of fawn, adult-female, or adult-male timestamp visits in piecewise follow-up intervals, with log interval duration as an offset. Elapsed time since the predator event, clock hour, and day of year were smooth terms. Camera and camera-session random-effect smooths were fitted; all focal wolf/jackal models retained both frailties. Wolf and jackal were analyzed in separate directional dyads, so the PAMM contained no joint wolf and jackal main effects or interaction. Its 0-6 h and 6-24 h summaries were contrasts against a late, taxon-specific data-supported reference interval, not against otherwise matched hours with no recent predator passage.

| Predator | Gazelle class | PAMM window | Hazard ratio (95% CI) | Randomization P | Global BH q |
|---|---|---|---:|---:|---:|
| Wolf | Fawn | 0-6 hours | 3.19 [1.84, 5.53] | 0.182 | 0.437 |
| Wolf | Fawn | 6-24 hours | 4.03 [2.63, 6.15] | 0.182 | 0.437 |
| Wolf | Female | 0-6 hours | 1.20 [0.90, 1.60] | 0.484 | 0.774 |
| Wolf | Female | 6-24 hours | 1.71 [1.37, 2.13] | 0.484 | 0.774 |
| Wolf | Male | 0-6 hours | 1.53 [1.14, 2.04] | 0.054 | 0.192 |
| Wolf | Male | 6-24 hours | 1.53 [1.18, 1.98] | 0.054 | 0.192 |
| Jackal | Fawn | 0-6 hours | 1.32 [1.04, 1.68] | 0.326 | 0.573 |
| Jackal | Fawn | 6-24 hours | 1.02 [0.85, 1.21] | 0.326 | 0.573 |
| Jackal | Female | 0-6 hours | 0.92 [0.81, 1.05] | 0.736 | 0.976 |
| Jackal | Female | 6-24 hours | 0.99 [0.92, 1.07] | 0.736 | 0.976 |
| Jackal | Male | 0-6 hours | 1.10 [0.97, 1.25] | 0.032 | 0.192 |
| Jackal | Male | 6-24 hours | 1.05 [0.96, 1.15] | 0.032 | 0.192 |

Although several model-based PAMM intervals excluded one, no pressure-to-gazelle association passed the camera-session circular-shift test after global false-discovery-rate adjustment. The randomization was the stronger PAMM check because it preserved site use, event totals, clock time, and within-taxon clustering.

## 2. Is this analysis redundant with the PAMM?

No. The PAMM already addressed the broad question of directional temporal clustering after predator records, but it did not estimate the requested matched hourly odds under explicit 0-6 h and >6-24 h independent-passage exposures. The new analysis uses 30-minute-independent predator events, binary hourly gazelle encounter histories, simultaneous additive wolf and jackal terms, and exact matching by camera, survey round, and apparent-solar hour. It should be viewed as a focused alternative estimand and robustness analysis, not as an unrelated second test.

## 3. Data structure and exposure counts

We analyzed 143,472 camera-hours from 5,978 known operational camera-days. Camera 45 round 2 was excluded because its 21 operational days occurred within 32 nominal dates but the 11 outage dates were unavailable. Excluding this round prevented unobserved outage hours from being coded as gazelle absences. Other rounds matched their start-inclusive, end-exclusive date spans exactly.
Predator exposures were derived from 852 independent wolf events and 10262 independent jackal events. A predator record after the start of a response hour could not classify that hour as exposed.

| Recent exposure category (preceding 24 h) | Camera-hours | Percent |
|---|---:|---:|
| Neither | 59,031 | 41.1% |
| Wolf only |  4,670 | 3.3% |
| Jackal only | 67,543 | 47.1% |
| Both | 12,228 | 8.5% |

| Class | Positive hours | Wolf 0-6 h | Wolf >6-24 h | Jackal 0-6 h | Jackal >6-24 h | Both predators within 24 h |
|---|---:|---:|---:|---:|---:|---:|
| Adult male | 2237 | 100 | 280 | 525 | 1166 | 249 |
| Juvenile | 176 | 11 | 25 | 48 | 98 | 26 |
| Adult female | 1449 | 52 | 162 | 326 | 718 | 143 |

Wolf-jackal recent-exposure interactions were not fitted. The prespecified adequacy threshold was at least 30 dual-exposure positive hours in every class; observed class counts were Juvenile = 26, Adult female = 143, Adult male = 249.

## 4. Model choice

For each demographic class, we fitted a conditional logistic case-crossover model with recent wolf and jackal exposure indicators in the same additive model. Strata were camera x survey round x apparent-solar hour. Consequently, each coefficient compared exposed and unexposed operational dates at the same camera, during the same one-hour diel bin and survey period. Cluster-robust standard errors accounted for repeated strata within cameras. This design more directly removes stable camera-level habitat and predator-index differences than a conventional random-intercept GLMM and controls the diel cycle without requiring an additional estimated smooth.

The response was at least one independent class encounter in the hour. Exposures were evaluated at the hour start, making the design conservative for predator detections occurring earlier within the same hour. Odds ratios are conditional associations and do not identify individual gazelles or demonstrate causal avoidance or attraction.
Reported P values are nominal and were not adjusted across classes and lag windows; inference therefore emphasizes confidence intervals, exposed outcomes, and the prespecified within-predator lag contrasts.

## 5. Results by demographic class

| Class | Recent predator exposure | Odds ratio (95% CI) | P | Exposed positive hours | Model positive hours | Assessment |
|---|---|---:|---:|---:|---:|---|
| Juvenile | Wolf: 0-6 h | 1.47 [0.63, 3.45] | 0.375 | 11 | 176 | Sparse exposed juvenile outcomes |
| Juvenile | Wolf: >6-24 h | 1.32 [0.80, 2.19] | 0.273 | 25 | 176 | Sparse exposed juvenile outcomes |
| Juvenile | Jackal: 0-6 h | 1.41 [0.76, 2.62] | 0.280 | 48 | 176 | Unsupported |
| Juvenile | Jackal: >6-24 h | 1.09 [0.73, 1.61] | 0.685 | 98 | 176 | Unsupported |
| Adult female | Wolf: 0-6 h | 0.88 [0.63, 1.23] | 0.440 | 52 | 1449 | Unsupported |
| Adult female | Wolf: >6-24 h | 1.06 [0.90, 1.26] | 0.496 | 162 | 1449 | Unsupported |
| Adult female | Jackal: 0-6 h | 1.12 [0.99, 1.27] | 0.077 | 326 | 1449 | Unsupported |
| Adult female | Jackal: >6-24 h | 1.18 [0.98, 1.41] | 0.081 | 718 | 1449 | Unsupported |
| Adult male | Wolf: 0-6 h | 1.05 [0.85, 1.28] | 0.658 | 100 | 2237 | Unsupported |
| Adult male | Wolf: >6-24 h | 1.14 [1.01, 1.30] | 0.041 | 280 | 2237 | Nominal support |
| Adult male | Jackal: 0-6 h | 1.06 [0.93, 1.22] | 0.383 | 525 | 2237 | Unsupported |
| Adult male | Jackal: >6-24 h | 1.19 [1.04, 1.35] | 0.011 | 1166 | 2237 | Nominal support |

Nominal primary associations were: Adult male Wolf: >6-24 h (OR = 1.14); Adult male Jackal: >6-24 h (OR = 1.19). These are unadjusted, prespecified coefficient-level tests.
None of the 0-6 h coefficients or within-predator lag contrasts had P < 0.05. The analysis therefore did not resolve an immediate reactive change in hourly gazelle detection following predator passage. The two adult-male >6-24 h estimates are modest nominal signals; because neither differed from its corresponding 0-6 h coefficient, they do not establish a delayed response.

Lag and between-predator contrasts compare coefficients within the same fitted model:

| Class | Contrast | Ratio of odds ratios (95% CI) | P |
|---|---|---:|---:|
| Juvenile | Wolf versus jackal: 0-6 h | 1.05 [0.44, 2.47] | 0.918 |
| Juvenile | Wolf versus jackal: >6-24 h | 1.22 [0.54, 2.78] | 0.635 |
| Juvenile | Jackal: 0-6 h versus >6-24 h | 1.30 [0.62, 2.72] | 0.493 |
| Juvenile | Wolf: 0-6 h versus >6-24 h | 1.11 [0.42, 2.91] | 0.830 |
| Adult female | Wolf versus jackal: 0-6 h | 0.78 [0.56, 1.10] | 0.153 |
| Adult female | Wolf versus jackal: >6-24 h | 0.90 [0.70, 1.16] | 0.419 |
| Adult female | Jackal: 0-6 h versus >6-24 h | 0.95 [0.77, 1.18] | 0.657 |
| Adult female | Wolf: 0-6 h versus >6-24 h | 0.83 [0.54, 1.26] | 0.379 |
| Adult male | Wolf versus jackal: 0-6 h | 0.99 [0.77, 1.26] | 0.907 |
| Adult male | Wolf versus jackal: >6-24 h | 0.96 [0.80, 1.16] | 0.697 |
| Adult male | Jackal: 0-6 h versus >6-24 h | 0.90 [0.74, 1.09] | 0.262 |
| Adult male | Wolf: 0-6 h versus >6-24 h | 0.92 [0.73, 1.14] | 0.438 |

Juvenile inference should be judged from exposed positive-hour counts and interval width, not from P values alone. An unsupported juvenile coefficient is not evidence that juveniles fail to respond when the relevant exposure-outcome combination is sparse.

## 6. Comparison with circular GAM and previous dynamic analysis

- **Main hierarchical circular GAM:** asks whether diel shape changes along continuous, site-level wolf and jackal activity gradients after accounting for hierarchical camera structure. Predator indices are spatially constant and do not represent a recent passage.
- **Previous camera x round x hour dynamic GAMM:** asks whether gazelle timing covaries with recurring local predator activity patterns averaged within camera-round-hour combinations. It represents repeated local temporal niches, not the occurrence of one immediately preceding predator event.
- **Existing PAMM:** asks how the recurrent gazelle timestamp-visit hazard changes smoothly after a focal predator timestamp in a same-camera session, with dyads fitted separately and a circular-shift null test.
- **New recent-exposure analysis:** asks whether binary gazelle detection odds differ after an actual independent local predator passage in the previous 0-6 h or >6-24 h, relative to matched dates at the same camera, round, and apparent-solar hour.

## 7. Manuscript recommendation

Treat this as a focused supplementary acute-response analysis. The modest adult-male >6-24 h associations should be labeled nominal and should not be interpreted as a delayed response because the direct lag contrasts were unsupported. The analysis should not replace the circular GAM's site-level diel inference or be described as causal avoidance or attraction.

The manuscript itself was not edited.

## Outputs

- `recent_predator_exposure_summary.csv`
- `recent_predator_exposure_results.md`
- `recent_predator_exposure_figure.png`
- `recent_predator_exposure_data_audit.csv`
