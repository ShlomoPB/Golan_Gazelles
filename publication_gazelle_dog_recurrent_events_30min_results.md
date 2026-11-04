# Directional gazelle-dog recurrent-event analysis

## Analysis

- Input: `../full-export_2026-04-15_2006.csv`.
- Sampling frame: `camelot_sessions`, from 2020-01-12 to 2020-12-15.
- Event definition: `gap`. Exact duplicate export rows were collapsed, then consecutive records within 30 min at one camera-session-taxon were grouped as one independent event.
- Variance family: `poisson`; Poisson piecewise-event rates were modeled with camera and camera-session frailties.
- Maximum follow-up after each primary event: 14 days, additionally censored by the next primary event or Camelot camera-session end.
- Competing-pressure censoring: `FALSE`; the default dyadic analysis does not discard follow-up when another pressure taxon is detected.
- Focal pressure set: `Dog`.
- Dog activity is the export category `Canis familiaris`; the available metadata do not distinguish livestock-guarding, accompanied, roaming, or feral dogs.
- Primary model: directional recurrent-event piecewise exponential additive mixed model (PAMM), with smooth effects of elapsed time, local time of day, and day of year, plus camera and camera-session frailties.
- Randomization test: 999 integer-day circular shifts of secondary events within camera sessions, preserving camera use, event totals, and clock time.
- The randomization effect compares the 0-1 day rate with the 1-7 day rate relative to neutral shifted schedules. Ratios below 1 indicate fewer short lags than expected; ratios above 1 indicate more short lags.
- Randomization p-values were adjusted by the Benjamini-Hochberg method across all 6 directional tests and, secondarily, within each 3-test directional family.

## Event audit

| Taxon | Independent events | Cameras | Sessions |
|---|---:|---:|---:|
| Fawn | 174 | 31 | 76 |
| Female | 1427 | 50 | 217 |
| Male | 2204 | 52 | 258 |
| Dog | 221 | 35 | 74 |

## Short-lag results

| Direction | Primary | Secondary | PAMM HR, 0-6 h vs late supported interval | 95% CI | Reference days | Randomization ratio | p | BH q, global | BH q, direction |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| Pressure to gazelle | Dog | Fawn | 0.28 | [0.04, 1.87] | 6.0-8.0 | 1.00 | 1.000 | 1.000 | 1.000 |
| Pressure to gazelle | Dog | Female | 0.45 | [0.20, 1.02] | 6.0-8.0 | 0.82 | 0.473 | 0.710 | 0.710 |
| Pressure to gazelle | Dog | Male | 0.77 | [0.40, 1.48] | 6.0-8.0 | 0.84 | 0.371 | 0.710 | 0.710 |
| Gazelle to pressure | Fawn | Dog | NA | [NA, NA] | NA-NA | 1.70 | 0.380 | 0.710 | 0.654 |
| Gazelle to pressure | Female | Dog | 0.95 | [0.44, 2.05] | 6.8-9.0 | 0.80 | 0.436 | 0.710 | 0.654 |
| Gazelle to pressure | Male | Dog | 0.74 | [0.36, 1.54] | 4.9-6.5 | 0.95 | 0.799 | 0.959 | 0.799 |

## Multiplicity-aware inference

No association met a 5% false-discovery-rate threshold across all 6 tests. Within the two prespecified directional families, the 5% threshold retained: none. Direction-family adjustment is secondary and should not replace the global result.

## Interpretation guardrails

The PAMM curves describe directional changes in rates of 30-min independent camera events after another taxon was recorded. The camera-session randomization retains spatial intensity, diel timing, and within-taxon temporal clustering. Neither analysis demonstrates individual-level following, predation, protection, or causal behavioral avoidance.

The simple median observed lag is conditional on a secondary event occurring before censoring and must not be interpreted as an unbiased waiting time. It is included only as a descriptive quantity; the recurrent-event PAMM is the inferential analysis.

## Key outputs

- `publication_gazelle_dog_recurrent_events_30min_gazelle_response_curves.png`
- `publication_gazelle_dog_recurrent_events_30min_pressure_response_curves.png`
- `publication_gazelle_dog_recurrent_events_30min_permutation_effects.png`
- `publication_gazelle_dog_recurrent_events_30min_lag_summary.csv`
- `publication_gazelle_dog_recurrent_events_30min_pamm_window_contrasts.csv`
- `publication_gazelle_dog_recurrent_events_30min_permutation_tests.csv`

## Method references

Ferry, N., Dupont, P., Bender, A., & Heurich, M. (2024). Introducing recurrent event analyses to assess species interactions based on camera-trap data: A comparison with time-to-first-event approaches. *Methods in Ecology and Evolution*, 15, 1233-1246. https://doi.org/10.1111/2041-210X.14338

Niedballa, J., et al. (2019). Assessing analytical methods for detecting spatiotemporal interactions between species from camera trapping data. *Remote Sensing in Ecology and Conservation*. https://doi.org/10.1002/rse2.107
