# Directional gazelle-pressure recurrent-event analysis

## Analysis

- Input: `../full-export_2026-04-15_2006.csv`.
- Independent-event definition: records of the same taxon, camera, and session separated by more than 30 minutes.
- Maximum follow-up after each primary event: 14 days, additionally censored by the next primary event or active camera-round end.
- Competing-pressure censoring: `FALSE`; the default dyadic analysis does not discard follow-up when another pressure taxon is detected.
- Human-associated activity includes people, soldiers, bicycles, and civilian, military, and construction vehicles, matching the earlier human-pressure covariate.
- Primary model: directional recurrent-event piecewise exponential additive mixed model (PAMM), with elapsed time and time of day as smooth terms, period as a fixed effect, and camera and camera-session frailties.
- Randomization test: 9 integer-day circular shifts of secondary events within camera sessions, preserving camera use, event totals, and clock time.
- The randomization effect compares the 0-1 day rate with the 1-7 day rate relative to neutral shifted schedules. Ratios below 1 indicate fewer short lags than expected; ratios above 1 indicate more short lags.

## Event audit

| Taxon | Independent events | Cameras | Sessions |
|---|---:|---:|---:|
| Fawn | 180 | 31 | 77 |
| Female | 1464 | 50 | 215 |
| Male | 2255 | 51 | 253 |
| Wolf | 836 | 56 | 216 |
| Jackal | 9758 | 59 | 314 |
| Boar | 2985 | 57 | 262 |
| Human | 5519 | 56 | 260 |

## Short-lag results

| Direction | Primary | Secondary | PAMM HR, 0-6 h vs late supported interval | 95% CI | Reference days | Randomization ratio | Randomization p |
|---|---|---|---:|---:|---:|---:|---:|
| Pressure to gazelle | Wolf | Fawn | 3.26 | [1.52, 7.00] | 10.5-14.0 | 1.23 | 0.200 |
| Pressure to gazelle | Wolf | Female | 1.50 | [1.06, 2.10] | 10.5-14.0 | 1.10 | 0.300 |
| Pressure to gazelle | Wolf | Male | 1.69 | [1.32, 2.16] | 10.5-14.0 | 1.27 | 0.100 |
| Pressure to gazelle | Jackal | Fawn | 1.71 | [1.18, 2.47] | 1.5-2.0 | 1.21 | 0.600 |
| Pressure to gazelle | Jackal | Female | 1.22 | [1.03, 1.43] | 1.5-2.0 | 1.11 | 0.200 |
| Pressure to gazelle | Jackal | Male | 1.33 | [1.17, 1.52] | 1.5-2.0 | 1.23 | 0.100 |
| Pressure to gazelle | Boar | Fawn | 1.43 | [0.79, 2.62] | 4.5-6.0 | 1.14 | 0.600 |
| Pressure to gazelle | Boar | Female | 1.34 | [1.09, 1.66] | 4.5-6.0 | 1.07 | 0.300 |
| Pressure to gazelle | Boar | Male | 1.54 | [1.31, 1.82] | 4.5-6.0 | 1.19 | 0.100 |
| Pressure to gazelle | Human | Fawn | 0.47 | [0.28, 0.79] | 2.6-3.5 | 0.62 | 0.100 |
| Pressure to gazelle | Human | Female | 0.86 | [0.72, 1.03] | 2.6-3.5 | 1.05 | 0.500 |
| Pressure to gazelle | Human | Male | 0.89 | [0.75, 1.06] | 2.6-3.5 | 1.12 | 0.200 |
| Gazelle to pressure | Fawn | Wolf | 2.06 | [1.05, 4.06] | 9.0-12.0 | 1.31 | 0.100 |
| Gazelle to pressure | Fawn | Jackal | 1.79 | [1.35, 2.37] | 9.0-12.0 | 1.27 | 0.100 |
| Gazelle to pressure | Fawn | Boar | 1.52 | [1.06, 2.20] | 9.0-12.0 | 1.05 | 0.800 |
| Gazelle to pressure | Fawn | Human | 1.43 | [1.06, 1.95] | 9.0-12.0 | 1.20 | 0.400 |
| Gazelle to pressure | Female | Wolf | 2.49 | [1.75, 3.54] | 7.5-10.0 | 1.23 | 0.100 |
| Gazelle to pressure | Female | Jackal | 2.21 | [1.91, 2.54] | 7.5-10.0 | 1.14 | 0.100 |
| Gazelle to pressure | Female | Boar | 2.07 | [1.65, 2.60] | 7.5-10.0 | 1.17 | 0.100 |
| Gazelle to pressure | Female | Human | 1.31 | [1.12, 1.53] | 7.5-10.0 | 1.15 | 0.100 |
| Gazelle to pressure | Male | Wolf | 2.28 | [1.74, 2.98] | 5.2-7.0 | 1.34 | 0.100 |
| Gazelle to pressure | Male | Jackal | 1.54 | [1.38, 1.71] | 5.2-7.0 | 1.18 | 0.100 |
| Gazelle to pressure | Male | Boar | 1.65 | [1.37, 1.99] | 5.2-7.0 | 1.32 | 0.100 |
| Gazelle to pressure | Male | Human | 1.12 | [0.98, 1.27] | 5.2-7.0 | 1.02 | 0.700 |

## Interpretation guardrails

The PAMM curves describe directional changes in camera visitation rates after another taxon was recorded. The camera-session randomization is the stronger test of fine-scale temporal association because it retains spatial intensity and diel activity. Neither analysis demonstrates predation, individual-level pursuit, or causal behavioral avoidance. Boar and human-associated activity are treated as perceived pressures, not taxonomic predators.

The simple median observed lag is conditional on a secondary event occurring before censoring and must not be interpreted as an unbiased waiting time. It is included only as a descriptive quantity; the recurrent-event PAMM is the inferential analysis.

## Key outputs

- `smoke_gazelle_pressure_recurrent_events_gazelle_response_curves.png`
- `smoke_gazelle_pressure_recurrent_events_pressure_response_curves.png`
- `smoke_gazelle_pressure_recurrent_events_permutation_effects.png`
- `smoke_gazelle_pressure_recurrent_events_lag_summary.csv`
- `smoke_gazelle_pressure_recurrent_events_pamm_window_contrasts.csv`
- `smoke_gazelle_pressure_recurrent_events_permutation_tests.csv`

## Method references

Ferry, N., Dupont, P., Bender, A., & Heurich, M. (2024). Introducing recurrent event analyses to assess species interactions based on camera-trap data: A comparison with time-to-first-event approaches. *Methods in Ecology and Evolution*, 15, 1233-1246. https://doi.org/10.1111/2041-210X.14338

Niedballa, J., et al. (2019). Assessing analytical methods for detecting spatiotemporal interactions between species from camera trapping data. *Remote Sensing in Ecology and Conservation*. https://doi.org/10.1002/rse2.107
