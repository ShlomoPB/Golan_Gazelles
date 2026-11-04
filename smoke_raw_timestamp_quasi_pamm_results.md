# Directional gazelle-pressure recurrent-event analysis

## Analysis

- Input: `../full-export_2026-04-15_2006.csv`.
- Sampling frame: `camelot_sessions`, from 2020-01-01 to 2020-12-15.
- Event definition: `timestamp`. Exact duplicate export rows at one camera-session-taxon timestamp were collapsed, but every distinct timestamp was retained; no independence interval was imposed.
- Variance family: `quasipoisson`; quasi-Poisson uncertainty accounts for extra-Poisson variation among timestamp-level counts.
- Maximum follow-up after each primary event: 14 days, additionally censored by the next primary event or Camelot camera-session end.
- Competing-pressure censoring: `FALSE`; the default dyadic analysis does not discard follow-up when another pressure taxon is detected.
- Human-associated activity includes people, soldiers, bicycles, and civilian, military, and construction vehicles, matching the earlier human-pressure covariate.
- Primary model: directional recurrent-event piecewise exponential additive mixed model (PAMM), with smooth effects of elapsed time, local time of day, and day of year, plus camera and camera-session frailties.
- Randomization test: 3 integer-day circular shifts of secondary events within camera sessions, preserving camera use, event totals, and clock time.
- The randomization effect compares the 0-1 day rate with the 1-7 day rate relative to neutral shifted schedules. Ratios below 1 indicate fewer short lags than expected; ratios above 1 indicate more short lags.

## Event audit

| Taxon | Distinct timestamp visits | Cameras | Sessions |
|---|---:|---:|---:|
| Fawn | 283 | 31 | 76 |
| Female | 3288 | 50 | 217 |
| Male | 3010 | 52 | 258 |
| Wolf | 1084 | 56 | 213 |
| Jackal | 13839 | 60 | 322 |
| Boar | 5499 | 58 | 270 |
| Human | 10808 | 57 | 263 |

## Short-lag results

| Direction | Primary | Secondary | PAMM HR, 0-6 h vs late supported interval | 95% CI | Reference days | Randomization ratio | Randomization p |
|---|---|---|---:|---:|---:|---:|---:|
| Pressure to gazelle | Wolf | Fawn | 3.19 | [1.84, 5.53] | 9.0-12.0 | 1.64 | 0.750 |
| Pressure to gazelle | Wolf | Female | 1.20 | [0.90, 1.60] | 9.0-12.0 | 1.06 | 0.750 |
| Pressure to gazelle | Wolf | Male | 1.53 | [1.14, 2.04] | 9.0-12.0 | 1.19 | 0.250 |
| Pressure to gazelle | Jackal | Fawn | 1.32 | [1.04, 1.68] | 1.1-1.5 | 1.25 | 0.500 |
| Pressure to gazelle | Jackal | Female | 0.92 | [0.81, 1.05] | 1.1-1.5 | 0.99 | 0.750 |
| Pressure to gazelle | Jackal | Male | 1.10 | [0.97, 1.25] | 1.1-1.5 | 1.13 | 0.250 |
| Pressure to gazelle | Boar | Fawn | 0.98 | [0.71, 1.34] | 3.4-4.5 | 1.16 | 0.500 |
| Pressure to gazelle | Boar | Female | 1.29 | [1.09, 1.52] | 3.4-4.5 | 1.08 | 0.500 |
| Pressure to gazelle | Boar | Male | 1.34 | [1.16, 1.55] | 3.4-4.5 | 1.19 | 0.250 |
| Pressure to gazelle | Human | Fawn | 0.22 | [0.16, 0.30] | 1.5-2.0 | 0.72 | 0.250 |
| Pressure to gazelle | Human | Female | 0.72 | [0.61, 0.85] | 1.5-2.0 | 0.96 | 0.500 |
| Pressure to gazelle | Human | Male | 0.81 | [0.70, 0.94] | 1.5-2.0 | 1.06 | 0.250 |
| Gazelle to pressure | Fawn | Wolf | 1.04 | [0.55, 1.94] | 7.5-10.0 | 1.09 | 0.750 |
| Gazelle to pressure | Fawn | Jackal | 1.56 | [1.19, 2.05] | 7.5-10.0 | 1.39 | 0.250 |
| Gazelle to pressure | Fawn | Boar | 1.39 | [0.96, 2.02] | 7.5-10.0 | 1.14 | 0.500 |
| Gazelle to pressure | Fawn | Human | 1.26 | [0.87, 1.83] | 7.5-10.0 | 0.95 | 0.500 |
| Gazelle to pressure | Female | Wolf | 1.65 | [1.24, 2.18] | 4.9-6.5 | 0.91 | 0.500 |
| Gazelle to pressure | Female | Jackal | 1.35 | [1.19, 1.53] | 4.9-6.5 | 1.12 | 0.250 |
| Gazelle to pressure | Female | Boar | 1.47 | [1.22, 1.77] | 4.9-6.5 | 1.07 | 0.500 |
| Gazelle to pressure | Female | Human | 1.25 | [1.07, 1.46] | 4.9-6.5 | 1.29 | 0.250 |
| Gazelle to pressure | Male | Wolf | 1.67 | [1.30, 2.14] | 4.5-6.0 | 1.29 | 0.500 |
| Gazelle to pressure | Male | Jackal | 1.19 | [1.08, 1.32] | 4.5-6.0 | 1.09 | 0.250 |
| Gazelle to pressure | Male | Boar | 1.30 | [1.10, 1.54] | 4.5-6.0 | 1.22 | 0.250 |
| Gazelle to pressure | Male | Human | 1.01 | [0.88, 1.17] | 4.5-6.0 | 0.95 | 0.750 |

## Interpretation guardrails

The PAMM curves describe directional changes in camera trigger rates after another taxon was recorded. Because consecutive timestamps may represent the same individual or group, this timestamp-level analysis estimates detection-sequence association rather than independent animal encounters. The camera-session randomization is the stronger test of fine-scale temporal association because it retains spatial intensity, diel timing, and within-taxon temporal clustering. Neither analysis demonstrates predation, individual-level pursuit, or causal behavioral avoidance. Boar and human-associated activity are treated as perceived pressures, not taxonomic predators.

The simple median observed lag is conditional on a secondary event occurring before censoring and must not be interpreted as an unbiased waiting time. It is included only as a descriptive quantity; the recurrent-event PAMM is the inferential analysis.

## Key outputs

- `smoke_raw_timestamp_quasi_pamm_gazelle_response_curves.png`
- `smoke_raw_timestamp_quasi_pamm_pressure_response_curves.png`
- `smoke_raw_timestamp_quasi_pamm_permutation_effects.png`
- `smoke_raw_timestamp_quasi_pamm_lag_summary.csv`
- `smoke_raw_timestamp_quasi_pamm_pamm_window_contrasts.csv`
- `smoke_raw_timestamp_quasi_pamm_permutation_tests.csv`

## Method references

Ferry, N., Dupont, P., Bender, A., & Heurich, M. (2024). Introducing recurrent event analyses to assess species interactions based on camera-trap data: A comparison with time-to-first-event approaches. *Methods in Ecology and Evolution*, 15, 1233-1246. https://doi.org/10.1111/2041-210X.14338

Niedballa, J., et al. (2019). Assessing analytical methods for detecting spatiotemporal interactions between species from camera trapping data. *Remote Sensing in Ecology and Conservation*. https://doi.org/10.1002/rse2.107
