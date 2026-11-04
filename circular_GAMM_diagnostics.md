# Circular GAMM diagnostics

## Data and effort

The frozen independent-event file contributed exactly 3925 class-events: 180 juvenile, 1,477 adult-female, and 2,268 adult-male events.
The analysis table contains 8,256 rows (344 positive-effort camera-rounds x 24 apparent-solar-hour bins).
Each hourly cell used `Days.Operational` as its exposure because exact within-day malfunction times were unavailable. Thus, the offset estimates events in a given solar-hour bin per operational camera-day.
Wolf and jackal indices were standardized across the 60 cameras (wolf mean = 0.579, SD = 0.202; jackal mean = 0.500, SD = 0.226).

## Primary specification

Separate negative-binomial `mgcv::bam` models used a log link, fREML, cyclic cubic bases over 0-24 h, shrinkage (`select = TRUE`), factor(round), and random intercept smooths for camera-round, camera, and polygon.
Because supported camera-specific timing heterogeneity materially improved adult fit, the final female and male models also contained shared-penalty camera-specific cyclic deviations (k = 4 per camera). The juvenile model retained random intercepts only because 180 events were insufficient to support 60 camera-level activity curves.
The ordinary wolf, jackal, and interaction terms describe overall event intensity; their time-varying cyclic smooths describe changes in diel shape.

## Class-specific diagnostics

### Juvenile

Formula: `n_events ~ wolf_z * jackal_z + round_factor + s(solar_hour, bs = "cc",      k = 10) + s(solar_hour, by = wolf_z, bs = "cc", k = 10) +      s(solar_hour, by = jackal_z, bs = "cc", k = 10) + s(solar_hour,      by = wolf_jackal_z, bs = "cc", k = 10) + s(camera_round_factor,      bs = "re") + s(camera_factor, bs = "re") + s(polygon_factor,      bs = "re") + offset(log(effort_days))`

The model used k = 10 with `random intercepts only`, theta = 4.107, and explained 46.6% of deviance. Convergence status: `full convergence`.
Observed zeros = 98.1%; mean model-implied zero probability = 98.0%; Pearson dispersion = 0.421.
Mean within-camera-round lag-1 Pearson residual correlation = 0.727. Maximum estimated ecological-smooth concurvity = 0.225.
Temporal-smooth k-index range = 0.889-0.889; minimum k-check P = 0.375.

| Time-varying term | EDF | Approximate P |
|---|---:|---:|
| Wolf time-varying effect | 1.02 | 0.086 |
| Jackal time-varying effect | 0.00 | 0.947 |
| Wolf x jackal time-varying effect | 2.02 | <0.001 |

Native `mgcv` residual, k-index, and concurvity diagnostics were used; `gratia` and `DHARMa` were not required or installed.

### Adult female

Formula: `n_events ~ wolf_z * jackal_z + round_factor + s(solar_hour, bs = "cc",      k = 8) + s(solar_hour, by = wolf_z, bs = "cc", k = 8) + s(solar_hour,      by = jackal_z, bs = "cc", k = 8) + s(solar_hour, by = wolf_jackal_z,      bs = "cc", k = 8) + s(solar_hour, by = camera_factor, bs = "cc",      k = 4, m = 1, id = 20) + s(camera_round_factor, bs = "re") +      s(camera_factor, bs = "re") + s(polygon_factor, bs = "re") +      offset(log(effort_days))`

The model used k = 8 with `camera-specific cyclic deviations`, theta = 3.700, and explained 52.8% of deviance. Convergence status: `full convergence`.
Observed zeros = 88.5%; mean model-implied zero probability = 88.0%; Pearson dispersion = 0.627.
Mean within-camera-round lag-1 Pearson residual correlation = 0.332. Maximum estimated ecological-smooth concurvity = 0.542.
Temporal-smooth k-index range = 0.873-0.873; minimum k-check P = 0.230.

| Time-varying term | EDF | Approximate P |
|---|---:|---:|
| Wolf time-varying effect | 3.83 | <0.001 |
| Jackal time-varying effect | 0.00 | 0.315 |
| Wolf x jackal time-varying effect | 0.08 | 0.062 |

Native `mgcv` residual, k-index, and concurvity diagnostics were used; `gratia` and `DHARMa` were not required or installed.

### Adult male

Formula: `n_events ~ wolf_z * jackal_z + round_factor + s(solar_hour, bs = "cc",      k = 8) + s(solar_hour, by = wolf_z, bs = "cc", k = 8) + s(solar_hour,      by = jackal_z, bs = "cc", k = 8) + s(solar_hour, by = wolf_jackal_z,      bs = "cc", k = 8) + s(solar_hour, by = camera_factor, bs = "cc",      k = 4, m = 1, id = 20) + s(camera_round_factor, bs = "re") +      s(camera_factor, bs = "re") + s(polygon_factor, bs = "re") +      offset(log(effort_days))`

The model used k = 8 with `camera-specific cyclic deviations`, theta = 6.464, and explained 44.4% of deviance. Convergence status: `full convergence`.
Observed zeros = 81.6%; mean model-implied zero probability = 81.0%; Pearson dispersion = 0.762.
Mean within-camera-round lag-1 Pearson residual correlation = 0.204. Maximum estimated ecological-smooth concurvity = 0.542.
Temporal-smooth k-index range = 0.904-0.904; minimum k-check P = 0.562.

| Time-varying term | EDF | Approximate P |
|---|---:|---:|
| Wolf time-varying effect | 0.00 | 0.064 |
| Jackal time-varying effect | 0.00 | 0.643 |
| Wolf x jackal time-varying effect | 3.01 | 0.002 |

Native `mgcv` residual, k-index, and concurvity diagnostics were used; `gratia` and `DHARMa` were not required or installed.

## Camera-specific curve sensitivity

A shared-penalty camera-specific cyclic deviation (k = 4 per camera) was fitted following the hierarchical timing logic of Iannarilli et al. (2025). It was promoted for adult classes only when it converged and improved AIC by more than 10; the juvenile version remained a sensitivity regardless of fit because only 180 events informed 60 cameras.

| Class | Convergence | Delta AIC (site shape - random-intercept candidate) | Site-shape EDF | Deviance explained | Selected final |
|---|---|---:|---:|---:|---|
| Juvenile | full convergence | 2.31 | 7.81 | 47.4% | FALSE |
| Adult female | full convergence | -209.64 | 63.93 | 52.8% | TRUE |
| Adult male | full convergence | -177.88 | 64.08 | 44.4% | TRUE |

AIC values from penalized fREML fits are treated descriptively rather than as formal likelihood-ratio tests. Primary inference remains based on the prespecified population-level smooths, shrinkage, diagnostics, and prediction simulations.

## Early/Late curve-shape sensitivity

The frozen survey midpoint was `2020-10-17`. A numeric Late-period cyclic difference smooth (k = 6) controlled broad period-specific changes in diel shape while factor(round) continued to control round-specific event intensity. This sensitivity did not fit predator-by-period interactions.

| Class | Effect | Primary EDF | Primary P | Season-adjusted EDF | Season-adjusted P | Delta AIC |
|---|---|---:|---:|---:|---:|---:|
| Juvenile | Wolf time-varying effect | 1.02 | 0.086 | 1.02 | 0.086 | 0.00 |
| Juvenile | Jackal time-varying effect | 0.00 | 0.947 | 0.00 | 0.947 | 0.00 |
| Juvenile | Wolf x jackal time-varying effect | 2.02 | <0.001 | 2.02 | <0.001 | 0.00 |
| Juvenile | Late versus Early curve-shape control | NA | NA | 0.00 | 0.484 | 0.00 |
| Adult female | Wolf time-varying effect | 3.83 | <0.001 | 3.87 | <0.001 | -41.63 |
| Adult female | Jackal time-varying effect | 0.00 | 0.315 | 0.00 | 0.280 | -41.63 |
| Adult female | Wolf x jackal time-varying effect | 0.08 | 0.062 | 0.02 | 0.065 | -41.63 |
| Adult female | Late versus Early curve-shape control | NA | NA | 2.97 | <0.001 | -41.63 |
| Adult male | Wolf time-varying effect | 0.00 | 0.064 | 0.15 | 0.052 | -19.86 |
| Adult male | Jackal time-varying effect | 0.00 | 0.643 | 0.00 | 0.630 | -19.86 |
| Adult male | Wolf x jackal time-varying effect | 3.01 | 0.002 | 2.94 | 0.002 | -19.86 |
| Adult male | Late versus Early curve-shape control | NA | NA | 2.92 | <0.001 | -19.86 |

The three supported predator-timing smooths retained the same inference after this adjustment. AIC differences are descriptive because the compared penalized fREML specifications differ in fixed-effect structure.

