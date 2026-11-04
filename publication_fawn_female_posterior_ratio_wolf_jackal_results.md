# Posterior fawn:female ratios across wolf-jackal contexts

## Estimand

Each cell is the posterior distribution of `100 * expected fawn abundance / expected female abundance`, calculated within the same posterior draw from a joint two-factor negative-binomial `lfMsNMix` model. Predictions exclude latent site factors (`include.w = FALSE`) and use continuous standardized predator covariates at `-1`, `0`, and `+1` SD.

Maximum abundance-slope Rhat: `1.12`; minimum abundance-slope ESS: `396`.

Maximum grid-cell log-ratio Rhat: `1.87`; minimum grid-cell log-ratio ESS: `39`.

## Ratio predictions

| Period | Wolf z | Jackal z | Median fawns/100 females | 95% CrI |
|---|---:|---:|---:|---:|
| `Overall` | `-1` | `-1` | `16.1` | `[3.7, 79.0]` |
| `Overall` | `0` | `-1` | `16.3` | `[3.9, 70.3]` |
| `Overall` | `1` | `-1` | `16.2` | `[3.4, 75.0]` |
| `Overall` | `-1` | `0` | `17.7` | `[4.5, 87.7]` |
| `Overall` | `0` | `0` | `15.8` | `[4.1, 71.2]` |
| `Overall` | `1` | `0` | `13.9` | `[3.3, 63.9]` |
| `Overall` | `-1` | `1` | `19.2` | `[3.8, 131.3]` |
| `Overall` | `0` | `1` | `15.1` | `[3.7, 83.0]` |
| `Overall` | `1` | `1` | `12.0` | `[2.6, 68.9]` |
| `Early (hotter)` | `-1` | `-1` | `3.9` | `[0.6, 41.6]` |
| `Early (hotter)` | `0` | `-1` | `5.2` | `[1.0, 38.1]` |
| `Early (hotter)` | `1` | `-1` | `7.0` | `[1.3, 45.7]` |
| `Early (hotter)` | `-1` | `0` | `7.1` | `[1.3, 50.4]` |
| `Early (hotter)` | `0` | `0` | `6.5` | `[1.2, 47.2]` |
| `Early (hotter)` | `1` | `0` | `5.9` | `[1.1, 50.1]` |
| `Early (hotter)` | `-1` | `1` | `13.3` | `[1.8, 104.4]` |
| `Early (hotter)` | `0` | `1` | `8.1` | `[1.3, 68.3]` |
| `Early (hotter)` | `1` | `1` | `5.0` | `[0.7, 72.9]` |
| `Late (mild-cold)` | `-1` | `-1` | `16.6` | `[2.8, 137.5]` |
| `Late (mild-cold)` | `0` | `-1` | `14.6` | `[2.7, 115.5]` |
| `Late (mild-cold)` | `1` | `-1` | `13.2` | `[2.5, 116.9]` |
| `Late (mild-cold)` | `-1` | `0` | `17.4` | `[3.1, 127.7]` |
| `Late (mild-cold)` | `0` | `0` | `15.1` | `[2.8, 109.0]` |
| `Late (mild-cold)` | `1` | `0` | `13.2` | `[2.4, 103.0]` |
| `Late (mild-cold)` | `-1` | `1` | `18.9` | `[2.7, 164.7]` |
| `Late (mild-cold)` | `0` | `1` | `15.9` | `[2.7, 118.5]` |
| `Late (mild-cold)` | `1` | `1` | `13.5` | `[2.0, 110.9]` |

## Differential predator coefficients

Contrasts are fawn minus female coefficients on the log expected-abundance scale.

| Period | Contrast | Median | 95% CrI | Pr(<0) |
|---|---|---:|---:|---:|
| `Overall` | `wolf` | `-0.13` | `[-0.54, 0.28]` | `0.744` |
| `Overall` | `jackal` | `-0.00` | `[-0.54, 0.47]` | `0.504` |
| `Overall` | `wolf_x_jackal` | `-0.12` | `[-0.62, 0.39]` | `0.684` |
| `Early (hotter)` | `wolf` | `-0.09` | `[-0.55, 0.43]` | `0.656` |
| `Early (hotter)` | `jackal` | `0.23` | `[-0.34, 0.81]` | `0.207` |
| `Early (hotter)` | `wolf_x_jackal` | `-0.38` | `[-1.01, 0.30]` | `0.892` |
| `Late (mild-cold)` | `wolf` | `-0.13` | `[-0.53, 0.25]` | `0.752` |
| `Late (mild-cold)` | `jackal` | `0.02` | `[-0.46, 0.48]` | `0.472` |
| `Late (mild-cold)` | `wolf_x_jackal` | `-0.03` | `[-0.51, 0.46]` | `0.546` |

## Posterior predictive checks

Site-grouped Freeman-Tukey discrepancies; values near 0.5 indicate better calibration.

| Period | Class | Bayesian p |
|---|---|---:|
| `Overall` | `fawn` | `0.461` |
| `Overall` | `male` | `0.428` |
| `Overall` | `female` | `0.470` |
| `Early (hotter)` | `fawn` | `0.398` |
| `Early (hotter)` | `male` | `0.531` |
| `Early (hotter)` | `female` | `0.437` |
| `Late (mild-cold)` | `fawn` | `0.344` |
| `Late (mild-cold)` | `male` | `0.511` |
| `Late (mild-cold)` | `female` | `0.443` |

## Interpretation constraint

This is a model-predicted relative-abundance ratio. It is not a direct estimate of birth rate, fawn survival, or realized recruitment. The figure should be treated as primary only if convergence and posterior predictive diagnostics are adequate; otherwise it remains exploratory and the beta-binomial composition model should carry the ratio inference.
