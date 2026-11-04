# Posterior fawn:female ratios across wolf-jackal contexts

## Estimand

Each cell is the posterior distribution of `100 * expected fawn abundance / expected female abundance`, calculated within the same posterior draw from a joint two-factor negative-binomial `lfMsNMix` model. Predictions exclude latent site factors (`include.w = FALSE`) and use continuous standardized predator covariates at `-1`, `0`, and `+1` SD.

Maximum abundance-slope Rhat: `15.58`; minimum abundance-slope ESS: `4`.

## Ratio predictions

| Period | Wolf z | Jackal z | Median fawns/100 females | 95% CrI |
|---|---:|---:|---:|---:|
| `Overall` | `-1` | `-1` | `54.2` | `[10.0, 3284.6]` |
| `Overall` | `0` | `-1` | `36.3` | `[12.3, 504.5]` |
| `Overall` | `1` | `-1` | `25.5` | `[9.1, 195.8]` |
| `Overall` | `-1` | `0` | `56.2` | `[11.3, 574.5]` |
| `Overall` | `0` | `0` | `43.8` | `[12.1, 513.6]` |
| `Overall` | `1` | `0` | `29.5` | `[10.3, 528.3]` |
| `Overall` | `-1` | `1` | `34.6` | `[9.0, 355.2]` |
| `Overall` | `0` | `1` | `45.8` | `[9.5, 566.7]` |
| `Overall` | `1` | `1` | `36.0` | `[6.3, 4286.9]` |
| `Early (hotter)` | `-1` | `-1` | `7.5` | `[2.3, 26.2]` |
| `Early (hotter)` | `0` | `-1` | `10.7` | `[4.5, 31.3]` |
| `Early (hotter)` | `1` | `-1` | `16.2` | `[6.4, 51.1]` |
| `Early (hotter)` | `-1` | `0` | `15.7` | `[4.4, 32.4]` |
| `Early (hotter)` | `0` | `0` | `13.8` | `[5.8, 33.0]` |
| `Early (hotter)` | `1` | `0` | `11.2` | `[6.5, 43.8]` |
| `Early (hotter)` | `-1` | `1` | `26.4` | `[4.5, 114.4]` |
| `Early (hotter)` | `0` | `1` | `17.4` | `[4.8, 55.7]` |
| `Early (hotter)` | `1` | `1` | `10.0` | `[3.1, 59.2]` |
| `Late (mild-cold)` | `-1` | `-1` | `73.6` | `[24.9, 1732.2]` |
| `Late (mild-cold)` | `0` | `-1` | `64.6` | `[26.6, 906.5]` |
| `Late (mild-cold)` | `1` | `-1` | `64.3` | `[22.7, 821.2]` |
| `Late (mild-cold)` | `-1` | `0` | `67.7` | `[36.0, 405.8]` |
| `Late (mild-cold)` | `0` | `0` | `57.8` | `[28.8, 206.5]` |
| `Late (mild-cold)` | `1` | `0` | `53.2` | `[18.9, 192.1]` |
| `Late (mild-cold)` | `-1` | `1` | `58.0` | `[20.2, 173.3]` |
| `Late (mild-cold)` | `0` | `1` | `46.4` | `[18.1, 105.7]` |
| `Late (mild-cold)` | `1` | `1` | `39.8` | `[9.9, 114.4]` |

## Differential predator coefficients

Contrasts are fawn minus female coefficients on the log expected-abundance scale.

| Period | Contrast | Median | 95% CrI | Pr(<0) |
|---|---|---:|---:|---:|
| `Overall` | `wolf` | `-0.16` | `[-0.94, 0.47]` | `0.680` |
| `Overall` | `jackal` | `0.01` | `[-0.80, 0.76]` | `0.493` |
| `Overall` | `wolf_x_jackal` | `0.19` | `[-0.56, 2.33]` | `0.320` |
| `Early (hotter)` | `wolf` | `0.00` | `[-0.62, 0.49]` | `0.493` |
| `Early (hotter)` | `jackal` | `0.25` | `[-0.34, 0.85]` | `0.213` |
| `Early (hotter)` | `wolf_x_jackal` | `-0.39` | `[-1.01, 0.14]` | `0.940` |
| `Late (mild-cold)` | `wolf` | `-0.13` | `[-0.78, 0.23]` | `0.800` |
| `Late (mild-cold)` | `jackal` | `-0.20` | `[-1.75, 0.33]` | `0.747` |
| `Late (mild-cold)` | `wolf_x_jackal` | `-0.06` | `[-0.49, 0.60]` | `0.600` |

## Posterior predictive checks

Site-grouped Freeman-Tukey discrepancies; values near 0.5 indicate better calibration.

| Period | Class | Bayesian p |
|---|---|---:|
| `Overall` | `fawn` | `0.387` |
| `Overall` | `male` | `0.387` |
| `Overall` | `female` | `0.480` |
| `Early (hotter)` | `fawn` | `0.453` |
| `Early (hotter)` | `male` | `0.347` |
| `Early (hotter)` | `female` | `0.427` |
| `Late (mild-cold)` | `fawn` | `0.387` |
| `Late (mild-cold)` | `male` | `0.440` |
| `Late (mild-cold)` | `female` | `0.560` |

## Interpretation constraint

This is a model-predicted relative-abundance ratio. It is not a direct estimate of birth rate, fawn survival, or realized recruitment. The figure should be treated as primary only if convergence and posterior predictive diagnostics are adequate; otherwise it remains exploratory and the beta-binomial composition model should carry the ratio inference.
