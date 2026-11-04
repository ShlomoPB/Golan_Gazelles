# Fixed-season sensitivity for the focused fawn:female model

## Rationale

Season has only two deliberately defined levels, so it was evaluated as a fixed Early-Late contrast rather than a random effect. All other fixed effects remained limited to standardized wolf activity, standardized jackal activity, their interaction, and standardized open-space cover. Nested polygon/camera random intercepts were unchanged.

## Fixed-season interaction model

| Term | Odds ratio [95% CI] | p |
|---|---:|---:|
| `wolf_z` | `1.01 [0.67, 1.51]` | `0.973` |
| `jackal_z` | `0.98 [0.60, 1.61]` | `0.942` |
| `open_z` | `1.57 [0.96, 2.57]` | `0.074` |
| `seasonLate (mild-cold)` | `1.51 [1.04, 2.18]` | `0.030` |
| `wolf_z:jackal_z` | `0.85 [0.52, 1.38]` | `0.513` |

## Nested model tests

| Comparison | LR chi-square | Delta AIC (larger - smaller) | p |
|---|---:|---:|---:|
| Add fixed Early-Late season effect | `4.56` | `-2.56` | `0.033` |
| Wolf x jackal interaction after fixed season adjustment | `0.39` | `1.61` | `0.531` |

## Interpretation

This sensitivity analysis tests whether a prespecified temporal shift confounds the focused ecological associations. It should not be selected merely because individual p-values cross 0.05. Effect sizes, uncertainty, the nested season test, and the ecological rationale should be considered together.
