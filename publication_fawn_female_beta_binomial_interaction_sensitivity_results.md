# Wolf x jackal sensitivity for fawn representation

## Model comparison

The primary season-adjusted additive beta-binomial model was compared with an otherwise identical model containing a wolf by jackal interaction. Landscape, human-activity, camera-design, season, and nested polygon/camera adjustments were unchanged.

Adding the interaction changed AIC from `441.08` to `442.72` (delta AIC `1.65`; LR chi-square `0.35`; p `0.552`).

## Interaction-model coefficients

| Term | Odds ratio [95% CI] | p |
|---|---:|---:|
| `wolf_z` | `0.89 [0.64, 1.22]` | `0.452` |
| `jackal_z` | `1.02 [0.65, 1.60]` | `0.941` |
| `seasonLate (mild-cold)` | `1.57 [1.01, 2.44]` | `0.044` |
| `open_z` | `1.58 [1.06, 2.34]` | `0.024` |
| `wolf_z:jackal_z` | `0.86 [0.53, 1.39]` | `0.538` |

## Conditional predator slopes

| Effect | Context | Odds ratio [95% CI] | p |
|---|---|---:|---:|
| Wolf slope | `Jackal z = -1` | `1.03 [0.67, 1.58]` | `0.893` |
| Wolf slope | `Jackal z = 0` | `0.89 [0.64, 1.22]` | `0.452` |
| Wolf slope | `Jackal z = 1` | `0.76 [0.38, 1.52]` | `0.441` |
| Jackal slope | `Wolf z = -1` | `1.18 [0.53, 2.66]` | `0.683` |
| Jackal slope | `Wolf z = 0` | `1.02 [0.65, 1.60]` | `0.941` |
| Jackal slope | `Wolf z = 1` | `0.87 [0.55, 1.40]` | `0.576` |

## Decision

The interaction model is retained as a sensitivity analysis. The additive model remains primary because the interaction did not improve fit, its interval included one, and all conditional predator-slope intervals included one. Open-space cover and the Late versus Early contrast remained supported in the interaction model.
