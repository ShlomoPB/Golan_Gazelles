# Focused beta-binomial fawn:female composition model

## Model

The camera-round response was `cbind(fawn, female)`. The only fixed-effect predictors were standardized wolf activity, standardized jackal activity, their interaction, and standardized open-space cover. Nested polygon/camera random intercepts were retained to account for repeated camera-rounds.

The analysis used `221` informative camera-rounds from `50` cameras in `12` polygons.

## Coefficients

| Term | Odds ratio [95% CI] | p |
|---|---:|---:|
| `wolf_z` | `1.01 [0.67, 1.51]` | `0.970` |
| `jackal_z` | `0.98 [0.60, 1.60]` | `0.934` |
| `open_z` | `1.54 [0.94, 2.53]` | `0.088` |
| `wolf_z:jackal_z` | `0.87 [0.53, 1.41]` | `0.562` |

## Interaction test

Compared with the focused additive model, adding the wolf x jackal interaction changed AIC from `437.22` to `438.91` (LR chi-square `0.31`, p `0.578`).

## Interpretation constraint

With an interaction in the model, the wolf coefficient is the wolf association when standardized jackal activity equals zero, and the jackal coefficient is the jackal association when standardized wolf activity equals zero. The interaction coefficient estimates how each predator slope changes with a 1-SD increase in the other predator. Unsupported intervals should not be interpreted as biological effects.

The response remains an apparent camera-based fawn representation index, not recruitment, survival, or a population-level demographic ratio.
