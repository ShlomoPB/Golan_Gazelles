# Pooled fawn:female beta-binomial analysis

## Analysis

Season was omitted to estimate a single full-survey predator association. The response remained at camera-round resolution rather than summing counts by camera, preserving repeated observations and their denominators. A nested polygon/camera random intercept accounted for repeated cameras, and the beta-binomial distribution accounted for extra-binomial variation.

Both models adjusted for standardized open-space, built-and-road, reserve, and human covariates plus path type, camera model, vegetation structure, and area behind the camera. Continuous predator effects are per 1-SD increase.

## Model comparison

| Comparison | Additive AIC | Interaction AIC | Delta AIC | LR chi-square | p |
|---|---:|---:|---:|---:|---:|
| Wolf x jackal interaction | `443.22` | `444.91` | `1.69` | `0.31` | `0.578` |

Adding the common Early-Late shift improved the additive model (pooled AIC `443.22`; season-adjusted AIC `441.08`; LR chi-square `4.15`, p `0.042`).

## Predator coefficients

| Model | Term | Odds ratio [95% CI] | p |
|---|---|---:|---:|
| Pooled interaction | `wolf_z` | `0.89 [0.65, 1.22]` | `0.471` |
| Pooled interaction | `jackal_z` | `1.04 [0.66, 1.64]` | `0.865` |
| Pooled interaction | `wolf_z:jackal_z` | `0.87 [0.53, 1.41]` | `0.565` |
| Pooled additive | `wolf_z` | `0.93 [0.70, 1.23]` | `0.617` |
| Pooled additive | `jackal_z` | `0.97 [0.65, 1.44]` | `0.874` |

## Interaction-model simple slopes

| Effect | Other predator z | Odds ratio [95% CI] | p |
|---|---:|---:|---:|
| Wolf slope | `-1` | `1.03 [0.67, 1.58]` | `0.906` |
| Wolf slope | `0` | `0.89 [0.65, 1.22]` | `0.471` |
| Wolf slope | `1` | `0.77 [0.39, 1.55]` | `0.466` |
| Jackal slope | `-1` | `1.20 [0.53, 2.70]` | `0.661` |
| Jackal slope | `0` | `1.04 [0.66, 1.64]` | `0.865` |
| Jackal slope | `1` | `0.90 [0.56, 1.44]` | `0.668` |

## Interpretation

The interaction model was not supported over the additive model. Across the full survey, neither predator had a clear association with fawn representation, whether the other predator was low, average, or high. The weak visual tendency toward lower fawn representation when both predators were high is uncertain and should not be interpreted as an additive or synergistic predation effect.

Removing season estimates a weighted average across periods. Because fawn representation increased later in the survey in the season-adjusted analysis, the pooled model answers a narrower spatial question but conceals that temporal shift. It is useful as a sensitivity analysis, not preferable to the season-adjusted model for primary demographic inference.

As in the seasonal analysis, the response is an apparent fawn:female camera composition index. Predator coefficients can reflect differences in fawn detections, female detections, or residual class-specific detectability; they do not directly estimate predation, recruitment, or survival.
