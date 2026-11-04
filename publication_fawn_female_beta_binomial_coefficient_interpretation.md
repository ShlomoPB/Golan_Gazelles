# Fawn:female beta-binomial coefficient interpretation

## Response scale

The response is the fawn proportion among classified fawn and adult-female counts. An odds ratio above 1 means greater fawn representation relative to females; it does not necessarily mean greater fawn abundance. Continuous covariates are standardized, so their odds ratios correspond to a 1-SD increase. Categorical effects are comparisons with the stated reference level.

## Model hierarchy

| Comparison | LR chi-square | df | p |
|---|---:|---:|---:|
| Overall wolf x jackal interaction | `0.31` | `1` | `0.578` |
| Common Early-Late shift | `4.19` | `1` | `0.041` |
| Common Early-Late shift in additive model | `4.15` | `1` | `0.042` |
| Any seasonal structure | `8.72` | `4` | `0.069` |
| Seasonal wolf x jackal terms | `3.91` | `2` | `0.142` |
| Seasonal change in predator response | `4.53` | `3` | `0.210` |
| Seasonal change in additive predator slopes | `0.97` | `2` | `0.615` |

## Preferred-model coefficients

The preferred model used additive wolf and jackal slopes plus a common Early-Late shift. Models with a wolf x jackal effect or season-specific predator slopes were not supported.

| Group | Term | Reference | Odds ratio [95% CI] | p |
|---|---|---|---:|---:|
| Predators and season | `wolf_z` | - | `0.93 [0.70, 1.23]` | `0.608` |
| Predators and season | `jackal_z` | - | `0.94 [0.63, 1.41]` | `0.773` |
| Predators and season | `seasonLate (mild-cold)` | Early (hotter) | `1.57 [1.01, 2.44]` | `0.045` |
| Landscape and human activity | `open_z` | - | `1.55 [1.04, 2.32]` | `0.031` |
| Landscape and human activity | `built_z` | - | `1.01 [0.50, 2.04]` | `0.988` |
| Landscape and human activity | `reserve_z` | - | `1.01 [0.64, 1.60]` | `0.967` |
| Landscape and human activity | `human_z` | - | `1.00 [0.71, 1.42]` | `0.990` |
| Path type | `path_typeAnimal path without vehicles` | 4x4 road with off-road vehicles | `1.51 [0.30, 7.71]` | `0.621` |
| Path type | `path_typeDirt path less than 1m wide` | 4x4 road with off-road vehicles | `3.09 [0.64, 14.91]` | `0.159` |
| Path type | `path_typeDirt path more than 1m wide` | 4x4 road with off-road vehicles | `1.19 [0.58, 2.43]` | `0.637` |
| Path type | `path_typePaved road` | 4x4 road with off-road vehicles | `0.88 [0.32, 2.46]` | `0.809` |
| Camera model | `camera_modelBrowning DARK OPS HD PRO X - Long Sensitivity` | Browning DARK OPS HD MAX | `2.68 [0.78, 9.18]` | `0.117` |
| Camera model | `camera_modelBrowning DARK OPS HD PRO X - Normal Sensitivity` | Browning DARK OPS HD MAX | `2.49 [0.83, 7.49]` | `0.105` |
| Camera model | `camera_modelBushnell TROPHY CAM HD AGGRESSOR NO-GLOW` | Browning DARK OPS HD MAX | `1.11 [0.37, 3.30]` | `0.855` |
| Vegetation structure | `vegetation_structureShrubs < 2 m` | Grassland | `11.70 [2.90, 47.25]` | `<0.001` |
| Vegetation structure | `vegetation_structureShrubs > 2 m` | Grassland | `1.95 [0.56, 6.77]` | `0.294` |
| Vegetation structure | `vegetation_structureTree Cover Overhead` | Grassland | `1.39 [0.53, 3.64]` | `0.499` |
| Area behind camera | `area_behind_cameraFence` | Clear | `2.15 [0.64, 7.28]` | `0.217` |
| Area behind camera | `area_behind_cameraImpassable due to topography` | Clear | `2.24 [0.73, 6.86]` | `0.159` |
| Area behind camera | `area_behind_cameraObstructed by vegetation` | Clear | `1.45 [0.60, 3.48]` | `0.408` |

## Conditional predator effects

Each predator slope is the change associated with a 1-SD increase while holding the other predator at the stated value.

| Season | Effect | Other predator z | Odds ratio [95% CI] | p |
|---|---|---:|---:|---:|
| Early (hotter) | Wolf slope | `-1` | `1.57 [0.85, 2.91]` | `0.152` |
| Early (hotter) | Wolf slope | `0` | `0.87 [0.54, 1.40]` | `0.568` |
| Early (hotter) | Wolf slope | `1` | `0.48 [0.18, 1.29]` | `0.146` |
| Early (hotter) | Jackal slope | `-1` | `2.24 [0.72, 6.91]` | `0.162` |
| Early (hotter) | Jackal slope | `0` | `1.24 [0.67, 2.30]` | `0.493` |
| Early (hotter) | Jackal slope | `1` | `0.69 [0.37, 1.27]` | `0.235` |
| Late (mild-cold) | Wolf slope | `-1` | `0.86 [0.55, 1.36]` | `0.515` |
| Late (mild-cold) | Wolf slope | `0` | `0.89 [0.62, 1.26]` | `0.503` |
| Late (mild-cold) | Wolf slope | `1` | `0.91 [0.43, 1.95]` | `0.814` |
| Late (mild-cold) | Jackal slope | `-1` | `0.94 [0.39, 2.30]` | `0.896` |
| Late (mild-cold) | Jackal slope | `0` | `0.97 [0.58, 1.63]` | `0.910` |
| Late (mild-cold) | Jackal slope | `1` | `1.00 [0.60, 1.68]` | `1.000` |

## Factor support

| Variable | Level | Camera-rounds | Cameras | Reference |
|---|---|---:|---:|---|
| season | Early (hotter) | `111` | `47` | Yes |
| season | Late (mild-cold) | `110` | `47` | No |
| path_type | 4x4 Road with off-road vehicles | `136` | `31` | Yes |
| path_type | Animal path without vehicles | `11` | `3` | No |
| path_type | Dirt path less than 1m wide | `6` | `1` | No |
| path_type | Dirt path more than 1m wide | `47` | `11` | No |
| path_type | Paved road | `21` | `5` | No |
| camera_model | Browning DARK OPS HD MAX | `49` | `16` | Yes |
| camera_model | Browning DARK OPS HD PRO X - Long Sensitivity | `43` | `19` | No |
| camera_model | Browning DARK OPS HD PRO X - Normal Sensitivity | `108` | `29` | No |
| camera_model | Bushnell TROPHY CAM HD AGGRESSOR NO-GLOW | `21` | `12` | No |
| vegetation_structure | Grassland | `161` | `35` | Yes |
| vegetation_structure | Shrubs < 2 m | `6` | `2` | No |
| vegetation_structure | Shrubs > 2 m | `28` | `6` | No |
| vegetation_structure | Tree Cover Overhead | `26` | `7` | No |
| area_behind_camera | Clear | `101` | `23` | Yes |
| area_behind_camera | Fence | `23` | `5` | No |
| area_behind_camera | Impassable due to topography | `28` | `7` | No |
| area_behind_camera | Obstructed by vegetation | `69` | `15` | No |

## Interpretation

- Open-space cover was a supported ecological covariate: a 1-SD increase was associated with approximately 55% higher odds of a classified detection being a fawn rather than an adult female.
- The preferred model supported approximately 57% higher odds of fawn representation in the Late period than the Early period, after adjustment for predator, landscape, human, and camera-design variables.
- No wolf or jackal effect was supported in the preferred model. In the full interaction model, their main coefficients are conditional Early-period slopes when the other predator equals its mean.
- The Early wolf x jackal coefficient and its Late-period change were individually suggestive but uncertain, and the global interaction tests were not significant. Simple slopes should therefore be treated as pattern descriptions, not evidence of predator effects.
- Short-shrub vegetation had a large positive coefficient relative to grassland, but only six informative camera-rounds from two cameras had this category. The estimate is sparse and should not be interpreted as a robust habitat effect.
- All other camera-design intervals included an odds ratio of 1, providing no clear evidence of differential fawn-versus-female representation after adjustment.
- The polygon random-intercept variance was effectively zero, while residual camera-level heterogeneity remained.
