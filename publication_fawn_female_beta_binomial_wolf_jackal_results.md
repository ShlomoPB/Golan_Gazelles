# Beta-binomial fawn:female composition across wolf-jackal context

## Model

Camera-round response: `cbind(fawn, female)`, excluding rounds with zero combined fawn and female detections. The model uses a beta-binomial logit link, continuous standardized wolf and jackal abundance, open-space, built-and-road, reserve, and human-activity controls, camera-design controls, and nested polygon/camera random intercepts.

Informative camera-rounds: `221`; cameras: `50`; independent deployment-midpoint season split: `2020-10-17`.

Predictions are population-level fixed-effect predictions standardized over the observed distribution of path type, camera model, vegetation structure, and area behind the camera. Landscape and human-activity controls are held at their standardized means.

## Model tests

| Comparison | LR chi-square | df | p |
|---|---:|---:|---:|
| `Overall wolf x jackal interaction` | `0.31` | `1` | `0.578` |
| `Common Early-Late shift` | `4.19` | `1` | `0.041` |
| `Common Early-Late shift in additive model` | `4.15` | `1` | `0.042` |
| `Any seasonal structure` | `8.72` | `4` | `0.069` |
| `Seasonal wolf x jackal terms` | `3.91` | `2` | `0.142` |
| `Seasonal change in predator response` | `4.53` | `3` | `0.210` |
| `Seasonal change in additive predator slopes` | `0.97` | `2` | `0.615` |

## Prediction grid

| Period | Wolf z | Jackal z | Fawns per 100 females | 95% interval |
|---|---:|---:|---:|---:|
| `Overall` | `-1` | `-1` | `8.9` | `[5.0, 19.6]` |
| `Overall` | `0` | `-1` | `9.1` | `[6.6, 15.3]` |
| `Overall` | `1` | `-1` | `9.4` | `[6.6, 16.0]` |
| `Overall` | `-1` | `0` | `10.6` | `[7.1, 19.0]` |
| `Overall` | `0` | `0` | `9.5` | `[7.6, 14.4]` |
| `Overall` | `1` | `0` | `8.5` | `[6.5, 13.6]` |
| `Overall` | `-1` | `1` | `12.5` | `[4.7, 40.2]` |
| `Overall` | `0` | `1` | `9.8` | `[5.9, 20.4]` |
| `Overall` | `1` | `1` | `7.7` | `[4.5, 16.7]` |
| `Early (hotter)` | `-1` | `-1` | `3.7` | `[1.4, 12.1]` |
| `Early (hotter)` | `0` | `-1` | `5.8` | `[3.4, 11.9]` |
| `Early (hotter)` | `1` | `-1` | `8.8` | `[5.8, 16.3]` |
| `Early (hotter)` | `-1` | `0` | `8.1` | `[4.5, 17.6]` |
| `Early (hotter)` | `0` | `0` | `7.1` | `[5.1, 12.3]` |
| `Early (hotter)` | `1` | `0` | `6.2` | `[3.9, 11.9]` |
| `Early (hotter)` | `-1` | `1` | `17.2` | `[4.9, 72.7]` |
| `Early (hotter)` | `0` | `1` | `8.7` | `[4.4, 21.5]` |
| `Early (hotter)` | `1` | `1` | `4.3` | `[1.8, 13.2]` |
| `Late (mild-cold)` | `-1` | `-1` | `13.4` | `[7.1, 28.5]` |
| `Late (mild-cold)` | `0` | `-1` | `11.6` | `[8.1, 19.7]` |
| `Late (mild-cold)` | `1` | `-1` | `10.1` | `[6.6, 18.3]` |
| `Late (mild-cold)` | `-1` | `0` | `12.7` | `[7.9, 24.5]` |
| `Late (mild-cold)` | `0` | `0` | `11.3` | `[8.6, 17.8]` |
| `Late (mild-cold)` | `1` | `0` | `10.1` | `[7.4, 16.2]` |
| `Late (mild-cold)` | `-1` | `1` | `12.0` | `[3.8, 43.7]` |
| `Late (mild-cold)` | `0` | `1` | `11.0` | `[5.9, 24.2]` |
| `Late (mild-cold)` | `1` | `1` | `10.1` | `[5.5, 21.7]` |

## Interpretation constraint

This model estimates fawn representation among classified fawn and adult-female camera counts. It accounts for repeated cameras, overdispersion, and measured camera-design covariates, but it does not fully estimate class-specific detection probability. It should be described as an apparent fawn:female composition index rather than realized recruitment, survival, or a population sex-age ratio.
