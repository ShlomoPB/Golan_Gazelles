# Juvenile:Adult-female Ratio Analysis

## Descriptive ratios

Across all camera-rounds with finite counts for both classes, the pooled juvenile:adult-female ratio was `0.102`, equivalent to `10.2` juveniles per 100 adult females.

| Round | Juveniles | Adult females | Ratio | Juveniles per 100 females |
|---|---:|---:|---:|---:|
| `X1` | `14` | `198` | `0.071` | `7.1` |
| `X2` | `38` | `434` | `0.088` | `8.8` |
| `X3` | `18` | `263` | `0.068` | `6.8` |
| `X4` | `29` | `379` | `0.077` | `7.7` |
| `X5` | `48` | `408` | `0.118` | `11.8` |
| `X6` | `38` | `230` | `0.165` | `16.5` |
| `X7` | `14` | `39` | `0.359` | `35.9` |

## Beta-binomial GLMM

Response: `cbind(juvenile, female)` at the camera x round level, using only rows with at least one identified juvenile or adult female.

Fixed effects: `Open.Spaces`, `Nature.Reserves`, `Built.and.Roads`, `Null.Wolf.Abundance`, `null.jackal.abundance`, `null.boar.abundance`, `human.abundance`, and `round`.

Random effects: `(1 | Polygon.Name / camera_id)`.

Beta-binomial dispersion parameter (`sigma`) = `69.347`.

| Term | Log-odds | Odds ratio [95% CI] | p |
|---|---:|---:|---:|
| `(Intercept)` | `-2.826` | `0.06 [0.03, 0.13]` | `<0.001` |
| `Open.Spaces.z` | `0.529` | `1.70 [1.02, 2.82]` | `0.042` |
| `Nature.Reserves.z` | `0.221` | `1.25 [0.73, 2.12]` | `0.413` |
| `Built.and.Roads.z` | `0.569` | `1.77 [0.97, 3.23]` | `0.064` |
| `Null.Wolf.Abundance.z` | `-0.049` | `0.95 [0.65, 1.39]` | `0.798` |
| `null.jackal.abundance.z` | `0.125` | `1.13 [0.70, 1.85]` | `0.616` |
| `null.boar.abundance.z` | `-0.190` | `0.83 [0.47, 1.45]` | `0.509` |
| `human.abundance.z` | `-0.129` | `0.88 [0.57, 1.35]` | `0.558` |
| `roundX2` | `0.155` | `1.17 [0.55, 2.46]` | `0.684` |
| `roundX3` | `-0.279` | `0.76 [0.32, 1.78]` | `0.521` |
| `roundX4` | `-0.064` | `0.94 [0.44, 2.01]` | `0.869` |
| `roundX5` | `0.344` | `1.41 [0.68, 2.92]` | `0.353` |
| `roundX6` | `0.587` | `1.80 [0.84, 3.83]` | `0.128` |
| `roundX7` | `0.841` | `2.32 [0.68, 7.93]` | `0.180` |

## Short interpretation

Terms with `p < 0.05` were:
- `Open.Spaces.z` was associated with `higher` juvenile representation among identified juveniles and adult females (odds ratio `1.70`, 95% CI `1.02 to 2.82`, `p = 0.042`).
