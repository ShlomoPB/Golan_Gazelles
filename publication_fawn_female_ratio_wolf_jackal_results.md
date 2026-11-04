# Predicted fawn:female ratios across wolf-jackal contexts

## Definition

Each cell is the ratio of the posterior mean predicted fawn abundance to the posterior mean predicted female abundance from the corresponding `spAbundance` wolf x jackal N-mixture models. Values in parentheses in the figures are the equivalent number of predicted fawns per 100 predicted females.

The Overall panel was calculated from the full seven-round model. It is not an arithmetic average of the two seasonal ratios. Early and Late values come from the independently fitted seasonal models used in the original seasonal heatmap.

These values describe model-based relative abundance, not observed demographic composition, recruitment, survival, or a population-level fawn:female ratio. The saved source grids contain posterior means but not joint posterior prediction draws, so uncertainty intervals cannot be reconstructed without refitting the models.

## Prediction grid

### Overall

| Wolf z | Jackal z | Predicted fawn | Predicted female | Fawn:female | Fawns per 100 females |
|---:|---:|---:|---:|---:|---:|
| 1 | -1 | 70.61 | 628.56 | 0.11 | 11.2 |
| 1 | 0 | 58.82 | 593.29 | 0.10 | 9.9 |
| 1 | 1 | 61.90 | 621.85 | 0.10 | 10.0 |
| 0 | -1 | 90.19 | 664.35 | 0.14 | 13.6 |
| 0 | 0 | 38.76 | 340.57 | 0.11 | 11.4 |
| 0 | 1 | 18.99 | 181.19 | 0.10 | 10.5 |
| -1 | -1 | 135.11 | 737.83 | 0.18 | 18.3 |
| -1 | 0 | 27.68 | 201.99 | 0.14 | 13.7 |
| -1 | 1 | 8.00 | 59.91 | 0.13 | 13.4 |

### Early (hotter)

| Wolf z | Jackal z | Predicted fawn | Predicted female | Fawn:female | Fawns per 100 females |
|---:|---:|---:|---:|---:|---:|
| 1 | -1 | 75.74 | 382.71 | 0.20 | 19.8 |
| 1 | 0 | 34.04 | 319.84 | 0.11 | 10.6 |
| 1 | 1 | 27.55 | 306.13 | 0.09 | 9.0 |
| 0 | -1 | 43.44 | 353.27 | 0.12 | 12.3 |
| 0 | 0 | 19.19 | 174.13 | 0.11 | 11.0 |
| 0 | 1 | 12.25 | 90.55 | 0.14 | 13.5 |
| -1 | -1 | 49.81 | 349.31 | 0.14 | 14.3 |
| -1 | 0 | 13.87 | 99.05 | 0.14 | 14.0 |
| -1 | 1 | 11.99 | 31.98 | 0.37 | 37.5 |

### Late (mild-cold)

| Wolf z | Jackal z | Predicted fawn | Predicted female | Fawn:female | Fawns per 100 females |
|---:|---:|---:|---:|---:|---:|
| 1 | -1 | 43.91 | 425.80 | 0.10 | 10.3 |
| 1 | 0 | 39.54 | 381.89 | 0.10 | 10.4 |
| 1 | 1 | 47.52 | 372.98 | 0.13 | 12.7 |
| 0 | -1 | 58.03 | 447.76 | 0.13 | 13.0 |
| 0 | 0 | 23.28 | 239.85 | 0.10 | 9.7 |
| 0 | 1 | 11.09 | 132.65 | 0.08 | 8.4 |
| -1 | -1 | 89.60 | 492.53 | 0.18 | 18.2 |
| -1 | 0 | 15.18 | 155.00 | 0.10 | 9.8 |
| -1 | 1 | 3.86 | 52.02 | 0.07 | 7.4 |

## Descriptive patterns

Across the full survey, predicted ratios ranged from `0.10` to `0.18` fawns per female.

The largest seasonal point estimate occurred in the Early period at wolf z = `-1` and jackal z = `1`: `0.37` fawns per female (`37.5` per 100 females). Both predicted abundances were low in that cell (fawn `11.99`; female `31.98`), so the high ratio reflects a stronger predicted female decline rather than high fawn abundance.

At that same wolf-jackal combination, the Late-period ratio was `0.07` (`7.4` per 100 females). Because posterior ratio intervals are unavailable, this seasonal contrast is descriptive and should not be presented as statistically supported.

## Outputs

- `publication_fawn_female_ratio_wolf_jackal_combined.png`
- `publication_fawn_female_ratio_wolf_jackal_overall.png`
- `publication_fawn_female_ratio_wolf_jackal_early_hotter.png`
- `publication_fawn_female_ratio_wolf_jackal_late_mild_cold.png`
- `publication_fawn_female_ratio_wolf_jackal_prediction_grid.csv`
