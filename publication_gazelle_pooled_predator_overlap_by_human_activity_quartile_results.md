# Gazelle overlap with pooled predator activity across human-activity quartiles

## Design

The top row represents the 15 cameras in Q1 (`human.abundance <= 0.268`) and the bottom row represents the 15 cameras in Q4 (`human.abundance >= 0.453`). The middle 30 cameras were excluded.

Wolf and jackal records were pooled and then declustered within camera session using a 30-min interval across taxa. This produced `10456` guild events from `10594` species-specific independent events. Jackals comprised `92.1%` of the source events, so the pooled curve is frequency weighted and primarily reflects jackal timing.

Overlap was estimated with `Dhat1` for fawns and `Dhat4` for adults. Q4-Q1 uncertainty used 500 camera-cluster bootstrap replicates, and null differences used 999 camera-label permutations. FDR correction was applied across eligible class comparisons.

## Results

| Class | Estimator | Q1 overlap [95% CI] | Q4 overlap [95% CI] | Q4-Q1 [95% CI] | Permutation P | FDR P |
|---|---|---:|---:|---:|---:|---:|
| Male | `Dhat4` | `0.644 [0.540, 0.753]` | `0.685 [0.610, 0.782]` | `0.041 [-0.084, 0.175]` | `0.726` | `0.964` |
| Female | `Dhat4` | `0.561 [0.410, 0.636]` | `0.556 [0.434, 0.701]` | `-0.005 [-0.124, 0.193]` | `0.964` | `0.964` |
| Fawn | `Dhat1` | `0.359 [NA, NA]` | `0.330 [NA, NA]` | `-0.029 [NA, NA]` | `NA` | `NA` |

Supported Q4-Q1 differences after FDR correction: `0` of `2` eligible comparisons.

The fawn comparison is descriptive only because its smallest stratum contained 13 fawn events.

Overlap describes temporal coincidence at sampled cameras. It does not isolate human activity effects from differences in predator composition, abundance, camera placement, or other site characteristics.

## Outputs

- Figure: `publication_gazelle_pooled_predator_overlap_by_human_activity_quartile.png` and `.pdf`.
- Results: `publication_gazelle_pooled_predator_overlap_by_human_activity_quartile_results.csv`.
- Density curves: `publication_gazelle_pooled_predator_overlap_by_human_activity_quartile_density_curves.csv`.
- Event audit: `publication_gazelle_pooled_predator_overlap_by_human_activity_quartile_event_counts.csv` and `_predator_guild_audit.csv`.
