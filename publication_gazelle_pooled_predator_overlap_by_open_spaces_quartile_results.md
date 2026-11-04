# Gazelle overlap with pooled predator activity across open-space quartiles

## Design

The top row represents the 15 cameras in Q1 (`Open.Spaces <= 0.652`) and the bottom row represents the 15 cameras in Q4 (`Open.Spaces >= 0.922`). The middle 30 cameras were excluded.

Wolf and jackal records were pooled and then declustered within camera session using a 30-min interval across taxa. This produced `10456` guild events from `10594` species-specific independent events. Jackals comprised `92.1%` of the source events, so the pooled curve is frequency weighted and primarily reflects jackal timing.

Overlap was estimated with `Dhat1` for fawns and `Dhat4` for adults. Q4-Q1 uncertainty used 500 camera-cluster bootstrap replicates, and null differences used 999 camera-label permutations. FDR correction was applied across eligible class comparisons.

## Results

| Class | Estimator | Q1 overlap [95% CI] | Q4 overlap [95% CI] | Q4-Q1 [95% CI] | Permutation P | FDR P |
|---|---|---:|---:|---:|---:|---:|
| Male | `Dhat4` | `0.684 [0.567, 0.773]` | `0.562 [0.488, 0.636]` | `-0.122 [-0.237, 0.016]` | `0.089` | `0.178` |
| Female | `Dhat4` | `0.578 [0.365, 0.741]` | `0.561 [0.373, 0.685]` | `-0.017 [-0.282, 0.202]` | `0.967` | `0.967` |
| Fawn | `Dhat1` | `0.551 [NA, NA]` | `0.353 [NA, NA]` | `-0.198 [NA, NA]` | `NA` | `NA` |

Supported Q4-Q1 differences after FDR correction: `0` of `2` eligible comparisons.

The fawn comparison is descriptive only because Q1 contained five fawn events from two cameras.

Overlap describes temporal coincidence at sampled cameras. It does not isolate habitat effects from differences in predator composition, abundance, camera placement, or other site characteristics.

## Outputs

- Figure: `publication_gazelle_pooled_predator_overlap_by_open_spaces_quartile.png` and `.pdf`.
- Results: `publication_gazelle_pooled_predator_overlap_by_open_spaces_quartile_results.csv`.
- Density curves: `publication_gazelle_pooled_predator_overlap_by_open_spaces_quartile_density_curves.csv`.
- Event audit: `publication_gazelle_pooled_predator_overlap_by_open_spaces_quartile_event_counts.csv` and `_predator_guild_audit.csv`.
