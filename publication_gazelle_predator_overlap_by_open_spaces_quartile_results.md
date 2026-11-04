# Gazelle-predator diel overlap across open-space quartiles

## Design and inference

Cameras were assigned to the bottom (`Q1 <= 0.652`) or top (`Q4 >= 0.922`) quartile of the site-level `Open.Spaces` cover covariate. The middle 50% of cameras was excluded. Q1 contained `15` cameras and Q4 contained `15` cameras.

For each gazelle class-predator pair, overlap was estimated separately in Q1 and Q4. A common estimator was used for both strata: `Dhat1` when the smallest of the four contributing event samples contained fewer than 75 events and `Dhat4` otherwise. Formal inference required at least 20 events in every contributing sample.

Uncertainty in each overlap and in `Q4 - Q1` was estimated with 500 camera-cluster bootstrap replicates, retaining all events from each resampled camera. The null distribution for the quartile difference was obtained from 999 permutations of the open-space-group labels among the 30 extreme cameras while preserving the 15/15 group sizes. Permutation p-values were adjusted across the eligible tests with the Benjamini-Hochberg false-discovery-rate procedure.

## Event counts

| Class | Predator | Open-space group | Gazelle events | Predator events | Gazelle cameras | Predator cameras |
|---|---|---|---:|---:|---:|---:|
| Fawn | Wolf | Q1: low open-space cover | 5 | 233 | 2 | 14 |
| Fawn | Wolf | Q4: high open-space cover | 90 | 227 | 11 | 15 |
| Fawn | Jackal | Q1: low open-space cover | 5 | 3372 | 2 | 14 |
| Fawn | Jackal | Q4: high open-space cover | 90 | 1627 | 11 | 15 |
| Female | Wolf | Q1: low open-space cover | 149 | 233 | 11 | 14 |
| Female | Wolf | Q4: high open-space cover | 731 | 227 | 15 | 15 |
| Female | Jackal | Q1: low open-space cover | 149 | 3372 | 11 | 14 |
| Female | Jackal | Q4: high open-space cover | 731 | 1627 | 15 | 15 |
| Male | Wolf | Q1: low open-space cover | 354 | 233 | 13 | 14 |
| Male | Wolf | Q4: high open-space cover | 910 | 227 | 15 | 15 |
| Male | Jackal | Q1: low open-space cover | 354 | 3372 | 13 | 14 |
| Male | Jackal | Q4: high open-space cover | 910 | 1627 | 15 | 15 |

## Results

| Class | Predator | Estimator | Q1 overlap [95% CI] | Q4 overlap [95% CI] | Q4-Q1 [95% CI] | Permutation P | FDR P | Direction |
|---|---|---|---:|---:|---:|---:|---:|---|
| Fawn | Jackal | `Dhat1` | `0.547 [NA, NA]` | `0.351 [NA, NA]` | `-0.196 [NA, NA]` | `NA` | `NA` | Data deficient |
| Female | Jackal | `Dhat4` | `0.577 [0.363, 0.744]` | `0.559 [0.374, 0.678]` | `-0.018 [-0.279, 0.203]` | `0.960` | `0.981` | Lower overlap at high-open-space cameras |
| Male | Jackal | `Dhat4` | `0.684 [0.568, 0.772]` | `0.560 [0.486, 0.638]` | `-0.123 [-0.235, 0.017]` | `0.089` | `0.356` | Lower overlap at high-open-space cameras |
| Fawn | Wolf | `Dhat1` | `0.559 [NA, NA]` | `0.369 [NA, NA]` | `-0.190 [NA, NA]` | `NA` | `NA` | Data deficient |
| Female | Wolf | `Dhat4` | `0.559 [0.309, 0.667]` | `0.575 [0.340, 0.712]` | `0.016 [-0.261, 0.280]` | `0.981` | `0.981` | Greater overlap at high-open-space cameras |
| Male | Wolf | `Dhat4` | `0.657 [0.505, 0.767]` | `0.580 [0.472, 0.653]` | `-0.077 [-0.226, 0.101]` | `0.390` | `0.780` | Lower overlap at high-open-space cameras |

## Publication assessment

Supported high-versus-low open-space differences after FDR correction: `0` of `4` eligible tests.

The analysis does not provide multiplicity-adjusted evidence that gazelle-predator temporal overlap differs between open-space extremes. The figure may be useful descriptively, but adding it as a main inferential result would overstate the evidence.

Fawn contrasts were retained as descriptive curves only because Q1 contained five fawn events from two cameras; they were excluded from bootstrap and permutation inference.

The quartile comparison discards the middle 50% of cameras and does not isolate open-space cover from predator activity, camera placement, or other site differences. Overlap coefficients quantify similarity in detection timing and do not demonstrate behavioral avoidance or attraction.

Overlap coefficients follow Ridout and Linkie (2009) and were calculated with the `overlap` R package.

## Outputs

- Main 2 x 3 figure: `publication_gazelle_predator_overlap_by_open_spaces_quartile.png` and `.pdf`.
- Results: `publication_gazelle_predator_overlap_by_open_spaces_quartile_results.csv`.
- Event and camera audit: `publication_gazelle_predator_overlap_by_open_spaces_quartile_event_counts.csv`, `_camera_summary.csv`, and `_camera_groups.csv`.
- Reproducibility draws: `publication_gazelle_predator_overlap_by_open_spaces_quartile_bootstrap_draws.csv` and `_permutation_draws.csv`.
