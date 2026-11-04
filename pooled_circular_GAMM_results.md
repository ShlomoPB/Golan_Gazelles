# Pooled species-level circular GAMM results

## Data audit

The full export contained 43,433 mountain-gazelle sighting rows. Collapsing successive records within 30 minutes inside each camera and Camelot session, without demographic splitting, produced exactly 3,688 species-level independent encounters.
All 3,688 encounters were assigned to one of 344 positive-effort camera-rounds after applying the documented session-specific clock corrections. The model used 8,256 camera-round-hour cells and did not sum demographic class-events.

## Requested pooled model

Wolf-associated change in diel timing: EDF = 3.44, approximate P = 0.345 (unsupported).
Jackal-associated change in diel timing: EDF = 1.54, approximate P = 0.044 (supported).
Wolf x jackal temporal interaction: EDF = 1.61, approximate P = 0.339 (unsupported).
The ordinary wolf, jackal, and interaction coefficients control differences in event intensity; the three cyclic varying-coefficient smooths above are the tests of diel timing.

## Predicted normalized curves

| Wolf | Jackal | Peak (solar h) | Proportion 06:00-18:00 | Mean h from solar noon |
|---|---|---:|---:|---:|
| -1 SD | -1 SD | 8.8 | 0.625 | 5.32 |
| -1 SD | +1 SD | 16.7 | 0.629 | 5.29 |
| +1 SD | -1 SD | 16.8 | 0.616 | 5.41 |
| +1 SD | +1 SD | 16.7 | 0.657 | 5.10 |

The four curves are population-level predictions from the requested model at all combinations of -1 and +1 SD wolf and jackal activity. Each curve integrates to one and therefore represents timing rather than abundance.

## Stability and diagnostics

| Model | Effect | EDF | Approximate P |
|---|---|---:|---:|
| Requested k = 8 | Wolf by time | 3.44 | 0.345 |
| Requested k = 8 | Jackal by time | 1.54 | 0.044 |
| Requested k = 8 | Wolf x jackal by time | 1.61 | 0.339 |
| Larger basis k = 12 | Wolf by time | 0.00 | 0.958 |
| Larger basis k = 12 | Jackal by time | 1.48 | 0.575 |
| Larger basis k = 12 | Wolf x jackal by time | 3.30 | 0.929 |
| Camera-round sensitivity | Wolf by time | 4.15 | 0.004 |
| Camera-round sensitivity | Jackal by time | 1.91 | 0.006 |
| Camera-round sensitivity | Wolf x jackal by time | 1.81 | 0.050 |

The requested model converged with theta = 1.324, explained 36.2% of deviance, and had Pearson dispersion = 0.925.
Observed and model-implied zero fractions were 0.748 and 0.747, respectively. Mean within-camera-round lag-1 Pearson residual correlation was 0.226.
The minimum ecological k-index was 0.780 (P = <0.001), and maximum ecological concurvity was 0.217.
The nominal jackal smooth was not retained with k = 12 (P = 0.575). Conversely, the camera-round sensitivity supported the wolf and jackal smooths and placed the interaction at the 0.05 boundary. The ecological k-index remained below one after increasing k. Thus, smooth-term inference was not stable across reasonable specifications, and no pooled predator-related timing effect should be treated as robust.
Point-estimate curves from the camera-round sensitivity were nearly identical to the requested curves (correlations 0.994-1.000), whereas the larger-basis curves were moderately similar (0.883-0.912). Full metrics are in `pooled_circular_GAMM_curve_stability.csv`.

## Comparison with demographic classes

| Analysis | Effect | EDF | Approximate P | Supported |
|---|---|---:|---:|---|
| Pooled species | Wolf by time | 3.44 | 0.345 | No |
| Pooled species | Jackal by time | 1.54 | 0.044 | Yes |
| Pooled species | Wolf x jackal by time | 1.61 | 0.339 | No |
| Juvenile | Wolf by time | 1.02 | 0.086 | No |
| Juvenile | Jackal by time | 0.00 | 0.947 | No |
| Juvenile | Wolf x jackal by time | 2.02 | <0.001 | Yes |
| Adult female | Wolf by time | 3.83 | <0.001 | Yes |
| Adult female | Jackal by time | 0.00 | 0.315 | No |
| Adult female | Wolf x jackal by time | 0.08 | 0.062 | No |
| Adult male | Wolf by time | 0.00 | 0.064 | No |
| Adult male | Jackal by time | 0.00 | 0.643 | No |
| Adult male | Wolf x jackal by time | 3.01 | 0.002 | Yes |

The pooled model should be compared with class-specific smooths by effect identity and curve shape, not by adding class-event sample sizes. A pooled encounter can contain multiple demographic classes, and the pooled dataset also retains gazelles without a focal demographic classification.

## Ecological interpretation

The requested fit suggests a possible pooled shift in diel timing along the jackal gradient, expressed mainly as a change in the balance between morning and late-afternoon activity. That result is exploratory because it was not robust to the larger-basis sensitivity. There was no stable pooled evidence for a wolf-associated timing shift or a wolf by jackal temporal interaction. These site-level associations do not identify responses to individual predator encounters and do not demonstrate causal avoidance.
The nominal pooled jackal effect was unsupported in every class-specific jackal smooth, whereas the retained class-specific models supported a female wolf effect and wolf by jackal interactions for juveniles and males. The pooled result therefore cannot be attributed to one demographic class. Its point curves are likely weighted toward adults because adults dominate classified encounters, but pooling also includes unclassified gazelles and encounters containing multiple classes and can mask opposing class-specific responses.

## Files

- `pooled_circular_GAMM_predicted_curves.png` and `.pdf`
- `pooled_circular_GAMM_predictions.csv`
- `pooled_circular_GAMM_curve_metrics.csv`
- `pooled_circular_GAMM_model_summary.csv`
- `pooled_circular_GAMM_diagnostics.csv`
- `pooled_circular_GAMM_species_level_independent_encounters.csv`
- `pooled_circular_GAMM_models.rds`
