# Camera-bootstrap context contrasts in gazelle-predator overlap

## Methods

For each gazelle class and predator, temporal-overlap coefficients were contrasted among all six pairs of median-split wolf-jackal contexts. Cameras were resampled with replacement within each context, and all gazelle and predator events from a sampled camera were retained together. This preserves within-camera dependence between the two taxa. Each contrast used 999 replicates.

A common overlap estimator was used on both sides of each contrast: `Dhat1` when any of the four contributing samples contained fewer than 75 events and `Dhat4` otherwise. Contrasts were excluded when either taxon had fewer than 5 events in either context.

Approximate two-sided p-values were calculated from the proportion of camera-bootstrap differences on either side of zero, with a +1 correction. P-values were adjusted across the context comparisons within each class-predator family and, secondarily, across all eligible comparisons. A contrast was labeled supported only when its percentile interval excluded zero and its within-family FDR p-value was below 0.05.

## Results

| Class | Predator | Context contrast | Estimator | Overlap | Difference [95% camera-bootstrap CI] | Within-family FDR p | Global FDR p | Status |
|---|---|---|---|---:|---:|---:|---:|---|
| Fawn | Jackal | High jackal / Low wolf minus High both | `Dhat1` | `NA` | `NA` | `NA` | `NA` | Fewer than 5 events in one context |
| Fawn | Jackal | High jackal / Low wolf minus Low jackal / High wolf | `Dhat1` | `NA` | `NA` | `NA` | `NA` | Fewer than 5 events in one context |
| Fawn | Jackal | Low both minus High both | `Dhat1` | `0.296 / 0.304` | `-0.009 [-0.122, 0.102]` | `0.916` | `0.916` | Not supported |
| Fawn | Jackal | Low both minus High jackal / Low wolf | `Dhat1` | `NA` | `NA` | `NA` | `NA` | Fewer than 5 events in one context |
| Fawn | Jackal | Low both minus Low jackal / High wolf | `Dhat1` | `0.296 / 0.342` | `-0.046 [-0.162, 0.050]` | `0.516` | `0.621` | Not supported |
| Fawn | Jackal | Low jackal / High wolf minus High both | `Dhat1` | `0.342 / 0.304` | `0.037 [-0.040, 0.130]` | `0.516` | `0.621` | Not supported |
| Fawn | Wolf | High jackal / Low wolf minus High both | `Dhat1` | `NA` | `NA` | `NA` | `NA` | Fewer than 5 events in one context |
| Fawn | Wolf | High jackal / Low wolf minus Low jackal / High wolf | `Dhat1` | `NA` | `NA` | `NA` | `NA` | Fewer than 5 events in one context |
| Fawn | Wolf | Low both minus High both | `Dhat1` | `0.324 / 0.346` | `-0.022 [-0.148, 0.134]` | `0.862` | `0.916` | Not supported |
| Fawn | Wolf | Low both minus High jackal / Low wolf | `Dhat1` | `NA` | `NA` | `NA` | `NA` | Fewer than 5 events in one context |
| Fawn | Wolf | Low both minus Low jackal / High wolf | `Dhat1` | `0.324 / 0.377` | `-0.053 [-0.194, 0.093]` | `0.834` | `0.714` | Not supported |
| Fawn | Wolf | Low jackal / High wolf minus High both | `Dhat1` | `0.377 / 0.346` | `0.031 [-0.092, 0.189]` | `0.834` | `0.758` | Not supported |
| Female | Jackal | High jackal / Low wolf minus High both | `Dhat1` | `0.545 / 0.432` | `0.113 [-0.044, 0.280]` | `0.352` | `0.510` | Not supported |
| Female | Jackal | High jackal / Low wolf minus Low jackal / High wolf | `Dhat1` | `0.545 / 0.622` | `-0.077 [-0.253, 0.153]` | `0.644` | `0.758` | Not supported |
| Female | Jackal | Low both minus High both | `Dhat4` | `0.462 / 0.429` | `0.033 [-0.093, 0.164]` | `0.644` | `0.812` | Not supported |
| Female | Jackal | Low both minus High jackal / Low wolf | `Dhat1` | `0.462 / 0.545` | `-0.083 [-0.261, 0.097]` | `0.528` | `0.621` | Not supported |
| Female | Jackal | Low both minus Low jackal / High wolf | `Dhat4` | `0.462 / 0.625` | `-0.164 [-0.315, 0.046]` | `0.352` | `0.510` | Not supported |
| Female | Jackal | Low jackal / High wolf minus High both | `Dhat4` | `0.625 / 0.429` | `0.197 [-0.002, 0.320]` | `0.324` | `0.435` | Not supported |
| Female | Wolf | High jackal / Low wolf minus High both | `Dhat1` | `0.518 / 0.474` | `0.045 [-0.170, 0.239]` | `0.908` | `0.916` | Not supported |
| Female | Wolf | High jackal / Low wolf minus Low jackal / High wolf | `Dhat1` | `0.518 / 0.654` | `-0.136 [-0.337, 0.101]` | `0.552` | `0.621` | Not supported |
| Female | Wolf | Low both minus High both | `Dhat4` | `0.482 / 0.473` | `0.010 [-0.175, 0.201]` | `0.908` | `0.916` | Not supported |
| Female | Wolf | Low both minus High jackal / Low wolf | `Dhat1` | `0.483 / 0.518` | `-0.036 [-0.234, 0.183]` | `0.908` | `0.916` | Not supported |
| Female | Wolf | Low both minus Low jackal / High wolf | `Dhat4` | `0.482 / 0.658` | `-0.175 [-0.369, 0.068]` | `0.552` | `0.510` | Not supported |
| Female | Wolf | Low jackal / High wolf minus High both | `Dhat4` | `0.658 / 0.473` | `0.185 [-0.050, 0.376]` | `0.552` | `0.510` | Not supported |
| Male | Jackal | High jackal / Low wolf minus High both | `Dhat4` | `0.741 / 0.539` | `0.201 [0.063, 0.296]` | `0.036` | `0.180` | Supported |
| Male | Jackal | High jackal / Low wolf minus Low jackal / High wolf | `Dhat4` | `0.741 / 0.602` | `0.139 [-0.014, 0.248]` | `0.156` | `0.468` | Not supported |
| Male | Jackal | Low both minus High both | `Dhat4` | `0.543 / 0.539` | `0.004 [-0.091, 0.124]` | `0.870` | `0.916` | Not supported |
| Male | Jackal | Low both minus High jackal / Low wolf | `Dhat4` | `0.543 / 0.741` | `-0.198 [-0.301, -0.046]` | `0.054` | `0.270` | Not supported |
| Male | Jackal | Low both minus Low jackal / High wolf | `Dhat4` | `0.543 / 0.602` | `-0.059 [-0.175, 0.066]` | `0.418` | `0.621` | Not supported |
| Male | Jackal | Low jackal / High wolf minus High both | `Dhat4` | `0.602 / 0.539` | `0.063 [-0.035, 0.183]` | `0.306` | `0.510` | Not supported |
| Male | Wolf | High jackal / Low wolf minus High both | `Dhat4` | `0.707 / 0.584` | `0.123 [-0.046, 0.251]` | `0.380` | `0.510` | Not supported |
| Male | Wolf | High jackal / Low wolf minus Low jackal / High wolf | `Dhat4` | `0.707 / 0.630` | `0.077 [-0.113, 0.204]` | `0.571` | `0.714` | Not supported |
| Male | Wolf | Low both minus High both | `Dhat4` | `0.552 / 0.584` | `-0.033 [-0.151, 0.104]` | `0.650` | `0.812` | Not supported |
| Male | Wolf | Low both minus High jackal / Low wolf | `Dhat4` | `0.552 / 0.707` | `-0.155 [-0.259, 0.005]` | `0.348` | `0.435` | Not supported |
| Male | Wolf | Low both minus Low jackal / High wolf | `Dhat4` | `0.552 / 0.630` | `-0.078 [-0.224, 0.048]` | `0.380` | `0.510` | Not supported |
| Male | Wolf | Low jackal / High wolf minus High both | `Dhat4` | `0.630 / 0.584` | `0.045 [-0.081, 0.205]` | `0.571` | `0.714` | Not supported |

## Multiplicity-aware summary

Supported within-family context contrasts: `1` of `30` eligible comparisons.
- **Male x Jackal:** High jackal / Low wolf minus High both, difference `0.201` (95% CI `0.063-0.296`; within-family FDR `p = 0.036`).

## Interpretation constraints

- These are camera-level contrasts in temporal overlap, not tests of individual behavioral avoidance or attraction.
- Contexts were defined by camera-level predator abundance and were not experimentally assigned. Habitat and other site differences may therefore confound overlap contrasts.
- Bootstrap p-values are approximate. The percentile confidence intervals and effect magnitudes should be emphasized over threshold-based significance.
- The fawn high-jackal/low-wolf sample contained only three events and was excluded from formal contrasts.

## Outputs

- Contrast table: `publication_gazelle_predator_overlap_context_contrasts_results.csv`.
- Figure: `publication_gazelle_predator_overlap_context_contrasts.png` and `.pdf`.
- Reproducible bootstrap draws: `publication_gazelle_predator_overlap_context_contrasts_bootstrap_draws.csv`.
