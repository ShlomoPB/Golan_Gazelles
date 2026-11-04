# Gazelle-predator diel activity grid

## Analysis

The figure splits the previous class-oriented diel figure into 12 panels: three gazelle classes by four median-split wolf-jackal camera contexts. Gazelle curves are the existing 30-min independent-event curves. Wolf (`Canis lupus`) and jackal (`Canis aureus`) records were extracted from `../full-export_2026-04-15_2006.csv`, restricted to the same documented active camera rounds, collapsed to 30-min independent events within camera session and taxon, and converted from Israel local civil time to apparent solar time.

Each curve is normalized independently to integrate to one. Curve height therefore represents relative timing within a taxon-context combination, not abundance or event rate.

## Predator event audit

| Processing stage | Wolf | Jackal |
|---|---:|---:|
| Raw full-export target records | 5653 | 65030 |
| Target records in active formal rounds | 5548 | 60688 |
| 30-min independent predator events | 836 | 9758 |

## Context event counts

| Gazelle class | Predator context | Cameras | Gazelle events | Wolf events | Jackal events |
|---|---|---:|---:|---:|---:|
| Fawn | Low both | 14 | 67 | 104 | 934 |
| Female | Low both | 14 | 531 | 104 | 934 |
| Male | Low both | 14 | 794 | 104 | 934 |
| Fawn | High jackal / Low wolf | 16 | 3 | 109 | 4613 |
| Female | High jackal / Low wolf | 16 | 56 | 109 | 4613 |
| Male | High jackal / Low wolf | 16 | 151 | 109 | 4613 |
| Fawn | Low jackal / High wolf | 16 | 72 | 305 | 1175 |
| Female | Low jackal / High wolf | 16 | 603 | 305 | 1175 |
| Male | Low jackal / High wolf | 16 | 771 | 305 | 1175 |
| Fawn | High both | 14 | 38 | 318 | 3036 |
| Female | High both | 14 | 276 | 318 | 3036 |
| Male | High both | 14 | 544 | 318 | 3036 |

## Descriptive gazelle-predator temporal overlap

These point estimates quantify curve similarity but are not formal tests of predator-linked activity shifts. `Dhat1` was used when the smaller sample was below 75 events and `Dhat4` otherwise.

| Class | Context | Predator | Events | Estimator | Overlap | Status |
|---|---|---|---:|---|---:|---|
| Fawn | Low both | Wolf | `67 / 104` | `Dhat1` | `0.324` | Descriptive estimate |
| Fawn | Low both | Jackal | `67 / 934` | `Dhat1` | `0.296` | Descriptive estimate |
| Fawn | High jackal / Low wolf | Wolf | `3 / 109` | `Dhat1` | `NA` | Fewer than 5 events |
| Fawn | High jackal / Low wolf | Jackal | `3 / 4613` | `Dhat1` | `NA` | Fewer than 5 events |
| Fawn | Low jackal / High wolf | Wolf | `72 / 305` | `Dhat1` | `0.377` | Descriptive estimate |
| Fawn | Low jackal / High wolf | Jackal | `72 / 1175` | `Dhat1` | `0.342` | Descriptive estimate |
| Fawn | High both | Wolf | `38 / 318` | `Dhat1` | `0.346` | Descriptive estimate |
| Fawn | High both | Jackal | `38 / 3036` | `Dhat1` | `0.304` | Descriptive estimate |
| Female | Low both | Wolf | `531 / 104` | `Dhat4` | `0.482` | Descriptive estimate |
| Female | Low both | Jackal | `531 / 934` | `Dhat4` | `0.462` | Descriptive estimate |
| Female | High jackal / Low wolf | Wolf | `56 / 109` | `Dhat1` | `0.518` | Descriptive estimate |
| Female | High jackal / Low wolf | Jackal | `56 / 4613` | `Dhat1` | `0.545` | Descriptive estimate |
| Female | Low jackal / High wolf | Wolf | `603 / 305` | `Dhat4` | `0.658` | Descriptive estimate |
| Female | Low jackal / High wolf | Jackal | `603 / 1175` | `Dhat4` | `0.625` | Descriptive estimate |
| Female | High both | Wolf | `276 / 318` | `Dhat4` | `0.473` | Descriptive estimate |
| Female | High both | Jackal | `276 / 3036` | `Dhat4` | `0.429` | Descriptive estimate |
| Male | Low both | Wolf | `794 / 104` | `Dhat4` | `0.552` | Descriptive estimate |
| Male | Low both | Jackal | `794 / 934` | `Dhat4` | `0.543` | Descriptive estimate |
| Male | High jackal / Low wolf | Wolf | `151 / 109` | `Dhat4` | `0.707` | Descriptive estimate |
| Male | High jackal / Low wolf | Jackal | `151 / 4613` | `Dhat4` | `0.741` | Descriptive estimate |
| Male | Low jackal / High wolf | Wolf | `771 / 305` | `Dhat4` | `0.630` | Descriptive estimate |
| Male | Low jackal / High wolf | Jackal | `771 / 1175` | `Dhat4` | `0.602` | Descriptive estimate |
| Male | High both | Wolf | `544 / 318` | `Dhat4` | `0.584` | Descriptive estimate |
| Male | High both | Jackal | `544 / 3036` | `Dhat4` | `0.539` | Descriptive estimate |

## Descriptive patterns

- Wolf and jackal activity was predominantly nocturnal in every context, whereas the gazelle curves were more concentrated in daylight. Fawn-predator overlap was consistently low where estimable (`0.296-0.377`).
- Female activity had its greatest predator overlap with `Wolf` in the `Low jackal / High wolf` context (`Dhat = 0.658`).
- Male activity had its greatest predator overlap with `Jackal` in the `High jackal / Low wolf` context (`Dhat = 0.741`). The relatively nocturnal male curve in that context visually aligns with both predator curves more than the male curves in the other contexts.
- The class-specific shifts do not show a single consistent rule in which all gazelle classes move toward or away from predator activity. The strongest visual evidence is therefore context- and class-specific, not a general predator-tracking or predator-avoidance pattern.

## Interpretation constraints

- Similar or opposing peaks are descriptive. Shared responses to daylight, temperature, habitat, prey behavior, and camera placement can produce apparent temporal alignment or separation.
- Predator curves are pooled across the cameras in each context and repeated across the three gazelle-class panels. They do not represent predator activity conditional on a particular gazelle event.
- The high-jackal/low-wolf fawn panel contains only three fawn events and is labeled data deficient rather than displaying a density curve.
- A formal test of linked temporal shifts would require a camera-aware model with predator activity as a time-varying predictor or a prespecified randomization analysis; visual curve displacement alone is not causal evidence.

## Outputs

- Figure: `publication_gazelle_predator_diel_by_context_median.png` and `.pdf`.
- Density curves: `publication_gazelle_predator_diel_by_context_median_density_curves.csv`.
- Panel counts: `publication_gazelle_predator_diel_by_context_median_panel_event_counts.csv`.
- Predator event audit: `publication_gazelle_predator_diel_by_context_median_predator_event_audit.csv`.
- Independent predator events: `publication_gazelle_predator_diel_by_context_median_predator_independent_events.csv`.
- Descriptive overlap: `publication_gazelle_predator_diel_by_context_median_gazelle_predator_overlap.csv`.
