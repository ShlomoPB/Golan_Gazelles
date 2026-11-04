# Gazelle nocturnality across median-split wolf-jackal contexts

## Methods

Following Procko et al. (2023), we calculated nocturnality for each independent detection as the absolute decimal-hour distance from apparent solar noon (range 0-12 h; 0 = solar noon and 12 = solar midnight). We fitted a separate Gaussian mixed model for each gazelle class, with the four-level wolf-jackal context as a fixed effect and camera as a random intercept. We tested the overall context effect by likelihood-ratio comparison with an intercept-only mixed model and report unadjusted pairwise contrasts from estimated marginal means. Contexts with fewer than 20 events or 5 cameras with detections were excluded from formal inference.

Procko, M., R. Naidoo, V. LeMay, and A. C. Burton. 2023. Human presence and infrastructure impact wildlife nocturnality differently across an assemblage of mammalian species. PLOS ONE 18:e0286131. https://doi.org/10.1371/journal.pone.0286131

## Context summaries

| Class | Predator context | Events | Cameras with events | Mean nocturnality (h) | Median (h) | Formal inference |
|---|---|---:|---:|---:|---:|---:|
| Fawn | Low both | 67 | 12 | 3.55 | 3.14 | Yes |
| Fawn | High jackal / Low wolf | 3 | 2 | 3.73 | 2.59 | No |
| Fawn | Low jackal / High wolf | 72 | 9 | 3.74 | 3.41 | Yes |
| Fawn | High both | 38 | 8 | 3.40 | 3.42 | Yes |
| Female | Low both | 531 | 13 | 4.66 | 4.21 | Yes |
| Female | High jackal / Low wolf | 69 | 10 | 4.49 | 4.57 | Yes |
| Female | Low jackal / High wolf | 601 | 14 | 5.77 | 5.50 | Yes |
| Female | High both | 276 | 13 | 4.50 | 4.24 | Yes |
| Male | Low both | 791 | 13 | 5.30 | 4.83 | Yes |
| Male | High jackal / Low wolf | 163 | 11 | 6.16 | 5.79 | Yes |
| Male | Low jackal / High wolf | 770 | 15 | 5.54 | 5.43 | Yes |
| Male | High both | 544 | 13 | 5.06 | 4.67 | Yes |

## Omnibus context tests

| Class | Likelihood-ratio chi-square | df | P |
|---|---:|---:|---:|
| Fawn | 0.55 | 2 | 0.760 |
| Female | 3.90 | 3 | 0.273 |
| Male | 6.00 | 3 | 0.111 |

## Unadjusted pairwise contrasts

Positive differences indicate that the first context was more nocturnal than the second.

| Class | Contrast | Difference [95% CI] (h) | P |
|---|---|---:|---:|
| Fawn | Low both - Low jackal / High wolf | -0.19 [-0.96, 0.59] | 0.640 |
| Fawn | Low both - High both | 0.15 [-0.78, 1.08] | 0.752 |
| Fawn | Low jackal / High wolf - High both | 0.34 [-0.58, 1.25] | 0.474 |
| Female | Low both - High jackal / Low wolf | -0.30 [-1.61, 1.01] | 0.649 |
| Female | Low both - Low jackal / High wolf | -0.87 [-1.89, 0.16] | 0.097 |
| Female | Low both - High both | 0.11 [-0.95, 1.17] | 0.839 |
| Female | High jackal / Low wolf - Low jackal / High wolf | -0.56 [-1.87, 0.75] | 0.399 |
| Female | High jackal / Low wolf - High both | 0.41 [-0.93, 1.76] | 0.545 |
| Female | Low jackal / High wolf - High both | 0.98 [-0.08, 2.04] | 0.071 |
| Male | Low both - High jackal / Low wolf | -0.65 [-1.67, 0.37] | 0.209 |
| Male | Low both - Low jackal / High wolf | -0.29 [-1.15, 0.56] | 0.505 |
| Male | Low both - High both | 0.56 [-0.33, 1.45] | 0.217 |
| Male | High jackal / Low wolf - Low jackal / High wolf | 0.36 [-0.64, 1.36] | 0.479 |
| Male | High jackal / Low wolf - High both | 1.21 [0.18, 2.25] | 0.021 |
| Male | Low jackal / High wolf - High both | 0.85 [-0.02, 1.73] | 0.056 |

## Pairwise contrasts with raw P < 0.05

- Male: High jackal / Low wolf - High both = 1.21 h (95% CI 0.18 to 2.25; P = 0.021)
Because these are unadjusted post hoc contrasts, they should be interpreted in relation to the class-specific omnibus test rather than as stand-alone evidence.

## Interpretation constraint

Nocturnality condenses each diel distribution to distance from solar noon. It can identify shifts toward night but cannot distinguish morning from evening redistribution. It should therefore be interpreted alongside the circular kernel curves and compareCkern tests. Predator contexts are observational site categories and do not establish predator-caused behavioral change.
