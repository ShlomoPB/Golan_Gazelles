# Exploratory wolf:jackal abundance-ratio analysis

## Ratio definition

For each of 60 cameras, the raw covariate was `Null.Wolf.Abundance / null.jackal.abundance`. Both indices were strictly positive, so no pseudocount was used. The raw ratio ranged from `0.0006` to `8.177` (median `1.183`). Analyses used its natural logarithm because the raw ratio was right-skewed; Spearman correlations are unchanged by this monotonic transformation.

A ratio above 1 indicates a wolf-dominated abundance index and a value below 1 a jackal-dominated index. This is a relative index, not a direct measure of predation risk or predator density.

## Primary correlations

| Outcome | Cameras | Spearman rho | 95% camera-bootstrap CI | Raw p | Primary FDR p |
|---|---:|---:|---:|---:|---:|
| Open-space cover | 60 | `0.145` | `[-0.138, 0.409]` | `0.271` | `0.541` |
| Gazelle abundance index | 60 | `0.369` | `[0.131, 0.571]` | `0.004` | `0.015` |
| Fawn share of fawn + female counts | 50 | `0.072` | `[-0.242, 0.368]` | `0.618` | `0.814` |
| All-gazelle daylight fraction | 48 | `-0.035` | `[-0.355, 0.259]` | `0.814` | `0.814` |

## Key interpretation

- The wolf:jackal ratio was not detectably associated with open-space cover (`rho = 0.145`, FDR `p = 0.541`), so the ratio was not simply an open-space gradient in this camera set.
- Cameras with more wolf-dominated predator indices had higher gazelle abundance (`rho = 0.369`, 95% CI `0.131-0.571`, FDR `p = 0.015`). The same pattern was stronger for effort-standardized observed counts.
- The positive count-rate association was clearest for adult females (`rho = 0.414`) and males (`rho = 0.480`), whereas the fawn count-rate trend was weaker.
- Fawn representation (`rho = 0.072`, FDR `p = 0.814`) and all-gazelle daylight use (`rho = -0.035`, FDR `p = 0.814`) showed essentially no ratio association.
- In the model separating predator balance from total predator activity, observed gazelle count rate increased with wolf share (`beta = 0.381`, FDR `p = 0.013`) but not with the total wolf-plus-jackal index (`beta = -0.016`, FDR `p = 0.891`). This is an association with relative predator composition, not evidence that wolves cause higher gazelle abundance.

## Secondary correlation screen

| Domain | Outcome | Cameras | Spearman rho | 95% camera-bootstrap CI | Secondary FDR p |
|---|---|---:|---:|---:|---:|
| Abundance | Observed identified count rate | 60 | `0.480` | `[0.264, 0.653]` | `<0.001` |
| Abundance | Fawn count rate | 60 | `0.284` | `[0.035, 0.500]` | `0.076` |
| Abundance | Female count rate | 60 | `0.414` | `[0.198, 0.603]` | `0.004` |
| Abundance | Male count rate | 60 | `0.480` | `[0.263, 0.651]` | `<0.001` |
| Activity | Fawn daylight fraction | 14 | `-0.033` | `[-0.629, 0.577]` | `0.952` |
| Activity | Female daylight fraction | 42 | `-0.241` | `[-0.535, 0.069]` | `0.272` |
| Activity | Male daylight fraction | 47 | `-0.037` | `[-0.334, 0.267]` | `0.952` |
| Other pressures | Human abundance | 60 | `-0.008` | `[-0.312, 0.286]` | `0.952` |
| Other pressures | Boar abundance | 60 | `-0.118` | `[-0.393, 0.167]` | `0.510` |
| Habitat | Built and roads | 60 | `0.161` | `[-0.121, 0.423]` | `0.402` |
| Habitat | Nature reserves | 60 | `0.119` | `[-0.159, 0.399]` | `0.510` |

## Denominator-aware and partial-model sensitivities

| Model | Cameras | Ratio effect | 95% CI | Odds ratio (quasibinomial only) | FDR p |
|---|---:|---:|---:|---:|---:|
| Quasibinomial fawn composition | 50 | `-0.080` | `[-0.440, 0.279]` | `0.923 [0.644, 1.322]` | `0.663` |
| Quasibinomial daylight use | 48 | `-0.059` | `[-0.226, 0.108]` | `0.943 [0.798, 1.114]` | `0.592` |
| Gazelle abundance adjusted for open space | 60 | `0.116` | `[-0.126, 0.359]` | `not applicable` | `0.527` |
| Observed count rate adjusted for open space | 60 | `0.187` | `[-0.050, 0.424]` | `not applicable` | `0.255` |
| Gazelle abundance using bounded wolf share | 60 | `0.284` | `[0.049, 0.519]` | `not applicable` | `0.064` |
| Observed count rate using bounded wolf share | 60 | `0.382` | `[0.160, 0.604]` | `not applicable` | `0.008` |

## Abundance-signal robustness checks

| Check | Outcome | Excluded camera | Cameras | Partial or Spearman rho | FDR p |
|---|---|---:|---:|---:|---:|
| Partial Spearman controlling Open.Spaces | gazelleabundance | `none` | 60 | `0.343` | `0.008` |
| Partial Spearman controlling Open.Spaces | observed_gazelle_rate_100_days | `none` | 60 | `0.469` | `<0.001` |
| Spearman excluding minimum-ratio camera | gazelleabundance | `57` | 59 | `0.409` | `0.002` |
| Spearman excluding minimum-ratio camera | observed_gazelle_rate_100_days | `57` | 59 | `0.510` | `<0.001` |

## Predator balance versus total predator index

These models include bounded wolf share, the standardized sum of wolf and jackal indices, and open-space cover simultaneously.

| Outcome | Component | Effect per 1 SD | 95% CI | FDR p |
|---|---|---:|---:|---:|
| Fawn composition | Predator balance | `-0.147` | `[-0.602, 0.307]` | `0.781` |
| Fawn composition | Total predator index | `-0.072` | `[-0.456, 0.313]` | `0.819` |
| All-gazelle daylight use | Predator balance | `-0.121` | `[-0.303, 0.061]` | `0.536` |
| All-gazelle daylight use | Total predator index | `0.075` | `[-0.090, 0.239]` | `0.758` |
| Gazelle abundance index | Predator balance | `0.279` | `[0.042, 0.516]` | `0.099` |
| Gazelle abundance index | Total predator index | `-0.067` | `[-0.305, 0.172]` | `0.781` |
| Observed identified count rate | Predator balance | `0.381` | `[0.156, 0.605]` | `0.013` |
| Observed identified count rate | Total predator index | `-0.016` | `[-0.242, 0.210]` | `0.891` |

## Construction diagnostics and cautions

- As expected mathematically, the ratio correlated positively with wolf abundance (`rho = 0.686`) and negatively with jackal abundance (`rho = -0.762`).
- Ratio associations can arise from changes in wolves, jackals, or both. Any promising result should therefore be checked against models containing the two predator indices separately.
- The ratio cannot distinguish low-wolf/low-jackal from high-wolf/high-jackal cameras with the same balance. The balance-plus-load models retain total predator activity as a separate covariate.
- Camera `57` had the minimum ratio because its wolf index was almost zero. Rank correlations were not driven by this camera, but ordinary linear log-ratio models were sensitive to its leverage; bounded wolf share and rank-based checks are more stable summaries.
- The fawn:female outcome is a camera-based composition index, not recruitment or survival. The quasibinomial sensitivity accounts for unequal fawn-plus-female totals across cameras.
- Daylight use is the fraction of 30-min independent events occurring at solar elevation >= 0 degrees. The quasibinomial sensitivity accounts for unequal event totals across cameras.
- These tests do not account for spatial autocorrelation, camera detectability, or other confounding beyond the stated open-space partial models. They are hypothesis-generating rather than confirmatory.

## Outputs

- Enriched site covariates: `SiteCovs_060324_with_wolf_jackal_ratio.csv`
- Analysis-ready camera table: `exploratory_wolf_jackal_ratio_camera_data.csv`
- Correlations: `exploratory_wolf_jackal_ratio_correlations.csv`
- Sensitivity models: `exploratory_wolf_jackal_ratio_sensitivity_models.csv`
- Robustness checks: `exploratory_wolf_jackal_ratio_robustness_checks.csv`
- Predator balance/load models: `exploratory_wolf_jackal_ratio_predator_balance_load_models.csv`
- Primary figure: `exploratory_wolf_jackal_ratio_primary_relationships.png`
- Correlation summary: `exploratory_wolf_jackal_ratio_correlation_summary.png`
