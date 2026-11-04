# Unified spatiotemporal gazelle analysis

## Scope and data treatment

- Timestamp interpretation: `local`.
- Formal deployment window: `2020-08-20` to `2020-12-15`.
- Independent early/late split: `2020-10-17`, the midpoint of the deployment window.
- Raw independent-event rows: `3932`.
- Event rows retained within active camera-rounds: `3906`.
- Event rows excluded outside active camera-rounds: `26`.
- Spatial use was modeled as independent-event rate per operational camera-day, not latent animal abundance.
- Juvenile:adult-female results are demographic composition indices, not direct recruitment or survival estimates.
- Temporal overlap intervals used camera-level cluster bootstrap resampling.

## Model comparisons

- Evidence that ecological slopes differ among classes in the spatial event-rate model: likelihood-ratio `p = 0.008`.
- Evidence for a wolf x jackal effect on fawn representation relative to females: likelihood-ratio `p = 0.079`.
- Evidence that daylight-use slopes differ among classes: likelihood-ratio `p = 0.777`.
- NB1 sensitivity model minus NB2 primary-model AIC: `98.3`; positive values favor NB2.

## Model diagnostics

| Model | AIC | Convergence | Positive-definite Hessian | Maximum gradient |
|---|---:|---:|---|---:|
| `space_common_nb` | `3685.5` | `0` | `TRUE` | `0.0007` |
| `space_class_specific_nb` | `3682.6` | `0` | `TRUE` | `0.0014` |
| `space_common_nb1` | `3782.0` | `0` | `TRUE` | `0.0018` |
| `space_class_specific_nb1` | `3780.9` | `0` | `TRUE` | `0.0005` |
| `space_class_specific_zinb` | `NA` | `0` | `FALSE` | `0.0003` |
| `space_encounter_binomial` | `866.5` | `0` | `TRUE` | `0.0001` |
| `space_positive_truncated_nb` | `2853.9` | `0` | `TRUE` | `0.0003` |
| `composition_additive` | `440.9` | `0` | `TRUE` | `0.0001` |
| `composition_wolf_jackal_interaction` | `439.8` | `0` | `TRUE` | `0.0000` |
| `diel_daylight_common` | `1649.4` | `0` | `TRUE` | `0.0000` |
| `diel_daylight_class_specific` | `1667.6` | `0` | `TRUE` | `0.0009` |

## Spatial model simulation checks

| Class | Statistic | Observed | Simulated mean | Simulated 95% interval | Upper-tail probability |
|---|---|---:|---:|---:|---:|
| `Fawn` | `zeros` | `266.00` | `267.36` | `246.97 to 288.02` | `0.565` |
| `Fawn` | `mean` | `0.52` | `0.62` | `0.30 to 1.17` | `0.590` |
| `Fawn` | `variance` | `1.91` | `7.42` | `0.80 to 33.94` | `0.745` |
| `Fawn` | `maximum` | `12.00` | `26.84` | `7.00 to 81.10` | `0.815` |
| `Female` | `zeros` | `128.00` | `142.19` | `120.97 to 167.05` | `0.845` |
| `Female` | `mean` | `4.27` | `5.90` | `3.28 to 10.35` | `0.830` |
| `Female` | `variance` | `86.78` | `477.80` | `37.72 to 1971.44` | `0.880` |
| `Female` | `maximum` | `126.00` | `193.90` | `52.90 to 577.03` | `0.645` |
| `Male` | `zeros` | `90.00` | `99.33` | `77.97 to 121.08` | `0.810` |
| `Male` | `mean` | `6.59` | `9.68` | `4.97 to 18.06` | `0.865` |
| `Male` | `variance` | `86.47` | `1231.74` | `117.14 to 4336.54` | `0.985` |
| `Male` | `maximum` | `62.00` | `294.34` | `80.97 to 870.05` | `1.000` |

## Class-specific spatial-use effects

Effects are changes in log independent-event rate per 1 SD increase, evaluated at mean values of the interacting predator covariate.

| Class | Effect | Estimate | 95% CI | Rate ratio |
|---|---|---:|---:|---:|
| `Fawn` | `wolf_z` | `0.34` | `-0.09 to 0.77` | `1.41` |
| `Fawn` | `jackal_z` | `-0.55` | `-1.04 to -0.06` | `0.58` |
| `Fawn` | `open_z` | `0.92` | `0.42 to 1.42` | `2.51` |
| `Fawn` | `built_z` | `-0.11` | `-0.74 to 0.52` | `0.90` |
| `Fawn` | `human_z` | `-0.13` | `-0.56 to 0.29` | `0.87` |
| `Fawn` | `wolf_z:jackal_z` | `0.25` | `-0.24 to 0.74` | `1.29` |
| `Female` | `wolf_z` | `0.60` | `0.24 to 0.97` | `1.83` |
| `Female` | `jackal_z` | `-0.81` | `-1.20 to -0.41` | `0.45` |
| `Female` | `open_z` | `0.60` | `0.19 to 1.00` | `1.82` |
| `Female` | `built_z` | `-0.69` | `-1.25 to -0.14` | `0.50` |
| `Female` | `human_z` | `-0.14` | `-0.48 to 0.20` | `0.87` |
| `Female` | `wolf_z:jackal_z` | `0.59` | `0.18 to 0.99` | `1.80` |
| `Male` | `wolf_z` | `0.46` | `0.11 to 0.82` | `1.59` |
| `Male` | `jackal_z` | `-0.78` | `-1.16 to -0.41` | `0.46` |
| `Male` | `open_z` | `0.34` | `-0.05 to 0.74` | `1.41` |
| `Male` | `built_z` | `-0.68` | `-1.22 to -0.14` | `0.51` |
| `Male` | `human_z` | `-0.03` | `-0.35 to 0.29` | `0.97` |
| `Male` | `wolf_z:jackal_z` | `0.48` | `0.09 to 0.86` | `1.61` |

## Direct tests of class differences in spatial slopes

| Comparison | Effect | Difference | 95% CI |
|---|---|---:|---:|
| `Fawn - Male` | `open_z` | `0.57` | `0.24 to 0.91` |
| `Female - Male` | `open_z` | `0.25` | `0.07 to 0.43` |
| `Fawn - Female` | `built_z` | `0.58` | `0.17 to 0.99` |
| `Fawn - Male` | `built_z` | `0.57` | `0.17 to 0.97` |

## Fawn:adult-female composition effects

| Effect | Log-odds estimate | 95% CI | Odds ratio |
|---|---:|---:|---:|
| `wolf_z` | `-0.18` | `-0.56 to 0.19` | `0.83` |
| `jackal_z` | `0.30` | `-0.17 to 0.77` | `1.34` |
| `open_z` | `0.59` | `0.12 to 1.06` | `1.80` |
| `built_z` | `0.80` | `0.20 to 1.41` | `2.23` |
| `human_z` | `-0.14` | `-0.54 to 0.27` | `0.87` |
| `reserve_z` | `0.07` | `-0.40 to 0.53` | `1.07` |
| `wolf_z:jackal_z` | `-0.44` | `-0.90 to 0.02` | `0.64` |

## Early-to-late changes in daylight use

| Class | Log-odds change | 95% CI | Odds ratio |
|---|---:|---:|---:|
| `Fawn` | `-0.65` | `-1.68 to 0.38` | `0.52` |
| `Female` | `0.47` | `0.15 to 0.79` | `1.60` |
| `Male` | `0.24` | `-0.01 to 0.49` | `1.27` |

## Camera-clustered temporal overlap

| Period | Comparison | Events | Cameras | Overlap [95% cluster-bootstrap CI] |
|---|---|---:|---:|---:|
| `All` | `Fawn vs Female` | `180 / 1466` | `31 / 50` | `0.789 [0.708, 0.856]` |
| `All` | `Fawn vs Male` | `180 / 2260` | `31 / 51` | `0.743 [0.663, 0.780]` |
| `All` | `Female vs Male` | `1466 / 2260` | `50 / 51` | `0.948 [0.855, 0.963]` |
| `Early` | `Fawn vs Female` | `64 / 708` | `18 / 47` | `0.712 [0.603, 0.800]` |
| `Early` | `Fawn vs Male` | `64 / 1076` | `18 / 50` | `0.720 [0.591, 0.762]` |
| `Early` | `Female vs Male` | `708 / 1076` | `47 / 50` | `0.954 [0.864, 0.952]` |
| `Late` | `Fawn vs Female` | `116 / 758` | `28 / 46` | `0.851 [0.748, 0.881]` |
| `Late` | `Fawn vs Male` | `116 / 1184` | `28 / 49` | `0.763 [0.649, 0.811]` |
| `Late` | `Female vs Male` | `758 / 1184` | `46 / 49` | `0.897 [0.833, 0.927]` |

## Interpretation guardrails

- Predator associations remain observational and cannot establish wolf-mediated suppression of jackals without temporally resolved predator data or a causal design.
- The zero-inflated spatial model was retained only as a sensitivity fit and should be rejected if its Hessian is not positive definite.
- NB1 reproduced the main coefficient directions and better matched some extreme-count statistics, but NB2 was retained when its information-criterion support was substantially stronger.
- The daylight model is conditional on successful demographic classification. Unclassified gazelle detections are needed to test whether fawns are harder to identify in infrared imagery.
- If camera timestamps were local rather than UTC, rerun with `--timestamp_mode=local` before interpreting solar time or daylight use.
- The early/late comparison is a survey-period contrast. It should not be described as a temperature effect until external or field temperature data are joined.
