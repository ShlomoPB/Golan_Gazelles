# Dynamic predator circular GAMM results

## Analysis scope

The analysis used 8,256 camera-round-hour cells, 3,925 independent gazelle class-events, 851 independent wolf events, and 10,323 independent jackal events.
Unlike the spatial-gradient GAMM, these models used predator rates that changed among camera-round-hour cells. Camera assignments to low or high site-level predator contexts were not used.

## Dynamic event-rate associations

| Class | Window | Effect | Beta | 95% CI | IRR | 95% CI | P |
|---|---|---|---:|---:|---:|---:|---:|
| Juvenile | 1h | Wolf activity | -0.028 | -0.222 to 0.165 | 0.972 | 0.801 to 1.179 | 0.773 |
| Juvenile | 1h | Jackal activity | -0.013 | -0.258 to 0.232 | 0.987 | 0.773 to 1.261 | 0.919 |
| Juvenile | 1h | Wolf x jackal interaction | -0.041 | -0.210 to 0.128 | 0.960 | 0.810 to 1.137 | 0.635 |
| Adult female | 1h | Wolf activity | 0.010 | -0.052 to 0.072 | 1.010 | 0.950 to 1.075 | 0.743 |
| Adult female | 1h | Jackal activity | -0.000 | -0.093 to 0.092 | 1.000 | 0.911 to 1.096 | 0.992 |
| Adult female | 1h | Wolf x jackal interaction | 0.006 | -0.055 to 0.066 | 1.006 | 0.946 to 1.068 | 0.858 |
| Adult male | 1h | Wolf activity | 0.050 | 0.005 to 0.095 | 1.052 | 1.005 to 1.100 | 0.029 |
| Adult male | 1h | Jackal activity | 0.051 | -0.013 to 0.115 | 1.052 | 0.987 to 1.122 | 0.118 |
| Adult male | 1h | Wolf x jackal interaction | -0.011 | -0.055 to 0.033 | 0.989 | 0.946 to 1.033 | 0.612 |
| Juvenile | 3h | Wolf activity | -0.064 | -0.295 to 0.168 | 0.938 | 0.744 to 1.183 | 0.591 |
| Juvenile | 3h | Jackal activity | -0.046 | -0.321 to 0.228 | 0.955 | 0.726 to 1.256 | 0.741 |
| Juvenile | 3h | Wolf x jackal interaction | -0.092 | -0.290 to 0.106 | 0.912 | 0.749 to 1.111 | 0.362 |
| Adult female | 3h | Wolf activity | 0.008 | -0.060 to 0.076 | 1.008 | 0.942 to 1.079 | 0.816 |
| Adult female | 3h | Jackal activity | -0.002 | -0.113 to 0.109 | 0.998 | 0.893 to 1.115 | 0.972 |
| Adult female | 3h | Wolf x jackal interaction | 0.006 | -0.066 to 0.079 | 1.006 | 0.936 to 1.082 | 0.861 |
| Adult male | 3h | Wolf activity | 0.029 | -0.023 to 0.082 | 1.030 | 0.977 to 1.086 | 0.276 |
| Adult male | 3h | Jackal activity | 0.113 | 0.036 to 0.191 | 1.120 | 1.036 to 1.211 | 0.004 |
| Adult male | 3h | Wolf x jackal interaction | -0.023 | -0.073 to 0.026 | 0.977 | 0.930 to 1.027 | 0.360 |

Adult-male coefficients for hourly wolf and jackal activity were positive after accounting for their camera means, round-level deviations, broad diel shape, survey round, exposure, and spatial clustering. For the 1-hour predictors, the wolf IRR was 1.052 (95% CI 1.005-1.100; P 0.029) and the jackal IRR was 1.052 (0.987-1.122; P 0.118).
Adult-male coefficients remained positive in the 3-hour sensitivity (wolf IRR 1.030, P 0.276; jackal IRR 1.120, P 0.004). However, the wolf interval excluded one only for the 1-hour definition, whereas the jackal interval excluded one only for the 3-hour definition; neither predator effect was independently supported in both windows. Adult-female hourly effects were unsupported in both windows.
Juvenile hourly wolf, jackal, and wolf x jackal effects were unsupported in both exposure definitions. The 1-hour interaction estimate was beta = -0.041 (P 0.635), and the 3-hour estimate was beta = -0.092 (P 0.362). Only seven 1-hour cells contained a juvenile event and both predator types, so the juvenile interaction remains weakly informed.

## Diel-shape modification

### Juvenile

Dynamic wolf-by-time: EDF = 0.00, P 0.668; dynamic jackal-by-time: EDF = 0.12, P 0.305; dynamic wolf x jackal by time: EDF = 0.00, P 0.508.

### Adult female

Dynamic wolf-by-time: EDF = 0.00, P 0.799; dynamic jackal-by-time: EDF = 0.00, P 0.820; dynamic wolf x jackal by time: EDF = 0.32, P 0.280.

### Adult male

Dynamic wolf-by-time: EDF = 0.71, P 0.208; dynamic jackal-by-time: EDF = 0.00, P 0.455; dynamic wolf x jackal by time: EDF = 0.00, P 0.424.

No predator-by-time smooth was supported consistently across the 1-hour and 3-hour analyses. The models therefore provided no robust evidence that short-term wolf or jackal activity altered the shape or timing of gazelle diel distributions. Adult males showed a positive but window-dependent event-rate signal, not a shift to a different part of the day.

## Interpretive limits

The hourly predictors are average event rates for a solar-hour bin across each camera round, not exact predator presence on the same date as each gazelle event. Associations therefore measure fine-scale temporal concordance or segregation within camera-rounds, not immediate behavioral responses to individual encounters.
All fits converged and the final baseline smooths used less than their maximum available degrees of freedom after increasing the basis sizes to 12 for juveniles and 16 for adults. However, full-data k-index tests remained low (all P < 0.01) and residual lag-1 correlation remained especially high for juveniles. The dynamic analysis should therefore be treated as exploratory rather than replacing the primary camera-context analysis.
Observed predator rates are estimated with error. The three-hour sensitivity evaluates whether conclusions survive modest temporal smoothing, but a fully joint multivariate model would be required to propagate predator-process uncertainty exactly.
Positive associations indicate temporal concordance and negative associations indicate segregation after controlling broad diel shape, survey round, exposure, and spatial clustering. Neither establishes causal attraction or avoidance.
