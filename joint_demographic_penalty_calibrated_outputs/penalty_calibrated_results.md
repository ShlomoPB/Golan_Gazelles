# Penalty-calibrated joint demographic diel test

The existing class-specific and original joint analyses were not modified.
This model replaced the numerical-rank Wald omnibus test with common
predator-time smooths plus symmetric sum-to-zero demographic deviation
smooths. Omnibus P values are the single-model penalty-aware smooth-term
tests implemented by mgcv.

## Omnibus tests

| Predator-time component | EDF | Statistic | Penalty-aware P | Original rank-Wald P |
|---|---:|---:|---:|---:|
| Wolf x time | 1.007 | 2.081 | 0.0778 | 1.0000 |
| Jackal x time | 0.000 | 0.000 | 0.7066 | 1.0000 |
| Wolf x jackal x time | 1.575 | 6.110 | 0.0060 | 0.9908 |

## Closed-testing pairwise follow-ups

| Predator-time component | Contrast | EDF | Statistic | Raw P | Holm P |
|---|---|---:|---:|---:|---:|
| Wolf x jackal x time | Juvenile - Adult female | 1.504 | 5.127 | 0.0060 | 0.0121 |
| Wolf x jackal x time | Juvenile - Adult male | 1.630 | 6.724 | 0.0037 | 0.0111 |
| Wolf x jackal x time | Adult female - Adult male | 0.000 | 0.000 | 0.2767 | 0.2767 |

## Diagnostics

The primary model was numerically converged (maximum absolute gradient = 6.90e-05; minimum Hessian eigenvalue = -1.41e-13; numerical PSD = TRUE).
It had 3008 coefficients at full rank, explained 56.4% of deviance, and had Pearson dispersion 0.683. All three pairwise models had maximum absolute gradients below 7.7e-05, numerical PSD Hessians, and no fitting warnings.

## Prespecified sensitivity audit

| Specification | EDF | Statistic | P | Max gradient |
|---|---:|---:|---:|---:|
| Primary apparent solar | 1.575 | 6.110 | 0.0060 | 6.90e-05 |
| k = 8 | 1.571 | 6.526 | 0.0048 | 6.42e-05 |
| k = 16 | 1.580 | 6.060 | 0.0061 | 7.74e-05 |
| Sun time | 1.509 | 5.268 | 0.0095 | 1.14e-04 |
| Alternate camera-only hierarchy | 4.175 | 28.866 | 0.0000 | 4.14e-04 |
| Independent 30-min histories | 1.484 | 5.274 | 0.0097 | 6.75e-05 |
| Sun time + camera-only hierarchy | 3.911 | 29.733 | 0.0000 | 7.19e-05 |

The wolf x jackal x time demographic deviation was supported in 7/7 prespecified specifications.

## Ordered H3 effect-size tests

Response strength was summarized as the RMS magnitude of each class's centered wolf x jackal temporal-interaction function. RMS measures the size of temporal redistribution but not its timing or direction.

| Class | RMS interaction effect | 95% simulation interval |
|---|---:|---:|
| Juvenile | 0.2737 | 0.1272 to 0.5350 |
| Adult female | 0.1113 | 0.0616 to 0.2224 |
| Adult male | 0.1113 | 0.0615 to 0.2228 |

| Directional contrast | Difference | 95% simulation interval | Delta-Wald P | Simulation tail area | Simulation Pr(>0) |
|---|---:|---:|---:|---:|---:|
| Juvenile > Adult female | 0.1623 | -0.0257 to 0.4200 | 0.1037 | 0.0456 | 0.954 |
| Juvenile > Adult male | 0.1623 | -0.0249 to 0.4198 | 0.1037 | 0.0458 | 0.954 |
| Adult female > Adult male | 0.0000 | -0.0018 to 0.0018 | 0.4995 | 0.4910 | 0.509 |
| Juvenile > mean of adult classes | 0.1623 | -0.0253 to 0.4199 | 0.1037 | 0.0456 | 0.954 |

| Ordered hypothesis | Delta-Wald IUT P | Simulation IUT tail area | Joint simulation probability | Delta support | Simulation tail area <0.05 |
|---|---:|---:|---:|---|---|
| Prespecified strict order: Juvenile > Adult female > Adult male | 0.4995 | 0.4910 | 0.484 | FALSE | FALSE |
| Post-hoc simplified order: Juvenile > both adult classes | 0.1037 | 0.0458 | 0.954 | FALSE | TRUE |

The strict juvenile > female > male order was prespecified. The juvenile > both adults hypothesis is a post-hoc simplification and must be reported as such. Intersection-union support requires every component directional contrast to be supported. Because RMS is a nonlinear nonnegative quantity near its boundary, delta-Wald and coefficient-simulation results are both shown; disagreement is treated as sensitivity dependence rather than definitive support.

## Astronomical-daylight direction test

Positive values indicate redistribution toward astronomical daylight. Daylight was averaged across the 5,978 operational camera-days with known dates and coordinates; camera 45 round 2 was excluded from this weighting because its exact outage dates were unavailable.

| Class | Predator component | Daylight shift (percentage points) | 95% simulation interval |
|---|---|---:|---:|
| Juvenile | Wolf gradient | 0.97 | -5.74 to 8.41 |
| Juvenile | Jackal gradient | 0.01 | -0.07 to 0.09 |
| Juvenile | Wolf x jackal interaction | 12.96 | -6.81 to 32.14 |
| Adult female | Wolf gradient | 1.63 | -6.41 to 10.33 |
| Adult female | Jackal gradient | 0.07 | -0.07 to 0.22 |
| Adult female | Wolf x jackal interaction | 13.77 | -0.88 to 28.70 |
| Adult male | Wolf gradient | 1.93 | -7.56 to 12.10 |
| Adult male | Jackal gradient | 0.08 | -0.04 to 0.22 |
| Adult male | Wolf x jackal interaction | 14.52 | -2.12 to 31.43 |

No class-specific daylight shift had a 95% interval excluding zero, and none of the nine direct demographic differences in daylight shift was supported. The push-to-diurnality component of H3 is therefore unsupported in the penalty-calibrated analysis.

## Result

The penalty-aware omnibus test supports demographic heterogeneity in the wolf x jackal temporal interaction (P = 0.006), but not in the wolf-only (P = 0.078) or jackal-only (P = 0.707) temporal responses.
Closed pairwise tests locate the interaction heterogeneity between juveniles and adult females (Holm P = 0.012) and between juveniles and adult males (Holm P = 0.011), not between the two adult classes (Holm P = 0.277).
The juvenile-adult contrast is greatest around 06:24 and 17:24, but no individual simultaneous-band interval excludes zero; inference is therefore global across the 24-hour function rather than hour-specific.
The strict juvenile > female > male magnitude ordering was unsupported because the two adult RMS interaction effects were effectively identical.
Juvenile RMS interaction magnitude exceeded each adult point estimate, but the post-hoc juvenile > both adults test was method-sensitive (simulation IUT tail area = 0.046; delta-Wald IUT P = 0.104) and should be described as suggestive rather than robust ordered support.
Predator activity did not produce a supported increase in astronomical daylight use in any demographic class, and demographic classes did not differ directly in their daylight shifts.

## Interpretation rule

The penalty-aware omnibus tests are primary. Pairwise tests are interpreted
only within an omnibus-supported predator family and use Holm correction.
Difference curves and simultaneous bands are supporting uncertainty displays.
