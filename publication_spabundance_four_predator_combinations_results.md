# Final spAbundance predictions at four predator combinations

Predictions use posterior medians and equal-tail 95% credible intervals from the 13,500 retained beta-coefficient draws in each final saved class-specific model. Values exclude the nonconverged abundance intercept and are relative to the mean-predator reference (wolf_z = 0, jackal_z = 0).

## Predictions

| Class | Predator combination | Median relative abundance | 95% CrI |
|---|---|---:|---:|
| Juvenile | Low wolf / Low jackal | 2.755 | 0.920-10.269 |
| Juvenile | Low wolf / High jackal | 0.216 | 0.046-0.772 |
| Juvenile | High wolf / Low jackal | 1.736 | 0.793-4.331 |
| Juvenile | High wolf / High jackal | 0.929 | 0.324-3.357 |
| Adult female | Low wolf / Low jackal | 1.997 | 1.068-4.039 |
| Adult female | Low wolf / High jackal | 0.168 | 0.077-0.367 |
| Adult female | High wolf / Low jackal | 2.068 | 1.103-4.227 |
| Adult female | High wolf / High jackal | 1.381 | 0.547-4.000 |
| Adult male | Low wolf / Low jackal | 1.897 | 1.090-3.514 |
| Adult male | Low wolf / High jackal | 0.257 | 0.135-0.505 |
| Adult male | High wolf / Low jackal | 1.744 | 1.108-2.939 |
| Adult male | High wolf / High jackal | 1.135 | 0.582-2.458 |

## Posterior ratios

| Class | Contrast (first / second) | Median ratio | 95% CrI | Pr(ratio > 1) |
|---|---|---:|---:|---:|
| Juvenile | Low wolf / Low jackal vs High wolf / High jackal | 2.929 | 0.456-18.670 | 0.882 |
| Juvenile | Low wolf / High jackal vs High wolf / High jackal | 0.230 | 0.021-1.597 | 0.069 |
| Juvenile | Low jackal vs High jackal at low wolf | 12.793 | 1.498-163.912 | 0.989 |
| Juvenile | Low jackal vs High jackal at high wolf | 1.861 | 0.330-10.023 | 0.769 |
| Juvenile | Low wolf vs High wolf at low jackal | 1.566 | 0.365-7.798 | 0.722 |
| Juvenile | Low wolf vs High wolf at high jackal | 0.230 | 0.021-1.597 | 0.069 |
| Adult female | Low wolf / Low jackal vs High wolf / High jackal | 1.452 | 0.356-5.727 | 0.702 |
| Adult female | Low wolf / High jackal vs High wolf / High jackal | 0.124 | 0.024-0.513 | 0.002 |
| Adult female | Low jackal vs High jackal at low wolf | 12.155 | 3.664-38.042 | 1.000 |
| Adult female | Low jackal vs High jackal at high wolf | 1.503 | 0.329-6.484 | 0.706 |
| Adult female | Low wolf vs High wolf at low jackal | 0.969 | 0.362-2.570 | 0.473 |
| Adult female | Low wolf vs High wolf at high jackal | 0.124 | 0.024-0.513 | 0.002 |
| Adult male | Low wolf / Low jackal vs High wolf / High jackal | 1.687 | 0.559-4.863 | 0.827 |
| Adult male | Low wolf / High jackal vs High wolf / High jackal | 0.228 | 0.065-0.705 | 0.004 |
| Adult male | Low jackal vs High jackal at low wolf | 7.474 | 2.610-20.562 | 1.000 |
| Adult male | Low jackal vs High jackal at high wolf | 1.528 | 0.546-4.149 | 0.805 |
| Adult male | Low wolf vs High wolf at low jackal | 1.089 | 0.484-2.483 | 0.580 |
| Adult male | Low wolf vs High wolf at high jackal | 0.228 | 0.065-0.705 | 0.004 |

The second and sixth requested contrasts are mathematically identical; both
are retained under the requested labels.

## Direct test of jackal-association weakening

| Class | Jackal-weakening ratio of ratios | 95% CrI | Pr(>1) |
|---|---:|---:|---:|
| Juvenile | 6.849 | 0.435-184.311 | 0.914 |
| Adult female | 7.971 | 1.171-59.046 | 0.982 |
| Adult male | 4.898 | 1.105-21.667 | 0.982 |

## Factual summary

The highest posterior-median cell was low wolf/low jackal for Juvenile and Adult male, but not for Adult female. No class was directly supported as the global maximum because at least one required 95% pairwise ratio interval included one.
High wolf/high jackal had a higher posterior median than low wolf/high jackal in all three classes. The draw-level ratio supported this direction for Adult female and Adult male but not for Juvenile because its 95% CrI included one.
The negative jackal association was weaker at high wolf activity in all three point estimates. The draw-level ratio-of-ratios 95% CrI directly supported this weakening for adult females and adult males, but not for juveniles.
Cell rankings alone do not establish pairwise differences. Direct support is determined from the draw-level ratio intervals in the contrast table.

Validation: maximum absolute difference from the final figure grid = 4.44e-16.
