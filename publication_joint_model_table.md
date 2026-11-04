# Joint-model Summary Table

Exploratory latent-factor negative-binomial N-mixture model with shared detection
`Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera + log_days_operational.z`.

## Species-level coefficients

| Class | Term | Posterior mean [95% CrI] | Rhat | ESS |
|---|---|---:|---:|---:|
| `fawn` | `Nature.Reserves.z` | `-0.11 [-0.67, 0.46]` | `1.02` | `395` |
| `male` | `Nature.Reserves.z` | `-0.60 [-1.09, -0.12]` | `1.04` | `626` |
| `female` | `Nature.Reserves.z` | `-0.44 [-1.00, 0.07]` | `1.03` | `427` |
| `fawn` | `Open.Spaces.z` | `0.89 [0.42, 1.37]` | `1.00` | `1070` |
| `male` | `Open.Spaces.z` | `0.49 [0.10, 0.88]` | `1.01` | `1690` |
| `female` | `Open.Spaces.z` | `0.66 [0.26, 1.08]` | `1.01` | `1562` |
| `fawn` | `null.jackal.abundance.z` | `-0.68 [-1.19, -0.18]` | `1.00` | `880` |
| `male` | `null.jackal.abundance.z` | `-0.79 [-1.21, -0.39]` | `1.00` | `1742` |
| `female` | `null.jackal.abundance.z` | `-0.80 [-1.24, -0.37]` | `1.00` | `1362` |
| `fawn` | `Null.Wolf.Abundance.z` | `0.35 [-0.05, 0.78]` | `1.05` | `1048` |
| `male` | `Null.Wolf.Abundance.z` | `0.47 [0.10, 0.84]` | `1.04` | `1701` |
| `female` | `Null.Wolf.Abundance.z` | `0.54 [0.16, 0.94]` | `1.04` | `1745` |

## Negative-binomial overdispersion

| Class | kappa | Rhat | ESS |
|---|---:|---:|---:|
| `fawn` | `2.73 [0.66, 8.89]` | `1.03` | `472` |
| `male` | `6.72 [2.61, 9.85]` | `1.01` | `3759` |
| `female` | `5.46 [1.83, 9.74]` | `1.01` | `1213` |

## Posterior predictive check

| Class | Fit(y) | Fit(yrep) | Bayesian p |
|---|---:|---:|---:|
| `fawn` | `338.18` | `106.55` | `0.025` |
| `male` | `159343.33` | `8214.08` | `0.000` |
| `female` | `301397.86` | `19792.39` | `0.000` |

## Residual cross-class correlations

| Pair | Posterior mean [95% CrI] |
|---|---:|
| `fawn-male` | `0.19 [-0.07, 0.51]` |
| `fawn-female` | `0.26 [-0.02, 0.53]` |
| `male-female` | `0.18 [-0.09, 0.51]` |
