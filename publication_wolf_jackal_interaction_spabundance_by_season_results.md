# Seasonal spAbundance wolf x jackal interaction models

## Analysis

The survey extended from `2020-08-20` to `2020-12-15` and was divided at its calendar midpoint (`2020-10-17`), independently of gazelle detections. Camera-rounds were assigned by deployment midpoint.

- Family: negative binomial.
- Abundance: `~ wolf_z * jackal_z`.
- Detection: `~ log_days_operational.z + cattle_sightings.z`.
- MCMC: `3` chains, `50000` iterations per chain, `5000` burn-in, thinning `10`.
- Predictions are relative to mean wolf and jackal activity within each period and exclude the abundance intercept.

## Fit support

| Class | Period | Sites | Nonmissing counts | Total detections |
|---|---|---:|---:|---:|
| `Fawn` | `Early` | `60` | `174` | `74` |
| `Fawn` | `Late` | `59` | `170` | `125` |
| `Female` | `Early` | `60` | `174` | `934` |
| `Female` | `Late` | `59` | `170` | `1017` |
| `Male` | `Early` | `60` | `174` | `1135` |
| `Male` | `Late` | `59` | `170` | `1220` |

## Predator coefficients

| Class | Period | Term | Mean [95% CrI] | Rhat | ESS |
|---|---|---|---:|---:|---:|
| `Fawn` | `Early` | `Null.Wolf.Abundance.z` | `0.46 [-0.39, 1.40]` | `1.001` | `6889` |
| `Fawn` | `Early` | `null.jackal.abundance.z` | `-0.40 [-1.49, 0.68]` | `1.001` | `6484` |
| `Fawn` | `Early` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.02 [-1.18, 1.23]` | `1.002` | `5653` |
| `Fawn` | `Late` | `Null.Wolf.Abundance.z` | `0.21 [-0.30, 0.76]` | `1.000` | `4797` |
| `Fawn` | `Late` | `null.jackal.abundance.z` | `-0.89 [-1.62, -0.24]` | `1.002` | `2902` |
| `Fawn` | `Late` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.65 [-0.04, 1.45]` | `1.002` | `3270` |
| `Female` | `Early` | `Null.Wolf.Abundance.z` | `0.68 [0.22, 1.15]` | `1.000` | `3537` |
| `Female` | `Early` | `null.jackal.abundance.z` | `-0.67 [-1.23, -0.14]` | `1.000` | `4821` |
| `Female` | `Early` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.46 [-0.15, 1.12]` | `1.001` | `2386` |
| `Female` | `Late` | `Null.Wolf.Abundance.z` | `0.40 [0.03, 0.80]` | `1.001` | `7444` |
| `Female` | `Late` | `null.jackal.abundance.z` | `-0.70 [-1.12, -0.30]` | `1.002` | `5997` |
| `Female` | `Late` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.51 [0.09, 0.94]` | `1.001` | `5436` |
| `Male` | `Early` | `Null.Wolf.Abundance.z` | `0.31 [-0.09, 0.72]` | `1.001` | `8367` |
| `Male` | `Early` | `null.jackal.abundance.z` | `-0.52 [-0.93, -0.10]` | `1.001` | `8817` |
| `Male` | `Early` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.25 [-0.18, 0.69]` | `1.002` | `7033` |
| `Male` | `Late` | `Null.Wolf.Abundance.z` | `0.34 [0.03, 0.69]` | `1.000` | `8455` |
| `Male` | `Late` | `null.jackal.abundance.z` | `-0.68 [-1.01, -0.36]` | `1.003` | `6411` |
| `Male` | `Late` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.52 [0.19, 0.87]` | `1.003` | `5439` |

## Late-minus-Early coefficient contrasts

These contrasts pair independent posterior draws from the separately fitted period models. They test differences in slopes under an approximation that does not model residual cross-period dependence at the same camera.

| Class | Term | Late - Early [95% CrI] | Pr(difference > 0) |
|---|---|---:|---:|
| `Fawn` | `Null.Wolf.Abundance.z` | `-0.26 [-1.32, 0.77]` | `0.309` |
| `Fawn` | `null.jackal.abundance.z` | `-0.49 [-1.78, 0.79]` | `0.221` |
| `Fawn` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.63 [-0.74, 2.06]` | `0.823` |
| `Female` | `Null.Wolf.Abundance.z` | `-0.27 [-0.88, 0.34]` | `0.188` |
| `Female` | `null.jackal.abundance.z` | `-0.03 [-0.71, 0.67]` | `0.463` |
| `Female` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.05 [-0.72, 0.81]` | `0.550` |
| `Male` | `Null.Wolf.Abundance.z` | `0.04 [-0.48, 0.56]` | `0.558` |
| `Male` | `null.jackal.abundance.z` | `-0.16 [-0.69, 0.36]` | `0.276` |
| `Male` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.27 [-0.28, 0.83]` | `0.834` |

## Convergence summary

Focal slope Rhat range: `1.000-1.003`; ESS range: `2386-8817`.
Intercept Rhat range: `1.107-2.989`. Absolute abundance should not be interpreted when intercept mixing is poor.

## Outputs

- Relative heatmap: `publication_wolf_jackal_interaction_spabundance_by_season_relative_heatmap.png` and `.pdf`.
- Relative curves: `publication_wolf_jackal_interaction_spabundance_by_season_relative_interaction_curves.png` and `.pdf`.
- Conditional wolf effects: `publication_wolf_jackal_interaction_spabundance_by_season_conditional_wolf_effects.png` and `.pdf`.
