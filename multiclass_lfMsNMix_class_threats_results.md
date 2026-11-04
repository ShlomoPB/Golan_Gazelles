# Multi-class lfMsNMix class-threat model

## Model

- Model type: `spAbundance::lfMsNMix`
- Family: negative binomial
- Latent factors: `2`
- Shared abundance formula used in the fitted model:
  `~ Nature.Reserves.z + Open.Spaces.z + null.boar.abundance.z + null.jackal.abundance.z + Null.Wolf.Abundance.z`
- Shared detection formula in the selected final run:
  `~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera + log_days_operational.z`

Note: `lfMsNMix` does not support truly different abundance formulas by class. The requested class-specific threat model was approximated by fitting the union of the requested covariates and then interpreting the class-specific posterior coefficients.

## Detection Comparison

The 2000-batch comparison favored keeping `log_days_operational.z` in the detection model. Relative to the design-only detection model, the `design_effort` fit improved community-level abundance mixing and female `kappa` mixing, while leaving the main class-specific slope directions unchanged.

## Final Class-specific Threat Coefficients

| Class | Term | Posterior mean | 95% CrI | Rhat | ESS |
|---|---|---:|---:|---:|---:|
| `fawn` | `Nature.Reserves.z` | `-0.11` | `-0.11 [-0.67, 0.46]` | `1.02` | `395` |
| `male` | `Nature.Reserves.z` | `-0.60` | `-0.60 [-1.09, -0.12]` | `1.04` | `626` |
| `female` | `Nature.Reserves.z` | `-0.44` | `-0.44 [-1.00, 0.07]` | `1.03` | `427` |
| `fawn` | `Open.Spaces.z` | `0.89` | `0.89 [0.42, 1.37]` | `1.00` | `1070` |
| `male` | `Open.Spaces.z` | `0.49` | `0.49 [0.10, 0.88]` | `1.01` | `1690` |
| `female` | `Open.Spaces.z` | `0.66` | `0.66 [0.26, 1.08]` | `1.01` | `1562` |
| `fawn` | `null.boar.abundance.z` | `0.31` | `0.31 [-0.18, 0.79]` | `1.02` | `914` |
| `male` | `null.boar.abundance.z` | `0.39` | `0.39 [-0.01, 0.79]` | `1.00` | `1817` |
| `female` | `null.boar.abundance.z` | `0.37` | `0.37 [-0.06, 0.80]` | `1.01` | `1280` |
| `fawn` | `null.jackal.abundance.z` | `-0.68` | `-0.68 [-1.19, -0.18]` | `1.00` | `880` |
| `male` | `null.jackal.abundance.z` | `-0.79` | `-0.79 [-1.21, -0.39]` | `1.00` | `1742` |
| `female` | `null.jackal.abundance.z` | `-0.80` | `-0.80 [-1.24, -0.37]` | `1.00` | `1362` |
| `fawn` | `Null.Wolf.Abundance.z` | `0.35` | `0.35 [-0.05, 0.78]` | `1.05` | `1048` |
| `male` | `Null.Wolf.Abundance.z` | `0.47` | `0.47 [0.10, 0.84]` | `1.04` | `1701` |
| `female` | `Null.Wolf.Abundance.z` | `0.54` | `0.54 [0.16, 0.94]` | `1.04` | `1745` |

## Kappa

| Model | Class | Mean | 95% CrI | Rhat | ESS |
|---|---|---:|---:|---:|---:|
| `baseline` | `fawn` | `2.50` | `2.50 [0.55, 8.55]` | `1.01` | `246` |
| `baseline` | `male` | `6.63` | `6.63 [2.40, 9.86]` | `1.00` | `900` |
| `baseline` | `female` | `3.71` | `3.71 [1.35, 9.23]` | `1.00` | `404` |
| `design_effort` | `fawn` | `2.68` | `2.68 [0.63, 8.87]` | `1.03` | `383` |
| `design_effort` | `male` | `6.71` | `6.71 [2.61, 9.86]` | `1.01` | `2106` |
| `design_effort` | `female` | `5.26` | `5.26 [1.73, 9.70]` | `1.00` | `652` |
| `final` | `fawn` | `2.73` | `2.73 [0.66, 8.89]` | `1.03` | `472` |
| `final` | `male` | `6.72` | `6.72 [2.61, 9.85]` | `1.01` | `3759` |
| `final` | `female` | `5.46` | `5.46 [1.83, 9.74]` | `1.01` | `1213` |

## Posterior Predictive Check

Manual site-grouped Freeman-Tukey PPC using every 100th posterior draw.

| Class | Fit(y) | Fit(yrep) | Bayesian p |
|---|---:|---:|---:|
| `fawn` | `338.18` | `106.55` | `0.025` |
| `male` | `159343.33` | `8214.08` | `0.000` |
| `female` | `301397.86` | `19792.39` | `0.000` |

## Latent Mean-abundance Correlation

These are posterior correlations among site-level log latent mean abundances (`log(mu)`) for each class pair. This is a directly comparable cross-class association metric across model variants.

| Model | Pair | Mean | 95% CrI |
|---|---|---:|---:|
| `baseline` | `fawn-male` | `0.24` | `0.24 [0.02, 0.46]` |
| `baseline` | `fawn-female` | `0.27` | `0.27 [0.07, 0.47]` |
| `baseline` | `male-female` | `0.22` | `0.22 [0.01, 0.43]` |
| `design_effort` | `fawn-male` | `0.39` | `0.39 [0.18, 0.59]` |
| `design_effort` | `fawn-female` | `0.38` | `0.38 [0.16, 0.60]` |
| `design_effort` | `male-female` | `0.34` | `0.34 [0.13, 0.55]` |
| `final` | `fawn-male` | `0.19` | `0.19 [-0.07, 0.51]` |
| `final` | `fawn-female` | `0.26` | `0.26 [-0.02, 0.53]` |
| `final` | `male-female` | `0.18` | `0.18 [-0.09, 0.51]` |
