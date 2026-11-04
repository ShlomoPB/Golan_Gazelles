# Multi-class latent-factor N-mixture model

## Model

- Model type: `spAbundance::lfMsNMix`
- Classes: `fawn`, `male`, `female`
- Family: negative binomial
- Latent factors: `2`
- Detection model: `log_days_operational.z + cattle_sightings.z`
- Shared abundance model: `Nature.Reserves.z + Open.Spaces.z + Built.and.Roads.z + null.jackal.abundance.z + Null.Wolf.Abundance.z`

Files:

- `gazelle_lfMsNMix_2f_NB_shared_env_fit.rds`
- `gazelle_lfMsNMix_2f_NB_shared_env_community_summary.txt`
- `gazelle_lfMsNMix_2f_NB_shared_env_species_summary.txt`
- `gazelle_lfMsNMix_2f_NB_shared_env_lambda_summary.csv`
- `gazelle_lfMsNMix_2f_NB_shared_env_class_correlations.csv`
- Longer-chain update:
- `gazelle_lfMsNMix_2f_NB_shared_env_1200b_25bl_5th_3ch_fit.rds`
- `gazelle_lfMsNMix_2f_NB_shared_env_1200b_25bl_5th_3ch_community_summary.txt`
- `gazelle_lfMsNMix_2f_NB_shared_env_1200b_25bl_5th_3ch_species_summary.txt`

## Cross-class association

Posterior mean latent abundance correlations from the two-factor model:

| Pair | Mean | 95% CrI |
|---|---:|---:|
| `fawn-male` | `0.60` | `[0.37, 0.76]` |
| `fawn-female` | `0.86` | `[0.56, 1.00]` |
| `male-female` | `0.89` | `[0.67, 1.00]` |

These are all positive, which supports the idea that the three classes co-vary strongly across sites even after accounting for the shared environmental covariates in the abundance model.

The longer-chain fit gave essentially the same correlation structure:

| Pair | Mean | 95% CrI |
|---|---:|---:|
| `fawn-male` | `0.60` | `[0.38, 0.76]` |
| `fawn-female` | `0.86` | `[0.53, 1.00]` |
| `male-female` | `0.88` | `[0.65, 1.00]` |

This stability across runs is the strongest evidence from the joint model.

The latent-factor loadings were:

| Loading | Mean | 95% CrI |
|---|---:|---:|
| `male-1` | `0.77` | `[0.39, 1.18]` |
| `female-1` | `1.11` | `[0.65, 1.59]` |
| `female-2` | `0.60` | `[0.07, 1.19]` |

The first factor primarily links all three classes positively, and the second factor adds extra structure for `female`.

## Ecological coefficients

Species-level abundance coefficients remained broadly consistent with the single-class non-spatial NB models:

- `Open.Spaces.z` was clearly positive for `fawn` and also positive for `female`.
- `Built.and.Roads.z` was negative for `male` and `female`, and weak for `fawn`.
- `null.jackal.abundance.z` was negative for all three classes.
- `Null.Wolf.Abundance.z` was positive for `male` and `female`, and weaker for `fawn`.

This means the multispecies model did not overturn the earlier ecological story, but it adds evidence that the classes are positively associated beyond those covariates.

The longer-chain fit improved mixing for the ecological slope terms and the NB overdispersion parameters:

- `kappa` improved to `Rhat` `1.01` for `fawn`, `1.00` for `male`, and `1.00` for `female`.
- Species-level slope terms mostly fell near `Rhat` `1.00` to `1.03`.
- The positive wolf effects for `male` and `female`, the negative jackal effects for all classes, the positive open-space effect for `fawn`, and the negative infrastructure effect for `male` all remained.

## Cautions

- The slope terms mixed reasonably, but the species-specific intercepts still mixed poorly.
- The longer run fixed most of the `kappa` mixing problem, but the species-specific abundance and detection intercepts still mixed poorly.
- This is a useful joint-structure model, but it is still not a true demographic interaction model. The positive correlations indicate shared latent abundance structure, not causal interactions by themselves.
