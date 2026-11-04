# Multispecies shared-environment lfMsNMix audit

Fit: `gazelle_lfMsNMix_2f_NB_shared_env_1200b_25bl_5th_3ch_fit.rds`.

## Residual latent-factor correlations

Correlations are derived within each posterior draw from the standardized latent-factor covariance matrix `Lambda %*% t(Lambda)`.

| Pair | Mean | 95% CrI |
|---|---:|---:|
| `fawn-male` | `0.60` | `[0.38, 0.76]` |
| `fawn-female` | `0.86` | `[0.53, 1.00]` |
| `male-female` | `0.88` | `[0.65, 1.00]` |

## Convergence

Free loading Rhat range: `1.002-1.007`; ecological abundance-slope Rhat range: `1.004-1.027`.
Abundance-intercept Rhat range: `1.157-2.133`; detection-intercept Rhat range: `1.182-2.305`. Absolute abundance and intercepts are not used for inference.

## Posterior predictive checks

| Class | Bayesian p |
|---|---:|
| `fawn` | `0.451` |
| `male` | `0.370` |
| `female` | `0.284` |
