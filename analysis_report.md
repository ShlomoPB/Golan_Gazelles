# Multispecies N-Mixture Analysis of Gazelle Abundance in Golan Heights

## Overview
This analysis examines variation in the abundance of male, female, and juvenile gazelles across 60 camera trap sites in the Golan Heights using Bayesian multispecies N-mixture models. The study evaluates associations with predator activity, focusing on wolf pack abundance as the primary wolf predictor, alongside jackals, boars, and environmental covariates (human presence, nature reserves, minefields, infrastructure, firing zones, water proximity).

## Data
- **Abundance data**: Counts from 60 cameras over 7 survey rounds (male, female, juvenile gazelles)
- **Site covariates**: Predator abundances, human activity, land use types, distance to water
- **Detection covariates**: Days operational, cattle sightings
- **Total observations**: 1,260 site-round combinations per species

## Model Structure
- **Distribution**: Poisson (primary) and Negative Binomial (comparison)
- **Abundance model**: Community-level random effects with species-specific intercepts
- **Detection model**: Species-specific logit-linear model
- **MCMC**: 12,000 posterior samples (3 chains x 4,000 samples)

## Key Results

### Community-Level Effects (Log Abundance Scale)
| Covariate | Mean | 95% CrI | Interpretation |
|-----------|------|---------|----------------|
| Minefields | -0.63 | (-1.84, 0.53) | Negative posterior mean, but the interval overlaps 0; no clear support for an effect |
| Proximity to Water | 0.28 | (-0.72, 1.30) | Positive posterior mean, but the interval overlaps 0; effect remains uncertain |
| Jackals | -0.38 | (-1.07, 0.34) | Negative posterior mean, but the interval overlaps 0; weak evidence only |
| Wolf Packs | 0.05 | (-0.73, 0.81) | Posterior centered near 0 with no clear support for an effect |
| Human Abundance | -0.15 | (-0.87, 0.55) | Negative posterior mean, but the interval overlaps 0; effect remains uncertain |

Across community-level abundance covariates, posterior means suggest several ecologically plausible directional patterns, but all reported 95% credible intervals overlap 0. These results are therefore better interpreted as weak or uncertain evidence than as strong support for directional effects.

### Species-Specific Patterns
- **Males**: Highest baseline abundance (intercept 3.35); posterior means suggest avoidance of built areas (-0.58) and a positive association with wolf packs (0.21), but these should only be described as species-level tendencies unless their credible intervals exclude 0.
- **Females**: Moderate abundance (3.08); posterior means indicate the strongest positive water association (0.55) and a negative human effect (-0.20), but both should be treated as uncertain unless supported by intervals excluding 0.
- **Juveniles**: Lowest abundance (1.07); posterior means indicate the strongest minefield avoidance (-1.07) and positive water (0.38) and cattle detection effects (0.18), but interval overlap with 0 should be checked before making directional claims.

These species-specific contrasts are useful for generating ecological hypotheses, but they should not be framed as definitive demographic differences unless the corresponding posterior intervals exclude 0 and model fit is adequate.

### Wolf Covariate Focus
- The current model uses `Wolf.pack.abundance.z` as the wolf predictor. The raw site data also contains `Null.Wolf.Abundance`, which is now included in the prepared covariates and can be tested directly as `Null.Wolf.Abundance.z`.
- Using `Null.Wolf.Abundance.z` may provide an alternate wolf signal with different ecological meaning than `Wolf.pack.abundance.z`.
- The wolf effect remains species-specific in posterior mean direction, but current intervals indicate substantial uncertainty; this is better described as suggestive heterogeneity rather than strong evidence for contrasting wolf responses among demographic groups.

### Correlated-Species Modeling with `lfMsNMix`
- `spAbundance` includes `lfMsNMix`, which is designed for latent factor multispecies N-mixture models with residual species correlation.
- The new workflow fits an abundance-based joint species distribution model with imperfect detection, using `Null.Wolf.Abundance.z` as the wolf covariate.
- `lfMsNMix` requires site coordinates in the data object for the latent factor spatial structure, and the prepared data now includes coordinates from `Camera Coordinates.xlsx`.
- `lfMsNMix` models species correlations via latent factors (`n.factors`), enabling residual association between male, female, and juvenile gazelle abundances after accounting for shared covariates.
- This is the appropriate next step for a correlated multispecies abundance model in `spAbundance`.

### Detection Probabilities
- Baseline detection: Males 0.20, Females 0.17, Juveniles 0.13
- Days operational: Posterior estimates suggest a positive effect on detection across species
- Cattle sightings: Posterior estimates suggest higher juvenile detectability

## Model Diagnostics
- **Convergence**: Excellent (Rhat < 1.05, ESS > 700)
- **Fit**: The Poisson model shows severe lack of fit (posterior predictive check p-values = 0), whereas the Negative Binomial model shows substantially improved fit (posterior predictive check p approximately 0.5)
- **Model Comparison**: The Negative Binomial model is decisively preferred based on posterior predictive checks (see `WAIC_COMPARISON_REPORT.md`)

Because substantive interpretation depends on adequate model fit, ecological inference should be based primarily on the Negative Binomial formulation rather than the poorly fitting Poisson model.

## Visualizations
- Effect plots for key covariates (wolf packs, jackals, minefields, water proximity)
- Spatial abundance maps showing predicted distributions

## Recommendations
1. **Adopt the Negative Binomial distribution** for primary inference and reporting because it provides substantially better fit (posterior predictive check p-values: 0.55 community, 0.48-0.65 species)
2. Consider spatial autocorrelation models for refined predictions
3. Investigate species x predator interactions
4. Address missing data in later survey rounds
5. See `WAIC_COMPARISON_REPORT.md` for full model comparison details

## Publication-Ready Outputs
- **Poisson Model**: [gazelle_msNMix_fit.rds](gazelle_msNMix_fit.rds) - 7.45 min runtime, 12,000 samples
- **NB Quick Fit**: [gazelle_msNMix_fit_nb.rds](gazelle_msNMix_fit_nb.rds) - 0.18 min runtime, 4,000 samples (validation only)
- **NB Full Fit**: [gazelle_msNMix_fit_nb_full.rds](gazelle_msNMix_fit_nb_full.rds) - Full publication model with 12,000 samples
- **Comparison Report**: [WAIC_COMPARISON_REPORT.md](WAIC_COMPARISON_REPORT.md)
