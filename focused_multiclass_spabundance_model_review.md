# Focused Multi-Class spAbundance Model Review

## Recommendation

The finalized separate-class negative-binomial N-mixture models should remain the primary analysis of fawn, adult-male, and adult-female responses to wolf and jackal activity. The multi-class model should answer one narrower question: whether the three latent relative-abundance or site-use processes retain residual spatial correlation after accounting for the focal habitat and predator gradients.

For that secondary question, retain the focused two-factor `spAbundance::lfMsNMix` model. Do not use its ecological slopes as replacements for the finalized separate-class estimates, and do not interpret its intercepts or latent abundance as absolute population size.

The focal ecological formula should be:

```r
~ Open.Spaces.z + Null.Wolf.Abundance.z * null.jackal.abundance.z
```

This specification retains both predator main effects required by the interaction hierarchy and controls for the principal habitat gradient before estimating residual class correlations. Nature-reserve and built-and-road cover are not included in the focused correlation model. The earlier broad additive fit remains a formula-sensitivity comparison.

The focused design is not collinear. Variance-inflation factors for open-space cover, wolf activity, jackal activity, and their interaction were 1.05-1.11 across the 60 cameras. Nature-reserve and built-and-road cover were moderately correlated (`r = 0.64`), so retaining all three land-use variables adds complexity that is not necessary for the focal hypothesis.

The two-factor `lfMsNMix` fit is therefore secondary rather than the sole or primary abundance analysis. The fitted latent abundance reached approximately 20,000-35,000 at some camera sites, and abundance and detection intercepts failed to converge. These features indicate weak abundance-detection identifiability and preclude inference about absolute abundance. The free loadings and resulting residual correlations did converge, so the model is retained for that limited estimand with transparent qualifications.

## Literature Basis

`lfMsNMix` is an abundance-based joint species distribution model that uses latent factors to estimate residual response correlations while accounting for imperfect detection. Its documentation states that factor dimension reduction is intended for a number of factors that is small relative to the number of species. With only three response classes, two factors are nearly saturated rather than a meaningful dimension reduction ([official `lfMsNMix` documentation](https://doserlab.com/files/spabundance-web/reference/lfmsnmix)).

The package vignette recommends comparing candidate structures with posterior predictive checks and WAIC, and warns that loading constraints make factor-model results sensitive to response ordering. It also recommends evaluating multiple posterior-predictive discrepancy statistics and grouping levels ([official N-mixture vignette](https://doserlab.com/files/spabundance-web/articles/nmixturemodels)).

The multi-class N-mixture framework is part of the Bayesian abundance-modeling system described by Doser et al. ([2024](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14332)).

## Focused N-Mixture Screening

We fitted three negative-binomial candidates with identical abundance and detection formulas and 50,000 iterations in each of three chains:

1. `msNMix`: no residual class-correlation structure.
2. `lfMsNMix`: one latent factor.
3. `lfMsNMix`: two latent factors.

The shared detection formula was:

```r
~ log_days_operational.z + cattle_sightings.z
```

The one-factor model was rejected. With one factor, residual correlations have absolute value one by construction, and the two free loadings switched signs among chains (`Rhat` = 22.68-23.35). It cannot provide an interpretable correlation comparison for three classes.

The two-factor focal slopes and free loadings mixed adequately (`Rhat` <= 1.011), but abundance and detection intercepts did not (`Rhat` = 1.13-2.77). Site-grouped Freeman-Tukey posterior predictive checks were satisfactory (Bayesian `P` = 0.427 for fawns, 0.402 for males, and 0.326 for females). Ungrouped checks indicated lack of fit, especially for males and females (`P` = 0.052, 0.0019, and 0.0023, respectively). Thus, the model reproduced camera-level totals better than the replicate-level count distribution.

The no-factor model had better ungrouped predictive checks but retained nonconverged abundance and detection intercepts (`Rhat` up to 1.73). Full WAIC comparison was not used for model selection because the weakly identified latent-abundance distributions required integration ceilings as high as 35,569 and made the comparison both computationally extreme and scientifically difficult to defend.

## Focused Two-Factor Estimates

Residual class correlations from the focused two-factor model were:

| Class pair | Mean | 95% CrI |
|---|---:|---:|
| Fawn-male | 0.66 | 0.43-0.81 |
| Fawn-female | 0.86 | 0.51-1.00 |
| Male-female | 0.91 | 0.73-1.00 |

These estimates are similar to those from the earlier broad additive model. The similarity suggests that the high residual correlations were not produced by including nature-reserve and built-and-road cover. However, they remain conditional residual associations and do not demonstrate demographic interactions.

The focused two-factor model supported positive open-space associations for all classes, positive wolf associations for all classes, and negative jackal associations for all classes. The wolf by jackal interaction was positive with a 95% credible interval excluding zero only for females.

| Class | Open spaces | Wolf | Jackal | Wolf x jackal |
|---|---:|---:|---:|---:|
| Fawn | 0.90 [0.42, 1.42] | 0.47 [0.05, 0.91] | -0.59 [-1.11, -0.10] | 0.38 [-0.12, 0.90] |
| Male | 0.42 [0.04, 0.80] | 0.53 [0.15, 0.93] | -0.59 [-0.99, -0.20] | 0.37 [-0.05, 0.81] |
| Female | 0.69 [0.28, 1.13] | 0.64 [0.23, 1.06] | -0.64 [-1.11, -0.21] | 0.50 [0.02, 0.99] |

The ecological coefficients above are retained as a model audit, not as the manuscript's primary class-specific estimates. The manuscript should use the finalized separate-class models for predator-response inference and the joint model only for the three residual correlations.

## Final Reporting Rules

1. Lead with the finalized separate-class N-mixture models for wolf, jackal, and wolf by jackal effects.
2. Present the focused two-factor `lfMsNMix` as a secondary joint analysis of residual class correlation.
3. Report only the free-loading diagnostics, derived correlations, and relevant posterior predictive checks from the joint model.
4. Do not report joint-model intercepts, absolute abundance, or its ecological slopes as primary results.
5. Describe the correlations as residual associations among latent relative local-abundance or site-use processes, not causal interactions.
6. Retain the earlier broad additive model as a sensitivity check because it produced similar correlations.

## Reproducible Files

- `fit_focused_multiclass_nmix.R`
- `evaluate_focused_multiclass_nmix.R`
- `publication_focused_multiclass_nmix_comparison_coefficients.csv`
- `publication_focused_multiclass_nmix_comparison_correlations.csv`
- `publication_focused_multiclass_nmix_comparison_diagnostics.csv`
- `publication_focused_multiclass_nmix_comparison_ppc.csv`
- `publication_focused_multiclass_nmix_comparison_common_nmax.csv`
- `publication_focused_lfmsnmix_class_correlations.csv`
- `publication_focused_lfmsnmix_diagnostics.csv`
- `publication_focused_lfmsnmix_ppc.csv`
