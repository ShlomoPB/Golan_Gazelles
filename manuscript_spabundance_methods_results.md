# Class-Specific Abundance Responses and Residual Demographic Correlations

## Methods

### Repeated-count response

To test whether fawns, adult females, and adult males differed in their responses to sympatric predators, we analyzed class-specific counts from 60 camera stations surveyed over seven rounds. After excluding camera-rounds without valid operation, each class contributed 344 nonmissing repeated counts. The count histories contained 199 fawn, 1,951 adult-female, and 2,355 adult-male detections.

N-mixture models separate repeated counts into a latent site-abundance process and an observation process conditional on latent abundance (Royle 2004). These models assume demographic closure across repeat visits and conditional independence of detections. Because demographic groups could enter or leave a camera's sampling area and unmarked camera detections may include repeated observations of individuals, we interpret model responses as relative local abundance or site use rather than absolute population size (Barker et al. 2018).

### Primary class-specific models

We fitted a separate negative-binomial N-mixture model to each demographic class with `spAbundance::NMix` 0.2.1 in R 4.5.3 (Doser et al. 2024; R Core Team 2026). For each class, the abundance formula was

```r
~ Null.Wolf.Abundance.z * null.jackal.abundance.z
```

The interaction specification included the wolf and jackal main effects. Because both activity covariates were standardized, each main effect describes a 1-SD increase when the other predator was at its mean, and the interaction tests whether one predator's association varied across the activity gradient of the other. The detection formula was

```r
~ log_days_operational.z + cattle_sightings.z
```

where camera operating days were log transformed and both detection predictors were standardized. We assigned normal priors with mean zero and variance 100 to abundance coefficients, normal priors with mean zero and variance 2.72 to detection coefficients, and a uniform prior from zero to 10 to the negative-binomial overdispersion parameter. We ran three chains for 50,000 iterations each, discarded 5,000 iterations per chain, and retained every tenth draw, yielding 13,500 posterior samples per class.

We assessed convergence with `Rhat` and effective sample size (ESS) and considered a coefficient supported when its 95% equal-tail credible interval (CrI) excluded zero. For supported interactions, we combined posterior main-effect and interaction draws to estimate the abundance ratio associated with a 1-SD increase in wolf activity when jackal activity was -1, 0, and +1 SD. We calculated corresponding conditional jackal effects across the wolf gradient. Because abundance and detection intercepts did not converge, response curves excluded the intercept and were normalized to the expected response at mean wolf and jackal activity.

As a framework sensitivity analysis, we fitted matched negative-binomial `pcount` models in `unmarked` 1.5.0 (Kellner et al. 2023). These models used the same wolf by jackal abundance formula and modeled detection as a function of path type, camera model, vegetation structure, and passability behind the camera. We treated agreement in coefficient direction as a robustness check rather than independent replication.

### Secondary latent-correlation model

To quantify residual spatial association among demographic classes, we jointly fitted fawn, adult-male, and adult-female counts in a two-factor negative-binomial N-mixture model with `spAbundance::lfMsNMix`. This model was used only to estimate residual class correlations; ecological inference remained based on the separate-class models.

The joint abundance formula was

```r
~ Open.Spaces.z +
  Null.Wolf.Abundance.z * null.jackal.abundance.z
```

Including open-space cover controlled for the principal habitat gradient before estimating residual class association. Nature-reserve and built-and-road cover were excluded from this focused model. The joint detection formula was the same effort-based formula used in the primary models:

```r
~ log_days_operational.z + cattle_sightings.z
```

The joint model estimated class-specific coefficients under community-level coefficient distributions. We used two factors because a one-factor model with three response classes constrains pairwise residual correlations to have absolute value one and failed to produce stable loadings in screening fits. We assigned mean-zero normal priors with variance 100 to community abundance means, mean-zero normal priors with variance 2.72 to community detection means, inverse-gamma priors with shape and scale 0.1 to among-class coefficient variances, and uniform priors from zero to 10 to class-specific negative-binomial overdispersion parameters.

We ran three chains for 50,000 iterations each, discarded 5,000 iterations per chain, and retained every tenth draw. Within each posterior draw, we calculated the latent covariance matrix as `Lambda %*% t(Lambda)` and standardized it to a correlation matrix. These correlations represent residual associations among class-specific latent relative-abundance or site-use processes after accounting for open-space cover and wolf-jackal context. They do not identify causal, behavioral, or demographic interactions.

We evaluated free-loading convergence and assessed calibration with marginal Freeman-Tukey posterior predictive checks grouped by camera site. We also inspected ungrouped and replicate-grouped discrepancies. Because class-specific abundance and detection intercepts did not converge and latent abundance was weakly identified, we did not interpret joint-model intercepts, absolute abundance, or ecological slopes.

## Results

### Class-specific wolf-jackal responses

All nine focal abundance slopes in the separate-class models mixed adequately (`Rhat` = 1.001-1.004; ESS = 3,549-9,538), as did the overdispersion parameters (`Rhat` = 1.001-1.005). Abundance and detection intercepts did not converge, so inference was restricted to relative predator effects.

Fawn relative abundance declined with jackal activity. At mean wolf activity, the jackal coefficient was -0.80 (95% CrI = -1.53 to -0.13), equivalent to an abundance ratio of 0.48 for a 1-SD increase in jackal activity. The wolf coefficient was uncertain (0.27, -0.27 to 0.85), as was the wolf by jackal interaction (0.50, -0.21 to 1.30). Thus, the fawn model supported a negative jackal gradient but not a wolf association or interaction (Fig. X).

Adult females showed a positive association with wolf activity (0.54, 95% CrI = 0.13-0.97), a negative association with jackal activity (-0.72, -1.19 to -0.27), and a positive wolf by jackal interaction (0.52, 0.04-1.02). At mean jackal activity, a 1-SD increase in wolf activity corresponded to a 1.76-fold change in expected relative abundance (95% CrI = 1.14-2.64). The wolf ratio increased from 1.05 (0.62-1.66) at low jackal activity to 3.13 (1.40-6.51) at high jackal activity. Conversely, the jackal ratio increased from 0.30 (0.16-0.52) at low wolf activity to 0.88 (0.39-1.74) at high wolf activity, indicating that the negative jackal association weakened as wolf activity increased.

Adult males showed the same conditional structure. Wolf activity was positively associated with relative abundance (0.35, 95% CrI = 0.01-0.71), jackal activity was negatively associated (-0.61, -0.96 to -0.26), and their interaction was positive (0.39, 0.03-0.77). The wolf abundance ratio increased from 0.98 (0.63-1.44) at low jackal activity to 1.45 (1.01-2.04) at mean jackal activity and 2.21 (1.19-3.91) at high jackal activity. The corresponding jackal ratio changed from 0.38 (0.22-0.62) at low wolf activity to 0.84 (0.49-1.35) at high wolf activity.

The `unmarked` sensitivity models recovered the same qualitative pattern. Jackal coefficients were negative for fawns (`P` = 0.011), females (`P` = 0.001), and males (`P` < 0.001). Wolf coefficients were positive for females (`P` = 0.027) and males (`P` = 0.020), and wolf by jackal interactions were positive for females (`P` = 0.017) and males (`P` = 0.006). The fawn interaction was weaker and unsupported (`P` = 0.092). Agreement across implementations supports a two-predator spatial context for adult gazelles, but these observational associations do not demonstrate that wolves suppress jackals or cause gazelle responses.

### Residual correlation among classes

The free latent-factor loadings in the focused two-factor model converged (`Rhat` = 1.001-1.015; ESS = 323-710). Site-grouped Freeman-Tukey posterior predictive checks were not extreme (Bayesian `P` = 0.408 for fawns, 0.402 for males, and 0.318 for females). Ungrouped checks indicated remaining replicate-level lack of fit, particularly for adult males and females (Bayesian `P` <= 0.002). Abundance-intercept `Rhat` values ranged from 1.35 to 1.58 and detection-intercept values ranged from 1.50 to 1.82; these parameters and absolute abundance were therefore excluded from inference.

Residual latent relative abundance was positively correlated among all demographic classes after accounting for open-space cover and wolf-jackal context (Table X). The posterior correlation was 0.66 (95% CrI = 0.44-0.81) between fawns and males, 0.86 (0.52-1.00) between fawns and females, and 0.91 (0.71-1.00) between males and females. These correlations indicate shared unexplained spatial structure, but they should be interpreted cautiously because the joint model did not reproduce replicate-level variation adequately.

| Class pair | Residual correlation | 95% CrI |
|---|---:|---:|
| Fawn-male | 0.66 | 0.44-0.81 |
| Fawn-female | 0.86 | 0.52-1.00 |
| Male-female | 0.91 | 0.71-1.00 |

## Table Caption

**Table X. Residual correlations among demographic-class latent relative-abundance processes.** Values are posterior means and 95% credible intervals from a two-factor negative-binomial `lfMsNMix` model jointly fitted to fawn, adult-male, and adult-female repeated counts. Correlations were derived from the standardized latent-factor covariance matrix and represent residual spatial association after accounting for open-space cover, wolf activity, jackal activity, and their interaction. Positive correlations indicate shared unexplained spatial structure, not causal demographic interactions.

Table file: `publication_focused_lfmsnmix_class_correlations.csv`.

## Figure Captions

**Figure X. Class-specific relative abundance across wolf-jackal activity contexts.** Curves show posterior medians and shaded 95% credible intervals from separate negative-binomial `spAbundance` N-mixture models for fawns, adult females, and adult males. Colors represent wolf activity at -1, 0, and +1 SD, and the horizontal axis shows standardized jackal activity. Relative abundance is normalized to the expected value at mean wolf and jackal activity and is displayed on a logarithmic scale. The interaction is supported for females and males but not fawns. Because abundance and detection intercepts did not converge, the figure excludes the intercept and should be interpreted as a comparison of relative gradients, not absolute population size.

Figure files: `publication_spabundance_wolf_jackal_relative_interaction_curves.png` and `.pdf`.

**Figure S1. Relative abundance across a 3 x 3 wolf-jackal interaction grid.** Cells show posterior median expected abundance at -1, 0, and +1 SD of wolf and jackal activity, divided by expected abundance at mean wolf and jackal activity. Values greater than one indicate higher and values below one indicate lower relative abundance than the mean-predator reference. The grid excludes the nonconverged abundance intercept and summarizes response-surface geometry rather than absolute abundance.

Figure files: `publication_spabundance_wolf_jackal_relative_heatmap.png` and `.pdf`.

## References

Barker, R. J., M. R. Schofield, W. A. Link, and J. R. Sauer. 2018. On the reliability of N-mixture models for count data. *Biometrics* 74:369-377. https://doi.org/10.1111/biom.12734

Doser, J. W., A. O. Finley, M. Kery, and E. F. Zipkin. 2024. spAbundance: An R package for single-species and multi-species spatially explicit abundance models. *Methods in Ecology and Evolution* 15:1024-1033. https://doi.org/10.1111/2041-210X.14332

Kellner, K. F., A. D. Smith, J. A. Royle, M. Kery, J. L. Belant, and R. B. Chandler. 2023. The unmarked R package: Twelve years of advances in occurrence and abundance modelling in ecology. *Methods in Ecology and Evolution* 14:1408-1415. https://doi.org/10.1111/2041-210X.14123

R Core Team. 2026. *R: A language and environment for statistical computing*. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Royle, J. A. 2004. N-mixture models for estimating population size from spatially replicated counts. *Biometrics* 60:108-115. https://doi.org/10.1111/j.0006-341X.2004.00142.x

Wickham, H. 2016. *ggplot2: Elegant graphics for data analysis*. Springer-Verlag, New York. https://doi.org/10.1007/978-3-319-24277-4
