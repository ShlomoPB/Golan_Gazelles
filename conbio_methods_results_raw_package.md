# Final Methods and Results Analysis Package

Generated from the retained corrected-data analyses on 2026-08-18. Early/Late and daylight diel models, seasonal N-mixture models, gazelle-predator temporal-overlap models, PAMMs, and generated maps are excluded.

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

# Fawn Representation Relative to Adult Females

## Methods

### Demographic-composition response

To test whether predator activity, habitat, and survey period were associated with the representation of fawns relative to adult females, we analyzed paired fawn and adult-female camera-trap counts at the camera-round level. We modeled the response as `cbind(fawn count, adult-female count)`, such that the fitted probability represented the probability that a classified fawn-or-adult-female detection was a fawn. We excluded camera-rounds in which both counts were zero because those observations contained no information about demographic composition. We interpret this response as an apparent, camera-based fawn representation index rather than as a direct estimate of recruitment, juvenile survival, or the population-level fawn:female ratio.

### Primary beta-binomial model

We fitted generalized linear mixed models with a beta-binomial distribution and logit link. The beta-binomial distribution allowed the variance in fawn representation to exceed that expected under an ordinary binomial model, a common feature of grouped ecological count data (Harrison 2015). We fitted models in R 4.5.3 (R Core Team 2026) with `glmmTMB` 1.1.13, which supports beta-binomial mixed models and flexible random-effects structures (Brooks et al. 2017).

Our primary ecological predictors were camera-level wolf and jackal activity indices, open-space cover, and survey period. We also adjusted for built-and-road cover, nature-reserve cover, human activity, path type, camera model, vegetation structure, and the passability of the area behind the camera. We standardized all continuous covariates to a mean of zero and standard deviation of one before model fitting. Consequently, each continuous-predictor odds ratio represents the change associated with a 1-SD increase while holding the other covariates constant; centering also makes additive effects interpretable at the average value of the other continuous predictors (Schielzeth 2010).

We divided the survey into Early and Late periods using the midpoint of the full deployment window (17 October 2020), determined independently of gazelle detections. We assigned camera-rounds according to their deployment midpoint and treated Early as the reference category. Because this variable represented one prespecified contrast with only two levels, we modeled it as a fixed effect rather than attempting to estimate a random-effect variance across seasons.

To account for repeated camera-rounds and spatial grouping, we included random intercepts for cameras nested within sampling polygons. The retained model was

```r
cbind(fawn, female) ~ wolf_z + jackal_z + season + open_z +
  built_z + reserve_z + human_z + path_type + camera_model +
  vegetation_structure + area_behind_camera +
  (1 | Polygon.Name / camera_id)
```

We selected this season-adjusted additive formulation from a defined model hierarchy. We used likelihood-ratio tests for nested comparisons to evaluate a common Early-Late shift, a wolf by jackal interaction, and season-specific predator slopes. We retained the common survey-period effect because it improved the additive model, whereas the predator interaction and season-specific predator slopes were unsupported. We verified optimizer convergence and a positive-definite Hessian for each fitted model.

We exponentiated logit-scale coefficients to obtain odds ratios (ORs) and calculated 95% Wald confidence intervals as the coefficient plus or minus 1.96 standard errors before exponentiation. An OR greater than one indicates greater fawn representation relative to adult females. We emphasized effect sizes and confidence intervals rather than classifying effects solely by whether `P` crossed 0.05.

### Predator-interaction sensitivity analysis

To test whether the association with one predator depended on the activity of the other, we fitted an otherwise identical sensitivity model that added `wolf_z:jackal_z`. Thus, season, landscape and human covariates, camera-design covariates, and the nested random-effects structure were held constant between the additive and interaction models. We compared these models by likelihood-ratio test and AIC. Because main effects in an interaction model are conditional, we also calculated wolf slopes when jackal activity was `-1`, `0`, and `+1` SD and jackal slopes at the corresponding wolf-activity values. We retained this model as a sensitivity analysis rather than the primary model when the interaction failed to improve fit.

We produced the coefficient figure with `ggplot2` 4.0.0 (Wickham 2016). We displayed the four focal primary-model terms as odds ratios with 95% confidence intervals on a logarithmic axis; adjustment covariates were omitted from the main figure but retained in the model.

## Results

We retained 221 informative camera-rounds from 50 cameras in 12 sampling polygons, comprising 199 fawn and 1,951 adult-female detections. The Early and Late periods contributed 111 and 110 camera-rounds, respectively. The primary model converged with a positive-definite Hessian.

Fawn representation was higher during the Late period than during the Early period (OR = 1.57, 95% CI = 1.01-2.44, `P` = 0.045; Fig. X). Adding the common Early-Late effect improved the additive model (likelihood-ratio chi-square = 4.15, df = 1, `P` = 0.042). A 1-SD increase in open-space cover was also associated with greater fawn representation relative to adult females (OR = 1.55, 95% CI = 1.04-2.32, `P` = 0.031).

We found no clear additive association between predator activity and fawn representation. The odds ratio for a 1-SD increase in wolf activity was 0.93 (95% CI = 0.70-1.23, `P` = 0.608), and that for a 1-SD increase in jackal activity was 0.94 (95% CI = 0.63-1.41, `P` = 0.773). Both confidence intervals included one. The remaining continuous landscape and human-activity adjustment terms were also unsupported. One vegetation-structure contrast, short shrubs relative to grassland, had a large coefficient, but that category was represented by only six informative camera-rounds from two cameras; we therefore treated it as a sparse adjustment estimate rather than a robust ecological result.

Adding the wolf by jackal interaction did not improve the season-adjusted model (likelihood-ratio chi-square = 0.35, df = 1, `P` = 0.552) and increased AIC from 441.08 to 442.72 (delta AIC = 1.65). The interaction odds ratio was 0.86 (95% CI = 0.53-1.39, `P` = 0.538). Conditional wolf slopes ranged from an OR of 1.03 (95% CI = 0.67-1.58) at low jackal activity to 0.76 (0.38-1.52) at high jackal activity. Conditional jackal slopes ranged from 1.18 (0.53-2.66) at low wolf activity to 0.87 (0.55-1.40) at high wolf activity. All conditional intervals included one. Open-space cover (OR = 1.58, 95% CI = 1.06-2.34) and the Late-period contrast (OR = 1.57, 95% CI = 1.01-2.44) remained supported in the interaction model. We therefore retained the additive model as primary and found no evidence that wolf activity, jackal activity, or their interaction was associated with fawn representation relative to adult females.

## Figure Caption

**Figure X. Predictors of fawn representation relative to adult females.** Points show odds ratios and horizontal lines show 95% Wald confidence intervals from the primary season-adjusted additive beta-binomial mixed model. The vertical line at an odds ratio of one denotes no association. Odds ratios greater than one indicate greater fawn representation among classified fawn and adult-female camera-trap counts. Continuous predictors were standardized, so their odds ratios represent a 1-SD increase. Estimates were adjusted for built-and-road cover, nature-reserve cover, human activity, path type, camera model, vegetation structure, area behind the camera, and nested polygon/camera random intercepts.

## References

Brooks, M. E., K. Kristensen, K. J. van Benthem, A. Magnusson, C. W. Berg, A. Nielsen, H. J. Skaug, M. Maechler, and B. M. Bolker. 2017. glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. *The R Journal* 9:378-400. https://doi.org/10.32614/RJ-2017-066

Harrison, X. A. 2015. A comparison of observation-level random effect and beta-binomial models for modelling overdispersion in binomial data in ecology and evolution. *PeerJ* 3:e1114. https://doi.org/10.7717/peerj.1114

R Core Team. 2026. *R: A language and environment for statistical computing*. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Schielzeth, H. 2010. Simple means to improve the interpretability of regression coefficients. *Methods in Ecology and Evolution* 1:103-113. https://doi.org/10.1111/j.2041-210X.2010.00012.x

Wickham, H. 2016. *ggplot2: Elegant graphics for data analysis*. Springer-Verlag, New York. https://doi.org/10.1007/978-3-319-24277-4

# Class-specific Diel Activity Across Wolf-Jackal Contexts

## Methods

### Event preparation and solar time

To test whether gazelle demographic classes differed in diel activity and whether their activity distributions varied across predator contexts, we analyzed 30-min independent camera-trap events for fawns, adult females, and adult males. We restricted detections to documented active camera-rounds within the formal survey period (20 August to 15 December 2020). All 3,925 classified independent events occurred during valid camera operation and were retained: 180 fawn, 1,477 adult-female, and 2,268 adult-male events.

Camera clocks were set to Israel local time at installation. We therefore parsed timestamps in the `Asia/Jerusalem` time zone and converted each event to apparent solar time from its date and camera longitude using an equation-of-time correction. This placed solar noon at approximately 12:00 and reduced seasonal clock-time displacement of activity peaks (Nouvellet et al. 2012). We treated event time as a circular variable on 0 to 2 pi radians.

### Overall class activity

We estimated overall circular activity densities for each demographic class and quantified pairwise overlap with the `overlap` 0.3.9 package (Meredith et al. 2024). The coefficient ranges from zero for nonoverlapping distributions to one for identical distributions (Ridout & Linkie 2009). We used `Dhat1` when the smaller sample contained fewer than 75 events and `Dhat4` otherwise. We obtained 95% intervals by resampling cameras with replacement 500 times and retaining all events from each sampled camera, thereby preserving within-camera dependence.

### Wolf-jackal predator contexts

To examine class-specific activity across a two-predator landscape while retaining all cameras, we divided the 60 cameras at the median site-level jackal activity index (0.495) and median wolf activity index (0.556). This produced four mutually exclusive contexts: low both, high jackal/low wolf, low jackal/high wolf, and high both. These labels describe relative camera-level activity indices and do not represent experimental treatments or fixed thresholds of predation risk.

We fitted circular activity densities separately for each class and context with `activity` 1.3.4 (Rowcliffe et al. 2014; Rowcliffe 2023). Densities were normalized within each class-context group and therefore describe the timing of detections, not event rates or abundance. We did not draw a density curve when a class-context group contained fewer than five independent events.

Within each gazelle class, we evaluated all six pairwise contrasts among predator contexts. We estimated `Dhat1` when the smaller context sample contained fewer than 75 events and `Dhat4` otherwise. We calculated 95% intervals with 500 camera-cluster bootstrap replicates and tested equality of the circular distributions with `activity::compareCkern` using 999 null replicates and a +1 Monte Carlo correction. We excluded comparisons with fewer than 20 events in either context from formal inference and adjusted `P` values across the six tests within each class with the Benjamini-Hochberg false-discovery-rate procedure (Benjamini & Hochberg 1995). Because `compareCkern` resamples event times rather than camera clusters, we interpret its adjusted `P` values together with the camera-bootstrap intervals and effect magnitudes.

## Results

### Overall class activity

Across the full survey, adult-female and adult-male activity distributions overlapped strongly (`Dhat4` = 0.948, 95% camera-bootstrap CI = 0.859-0.960). Fawn activity was more distinct from adult activity, with overlap of 0.787 (0.695-0.855) with females and 0.742 (0.650-0.782) with males.

### Activity across wolf-jackal contexts

The median split assigned 14 cameras to low both, 16 to high jackal/low wolf, 16 to low jackal/high wolf, and 14 to high both. Event counts in these contexts were 67, 3, 72, and 38 for fawns; 531, 69, 601, and 276 for females; and 791, 163, 770, and 544 for males, respectively. The high-jackal/low-wolf fawn context was therefore data deficient and was not assigned a density curve or included in formal comparisons (Fig. Y).

Female diel distributions differed among four of the six predator-context contrasts. The lowest overlap occurred between low-jackal/high-wolf and high-both cameras (`Dhat4` = 0.788, 95% camera-bootstrap CI = 0.652-0.883; within-class FDR `P` = 0.002). High-jackal/low-wolf activity differed from low-jackal/high-wolf activity (`Dhat1` = 0.794, 0.585-0.859; FDR `P` = 0.011), and low-both activity differed from both low-jackal/high-wolf (`Dhat4` = 0.799, 0.682-0.874; FDR `P` = 0.002) and high-both (`Dhat4` = 0.814, 0.632-0.896; FDR `P` = 0.002). The two contrasts involving high-jackal/low-wolf cameras and either low-both or high-both cameras were unsupported after false-discovery-rate correction (`P` >= 0.134).

Male context differences were statistically detectable but generally smaller in magnitude. The lowest supported overlap was between high-jackal/low-wolf and high-both cameras (`Dhat4` = 0.838, 95% camera-bootstrap CI = 0.652-0.894; within-class FDR `P` = 0.012). Low-both activity differed from high-both (`Dhat4` = 0.902; FDR `P` = 0.024) and low-jackal/high-wolf activity (`Dhat4` = 0.912; FDR `P` = 0.024), and low-jackal/high-wolf activity differed from high-both (`Dhat4` = 0.906; FDR `P` = 0.027). Because these latter overlap coefficients exceeded 0.90, they indicate modest redistribution detectable with large samples rather than strong temporal separation.

None of the eligible fawn context comparisons was supported after within-class false-discovery-rate correction. The lowest eligible fawn overlap was between low-jackal/high-wolf and high-both cameras (`Dhat1` = 0.796, 95% camera-bootstrap CI = 0.477-0.854; FDR `P` = 0.268). Comparisons involving high-jackal/low-wolf cameras were not formally tested because that context contained only three fawn events. The low count is an important spatial-demographic observation, but it cannot establish a diel response.

Gazelle detections remained predominantly diurnal or crepuscular across all four predator contexts (Fig. Y). The supported female and male differences therefore represented redistribution mainly within daylight and twilight rather than shifts toward nocturnal activity. Females in high-jackal/low-wolf cameras were broadly active around midday, whereas females and males in high-wolf contexts showed stronger dusk activity. Fawns showed morning-centered activity in low-both and high-both contexts and a broader morning-to-dusk distribution in low-jackal/high-wolf cameras. These associations do not establish that wolf or jackal activity caused the observed temporal differences because habitat and other site characteristics may also differ among contexts.

## Figure Captions

**Figure Y. Class-specific gazelle diel activity across wolf-jackal contexts.** Panels show normalized activity distributions for fawns, adult females, and adult males, with curves representing the four camera contexts defined by median wolf and jackal activity indices. Curves are based on 30-min independent gazelle events converted from Israel local civil time to apparent solar time. Shading denotes approximate night, civil twilight, and daylight. Each curve integrates to one and therefore represents timing rather than event rate or abundance. The high-jackal/low-wolf fawn group contained only three events, which are shown as event ticks without an estimated density curve.

Figure files: `publication_predator_context_diel_by_class_median_combined.png` and `.pdf`.

**Figure Z. Within-class temporal overlap among wolf-jackal contexts.** Points show pairwise coefficients of overlap and horizontal lines show 95% camera-cluster bootstrap intervals. Lower overlap indicates a larger difference in diel timing between contexts. `Dhat1` was used when the smaller sample contained fewer than 75 events and `Dhat4` otherwise. Statistical comparisons used `compareCkern` with 999 null replicates and Benjamini-Hochberg correction within each class. Comparisons with fewer than 20 events in either context are descriptive only.

Figure files: `publication_predator_context_diel_pairwise_overlap_median.png` and `.pdf`.

## References

Benjamini, Y., and Y. Hochberg. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B* 57:289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Meredith, M., M. Ridout, and L. A. D. Campbell. 2024. *overlap: Estimates of coefficient of overlapping for animal activity patterns*. R package version 0.3.9. https://doi.org/10.32614/CRAN.package.overlap

Nouvellet, P., G. S. A. Rasmussen, D. W. Macdonald, and F. Courchamp. 2012. Noisy clocks and silent sunrises: measurement methods of daily activity pattern. *Journal of Zoology* 286:179-184. https://doi.org/10.1111/j.1469-7998.2011.00864.x

R Core Team. 2026. *R: A language and environment for statistical computing*. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Ridout, M. S., and M. Linkie. 2009. Estimating overlap of daily activity patterns from camera trap data. *Journal of Agricultural, Biological, and Environmental Statistics* 14:322-337. https://doi.org/10.1198/jabes.2009.08038

Rowcliffe, J. M. 2023. *activity: Animal activity statistics*. R package version 1.3.4. https://doi.org/10.32614/CRAN.package.activity

Rowcliffe, J. M., R. Kays, B. Kranstauber, C. Carbone, and P. A. Jansen. 2014. Quantifying levels of animal activity using camera trap data. *Methods in Ecology and Evolution* 5:1170-1179. https://doi.org/10.1111/2041-210X.12278

