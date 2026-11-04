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
