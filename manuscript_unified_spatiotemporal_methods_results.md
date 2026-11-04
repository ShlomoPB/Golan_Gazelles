# Unified spatiotemporal analysis of demographic-class responses

## Statistical methods

We analyzed 30-min independent camera-trap events for fawns, adult females, and adult males using complementary latent-abundance, spatial-use, demographic-composition, and diel-activity models. Analyses were restricted to active camera-rounds within the formal survey window (20 August to 15 December 2020). Of 3,932 classified independent events, 3,906 occurred during valid camera operation and were retained. Continuous site covariates were standardized to a mean of zero and standard deviation of one. The survey was divided at its calendar midpoint, 17 October 2020, to provide a detection-independent early versus late comparison. Because temperature measurements were unavailable, this contrast was interpreted as survey period rather than as a direct thermal effect.

### Latent abundance

We used single-class negative-binomial N-mixture models as the primary spatial abundance analysis, with parallel frequentist (`unmarked::pcount`) and Bayesian (`spAbundance::NMix`) implementations. Detection was modeled with the shared camera-placement and camera-performance structure selected from the species-level detection analyses: path type, camera model, vegetation structure, and area behind the camera. Using the same detection model across demographic classes reduced the risk that class differences were produced by class-specific observation models.

Class-specific abundance models were used to estimate each class's dominant habitat or pressure association, and predator-interaction models tested whether wolf and jackal activity had additive or context-dependent associations with latent abundance. Because N-mixture models depend on closure, repeated-count, and detection-identifiability assumptions that are imperfect for camera-trap class counts, we treated agreement between `unmarked`, `spAbundance`, and independent-event models as the strongest evidence, rather than relying on any single model family in isolation.

### Spatial use

As a robustness and interpretation layer, we modeled the number of independent events for each class in each camera-round using negative-binomial generalized linear mixed models in `glmmTMB`. The logarithm of operational camera-days was included as an offset, so the response represents independent-event rate per operational day rather than latent animal abundance. The primary NB2 model included demographic class interactions with standardized wolf activity, jackal activity, their interaction, open-space cover, built and road cover, and human activity. Nature-reserve cover and survey period were included as additive covariates. Nested random intercepts for camera within polygon accounted for repeated measurements and spatial grouping, and class-specific dispersion parameters accommodated differences in count variability among classes.

We compared the class-specific model with a model in which ecological slopes were shared among classes using a likelihood-ratio test. Direct linear contrasts tested differences between class-specific slopes. To determine whether covariate effects primarily reflected camera-round occurrence or the number of events after occurrence, we also fit a binomial encounter model and a zero-truncated NB2 positive-intensity model. NB1 and zero-inflated NB2 formulations were fit as sensitivity analyses. Model adequacy was assessed from convergence codes, Hessian diagnostics, gradients, information criteria, and simulation-based checks of zeros, means, variances, and maxima.

### Demographic composition

We quantified fawn representation relative to adult females with a beta-binomial mixed model of the paired fawn and female counts. The model included wolf activity, jackal activity, their interaction, open-space cover, built and road cover, human activity, nature-reserve cover, survey round, and random intercepts for camera and polygon. This response is a camera-based demographic-composition index and should not be interpreted as recruitment, survival, or a population-level fawn:female ratio.

### Diel activity

Camera clocks were set to Israel local time at installation, so event timestamps were parsed as local clock time in `Asia/Jerusalem` before conversion to apparent solar time using event date and camera coordinates. Daylight was classified from calculated solar elevation. We modeled the number of daylight versus non-daylight events with beta-binomial mixed models containing class, survey period, predator and landscape covariates, and camera and polygon random intercepts. We compared common and class-specific ecological slopes by likelihood-ratio test.

We additionally estimated class- and period-specific diel curves using a cyclic negative-binomial generalized additive model. Independent events were aggregated into 2-h solar-time bins, operational exposure was included as an offset, separate cyclic smooths were estimated for each class-period combination, and camera and survey round were included as random-effect smooths. Pairwise activity overlap was estimated with the `overlap` package, with uncertainty obtained by resampling cameras rather than individual events to preserve within-camera dependence.

## Results

### Latent abundance models

The single-class N-mixture models and the independent-event models produced the same broad spatial conclusions. In the class-specific `spAbundance` models, fawn abundance was most strongly associated with open habitat (posterior mean 1.05, 95% CrI: 0.53 to 1.61), while boar, human, and wolf effects were weaker. Male abundance was negatively associated with built and road cover (-0.64, -1.04 to -0.24) and jackal activity (-0.73, -1.11 to -0.36). Female abundance was most clearly associated with predator covariates, declining with jackal activity (-0.81, -1.34 to -0.30) and increasing with wolf activity (0.42, 0.01 to 0.84), while the nature-reserve effect was weak (-0.10, -0.55 to 0.35).

Predator-interaction N-mixture models showed consistent signs across `unmarked` and `spAbundance`, making the wolf by jackal interaction one of the strongest cross-framework results. For fawns, jackal activity was negative in both frameworks (`unmarked` estimate -0.742, \(P = 0.011\); `spAbundance` posterior mean -0.805, 95% CrI: -1.55 to -0.12), whereas wolf and wolf by jackal terms were positive but less certain. For females, wolf activity was positive, jackal activity was negative, and the wolf by jackal interaction was positive in both frameworks (`unmarked` interaction 0.444, \(P = 0.017\); `spAbundance`: wolf 0.536, 0.20 to 0.91; jackal -0.656, -1.03 to -0.28; interaction 0.595, 0.16 to 1.05). Males showed the same predator-context pattern (`unmarked` interaction 0.423, \(P = 0.006\); `spAbundance`: wolf 0.341, 0.03 to 0.66; jackal -0.648, -0.96 to -0.32; interaction 0.402, 0.04 to 0.78). Thus, the N-mixture analyses support a consistent negative jackal association across classes and a stronger positive wolf association for adults than for fawns.

Ecologically, the positive interaction means that the negative association with jackals was weakest where wolf activity was high, or equivalently that the positive wolf association was strongest where jackal activity was also high. This pattern is consistent with gazelles responding to a two-predator landscape rather than to each predator independently. However, it should be interpreted as a predator-context association, not as direct evidence that wolves suppress jackals, unless temporally resolved predator data are used to test that mechanism.

### Independent-event spatial use

Ecological slopes differed among demographic classes (likelihood-ratio test, \(P = 0.008\)). Jackal activity was negatively associated with event rate for all three classes. A 1-SD increase in jackal activity corresponded to estimated rate ratios of 0.58 for fawns (95% CI: 0.35-0.94), 0.45 for females (0.30-0.66), and 0.46 for males (0.31-0.66). Wolf activity was positively associated with adults, with rate ratios of 1.83 for females (1.27-2.63) and 1.59 for males (1.12-2.27), but its fawn association was less precise (1.41, 0.91-2.16).

Open-space cover was most strongly associated with fawn event rate (rate ratio 2.51, 95% CI: 1.52-4.13), followed by females (1.82, 1.21-2.72); the male effect was weaker and included zero on the log scale (1.41, 0.95-2.10). Direct contrasts confirmed that the open-space slope was larger for fawns than males (difference on the log scale 0.58, 95% CI: 0.24-0.91) and larger for females than males (0.25, 0.07-0.43).

Built and road cover was negatively associated with female and male event rates, which were approximately halved per 1-SD increase (female rate ratio 0.50, 95% CI: 0.29-0.87; male 0.51, 0.30-0.87). The fawn association was weak (0.90, 0.48-1.69). Consequently, fawns were less negatively associated with infrastructure than either females (slope difference 0.58, 95% CI: 0.17-0.99) or males (0.57, 0.17-0.97). Human activity had no independently supported association after accounting for predator and landscape covariates.

The positive wolf by jackal interaction was supported for females and males but not fawns, matching the N-mixture interaction models. This cross-model agreement suggests that adult gazelle spatial use is shaped by the joint predator context, whereas fawn responses are more strongly dominated by open habitat and jackal activity.

### Encounter and positive-event intensity

The hurdle decomposition indicated that much of the spatial pattern arose from whether a class was encountered during a camera-round. Fawn encounter probability increased with wolf activity and open-space cover. Female encounters increased with wolf activity and open-space cover and declined with jackal activity, while male encounters increased with wolf activity and declined with jackal activity and built and road cover. Positive-event intensity after encounter retained positive wolf and open-space associations and a negative jackal association for females; male positive intensity declined with jackal activity. Thus, the principal patterns were not solely caused by a few high-count camera-rounds.

### Fawn representation relative to females

Fawn representation increased with open-space cover (odds ratio 1.80, 95% CI: 1.13-2.89) and built and road cover (2.23, 1.22-4.10). The infrastructure result does not indicate that fawns selected built areas. Instead, the event-rate models showed little fawn response but a strong female decline, causing fawns to comprise a greater fraction of classified fawn-plus-female observations where female use was lower.

The wolf by jackal composition interaction was negative but uncertain (odds ratio 0.64, 95% CI: 0.41-1.02; likelihood-ratio \(P = 0.079\)). This provides suggestive, not conclusive, evidence that fawn representation may be lowest where both predator indices are high. Neither nature-reserve cover nor human activity had a supported independent association with composition.

### Temporal activity and survey period

Adult female and male activity distributions overlapped strongly across the full survey using the confirmed local-clock interpretation (overlap 0.948, 95% camera-clustered bootstrap CI: 0.856-0.958). Fawn activity was more distinct from both females (0.789, 0.707-0.860) and males (0.743, 0.661-0.783). Fawn-female overlap increased from 0.712 in the early period to 0.851 in the late period, whereas female-male overlap declined from 0.954 to 0.897. These values describe distributional overlap and are not formal tests of equality.

The cyclic GAM detected non-uniform diel activity for fawns and females in the early period and for all classes in the late period; early male activity was comparatively weakly structured. The fitted curves indicate a stronger daytime concentration for fawns and greater adult similarity, with more pronounced late-period diel structure.

Early-to-late daylight use differed most clearly for adult females. Females showed greater daylight use in the late period (odds ratio 1.60, 95% CI: 1.16-2.20), males showed a weaker positive shift (1.27, 0.99-1.64), and fawns remained uncertain (0.52, 0.19-1.47). Because the early-late split is a survey-period contrast rather than a measured thermal covariate, this should be interpreted as a seasonal activity shift rather than direct evidence of temperature-driven behavior.

### Model assessment

All primary NB2, hurdle, composition, and daylight models converged with positive-definite Hessians. The zero-inflated model did not have a positive-definite Hessian and was rejected. NB2 was favored over NB1 by 98.3 AIC units and retained the same principal coefficient directions. Simulation checks showed that NB2 reproduced observed zero frequencies but tended to overpredict adult count variance and extreme values, especially for males. NB1 better reproduced some extremes but overpredicted adult zeros and had substantially poorer information-criterion support. Inference should therefore emphasize effect directions and contrasts rather than predictions in the extreme upper tail.

The current R environment reports that `glmmTMB` 1.1.13 was built against `TMB` 1.9.18 while `TMB` 1.9.17 is installed. This did not prevent convergence, but `glmmTMB` should be rebuilt against the installed `TMB` version before the archived publication run.

## Ecological synthesis

The integrated evidence supports three main conclusions. First, gazelles responded to the predator community as a joint wolf-jackal context: jackal activity was consistently negative, wolf activity was more positive for adults, and both `unmarked` and `spAbundance` supported positive wolf by jackal interactions for females and males. Second, open habitat was particularly important for fawns, while adults were more negatively associated with infrastructure. Third, fawns differed from adults temporally as well as spatially, and adult females showed increased daylight activity later in the survey.

The strongest next extension would replace static, survey-wide predator indices with predator activity matched to camera-round or shorter time windows. A dynamic model could then test whether gazelle event rate or timing changes after recent wolf or jackal detections while controlling for camera effort, habitat, season, and repeated site use. That design would distinguish immediate behavioral responses from broad spatial covariance more effectively than additional categorization of cameras into predator quantiles.
