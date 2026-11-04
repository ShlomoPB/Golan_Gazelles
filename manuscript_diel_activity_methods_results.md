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

Within each gazelle class, we evaluated all six pairwise contrasts among predator contexts. We estimated `Dhat1` when the smaller context sample contained fewer than 75 events and `Dhat4` otherwise. We calculated 95% intervals with 500 camera-cluster bootstrap replicates and tested equality of the circular distributions with `activity::compareCkern` using 999 null replicates and a +1 Monte Carlo correction. We excluded comparisons with fewer than 20 events in either context from formal inference and report unadjusted `P` values. Because `compareCkern` resamples event times rather than cameras, we interpreted its `P` values together with camera-bootstrap intervals and effect magnitudes.

### Nocturnality

To distinguish redistribution within the day from shifts toward nighttime, we calculated a continuous nocturnality value for every independent detection following Procko et al. (2023). This metric was the absolute decimal-hour distance from apparent solar noon and ranged from 0 at solar noon to 12 at solar midnight. For each demographic class, we fitted a Gaussian mixed model in `glmmTMB` 1.1.13 (Brooks et al. 2017), with predator context as a fixed effect and camera as a random intercept. We tested the overall context effect by likelihood-ratio comparison with an intercept-only mixed model and estimated pairwise context differences in hours with `emmeans` 1.11.2-8 (Lenth 2026). We report unadjusted two-sided `P` values. Contexts represented by fewer than 20 events or fewer than five cameras with detections were excluded from formal nocturnality inference.

## Results

### Overall class activity

Across the full survey, adult-female and adult-male activity distributions overlapped strongly (`Dhat4` = 0.948, 95% camera-bootstrap CI = 0.859-0.960). Fawn activity was more distinct from adult activity, with overlap of 0.787 (0.695-0.855) with females and 0.742 (0.650-0.782) with males.

### Activity across wolf-jackal contexts

The median split assigned 14 cameras to low both, 16 to high jackal/low wolf, 16 to low jackal/high wolf, and 14 to high both. Event counts in these contexts were 67, 3, 72, and 38 for fawns; 531, 69, 601, and 276 for females; and 791, 163, 770, and 544 for males, respectively. The high-jackal/low-wolf fawn context was therefore data deficient and was not assigned a density curve or included in formal comparisons (Fig. Y).

Female diel distributions differed among four of the six predator-context contrasts. The lowest overlap occurred between low-jackal/high-wolf and high-both cameras (`Dhat4` = 0.788, 95% camera-bootstrap CI = 0.652-0.883; `P` = 0.001). High-jackal/low-wolf activity differed from low-jackal/high-wolf activity (`Dhat1` = 0.794, 0.585-0.859; `P` = 0.007), and low-both activity differed from both low-jackal/high-wolf (`Dhat4` = 0.799, 0.682-0.874; `P` = 0.001) and high-both (`Dhat4` = 0.814, 0.632-0.896; `P` = 0.001). The two contrasts involving high-jackal/low-wolf cameras and either low-both or high-both cameras were unsupported (`P` = 0.179 and 0.112, respectively).

Male context differences were statistically detectable but generally smaller in magnitude. The lowest supported overlap was between high-jackal/low-wolf and high-both cameras (`Dhat4` = 0.838, 95% camera-bootstrap CI = 0.652-0.894; `P` = 0.002). Low-both activity differed from high-both (`Dhat4` = 0.902; `P` = 0.009) and low-jackal/high-wolf activity (`Dhat4` = 0.912; `P` = 0.012), and low-jackal/high-wolf activity differed from high-both (`Dhat4` = 0.906; `P` = 0.018). Because these latter overlap coefficients exceeded 0.90, they indicate modest redistribution detectable with large samples rather than strong temporal separation.

None of the eligible fawn context comparisons was supported. The lowest eligible fawn overlap was between low-jackal/high-wolf and high-both cameras (`Dhat1` = 0.796, 95% camera-bootstrap CI = 0.477-0.854; `P` = 0.179). Comparisons involving high-jackal/low-wolf cameras were not formally tested because that context contained only three fawn events. The low count is an important spatial-demographic observation, but it cannot establish a diel response.

Gazelle detections remained predominantly diurnal or crepuscular across all four predator contexts (Fig. Y). The supported female and male differences therefore represented redistribution mainly within daylight and twilight rather than shifts toward nocturnal activity. Females in high-jackal/low-wolf cameras were broadly active around midday, whereas females and males in high-wolf contexts showed stronger dusk activity. Fawns showed morning-centered activity in low-both and high-both contexts and a broader morning-to-dusk distribution in low-jackal/high-wolf cameras. These associations do not establish that wolf or jackal activity caused the observed temporal differences because habitat and other site characteristics may also differ among contexts.

### Nocturnality across predator contexts

Mean distance from solar noon did not differ overall among predator contexts for fawns (likelihood-ratio chi-square = 0.55, df = 2, `P` = 0.760), females (chi-square = 3.90, df = 3, `P` = 0.273), or males (chi-square = 6.00, df = 3, `P` = 0.111; Fig. Z). Context-specific raw means ranged from 3.40 to 3.74 h for fawns, 4.49 to 5.77 h for females, and 5.06 to 6.16 h for males. The high-jackal/low-wolf fawn context was excluded from its model because it contained three events from two cameras.

One unadjusted male contrast indicated greater nocturnality at high-jackal/low-wolf than at high-both cameras (estimated difference = 1.21 h, 95% CI = 0.18-2.25, `P` = 0.021). However, the overall male context test was unsupported and multiple raw pairwise contrasts were examined, so we treated this individual contrast as exploratory rather than evidence of a general nocturnality response. Taken together, the nocturnality models did not support systematic movement toward or away from nighttime. Thus, the differences detected by the circular analyses primarily reflected changes in the shape and timing of activity within the diurnal-crepuscular period, which the scalar nocturnality metric intentionally cannot distinguish.

## Figure Captions

**Figure Y. Class-specific gazelle diel activity across wolf-jackal contexts.** Panels show normalized activity distributions for fawns, adult females, and adult males, with curves representing the four camera contexts defined by median wolf and jackal activity indices. Curves are based on 30-min independent gazelle events converted from Israel local civil time to apparent solar time. Shading denotes approximate night, civil twilight, and daylight. Each curve integrates to one and therefore represents timing rather than event rate or abundance. The high-jackal/low-wolf fawn group contained only three events, which are shown as event ticks without an estimated density curve.

Figure files: `publication_predator_context_diel_by_class_median_combined.png` and `.pdf`.

**Figure Z. Gazelle nocturnality across wolf-jackal contexts.** Points and horizontal lines show estimated mean hours from apparent solar noon and 95% confidence intervals from class-specific Gaussian mixed models with camera random intercepts. Higher values indicate activity closer to solar midnight. The high-jackal/low-wolf fawn group is shown descriptively and was excluded from formal inference because it contained only three events from two cameras.

Figure files: `publication_predator_context_nocturnality_median.png` and `.pdf`.

**Figure S1. Within-class temporal overlap among wolf-jackal contexts.** Points show pairwise coefficients of overlap and horizontal lines show 95% camera-cluster bootstrap intervals. Lower overlap indicates a larger difference in diel timing between contexts. `Dhat1` was used when the smaller sample contained fewer than 75 events and `Dhat4` otherwise. Statistical comparisons used `compareCkern` with 999 null replicates and unadjusted `P` values. Comparisons with fewer than 20 events in either context are descriptive only.

Figure files: `publication_predator_context_diel_pairwise_overlap_median.png` and `.pdf`.

## References

Brooks, M. E., K. Kristensen, K. J. van Benthem, A. Magnusson, C. W. Berg, A. Nielsen, H. J. Skaug, M. Maechler, and B. M. Bolker. 2017. `glmmTMB` balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. *The R Journal* 9:378-400. https://doi.org/10.32614/RJ-2017-066

Lenth, R. V. 2026. *emmeans: Estimated marginal means, aka least-squares means*. R package version 1.11.2-8. https://CRAN.R-project.org/package=emmeans

Meredith, M., M. Ridout, and L. A. D. Campbell. 2024. *overlap: Estimates of coefficient of overlapping for animal activity patterns*. R package version 0.3.9. https://doi.org/10.32614/CRAN.package.overlap

Nouvellet, P., G. S. A. Rasmussen, D. W. Macdonald, and F. Courchamp. 2012. Noisy clocks and silent sunrises: measurement methods of daily activity pattern. *Journal of Zoology* 286:179-184. https://doi.org/10.1111/j.1469-7998.2011.00864.x

R Core Team. 2026. *R: A language and environment for statistical computing*. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Ridout, M. S., and M. Linkie. 2009. Estimating overlap of daily activity patterns from camera trap data. *Journal of Agricultural, Biological, and Environmental Statistics* 14:322-337. https://doi.org/10.1198/jabes.2009.08038

Procko, M., R. Naidoo, V. LeMay, and A. C. Burton. 2023. Human presence and infrastructure impact wildlife nocturnality differently across an assemblage of mammalian species. *PLOS ONE* 18:e0286131. https://doi.org/10.1371/journal.pone.0286131

Rowcliffe, J. M. 2023. *activity: Animal activity statistics*. R package version 1.3.4. https://doi.org/10.32614/CRAN.package.activity

Rowcliffe, J. M., R. Kays, B. Kranstauber, C. Carbone, and P. A. Jansen. 2014. Quantifying levels of animal activity using camera trap data. *Methods in Ecology and Evolution* 5:1170-1179. https://doi.org/10.1111/2041-210X.12278
