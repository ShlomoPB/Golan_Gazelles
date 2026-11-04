# Hierarchical Diel Activity Analysis

## Methods

### Event histories and solar-time standardization

To test whether mountain gazelle diel activity varied continuously along wolf and jackal activity gradients, we fitted hierarchical circular generalized additive models (GAMs) separately for juveniles, adult females, and adult males. We also fitted a pooled species-level model as a secondary analysis. The pooled response included all gazelle detections regardless of whether sex or age could be classified; it was not calculated by summing the three demographic classes. We classified detections as juvenile when life stage was recorded as juvenile, as adult female when sex was female and life stage was adult, and as adult male when sex was male and life stage was adult. These definitions made the three class-specific responses mutually exclusive.

We interpreted camera timestamps as Israel local civil time and applied the documented corrections for cameras with incorrectly set calendar dates before assigning detections to camera and survey round. We transformed each corrected timestamp to sun time using the detection date and camera coordinates with `overlap::sunTime()` in `overlap` 0.3.9. This transformation mapped local sunrise to 06:00 and sunset to 18:00 while preserving the relative duration of daytime and nighttime, allowing us to compare diel timing across dates despite seasonal changes in photoperiod (Nouvellet et al. 2012). We treated time as circular over 0-24 h, with cyclic knots at 0 and 24 h.

For the primary analysis, we collapsed the corrected raw records into binary camera-day-hour outcomes. For each of the 24 sun-time hours within every operational camera-round, the number of successes was the number of operational dates with at least one detection and the number of failures was the reported number of operational days minus the number of success dates. Consequently, multiple photographs or detections within the same camera-date-hour could not increase the response above one. We retained all zero-detection hours. The resulting grouped-binomial histories represented 344 positive-effort camera-rounds from 60 cameras and 5,999 camera-days, equivalent to 143,976 hourly exposure trials. As a sensitivity analysis, we reconstructed histories from detections separated by more than 30 min. This reconstruction returned the authoritative total of 3,688 pooled independent gazelle encounters.

### Hierarchical circular GAMs

We fitted models in `mgcv` 1.9-4 using a binomial distribution with a logit link and the grouped response `cbind(success, failure)` (Wood 2017). We adapted the hierarchical diel-activity framework of Iannarilli et al. (2025), which allows the population-level diel curve and site-level deviations to be estimated in a single model rather than treating individual detections as independent observations. Every candidate model included survey round as a categorical fixed effect, a cyclic cubic regression spline for baseline diel activity, and a camera random intercept. We compared additional hierarchical structures containing a camera-round random intercept and, for responses with sufficient data, camera-specific cyclic deviations. Camera-specific deviations had separate coefficients but shared a smoothing parameter through the `mgcv` `id` argument. We selected the least complex hierarchical structure within 2 conditional AIC units of the minimum. The retained structure included camera-round and camera random intercepts for juveniles and additionally included camera-specific cyclic deviations for adult females, adult males, and pooled gazelles.

We standardized the empirical-Bayes wolf and jackal site-level activity indices to means of zero and standard deviations of one and retained them as continuous camera-level predictors. They therefore described spatial gradients among cameras rather than time-varying predator detections within a camera. We multiplied the standardized indices to form the wolf x jackal interaction. To separate changes in overall encounter frequency from changes in diel timing, each predator included in a candidate model entered both as a conventional parametric term and as a varying-coefficient cyclic smooth over sun time. This continuous varying-coefficient formulation followed the general approach used to evaluate predator-associated changes in diel activity by Cunningham et al. (2019). Our full model was

```text
cbind(success, failure) ~
  survey round + wolf_z * jackal_z +
  s(sun time, cyclic) +
  s(sun time, by = wolf_z, cyclic) +
  s(sun time, by = jackal_z, cyclic) +
  s(sun time, by = wolf_z * jackal_z, cyclic) +
  hierarchical camera terms.
```

The parametric predator terms controlled for differences in overall encounter frequency, whereas the predator-by-time smooths tested whether the shape of the 24-h activity distribution changed along each predator gradient. We considered five ecological structures: a baseline diel model, a wolf-only model containing wolf main and wolf-by-time effects, a jackal-only model containing jackal main and jackal-by-time effects, an additive wolf + jackal model, and the full wolf x jackal model. We used cyclic cubic splines with a maximum basis dimension of 12, shrinkage selection (`select = TRUE`), and discrete fast restricted maximum likelihood in `mgcv::bam()`. We selected the most parsimonious ecological structure within 2 conditional AIC units of the minimum. This procedure selected the wolf x jackal structure for juveniles and adult females, the additive wolf + jackal structure for adult males, and the wolf-only structure for pooled gazelles.

### Inference, predictions, and sensitivity analyses

We evaluated the predator-by-time terms with the global smooth tests returned by `mgcv`. A small P value for one of these terms indicated that the diel distribution changed along that predator gradient; it did not identify a single hour at which activity differed. We classified an effect as robust only if its global smooth test had P < 0.05 in the primary full model and in all five prespecified sensitivity analyses and if the term was included in the selected ecological structure. The sensitivity analyses used basis dimensions of 8 and 16, apparent-solar-clock time instead of sunrise-sunset standardized time, an alternative random-effect hierarchy, and the 30-min independent encounter histories. We examined convergence warnings, effective degrees of freedom, basis-dimension diagnostics, concurvity, deviance explained, Pearson dispersion, fitted probabilities, and lag-1 Pearson residual correlation within camera-rounds.

For interpretation, we predicted population-level activity at the four combinations of wolf and jackal activity at -1 and +1 SD while holding survey round at round 4 and excluding camera and camera-round random effects. We converted predicted encounter probabilities to normalized 24-h activity densities that integrated to one within each predator combination. We obtained pointwise 95% uncertainty intervals from 1,000 multivariate-normal simulations of the fitted coefficient vector using the unconditional covariance matrix. We also calculated high-minus-low predator difference curves; for the interaction, we calculated a difference-in-differences contrast across the four wolf-jackal combinations. Because these data were observational and predator indices represented spatial site-level gradients, we interpreted supported terms as associations with diel timing rather than causal behavioral responses.

## References

Cunningham, C. X., V. Scoleri, C. N. Johnson, L. A. Barmuta, and M. E. Jones. 2019. Temporal partitioning of activity: rising and falling top-predator abundance triggers community-wide shifts in diel activity. *Ecography* 42:2157-2168. https://doi.org/10.1111/ecog.04485

Iannarilli, F., B. D. Gerber, J. Erb, and J. R. Fieberg. 2025. A how-to guide for estimating animal diel activity using hierarchical models. *Journal of Animal Ecology* 94:182-194. https://doi.org/10.1111/1365-2656.14213

Meredith, M., M. Ridout, and L. A. D. Campbell. 2024. *overlap: Estimates of coefficient of overlapping for animal activity patterns*. R package version 0.3.9. https://doi.org/10.32614/CRAN.package.overlap

Nouvellet, P., G. S. A. Rasmussen, D. W. Macdonald, and F. Courchamp. 2012. Noisy clocks and silent sunrises: measurement methods of daily activity pattern. *Journal of Zoology* 286:179-184. https://doi.org/10.1111/j.1469-7998.2011.00864.x

Wood, S. N. 2017. *Generalized Additive Models: An Introduction with R*. Second edition. Chapman and Hall/CRC, Boca Raton, Florida, USA.
