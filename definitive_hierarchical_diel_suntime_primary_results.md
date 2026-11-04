# Definitive hierarchical diel GAM results

## 1. Time-transformation audit

The previous scripts used apparent solar clock time: UTC hour + longitude/15 + the equation-of-time correction. This aligns solar noon but does not map sunrise and sunset to fixed positions; 06:00 and 18:00 were therefore six apparent-solar hours before and after noon, not observed sunrise and sunset.
The primary analysis instead used `overlap::sunTime()`, which maps date- and location-specific sunrise to 06:00 and sunset to 18:00. Camera-clock corrections were applied first. The former apparent-solar-clock analysis was retained as a sensitivity check.

## 2. Data-structure audit

The raw export contained 43,433 mountain-gazelle sighting rows. Binary camera-day-hour histories were created directly from corrected raw timestamps using mutually exclusive juvenile, adult-female, and adult-male classifications. The pooled response used any gazelle sighting, regardless of demographic classification. The pooled 30-minute reconstruction returned exactly 3,688 encounters.
Each positive-effort camera-round contributed `Days.Operational` trials to every sun-hour. Success was the number of operational dates with one or more encounters and failure was `Days.Operational - success`; zero-success hours were retained. Exact within-round malfunction dates were unavailable, so the authoritative round-level effort totals were used without inventing downtime dates.
Raw sun-time histories contained 181 juvenile, 1,572 adult-female, 2,327 adult-male, and 3,878 pooled positive camera-day-hours. The independent-history sensitivity contained 178, 1,456, 2,248, and 3,639 positive camera-day-hours, respectively, and was fitted with the same full-interaction hierarchy.

## 3. Selected model structures

- **Juvenile:** hierarchy = `camera_round`; selected ecological structure = wolf x jackal.
- **Adult female:** hierarchy = `site_shape`; selected ecological structure = wolf x jackal.
- **Adult male:** hierarchy = `site_shape`; selected ecological structure = additive wolf + jackal.
- **Pooled species:** hierarchy = `site_shape`; selected ecological structure = wolf-only.
For adult and pooled site-shape candidates, camera-specific cyclic deviations followed the Iannarilli factor-by-site construction. Because 60 separate smoothing parameters were computationally impractical, the documented `mgcv` `id` argument linked their penalties; cameras retained separate deviation coefficients but shared one smoothing parameter.
The candidate ecological structures were baseline-only, wolf-only, jackal-only, additive wolf + jackal, and wolf x jackal. Conditional AIC from the discrete fREML fits, with effective degrees of freedom and the `mgcv` smoothing-uncertainty correction, selected the most parsimonious structure within 2 AIC units of the minimum. For juveniles, this ranking was checked against exact ML and the full-interaction AIC differed by only 0.05.

## 4. Predator-by-time results and stability

- **Juvenile:** full-model smooth P values: wolf 0.085, jackal 0.966, wolf x jackal <0.001. Supported in all six specifications: interaction.
- **Adult female:** full-model smooth P values: wolf <0.001, jackal 0.508, wolf x jackal 0.116. Supported in all six specifications: wolf.
- **Adult male:** full-model smooth P values: wolf 0.056, jackal 0.413, wolf x jackal 0.002. Supported in all six specifications: none.
- **Pooled species:** full-model smooth P values: wolf <0.001, jackal 0.537, wolf x jackal 0.001. Supported in all six specifications: wolf.

A result was called robust only when its smooth was supported in the primary sun-time full model and in the k = 8, k = 16, apparent-solar-time, alternate-hierarchy, and independent-30-minute-history sensitivities, and when the corresponding effect occurred in the selected ecological structure. This prevents significance under one specification, response history, or unselected model from driving the ecological conclusion.

## 5. Predicted curve interpretation

Juveniles showed an interaction rather than additive main effects. Predicted activity peaked near 08:00 when predator indices were concordant (low-low or high-high) and near 17:15 in the two discordant contexts. The interaction difference-in-differences was positive from approximately 06:51 to 07:51 and negative from 15:27 to 18:03 using pointwise 95% simulation intervals.
For adult females, increasing wolf activity shifted the predicted peak from approximately 12:15 to 17:33. The wolf smooth was globally supported, but the pointwise high-minus-low difference interval did not exclude zero at any single hour; hour-specific claims should therefore remain descriptive. Male curves were effectively identical and peaked near 17:45 in all retained-model predictions.
For pooled gazelles, increasing wolf activity shifted the predicted daylight share from 0.584 to 0.610 and reduced mean distance from solar noon from 5.61 to 5.44 hours, while the peak remained near sunset (17:45 versus 17:39). This pooled wolf gradient was supported under both raw and independent binary histories.

## 6. Comparison with previous analyses

The earlier negative-binomial class-count GAMM agreed on the juvenile wolf x jackal interaction and female wolf gradient. It also gave a nominal male interaction, but the definitive binomial analysis selected the additive wolf + jackal structure and no male timing smooth; the male interaction is therefore sensitivity-only. Median-split kernels found no supported juvenile context contrast but several female and male contrasts, a broader categorical pattern that does not isolate continuous predator effects.
The older pooled negative-binomial model found no stable predator timing effect, whereas the definitive pooled binomial model retained a stable wolf gradient. The pooled wolf result is therefore robust within the hourly-occurrence framework but sensitive to modeling event counts versus hourly encounter opportunities. It should be presented as a secondary species-level benchmark rather than the main ecological headline.

## 7. Recommendation

Use the hierarchical binomial GAM as the primary diel analysis. It supports a robust wolf x jackal timing interaction for juveniles, a robust wolf-associated timing gradient for adult females, no retained predator-associated timing change for adult males, and a response-formulation-sensitive wolf gradient in the pooled species benchmark. The pooled effect is most consistent with the female response among classified animals, but it is not a formal demographic decomposition because the pooled response also includes unclassified and subadult gazelles.
Interpret the juvenile pattern as a nonlinear association with two-predator context: concordant low-low and high-high sites were morning-dominant, whereas discordant predator contexts shifted predicted activity later in the day. Interpret the female and pooled wolf gradients as redistribution of detections across the diel cycle, not evidence of causal predator avoidance. Do not retain the earlier count GAMM or median-split kernel analysis as co-primary inference.
Residual dependence remains, especially for juveniles (mean within-camera-round lag-1 Pearson correlation 0.71; k-index 0.86, P = 0.047), and all retained fits were underdispersed. The juvenile effect was nevertheless stable across basis size, hierarchy, time transformation, and raw versus independent histories. These limitations should be reported and causal language avoided.

## References

Iannarilli, F., Gerber, B. D., Erb, J., & Fieberg, J. R. (2025). A 'how-to' guide for estimating animal diel activity using hierarchical models. *Journal of Animal Ecology*, 94, 182-194. https://doi.org/10.1111/1365-2656.14213

Cunningham, C. X., Scoleri, V., Johnson, C. N., Barmuta, L. A., & Jones, M. E. (2019). Temporal partitioning of activity: rising and falling top-predator abundance triggers community-wide shifts in diel activity. *Ecography*, 42, 2157-2168. https://doi.org/10.1111/ecog.04485
