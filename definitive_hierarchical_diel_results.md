# Definitive hierarchical diel GAM results

## 1. Primary time scale and effort

The primary analysis used equal-duration apparent-solar-clock hours, calculated as UTC time plus longitude/15 and the date-specific equation-of-time correction. This additive transformation aligned solar noon at 12:00 while preserving true one-hour bin durations. The sunrise-sunset-standardized `overlap::sunTime()` analysis was retained as the prespecified time-scale sensitivity.

The raw export contained 43,433 mountain-gazelle sighting rows. Binary camera-day-hour histories were created directly from corrected raw timestamps using mutually exclusive juvenile, adult-female, and adult-male classifications. The pooled response used any gazelle sighting, regardless of demographic classification. The pooled 30-minute reconstruction returned exactly 3,688 encounters.
Each positive-effort camera-round contributed the recorded `Days.Operational` trials to every apparent-solar hour. Success was the number of operational dates with one or more encounters and failure was `Days.Operational - success`; zero-success hours were retained. Camera 45 round 2 contributed its recorded 21 operational days across a 32-day span. Exact missing dates and hourly outage information were unavailable and were not invented.
The primary raw histories contained 183 juvenile, 1576 adult female, 2307 adult male, 3857 pooled species positive camera-day-hours. The apparent-solar independent-history sensitivity contained 177, 1451, 2245, 3625 positive camera-day-hours in the same class order.

## 2. Model selection

The final column in each comparison reports competitor cAIC minus selected-model cAIC; positive values favor the selected model. A negative value indicates that the simpler selected structure was retained under the prespecified within-2-cAIC parsimony rule.

| Response | Selected hierarchy | Nearest hierarchy competitor (Delta cAIC) | Selected ecology | Nearest ecological competitor (Delta cAIC) |
|---|---|---|---|---|
| Juvenile | camera + camera-round random intercepts | camera random intercept (51.339) | wolf x jackal | additive wolf + jackal (9.666) |
| Adult female | camera-specific cyclic deviation + camera and camera-round random intercepts | camera + camera-round random intercepts (463.342) | wolf x jackal | additive wolf + jackal (2.495) |
| Adult male | camera-specific cyclic deviation + camera and camera-round random intercepts | camera + camera-round random intercepts (276.134) | wolf x jackal | additive wolf + jackal (3.262) |
| Pooled species | camera-specific cyclic deviation + camera and camera-round random intercepts | camera + camera-round random intercepts (570.719) | wolf-only | additive wolf + jackal (0.196) |

For adult and pooled site-shape candidates, camera-specific cyclic deviations followed the Iannarilli factor-by-site construction. Because 60 separate smoothing parameters were computationally impractical, the documented `mgcv` `id` argument linked their penalties; cameras retained separate deviation coefficients but shared one smoothing parameter.
The candidate ecological structures were baseline-only, wolf-only, jackal-only, additive wolf + jackal, and wolf x jackal. Conditional AIC was calculated with `AIC()` from the discrete fREML fits, using effective degrees of freedom and the `mgcv` smoothing-uncertainty correction. We selected the least complex structure within 2 cAIC units of the minimum.

## 3. Predator-by-time tests and robustness

The table reports all predator-by-time smooths from the primary full interaction model. Robustness required support in all seven specifications and retention of the corresponding effect in the selected ecological structure.

| Response | Predator-by-time effect | EDF | P | Retained | Supported specifications | Classification |
|---|---|---:|---:|---|---:|---|
| Juvenile | Wolf x time | 1.188 | 0.059 | yes | 0/7 | unsupported |
| Juvenile | Jackal x time | 0.000 | 0.900 | yes | 0/7 | unsupported |
| Juvenile | Wolf x jackal x time | 2.067 | <0.001 | yes | 7/7 | robust |
| Adult female | Wolf x time | 3.916 | <0.001 | yes | 7/7 | robust |
| Adult female | Jackal x time | 0.000 | 0.432 | yes | 2/7 | unsupported |
| Adult female | Wolf x jackal x time | 0.001 | 0.031 | yes | 6/7 | sensitivity-dependent |
| Adult male | Wolf x time | 0.000 | 0.037 | yes | 4/7 | sensitivity-dependent |
| Adult male | Jackal x time | 0.000 | 0.521 | yes | 1/7 | unsupported |
| Adult male | Wolf x jackal x time | 2.991 | <0.001 | yes | 6/7 | sensitivity-dependent |
| Pooled species | Wolf x time | 3.041 | <0.001 | yes | 7/7 | robust |
| Pooled species | Jackal x time | 0.000 | 0.606 | no | 2/7 | not retained |
| Pooled species | Wolf x jackal x time | 2.595 | <0.001 | no | 7/7 | not retained |

The seven specifications were the primary apparent-solar full model, k = 8, k = 16, sun-time, alternate-hierarchy, apparent-solar independent-30-minute-history, and the cached sun-time plus camera-only alternate-hierarchy fit. This prevents significance under one basis size, hierarchy, time transformation, response history, unselected model, or the prespecified combined time-scale and hierarchy sensitivity from driving inference.

## 4. Predicted curves

- **Juvenile:** predicted apparent-solar peak times across the four wolf-jackal combinations were 8.85, 16.75, 16.75, 8.55 h; robust predator-by-time effect(s): interaction.
- **Adult female:** predicted apparent-solar peak times across the four wolf-jackal combinations were 17.35, 17.35, 17.25, 17.25 h; robust predator-by-time effect(s): wolf.
- **Adult male:** predicted apparent-solar peak times across the four wolf-jackal combinations were 17.45, 17.55, 17.55, 17.45 h; robust predator-by-time effect(s): none.
- **Pooled species:** predicted apparent-solar peak times across the four wolf-jackal combinations were 17.55, 17.55, 17.45, 17.45 h; robust predator-by-time effect(s): wolf.

Juvenile activity retained the strong two-predator interaction: concordant low-low and high-high contexts peaked in the morning (08:33-08:51), whereas both discordant contexts peaked in late afternoon (16:45). The interaction difference-in-differences interval was positive from 06:57 to 10:03 and negative from 14:45 to 19:09.
Adult-female peak timing changed little (17:15-17:21), but the robust wolf-by-time smooth described lower midday and greater late-afternoon activity as wolf activity increased. The high-minus-low wolf contrast did not exclude zero at any individual hour, so hour-specific differences remain descriptive.
Adult males showed a global wolf-jackal-by-time interaction in the primary apparent-solar model, but this effect was sensitivity-dependent rather than robust because the cached sun-time plus camera-only sensitivity was unsupported (EDF < 0.001, P = 0.478). Concordant predator contexts had approximately 64.7% of normalized activity between apparent-solar 06:00 and 18:00, compared with approximately 59.6% in discordant contexts, and mean distance from solar noon differed by approximately 0.40 h. Peak time changed by only 0.10 h, and no pointwise male interaction interval excluded zero; these primary-model patterns are descriptive and do not constitute a robust male predator-associated timing effect.
The pooled wolf-only model retained a robust wolf-by-time smooth and a nearly unchanged late-afternoon peak (17:33 versus 17:27). Because the pooled response includes unclassified and subadult gazelles, it remains a secondary species-level benchmark rather than a demographic decomposition.

The curves are normalized within predator context and therefore describe redistribution across the 24-h cycle rather than absolute encounter probability. Apparent solar time aligns solar noon but does not set sunrise and sunset to 06:00 and 18:00.

## 5. Comparison with the previous sun-time-primary analysis

| Response | Previous hierarchy | Apparent-solar hierarchy | Previous ecology | Apparent-solar ecology | Previous robust effects | Apparent-solar robust effects |
|---|---|---|---|---|---|---|
| Juvenile | camera + camera-round random intercepts | camera + camera-round random intercepts | wolf x jackal | wolf x jackal | interaction | interaction |
| Adult female | camera-specific cyclic deviation + camera and camera-round random intercepts | camera-specific cyclic deviation + camera and camera-round random intercepts | wolf x jackal | wolf x jackal | wolf | wolf |
| Adult male | camera-specific cyclic deviation + camera and camera-round random intercepts | camera-specific cyclic deviation + camera and camera-round random intercepts | additive wolf + jackal | wolf x jackal | none | none |
| Pooled species | camera-specific cyclic deviation + camera and camera-round random intercepts | camera-specific cyclic deviation + camera and camera-round random intercepts | wolf-only | wolf-only | wolf | wolf |

The juvenile, adult-female, adult-male, and pooled robustness conclusions were unchanged from the sun-time-primary audit. Although the adult-male full interaction model was favored by 3.262 cAIC under apparent solar time and its interaction smooth was supported in all six one-at-a-time specifications, the cached sun-time plus camera-only alternate-hierarchy sensitivity was unsupported (P = 0.478). Under the prespecified requirement for support across all sensitivities, the male interaction is therefore sensitivity-dependent rather than robust.
Changes in hierarchy or selected candidate structure do not by themselves constitute ecological effects; focal inference follows the prespecified robustness rule and the predicted effect magnitude.

## 6. Recommendation

Use the equal-duration apparent-solar hierarchical binomial GAM as the definitive diel framework. The focal conclusions remain a robust juvenile wolf-jackal temporal interaction, a robust female wolf-associated timing gradient, and no robust predator-associated temporal effect for adult males. The primary adult-male interaction may be reported only as sensitivity-dependent because the sun-time plus camera-only sensitivity was unsupported. Interpret all predator-by-time terms as observational associations with redistribution of detections across the diel cycle, not causal predator avoidance. Retain sun time as a prespecified sensitivity rather than as the primary response.

## References

Iannarilli, F., Gerber, B. D., Erb, J., & Fieberg, J. R. (2025). A 'how-to' guide for estimating animal diel activity using hierarchical models. *Journal of Animal Ecology*, 94, 182-194. https://doi.org/10.1111/1365-2656.14213

Cunningham, C. X., Scoleri, V., Johnson, C. N., Barmuta, L. A., & Jones, M. E. (2019). Temporal partitioning of activity: rising and falling top-predator abundance triggers community-wide shifts in diel activity. *Ecography*, 42, 2157-2168. https://doi.org/10.1111/ecog.04485
