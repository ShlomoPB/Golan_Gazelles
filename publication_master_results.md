# Publication Master Results

**Audit date:** 2026-08-18  
**Scope:** Authoritative results for manuscript writing. Values are traced to saved data, model objects, and consolidated tables. No gazelle-predator temporal-overlap result is included.

## Survey Data

| Quantity | Authoritative or audited value | Interpretation |
|---|---:|---|
| Cameras | 60 | Unique site IDs |
| Survey rounds | 7 | Nominal repeated surveys |
| Camera-round combinations | 420 | 60 x 7 |
| Positive-effort camera-rounds | 344 | One lacks dates |
| Missing-effort camera-rounds | 76 | Retained as `NA`, not zero |
| Recorded operational camera-days | 5,999 | Includes 29 days for camera 45, round 4, whose dates are missing |
| Date-resolved operational camera-days | 5,970 | 343 rounds used for event-to-active-round diel filtering |
| Prepared 30-min gazelle event rows | 3,932 | Before active-round filtering |
| Active-round gazelle event rows | 3,906 | 26 events excluded because they fell outside date-resolved active rounds |
| Final fawn count history | 199 | 31 cameras; summed classified animals |
| Final adult-female count history | 1,951 | 50 cameras; summed classified animals |
| Final adult-male count history | 2,355 | 52 cameras; summed classified animals |
| Beta-binomial sample | 221 camera-rounds | 50 cameras in 12 polygons |

The count histories and independent-event rows are different analysis units. For diel analyses, event-row sample sizes were 180 fawn events, 1,466 adult-female events, and 2,260 adult-male events after active-round filtering. Their corresponding summed classified-animal counts were 199, 1,931, and 2,346.

Subadult-related records were excluded from all focal classes. In the full export, 48 collapsed 30-min events met either the `Not Adult Male` or `Sub-adult` definition and summed to 62 event-level maximum individuals. The older value of 51 not-adult males was not reproducible from either current full-export file.

A single definitive count of unidentifiable gazelles cannot be reported until the curation rule is fixed. Strictly filtering the raw export to `Sex == "Unidentifiable"` gives 381 collapsed category-specific 30-min events and 443 event-level maximum individuals, but 77 of those events were labeled juvenile by life stage and may contribute to curated fawn classifications. Missing-sex, too-far, and unidentifiable-life-stage categories add further ambiguity.

## Spatial Models

Coefficients are posterior means on the log expected-abundance scale from separate negative-binomial `spAbundance::NMix` models. Intervals are 95% equal-tail credible intervals. Main effects are conditional on the other standardized predator index being zero.

| Class | Wolf coefficient | Jackal coefficient | Wolf x jackal coefficient |
|---|---:|---:|---:|
| Fawn | 0.27 [-0.27, 0.85] | -0.80 [-1.53, -0.13] | 0.50 [-0.21, 1.30] |
| Adult female | 0.54 [0.13, 0.97] | -0.72 [-1.16, -0.27] | 0.52 [0.03, 1.02] |
| Adult male | 0.36 [0.01, 0.72] | -0.61 [-0.95, -0.25] | 0.40 [0.03, 0.79] |

Fawn relative abundance declined along the jackal gradient at mean wolf index, while the fawn wolf effect and interaction were uncertain. Adult females and males showed positive wolf effects, negative jackal effects, and positive interactions. The positive interaction means that the negative jackal gradient weakened as wolf index increased and, equivalently, the positive wolf gradient strengthened as jackal index increased.

Paired posterior conditional abundance ratios make that interaction explicit. For females, the ratio for a 1-SD wolf increase changed from 1.02 [0.63, 1.67] at low jackal index to 1.71 [1.14, 2.63] at mean and 2.85 [1.37, 6.49] at high jackal index. For males, the corresponding ratios were 0.96 [0.63, 1.44], 1.42 [1.01, 2.06], and 2.11 [1.19, 4.02]. For fawns, all conditional wolf intervals included one. Full conditional results are in `publication_spabundance_wolf_jackal_relative_conditional_effect_ratios.csv`.

All nine focal slopes mixed adequately (`Rhat` 1.000-1.004; ESS 3,549-9,547). The abundance and detection intercepts did not converge (`Rhat` 1.21-2.93 and 1.23-3.49, respectively), so these results describe relative gradients rather than reliable absolute abundance. The matched `unmarked` sensitivity analysis recovered the same qualitative class pattern.

## Residual Class Correlations

The focused two-factor `lfMsNMix` model estimated positive residual latent relative-abundance correlations after accounting for open-space cover and wolf, jackal, and wolf x jackal terms.

| Class pair | Posterior correlation | 95% CrI |
|---|---:|---:|
| Fawn-adult male | 0.66 | 0.43-0.81 |
| Fawn-adult female | 0.86 | 0.51-1.00 |
| Adult male-adult female | 0.91 | 0.73-1.00 |

Free factor loadings had `Rhat` 1.006-1.009 and ESS 310-612. Site-grouped Freeman-Tukey posterior predictive Bayesian P values were 0.429, 0.397, and 0.330 for fawns, males, and females. Ungrouped checks indicated replicate-level lack of fit, and abundance/detection intercepts did not converge. Consequently, only the residual correlations are retained; they indicate shared unexplained spatial structure, not causal demographic interaction.

## Fawn Representation

The pooled camera-based fawn:adult-female ratio from the final count histories was `199 / 1,951 = 0.102`, equivalent to 10.2 fawns per 100 adult females. This is a descriptive demographic-composition index, not recruitment or survival.

The primary beta-binomial model used 221 informative camera-rounds from 50 cameras in 12 polygons and converged with a positive-definite Hessian.

| Focal predictor | Odds ratio | 95% CI | P | Supported at 95%? |
|---|---:|---:|---:|---|
| Late vs Early survey period | 1.57 | 1.01-2.44 | 0.045 | Yes |
| Open-space cover, 1 SD | 1.55 | 1.04-2.32 | 0.031 | Yes |
| Wolf index, 1 SD | 0.93 | 0.70-1.23 | 0.608 | No |
| Jackal index, 1 SD | 0.94 | 0.63-1.41 | 0.773 | No |

Fawn representation was higher in the Late survey period and at sites with greater open-space cover. Additive wolf and jackal effects were unsupported. Adding the wolf x jackal interaction did not improve fit (likelihood-ratio chi-square = 0.35, df = 1, P = 0.552; delta AIC = 1.65), and the interaction odds ratio was 0.86 [0.53, 1.39]. The additive model therefore remains primary.

## Diel Activity

Using solar elevation of at least zero degrees to define daylight, 149 of 180 fawn events (82.8%), 901 of 1,466 adult-female events (61.5%), and 1,286 of 2,260 adult-male events (56.9%) occurred during daylight. These are proportions of independent event rows, not abundance estimates.

Overall class-class diel overlap was highest between adult females and males and lower in comparisons involving fawns.

| Class pair | Overlap | Camera-bootstrap 95% interval |
|---|---:|---:|
| Fawn-adult female | 0.789 | 0.707-0.860 |
| Fawn-adult male | 0.743 | 0.661-0.783 |
| Adult female-adult male | 0.948 | 0.856-0.958 |

The median site indices were 0.495 for jackals and 0.556 for wolves. The median split assigned 14 cameras to low both, 16 to high jackal/low wolf, 16 to low jackal/high wolf, and 14 to high both.

| Class | Low both | High jackal / low wolf | Low jackal / high wolf | High both |
|---|---:|---:|---:|---:|
| Fawn | 67 | 3 | 72 | 38 |
| Adult female | 531 | 56 | 603 | 276 |
| Adult male | 794 | 151 | 771 | 544 |

No eligible fawn contrast was supported after within-class FDR adjustment; the high-jackal/low-wolf group was data deficient and excluded from formal tests. Five of six adult-female contrasts and four of six adult-male contrasts were supported. The smallest supported overlaps included adult females in high-jackal/low-wolf versus low-jackal/high-wolf contexts (`Dhat1` = 0.766, camera-bootstrap interval 0.566-0.827, FDR P = 0.003) and adult males in high-jackal/low-wolf versus high-both contexts (`Dhat4` = 0.820, 0.635-0.883, FDR P = 0.006).

These tests support context-associated redistribution of adult gazelle activity distributions. They do not test direct temporal avoidance of wolves or jackals, and the camera contexts are observational relative-index categories rather than predation-risk treatments.

## Excluded Results

Gazelle-wolf overlap, gazelle-jackal overlap, predator-overlay figures, domestic-dog recurrent-event models, wolf:jackal ratio explorations, occupancy models, season-specific N-mixture models, and failed model-derived fawn:female surfaces are excluded from this manuscript and Supporting Information.
