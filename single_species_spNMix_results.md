# Single-species spatial N-mixture models with `spAbundance`

## Model setup

- Model type: `spAbundance::spNMix`
- Family: Negative binomial
- Spatial structure: exponential NNGP with 10 neighbors
- Detection model: `Path.Type + Camera.Sensitivity + Area.Behind.Camera + Vegetation.Structure`
- Wolf covariate: `Null.Wolf.Abundance.z`
- Habitat covariate chosen from the occupancy analysis:
  - Fawn: `Open.Spaces.z`
  - Male: `Built.and.Roads.z`
  - Female: `Nature.Reserves.z`

## Abundance results

### Fawn

File:

- `fawn_spNMix_Open_Spaces_Null_Wolf_Abundance_det_design_NB_summary.txt`

Abundance coefficients on the log scale:

| Term | Mean | 95% CrI | Rhat | ESS |
|---|---:|---:|---:|---:|
| `Open.Spaces.z` | 1.043 | [0.457, 1.786] | 1.004 | 269 |
| `null.jackal.abundance.z` | -0.629 | [-1.291, -0.023] | 1.037 | 281 |
| `Null.Wolf.Abundance.z` | 0.319 | [-0.178, 0.821] | 1.005 | 277 |

Interpretation:

- Fawn abundance increases with open space.
- Jackal abundance shows a negative association.
- The wolf effect is uncertain.

WAIC:

- 629.32

### Male

File:

- `male_spNMix_Built_and_Roads_Null_Wolf_Abundance_det_design_NB_summary.txt`

Abundance coefficients on the log scale:

| Term | Mean | 95% CrI | Rhat | ESS |
|---|---:|---:|---:|---:|
| `Built.and.Roads.z` | -0.849 | [-1.299, -0.400] | 1.138 | 122 |
| `null.jackal.abundance.z` | -0.844 | [-1.221, -0.441] | 1.104 | 273 |
| `Null.Wolf.Abundance.z` | 0.366 | [0.044, 0.697] | 1.002 | 331 |

Interpretation:

- Male abundance decreases with infrastructure.
- Jackal abundance is also negatively associated with male abundance.
- The wolf effect is positive in this fit.

### Female

File:

- `female_spNMix_Nature_Reserves_Null_Wolf_Abundance_det_design_NB_summary.txt`

Abundance coefficients on the log scale:

| Term | Mean | 95% CrI | Rhat | ESS |
|---|---:|---:|---:|---:|
| `Nature.Reserves.z` | -0.289 | [-0.960, 0.367] | 1.002 | 117 |
| `null.jackal.abundance.z` | -0.637 | [-1.125, -0.149] | 1.177 | 185 |
| `Null.Wolf.Abundance.z` | 0.478 | [0.090, 0.875] | 1.039 | 266 |

Interpretation:

- Female abundance shows a negative jackal association.
- The reserve effect is not clearly different from 0 in this N-mixture fit.
- The wolf effect is positive in this fit.

## Comparison with occupancy

- The strongest cross-model consistency is the negative jackal association: it remains negative for fawns, males, and females.
- Fawn open-space use and male infrastructure avoidance are also recovered in the abundance models.
- The female reserve effect is weaker in the N-mixture model than it was in the occupancy model.
- Wolf effects are more positive in the abundance models than in the occupancy models, especially for adults.

## Important cautions

- These are spatial negative binomial N-mixture models, which are difficult to fit because both spatial random effects and NB overdispersion can absorb extra variation.
- Intercept mixing was poor in all three fits (`Rhat` for the intercept is well above 1.1).
- Male and female abundance-side convergence is only moderate for some coefficients, especially the intercept and jackal effect.
- `waicAbund()` completed for the fawn model but did not finish within a long timeout for the male and female models.
- These results are usable as exploratory ecological evidence, but I would not treat them as publication-grade final abundance inference without either:
  - longer chains and stronger convergence checks, or
  - simpler models such as spatial Poisson N-mixture fits for comparison.

## Simpler detection model comparison

I re-fit the same three spatial NB N-mixture models with a simpler detection model:

- `log_days_operational.z + cattle_sightings.z`

Summary files:

- `fawn_spNMix_Open_Spaces_Null_Wolf_Abundance_det_effort_NB_summary.txt`
- `male_spNMix_Built_and_Roads_Null_Wolf_Abundance_det_effort_NB_summary.txt`
- `female_spNMix_Nature_Reserves_Null_Wolf_Abundance_det_effort_NB_summary.txt`

### Convergence comparison

- Fawn: abundance coefficients were similar under both detection models. The simpler detection model improved mixing for `kappa` (`Rhat` 1.08 to 1.04) and `phi` (`ESS` 801 to 1146), but the intercept mixed even worse (`Rhat` 1.18 to 1.72).
- Male: the simpler detection model improved convergence clearly. `Built.and.Roads.z` improved from `Rhat` 1.14 to 1.06, the intercept from 1.53 to 1.16, `sigma.sq` from 1.10 to 1.03, and `phi` from 1.37 to 1.03.
- Female: the simpler detection model improved jackal convergence (`Rhat` 1.18 to 1.01) and `sigma.sq` (`Rhat` 1.10 to 1.04), but the intercept remained poor and `phi` mixed worse (`Rhat` 1.05 to 1.23).

### Ecological stability

- The main abundance-side conclusions did not change materially.
- Fawns still showed a positive `Open.Spaces` effect and a negative jackal effect.
- Males still showed a negative `Built.and.Roads` effect, a negative jackal effect, and a positive wolf effect.
- Females still showed a negative jackal effect and a weak/uncertain reserve effect.

### Practical recommendation

- For the male model, the simpler effort-based detection model is preferable because it improves convergence without changing the ecological interpretation.
- For fawns and females, simplifying detection helps some parameters but does not fully solve mixing problems. If these models are to be used for inference, the next better simplification is likely on the abundance side or by switching from NB to Poisson for comparison.

## Spatial NB versus spatial Poisson

I re-fit the same three models with the simpler effort-based detection model and a spatial Poisson observation model instead of spatial negative binomial:

- `fawn_spNMix_Open_Spaces_Null_Wolf_Abundance_det_effort_Poisson_summary.txt`
- `male_spNMix_Built_and_Roads_Null_Wolf_Abundance_det_effort_Poisson_summary.txt`
- `female_spNMix_Nature_Reserves_Null_Wolf_Abundance_det_effort_Poisson_summary.txt`

### Convergence comparison

- Fawn: Poisson did not improve the fit overall. The abundance intercept improved somewhat (`Rhat` 1.72 to 1.32), and `phi` mixed well in both models, but the core ecological coefficients got worse under Poisson: `Open.Spaces.z` worsened from `Rhat` 1.01 to 1.19, jackal from 1.02 to 1.34, and wolf from 1.01 to 1.34.
- Male: Poisson was dramatically worse. The abundance intercept jumped from `Rhat` 1.16 to 9.86, `Built.and.Roads.z` from 1.06 to 2.49, `sigma.sq` from 1.03 to 4.33, and `phi` from 1.03 to 5.93.
- Female: Poisson was also dramatically worse. The abundance intercept jumped from `Rhat` 1.64 to 10.83, `Nature.Reserves.z` from 1.04 to 10.88, jackal from 1.01 to 2.50, `sigma.sq` from 1.04 to 5.96, and `phi` from 1.23 to 6.08.

### Ecological stability under Poisson

- Fawn: the positive open-space effect remained, but jackal and wolf effects became much less trustworthy because the chains mixed poorly.
- Male: the signs of the abundance effects were similar to the NB fit, but the Poisson model did not converge well enough to support inference.
- Female: the Poisson fit produced the same general negative jackal / positive wolf direction, but convergence was too poor for defensible interpretation.

### Practical recommendation after Poisson comparison

- Spatial Poisson does not solve the convergence problem here, and for males and females it is clearly worse than spatial NB.
- The best current abundance model among these fits is still the male spatial NB model with the simpler effort-based detection submodel.
- For fawns and females, the remaining instability appears to come from the abundance-side model structure rather than from detection-model complexity alone.
- For publication-grade analysis, the next better comparison would be either a non-spatial N-mixture model or a spatial model with stronger prior regularization and longer chains.

## Spatial versus non-spatial NB N-mixture

I then fit matched non-spatial negative-binomial N-mixture models with `spAbundance::NMix`, keeping the same class-specific habitat covariates and the same simpler detection model:

- `fawn_NMix_Open_Spaces_Null_Wolf_Abundance_det_effort_NB_summary.txt`
- `male_NMix_Built_and_Roads_Null_Wolf_Abundance_det_effort_NB_summary.txt`
- `female_NMix_Nature_Reserves_Null_Wolf_Abundance_det_effort_NB_summary.txt`

### Convergence comparison

- Fawn: removing the spatial process substantially improved the ecological coefficients. `Open.Spaces.z` improved from `Rhat` 1.01 to 1.01 with `ESS` 284 to 911, jackal from `Rhat` 1.02 / `ESS` 345 to `Rhat` 1.01 / `ESS` 869, and wolf from `Rhat` 1.01 / `ESS` 423 to `Rhat` 1.01 / `ESS` 1426. The abundance intercept remained poor in both models.
- Male: the non-spatial model mixed much better for the abundance coefficients. `Built.and.Roads.z` improved from `Rhat` 1.06 / `ESS` 216 to `Rhat` 1.00 / `ESS` 2430, jackal from `Rhat` 1.02 / `ESS` 270 to `Rhat` 1.00 / `ESS` 3163, and wolf from `Rhat` 1.01 / `ESS` 330 to `Rhat` 1.00 / `ESS` 3513. The abundance intercept was still not good, but the core slope terms were much more stable.
- Female: the non-spatial model also mixed much better for the abundance coefficients. `Nature.Reserves.z` improved from `Rhat` 1.04 / `ESS` 163 to `Rhat` 1.00 / `ESS` 2859, jackal from `Rhat` 1.01 / `ESS` 238 to `Rhat` 1.00 / `ESS` 2435, and wolf from `Rhat` 1.02 / `ESS` 354 to `Rhat` 1.00 / `ESS` 2660. As with the other classes, the abundance intercept remained weakly mixed.

### Ecological stability

- Fawn: the non-spatial model kept the same story as the spatial NB fit. `Open.Spaces.z` remained positive (`0.93`, 95% CrI `[0.43, 1.47]`), jackal remained negative (`-0.59`, `[-1.21, -0.03]`), and the wolf effect remained uncertain (`0.24`, `[-0.18, 0.66]`).
- Male: the non-spatial model also kept the same directional pattern. `Built.and.Roads.z` remained negative (`-0.59`, `[-0.99, -0.18]`), jackal remained negative (`-0.75`, `[-1.13, -0.37]`), and the wolf effect stayed positive but weaker (`0.25`, `[-0.01, 0.51]`).
- Female: the non-spatial model changed the female reserve effect from weakly negative to essentially null (`0.01`, `[-0.39, 0.41]`), while the jackal effect remained negative (`-0.81`, `[-1.30, -0.31]`) and the wolf effect remained borderline positive (`0.39`, `[-0.02, 0.80]`).

### Practical recommendation after non-spatial comparison

- If the goal is stable inference on habitat and predator effects, the non-spatial NB models are currently more defensible than the spatial NB or spatial Poisson versions.
- The main ecological conclusions that survive this comparison are the positive fawn association with open space, the negative male association with infrastructure, and the negative jackal association across all three classes.
- The female reserve effect does not survive strongly once the spatial random effect is removed, so that result should be treated cautiously.
- For research-grade reporting, I would present the non-spatial NB models as the more stable abundance analysis and treat the spatial N-mixture fits as sensitivity analyses rather than the primary inference.

## Expanded non-spatial NB models with boar and human abundance

I then expanded the non-spatial NB models to include two additional site covariates in the abundance component:

- `null.boar.abundance.z`
- `human.abundance.z`

The abundance formula for each class was therefore:

- class-specific habitat term + `null.jackal.abundance.z` + `null.boar.abundance.z` + `human.abundance.z` + `Null.Wolf.Abundance.z`

Summary files:

- `fawn_NMix_Open_Spaces_Null_Wolf_Abundance_det_effort_NB_summary.txt`
- `male_NMix_Built_and_Roads_Null_Wolf_Abundance_det_effort_NB_summary.txt`
- `female_NMix_Nature_Reserves_Null_Wolf_Abundance_det_effort_NB_summary.txt`

### Ecological results

- Fawn: `Open.Spaces.z` remained clearly positive (`1.05`, 95% CrI `[0.52, 1.61]`). The jackal effect weakened slightly and now overlaps 0 (`-0.54`, `[-1.18, 0.05]`). Boar was positive but uncertain (`0.27`, `[-0.33, 0.86]`), human abundance was negative but uncertain (`-0.33`, `[-0.90, 0.26]`), and the wolf effect remained uncertain (`0.25`, `[-0.20, 0.71]`).
- Male: the earlier pattern remained very stable. `Built.and.Roads.z` stayed clearly negative (`-0.64`, `[-1.03, -0.24]`), jackal stayed clearly negative (`-0.73`, `[-1.12, -0.35]`), boar was positive but uncertain (`0.22`, `[-0.15, 0.59]`), human abundance was near 0 (`-0.13`, `[-0.45, 0.22]`), and the wolf effect was borderline positive (`0.26`, `[-0.00, 0.53]`).
- Female: the reserve effect remained weak and uncertain (`-0.11`, `[-0.57, 0.35]`), jackal remained clearly negative (`-0.82`, `[-1.36, -0.30]`), boar was positive but uncertain (`0.29`, `[-0.28, 0.86]`), human abundance was near 0 (`-0.12`, `[-0.60, 0.40]`), and the wolf effect became clearly positive (`0.42`, `[0.01, 0.84]`).

### Convergence

- The abundance-side slope terms remained well mixed despite the extra covariates.
- Fawn: all abundance slopes had `Rhat` between `1.00` and `1.03`, though the intercept remained weakly mixed (`Rhat` `1.45`).
- Male: all abundance slopes had `Rhat` between `1.00` and `1.00`, with the intercept improved relative to earlier fits but still not ideal (`Rhat` `1.09`).
- Female: all abundance slopes had `Rhat` between `1.00` and `1.01`, with the intercept still weaker (`Rhat` `1.24`).

### Interpretation

- Adding boar and human abundance does not overturn the main non-spatial NB results.
- The strongest stable signal remains the negative jackal association across all three classes.
- Fawn abundance still tracks open space and male abundance still declines with infrastructure.
- Neither boar nor human abundance shows a clearly supported independent effect in these models.

## Focused human and boar candidate set

To isolate the effects of human and boar activity on each gazelle class, I fit a reduced non-spatial NB candidate set for each class using only:

- habitat only
- habitat + `human.abundance.z`
- habitat + `null.boar.abundance.z`
- habitat + `human.abundance.z` + `null.boar.abundance.z`

Jackal and wolf covariates were excluded from this candidate set so the human and boar terms were not competing with predator covariates for the same site-level variation.

### Main result

- Fawn: `human.abundance.z` was consistently negative but uncertain (`-0.39`, 95% CrI about `[-0.95, 0.15]`), and `null.boar.abundance.z` was consistently positive but uncertain (`0.26` to `0.27`, CrI overlapping 0). `Open.Spaces.z` remained clearly positive in every candidate model (`0.92` to `1.01`).
- Male: `human.abundance.z` was again weakly negative but uncertain (`-0.19` to `-0.20`, CrI overlapping 0), and `null.boar.abundance.z` was weakly positive but uncertain (`0.12` to `0.15`, CrI overlapping 0). `Built.and.Roads.z` remained negative in every candidate model, but the interval slightly overlapped 0 in this reduced set.
- Female: `human.abundance.z` remained weakly negative but uncertain (`-0.20`, CrI overlapping 0), and `null.boar.abundance.z` remained weakly positive but uncertain (`0.14`, CrI overlapping 0). `Nature.Reserves.z` remained positive but uncertain in every candidate model.

### Interpretation

- Across all three classes, there is no clear evidence from these N-mixture models that human abundance or boar abundance has an independent effect on gazelle abundance.
- The directions are consistent across classes: human tends to be negative and boar tends to be positive, but both are weak and imprecise.
- The strongest class-specific habitat terms still dominate the signal relative to these human and boar covariates.

Comparison table files:

- `human_boar_class_comparison_table.md`
- `human_boar_class_comparison_table.csv`
