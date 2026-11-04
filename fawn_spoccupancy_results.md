# Fawn occupancy with `spOccupancy`

## Response and design

- Response: juvenile gazelle counts converted to detection/non-detection (`count > 0`) for 60 camera sites across 7 survey rounds.
- Occupancy covariates: `Nature.Reserves.z`, `null.jackal.abundance.z`, and one wolf covariate.
- Detection covariates: `log_days_operational.z`, `cattle_sightings.z`.
- Spatial structure: `spOccupancy::spPGOcc` with an exponential covariance function and projected camera coordinates transformed to EPSG:2039 and stored in kilometers.

## Fitted models

Two spatial occupancy models were fit with identical settings except for the wolf predictor:

1. `Null.Wolf.Abundance.z`
2. `Wolf.pack.abundance.z`

Common MCMC settings:

- 3 chains
- 400 batches x 25 iterations per batch
- 10,000 samples per chain
- burn-in = 1,000
- thinning = 5
- total posterior samples retained = 5,400

## Posterior summaries

### Model using `Null.Wolf.Abundance.z`

| Term | Mean | 95% CrI | Interpretation |
|---|---:|---:|---|
| `Nature.Reserves.z` | 0.064 | [-0.999, 1.290] | No clear reserve effect |
| `null.jackal.abundance.z` | -1.318 | [-2.806, -0.252] | Clear negative association with fawn occupancy |
| `Null.Wolf.Abundance.z` | 0.482 | [-0.425, 1.599] | Positive posterior mean, but interval overlaps 0 |

Detection:

| Term | Mean | 95% CrI | Interpretation |
|---|---:|---:|---|
| `log_days_operational.z` | 0.283 | [-0.028, 0.602] | Weak positive detection effect |
| `cattle_sightings.z` | -0.101 | [-0.389, 0.163] | No clear detection effect |

Fit statistic:

- WAIC = 329.06

### Model using `Wolf.pack.abundance.z`

| Term | Mean | 95% CrI | Interpretation |
|---|---:|---:|---|
| `Nature.Reserves.z` | 0.110 | [-0.974, 1.386] | No clear reserve effect |
| `null.jackal.abundance.z` | -1.346 | [-2.951, -0.247] | Clear negative association with fawn occupancy |
| `Wolf.pack.abundance.z` | 0.052 | [-0.965, 1.093] | No clear wolf-pack effect |

Detection:

| Term | Mean | 95% CrI | Interpretation |
|---|---:|---:|---|
| `log_days_operational.z` | 0.284 | [-0.027, 0.607] | Weak positive detection effect |
| `cattle_sightings.z` | -0.103 | [-0.381, 0.162] | No clear detection effect |

Fit statistic:

- WAIC = 329.41

## Ecological reading

- Across both wolf formulations, the most consistent signal is a negative association between jackal abundance and fawn occupancy.
- Nature reserves do not show a clear positive or negative relationship with fawn occupancy in this model.
- The `Null.Wolf.Abundance` version has a slightly lower WAIC than the `Wolf.pack.abundance` version, but the difference is small and both wolf effects remain uncertain because their 95% credible intervals overlap 0.

## Cautions

- The juvenile detection-history missingness does not perfectly match `Days.Operational` missingness; the script warns about this and proceeds.
- Spatial variance mixed less efficiently than the regression coefficients (`sigma.sq` ESS is comparatively low), so a publication-grade final analysis should consider longer chains and posterior predictive checks before drawing strong inference about the spatial random effect.

## Alternate detection model

I also fit an alternate spatial occupancy model using `Null.Wolf.Abundance.z` and replacing the effort-based detection model with:

- `Path.Type`
- `Camera.Sensitivity`
- `Area.Behind.Camera`
- `Vegetation.Structure`

Saved outputs:

- `fawn_spPGOcc_Null_Wolf_Abundance_det_design_fit.rds`
- `fawn_spPGOcc_Null_Wolf_Abundance_det_design_summary.txt`

### Occupancy under the alternate detection model

| Term | Mean | 95% CrI | Interpretation |
|---|---:|---:|---|
| `Nature.Reserves.z` | 0.095 | [-0.944, 1.280] | No clear reserve effect |
| `null.jackal.abundance.z` | -1.296 | [-2.718, -0.204] | Clear negative association with fawn occupancy |
| `Null.Wolf.Abundance.z` | 0.507 | [-0.409, 1.633] | Positive posterior mean, but still uncertain |

### Detection under the alternate detection model

- Detection was higher on `Dirt path less than 1m wide` and `Dirt path more than 1m wide`.
- Detection was lower when the area behind the camera was `Impassable due to topography`.
- Detection was higher in `Shrubs < 2 m`.
- Detection was lower in `Shrubs > 2 m`, although that interval slightly overlaps 0.

### Model comparison

- Effort-based detection model WAIC: 329.06
- Alternate detection model WAIC: 317.47

The richer detection model improves predictive fit substantially, but the main occupancy conclusion does not change: fawn occupancy remains most consistently negatively associated with jackal abundance, not clearly associated with nature reserves, and only uncertainly associated with the wolf covariate.

## Same model for males and females

I fit the same spatial occupancy model for males and females using:

- occupancy: `Nature.Reserves.z + null.jackal.abundance.z + Null.Wolf.Abundance.z`
- detection: `Path.Type + Camera.Sensitivity + Area.Behind.Camera + Vegetation.Structure`

Saved outputs:

- `male_spPGOcc_Null_Wolf_Abundance_det_design_fit.rds`
- `male_spPGOcc_Null_Wolf_Abundance_det_design_summary.txt`
- `female_spPGOcc_Null_Wolf_Abundance_det_design_fit.rds`
- `female_spPGOcc_Null_Wolf_Abundance_det_design_summary.txt`

### Occupancy comparison across demographic classes

| Class | Nature reserves | Jackals | Wolves | Interpretation |
|---|---:|---:|---:|---|
| Male | -1.183 [-2.316, -0.124] | -1.190 [-2.326, -0.138] | 1.147 [0.056, 2.353] | Males avoid reserves and jackals, and show a positive wolf association |
| Female | -1.491 [-2.632, -0.481] | -1.405 [-2.562, -0.358] | 1.079 [-0.010, 2.284] | Females avoid reserves and jackals; wolf effect is borderline positive but still uncertain |
| Fawn | 0.095 [-0.944, 1.280] | -1.296 [-2.718, -0.204] | 0.507 [-0.409, 1.633] | Fawns avoid jackals, but show no clear reserve effect and no clear wolf effect |

### WAIC by class

- Male: 284.59
- Female: 376.11
- Fawn: 317.47

### Ecological contrast

- Jackal abundance is the most consistent negative predictor across all three demographic classes.
- Reserve use differs by class: males and females show lower occupancy in reserves, while fawns do not show a clear reserve response.
- The wolf signal is strongest for males, weaker and borderline for females, and uncertain for fawns.

## Swapping the habitat covariate

Using the same design-based detection model and `Null.Wolf.Abundance.z`, I replaced `Nature.Reserves.z` with `Open.Spaces.z` or `Built.and.Roads.z`.

### Fawns

- `Nature.Reserves.z`: 0.095 [-0.944, 1.280], WAIC 317.47
- `Open.Spaces.z`: 2.141 [0.865, 3.707], WAIC 304.82
- `Built.and.Roads.z`: -0.075 [-1.344, 1.288], WAIC 317.08

Interpretation:

- Fawn occupancy is much better explained by `Open.Spaces` than by reserves or infrastructure.
- Under the `Open.Spaces` model, the jackal effect weakens slightly and its 95% CrI just overlaps 0.

### Males

- `Nature.Reserves.z`: -1.183 [-2.316, -0.124], WAIC 284.59
- `Open.Spaces.z`: 0.338 [-0.790, 1.507], WAIC 288.82
- `Built.and.Roads.z`: -1.818 [-3.063, -0.703], WAIC 276.32

Interpretation:

- For males, `Built.and.Roads` performs best and shows a strong negative association.
- `Open.Spaces` does not explain male occupancy as well as reserves or infrastructure.

### Females

- `Nature.Reserves.z`: -1.491 [-2.632, -0.481], WAIC 376.11
- `Open.Spaces.z`: 0.826 [-0.230, 1.971], WAIC 382.51
- `Built.and.Roads.z`: -1.150 [-2.270, -0.142], WAIC 379.81

Interpretation:

- For females, `Nature.Reserves` remains the best of the three habitat terms.
- `Built.and.Roads` is also negative and credible, but fits worse than reserves.
- `Open.Spaces` is not clearly associated with female occupancy.

## Partitioning reserve effects

I then fit paired habitat models to test whether the reserve effect survives once correlated landscape structure is included explicitly:

- `Nature.Reserves.z + Built.and.Roads.z`
- `Nature.Reserves.z + Open.Spaces.z`

### Nature reserves plus infrastructure

| Class | Nature reserves | Built and roads | Interpretation |
|---|---:|---:|---|
| Fawn | 0.188 [-1.070, 1.642] | -0.129 [-1.644, 1.360] | Neither effect is clear |
| Male | -0.282 [-1.652, 1.131] | -1.716 [-3.184, -0.388] | Reserve effect disappears; infrastructure remains strongly negative |
| Female | -1.182 [-2.456, 0.009] | -0.604 [-1.842, 0.600] | Reserve effect weakens substantially and becomes borderline |

### Nature reserves plus open spaces

| Class | Nature reserves | Open spaces | Interpretation |
|---|---:|---:|---|
| Fawn | -0.328 [-1.537, 0.785] | 2.306 [0.921, 3.976] | Open space carries the signal; reserve effect disappears |
| Male | -1.225 [-2.346, -0.164] | 0.423 [-0.784, 1.682] | Reserve effect remains negative after controlling for open space |
| Female | -1.636 [-2.798, -0.591] | 1.178 [-0.040, 2.490] | Reserve effect remains strongly negative; open-space effect is uncertain |

### Takeaway

- For fawns, the reserve effect is not robust. Open space is the clearer explanatory variable.
- For males, the reserve effect is not robust to infrastructure, but it does persist when controlling only for open space.
- For females, the reserve effect is the most stable of the three classes and remains negative when controlling for open space, though it weakens when infrastructure is added.
