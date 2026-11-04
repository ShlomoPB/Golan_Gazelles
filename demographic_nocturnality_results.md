# Astronomical demographic nocturnality analysis

Astronomical sunset-to-sunrise classification yielded 17.2% nighttime detections for juveniles, 38.9% for adult females, and 43.2% for adult males. The corresponding supplied flash-based values were 4.1%, 21.0%, and 37.6%. Overall astronomical daylight among retained class-event rows was 59.6%, compared with the supplied flash-based value of 79.1%.

## Class-specific estimates

| Class | Total | Night | Daylight | Nocturnal % (camera-bootstrap 95% CI) | Diurnal % | Old flash nocturnal % |
|---|---:|---:|---:|---:|---:|---:|
| Juvenile | 180 | 31 | 149 | 17.2 (11.8-21.6) | 82.8 | 4.1 |
| Adult female | 1477 | 574 | 903 | 38.9 (30.8-45.7) | 61.1 | 21.0 |
| Adult male | 2268 | 980 | 1288 | 43.2 (38.4-48.4) | 56.8 | 37.6 |

## Astronomical implementation audit

Night was defined as solar elevation <0 degrees using each corrected local camera timestamp and camera coordinate. Across all 60 cameras and survey dates (2020-08-20 to 2020-12-15), sunrise ranged from 05:52 to 06:53 and sunset from 16:25 to 19:14.
The independent-event export recorded camera 51 at latitude 33.997700, whereas the authoritative camera map records latitude 32.997700 (a one-degree discrepancy); their longitudes agree, and all other camera coordinates agree. Solar calculations used the authoritative map coordinate while retaining the event-export coordinate as an audit field.

| Camera | Date | Sunrise | Elevation 5 min before/after | Sunset | Elevation 5 min before/after |
|---:|---|---|---:|---|---:|
| 1 | 2020-08-20 | 2020-08-20 06:07 | -1.01 / 1.02 | 2020-08-20 19:12 | 1.01 / -1.01 |
| 30 | 2020-10-17 | 2020-10-17 06:46 | -1.03 / 1.03 | 2020-10-17 17:56 | 1.02 / -1.03 |
| 1 | 2020-12-15 | 2020-12-15 06:36 | -0.93 / 0.92 | 2020-12-15 16:27 | 0.92 / -0.93 |

The sign changes around each manually inspected crossing verify the geometric-horizon implementation.

## Camera-cluster bootstrap

We resampled the 60 positive-effort study cameras 10,000 times and retained all class-specific detections from each sampled camera; cameras with no retained events contributed empty clusters.

| Contrast | Difference in nocturnal proportion | Camera-bootstrap 95% CI |
|---|---:|---:|
| Juvenile - Adult female | -0.216 | -0.309 to -0.118 |
| Juvenile - Adult male | -0.260 | -0.341 to -0.191 |
| Adult female - Adult male | -0.043 | -0.127 to 0.026 |

## Formal class comparison

The retained model was: Camera + camera-round random intercepts with demographic class and survey round as fixed effects.
The camera-only model converged (code 0; positive-definite Hessian = TRUE; singular = FALSE; Pearson dispersion = 0.991; maximum absolute gradient = 9.26e-04).
Adding camera-round produced variance 0.3664, changed AIC by -71.82 (candidate minus simple), and changed planned log odds ratios by at most 0.033. Candidate convergence code = 0; positive-definite Hessian = TRUE; singular = FALSE; Pearson dispersion = 0.927; maximum absolute gradient = 8.27e-04.

| Contrast | Adjusted odds ratio | 95% CI | Holm P |
|---|---:|---:|---:|
| Juvenile / Adult female | 0.373 | 0.243-0.572 | <0.001 |
| Juvenile / Adult male | 0.256 | 0.167-0.391 | <0.001 |
| Adult female / Adult male | 0.686 | 0.585-0.804 | <0.001 |

## Interpretation

The astronomical definition retained the previous qualitative ordering: juveniles were most diurnal, adult females intermediate, and adult males most nocturnal.
Astronomical nighttime fractions were higher than the supplied flash-based fractions, especially for juveniles and females. This is expected because the geometric definition begins exactly at sunset and includes twilight, whereas flash activation depends on ambient-light thresholds and camera settings. The old flash percentages could not be recomputed event-by-event because the authoritative independent-event file contains no flash field; they are therefore treated as supplied benchmarks rather than paired values.
These proportions describe when detections occurred, not absolute nighttime activity: astronomical night and daylight had unequal, seasonally changing durations. This analysis does not replace the circular GAM analysis of the complete 24-hour activity distribution.

## Manuscript recommendation

B. Report percentages, camera-bootstrap confidence intervals, and the formal class comparison because the retained round-adjusted model supports the adult-female versus adult-male contrast whereas its unadjusted camera-bootstrap interval includes zero.
