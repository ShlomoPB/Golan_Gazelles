# Gazelle diel activity across predator-activity quartiles

The analysis used 3925 validated class-specific 30-minute gazelle events. For each predator, Q1 and Q4 each contained 15 of the 60 cameras.

Nocturnality was calculated as absolute apparent-solar hours from solar noon (0 = noon; 12 = midnight). Negative Q4-Q1 changes indicate greater diurnality at high predator-activity sites.

## Quartile thresholds

# A tibble: 4 × 7
  predator quartile          cameras activity_min activity_max q25_threshold
  <fct>    <fct>               <int>        <dbl>        <dbl>         <dbl>
1 Wolf     Q1: low activity       15     0.000194        0.450         0.454
2 Wolf     Q4: high activity      15     0.710           1             0.454
3 Jackal   Q1: low activity       15     0.0780          0.311         0.331
4 Jackal   Q4: high activity      15     0.668           1             0.331
# ℹ 1 more variable: q75_threshold <dbl>

## Results

| Predator | Class | Events Q1 / Q4 | Q1-Q4 overlap [95% CI] | compareCkern P | Mean nocturnality Q1 / Q4 (h) | Delta nocturnality [95% CI] (h) | Nocturnality P |
|---|---|---:|---:|---:|---:|---:|---:|
| Wolf | Fawn | 30 / 73 | Dhat1 = 0.91 [0.53, 0.89] | 0.977 (Not supported) | 3.70 / 3.66 | -0.04 [-0.70, 1.66] | 0.915 (Not supported) |
| Wolf | Female | 216 / 609 | Dhat4 = 0.87 [0.67, 0.87] | 0.013 (P < 0.05) | 5.49 / 5.47 | -0.02 [-1.65, 1.69] | 0.995 (Not supported) |
| Wolf | Male | 369 / 852 | Dhat4 = 0.90 [0.80, 0.91] | 0.018 (P < 0.05) | 6.10 / 5.40 | -0.69 [-1.18, -0.15] | 0.072 (Not supported) |
| Jackal | Fawn | 54 / 14 | Dhat1 = 0.45 [0.30, 0.59] | Data deficient | 3.72 / 3.68 | -0.04 [-1.17, 1.20] | Data deficient |
| Jackal | Female | 463 / 107 | Dhat4 = 0.87 [0.69, 0.89] | 0.424 (Not supported) | 5.25 / 4.61 | -0.64 [-1.96, 0.62] | 0.404 (Not supported) |
| Jackal | Male | 855 / 225 | Dhat4 = 0.87 [0.74, 0.89] | 0.017 (P < 0.05) | 5.72 / 5.80 | +0.09 [-1.27, 1.27] | 0.854 (Not supported) |

Using unadjusted compareCkern P values, whole-distribution differences had P < 0.05 for female at wolf-quartile sites, male at wolf-quartile sites, male at jackal-quartile sites.
Using unadjusted two-sided camera-label permutation tests, nocturnality changes had P < 0.05 for no comparison.
The fawn contrast across jackal quartiles was data deficient because Q4 contained 14 events from four cameras.

## Interpretation constraints

Quartiles are descriptive site categories. The compareCkern test resamples event times rather than camera clusters, so its unadjusted P values should be interpreted alongside camera-bootstrap intervals and sample sizes. Nocturnality P values use camera-label permutations. Predator indices are relative camera-level activity measures, not direct predation-risk estimates.
