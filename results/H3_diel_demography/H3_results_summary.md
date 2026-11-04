# Revised H3 result summary

We analyzed 3,925 class-specific 30-minute events from 53 cameras in 12 polygons, representing 3,507 reconstructed encounters. The full class by wolf by jackal binomial mixed model converged with polygon, camera, and encounter random intercepts. Fawns were most diurnal (82.8% daylight), followed by females (61.1%) and males (56.8%). Model-adjusted daylight probabilities at mean predator and human activity were fawn 0.86 [0.77, 0.92]; female 0.68 [0.61, 0.74]; male 0.57 [0.50, 0.64], and all pairwise class contrasts were supported. Neither wolf activity (global P = 0.762) nor jackal activity (P = 0.559) showed a supported directional association with daylight. At the other predator's mean, wolf odds ratios ranged from 1.12 to 1.58 and jackal odds ratios from 1.10 to 1.17; every interval included one and no conditional predator slope remained supported after FDR correction. The wolf-jackal interaction (P = 0.438) and its class-specific component (P = 0.453) were unsupported. Higher human activity was associated with lower daylight odds (OR 0.75 [0.58, 0.97]), whereas Late-period odds were higher (OR 1.35 [1.12, 1.62]). A continuous-date sensitivity lowered AIC by 11.4 but left predator slopes nearly unchanged. The nocturnality-hours model showed directionally greater diurnality with increasing predators, but none was supported. Circular tests retained adult whole-distribution differences, indicating temporal redistribution without a consistent day-night shift. Predator indices represent relative camera-level activity, not causal predation risk.

## Primary formula

`daylight ~ class * wolf_z * jackal_z + human_z + period + (1 |      polygon) + (1 | camera_factor) + (1 | encounter_factor)`

## Direct answers to H3 questions

1. Fawns were most diurnal: 82.8% daylight versus 61.1% for females and 56.8% for males; all model-adjusted pairwise class contrasts were supported.
2. Wolf activity was not associated with a general shift toward daylight (global LRT P = 0.762).
3. Jackal activity was not associated with a general shift toward daylight (global LRT P = 0.559).
4. Wolf slopes did not differ among classes (LRT P = 0.800).
5. Jackal slopes did not differ among classes (LRT P = 0.667).
6. The joint wolf-jackal interaction was unsupported (LRT P = 0.438).
7. The class-specific three-way interaction was unsupported (LRT P = 0.453).
8. Fawns had high baseline daylight use and no supported additional predator-associated shift; sparse joint-gradient cells make invariance uncertain rather than demonstrated.
9. The directional models found no supported daylight shift, whereas circular tests retained 4 female and 4 male whole-distribution differences. These results indicate within-day or crepuscular redistribution rather than a consistent day-night displacement.
10. Baseline demographic differences are well supported; predator-associated directional effects are uncertain, especially for fawns.

## Central global tests

# A tibble: 7 × 7
  test                       method statistic    df p_value reduced_AIC full_AIC
  <chr>                      <chr>      <dbl> <dbl>   <dbl>       <dbl>    <dbl>
1 Baseline demographic-clas… Wald …     37.0      2 9.34e-9         NA     5031.
2 Any wolf association acro… Likel…      3.37     6 7.62e-1       5022.    5031.
3 Any jackal association ac… Likel…      4.88     6 5.59e-1       5024.    5031.
4 Wolf slopes differ among … Likel…      1.65     4 8.00e-1       5024.    5031.
5 Jackal slopes differ amon… Likel…      2.38     4 6.67e-1       5025.    5031.
6 Wolf by jackal interactio… Likel…      2.71     3 4.38e-1       5027.    5031.
7 Class by wolf by jackal i… Likel…      1.58     2 4.53e-1       5028.    5031.

## Daylight proportions

# A tibble: 3 × 4
  class  events daylight_events daylight_proportion
  <fct>   <int>           <int>               <dbl>
1 Fawn      180             149               0.828
2 Female   1477             903               0.611
3 Male     2268            1288               0.568
