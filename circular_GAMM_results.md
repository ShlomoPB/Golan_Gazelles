# Circular GAMM results

## Analysis summary

We modeled independent class-event counts in 8,256 camera-round-hour cells with negative-binomial hierarchical cyclic GAMMs. Predator covariates were continuous standardized site-level activity indices, not predator density.
Population-level predictions excluded all random intercepts, used round 4 and one operational camera-day, and were normalized within each predator combination for timing comparisons.

## Time-varying predator effects

### Juvenile

The wolf-by-time smooth had EDF = 1.02 and approximate P = 0.086; the jackal-by-time smooth had EDF = 0.00 and P = 0.947; and the wolf-by-jackal-by-time smooth had EDF = 2.02 and P = <0.001.
Across the four -1/+1 SD combinations, the modeled proportion of normalized activity between 06:00 and 18:00 apparent solar time ranged from 0.84 to 0.88. Principal point-estimate peaks ranged from 8.6 to 16.2 h.
Simulation intervals excluded zero for 8.0 h at most among the four conditional predator difference curves.

### Adult female

The wolf-by-time smooth had EDF = 3.83 and approximate P = <0.001; the jackal-by-time smooth had EDF = 0.00 and P = 0.315; and the wolf-by-jackal-by-time smooth had EDF = 0.08 and P = 0.062.
Across the four -1/+1 SD combinations, the modeled proportion of normalized activity between 06:00 and 18:00 apparent solar time ranged from 0.73 to 0.73. Principal point-estimate peaks ranged from 9.7 to 16.7 h.
Simulation intervals excluded zero for 0.0 h at most among the four conditional predator difference curves.

## Early/Late curve-shape sensitivity

We added a cyclic Late-versus-Early difference smooth using the frozen survey midpoint (`2020-10-17`) while retaining round-specific intensity effects. The period smooth was unsupported for juveniles (EDF = 0.00, P = 0.484) but supported for adult females (EDF = 2.97, P < 0.001) and adult males (EDF = 2.92, P < 0.001).

The supported predator effects were stable: the juvenile wolf-by-jackal-by-time effect remained supported (P < 0.001), as did the adult-female wolf-by-time effect (P < 0.001) and adult-male wolf-by-jackal-by-time effect (P = 0.002). Thus, broad seasonal changes in diel shape did not account for the predator-associated smooths.

### Adult male

The wolf-by-time smooth had EDF = 0.00 and approximate P = 0.064; the jackal-by-time smooth had EDF = 0.00 and P = 0.643; and the wolf-by-jackal-by-time smooth had EDF = 3.01 and P = 0.002.
Across the four -1/+1 SD combinations, the modeled proportion of normalized activity between 06:00 and 18:00 apparent solar time ranged from 0.59 to 0.64. Principal point-estimate peaks ranged from 16.6 to 16.9 h.
Simulation intervals excluded zero for 0.0 h at most among the four conditional predator difference curves.

## Relationship to the categorical kernel analysis

The existing kernel analysis supported 0 of 3 eligible context comparisons for juveniles, 4 of 6 for adult females, and 4 of 6 for adult males. The GAMM should be judged by stability, continuous-gradient inference, and the coherence of predicted temporal changes rather than exact agreement with those categorical tests.

## Limitations

`Days.Operational` was the best available exposure measure but assumes approximately equal exposure across the 24 solar-hour bins within a camera-round.
Predator indices are camera-level relative activity gradients, so associations do not demonstrate direct predator avoidance or causal behavioral displacement.
The juvenile sample is much smaller than the adult samples; unsupported juvenile effects should not be interpreted as proof of no response.
The juvenile model retained substantial within-camera-round residual temporal correlation (mean lag-1 Pearson correlation = 0.727), so its interaction should be treated as exploratory despite its smooth-term test.
Adult global smooth tests were supported, but their normalized pointwise difference intervals did not exclude zero at individual times; inference is therefore distribution-wide rather than evidence for sharply bounded hourly windows.
Approximate smooth-term P values are interpreted with shrinkage, diagnostics, and simulation-based effect curves rather than as stand-alone evidence.

## Recommendation for manuscript

**PRIMARY DIEL ANALYSIS RECOMMENDATION: GAMM**

The circular GAMMs converged, passed the prespecified basis and concurvity checks, and yielded temporally coherent predator-associated differences. They provide stronger primary inference by retaining the continuous predator gradients, repeated camera-round structure, and effort. Supported predator smooths were stable after controlling for Early/Late changes in diel curve shape.

### Manuscript-ready Methods

To test whether gazelle diel timing varied continuously across predator-activity gradients, we fitted separate hierarchical cyclic generalized additive models for juvenile, adult-female, and adult-male independent events following Iannarilli et al. (2025). We aggregated 30-min independent class-events into 24 one-hour apparent-solar-time bins for each positive-effort camera-round and retained zero-event cells. Models used a negative-binomial distribution with a log link and the logarithm of operational camera-days as an offset. Each model included ordinary standardized wolf and jackal activity terms and their interaction, factor(round), a baseline cyclic cubic smooth of apparent solar time, and cyclic varying-coefficient smooths for wolf activity, jackal activity, and their interaction. Random-intercept smooths accounted for camera-round, camera, and sampling-polygon clustering. Supported shared-penalty, camera-specific cyclic deviations accounted for site-to-site variation in adult curve shape; we omitted these deviations from the juvenile model because 180 events were insufficient to estimate 60 camera-specific curves reliably. We fitted models with `mgcv::bam` using fREML, shrinkage selection, and conservative basis dimensions (k = 8 for adults). A diagnostic expansion from k = 8 to k = 10 for juveniles resolved a marginal basis-adequacy warning without changing the ecological inference. We evaluated convergence, basis adequacy, concurvity, residual patterns, zero frequencies, and agreement between fitted and observed hourly rates. As a seasonal sensitivity, we added a cyclic Late-versus-Early difference smooth based on the frozen 17 October survey midpoint; this controlled period-specific curve shape without fitting predator-by-period interactions. We generated population-level curves at all combinations of -1 and +1 SD wolf and jackal activity and obtained 95% intervals from 1,500 simulations of the model coefficients and covariance matrix.

### Manuscript-ready Results

The hierarchical circular analysis included 8,256 camera-round-hour cells and 3,925 independent demographic-class events. Juvenile diel shape varied with the wolf by jackal gradient (EDF = 2.02, P < 0.001), whereas its wolf-only and jackal-only timing smooths were unsupported. Predicted juvenile activity peaked in the morning (8.6 h) when predator indices were concordant (both low or both high) and in late afternoon (16.1 h) when indices were discordant. Given only 180 juvenile events and remaining residual temporal dependence, we treated this interaction as exploratory.

Adult-female diel shape varied with wolf activity (EDF = 3.83, P < 0.001), but not with jackal activity or their interaction. The point-estimate principal peak shifted from 9.7 h at low wolf activity to 16.7 h at high wolf activity, while the predicted 06:00-18:00 proportion remained nearly constant (0.73), indicating redistribution within the diurnal-crepuscular period rather than increased nocturnality.

Adult-male diel shape varied with the wolf by jackal gradient (EDF = 3.01, P = 0.002); wolf-only and jackal-only timing smooths were unsupported. Modeled daytime activity was higher when predator indices were concordant (0.64) than discordant (0.59), although the principal peak remained near 16.6-16.9 h. The juvenile interaction, female wolf effect, and male interaction remained supported after controlling for Early/Late curve shape. For both adult classes, fully propagated pointwise intervals overlapped zero throughout the day, so these results indicate global changes in curve shape rather than precisely delimited hourly shifts.

### Proposed figure caption

**Figure X. Modeled diel activity of mountain gazelle demographic classes across continuous wolf and jackal activity gradients.** Curves show population-level normalized predictions from class-specific negative-binomial hierarchical cyclic GAMMs at -1 and +1 SD of each predator index; shading shows simulation-based 95% intervals and rugs show observed independent class-events. Curves integrate to one within each predator combination and therefore describe event timing rather than event abundance.

### Kernel sensitivity placement

Move the four median-split kernel curves, pairwise overlap coefficients, camera-cluster bootstrap intervals, and `compareCkern` tests to Supporting Information. In the main text, state: 'Median-split kernel analyses provided a categorical sensitivity check on the continuous hierarchical GAMM results.'
