# Confirmatory joint demographic diel GAM results

## Authoritative source and test interpretation

These results were extracted without refitting from the penalty-calibrated joint
GAM outputs in `joint_demographic_penalty_calibrated_outputs`. The confirmatory
tests are the `mgcv` penalty-aware boundary tests returned by `reTest()` for the
symmetric sum-to-zero demographic deviation smooths. The reported statistic is
therefore not labeled as a conventional chi-square statistic, and
`reference_rank` is the rank of the weighted reference distribution rather than
ordinary chi-square degrees of freedom. EDF is the estimated effective degrees
of freedom of the tested deviation smooth.

## Omnibus demographic-heterogeneity tests

| Predator-associated component | EDF | Reference rank | Statistic | P | Primary conclusion | Seven-specification sensitivity |
|---|---:|---:|---:|---:|---|---|
| Wolf x time | 1.007 | 30 | 2.081 | 0.0778 | Unsupported | Not reopened because the primary omnibus test was unsupported |
| Jackal x time | <0.001 | 30 | <0.001 | 0.7066 | Unsupported | Not reopened because the primary omnibus test was unsupported |
| Wolf x jackal x time | 1.575 | 30 | 6.110 | 0.0060 | Supported | Supported in 7/7 specifications |

The wolf x jackal demographic-deviation test remained supported under the
primary apparent-solar model, `k = 8`, `k = 16`, sun time, the alternate
camera-only hierarchy, independent 30-min histories, and sun time plus the
camera-only hierarchy (all P <= 0.0097). The unsupported wolf-only and
jackal-only families were not reopened in the prespecified closed-testing
sensitivity audit.

## Planned pairwise contrasts

Pairwise models were fitted only for the omnibus-supported wolf x jackal x time
family.

| Contrast | EDF | Reference rank | Statistic | Raw P | Holm P | Conclusion |
|---|---:|---:|---:|---:|---:|---|
| Juvenile - Adult female | 1.504 | 20 | 5.127 | 0.0060 | 0.0121 | Responses differ |
| Juvenile - Adult male | 1.630 | 20 | 6.724 | 0.0037 | 0.0111 | Responses differ |
| Adult female - Adult male | <0.001 | 20 | <0.001 | 0.2767 | 0.2767 | No supported difference |

These direct pairwise results apply to the primary apparent-solar joint model.
The seven-specification audit refitted the all-class omnibus test, not the three
pairwise models within every sensitivity specification. Consequently, the
omnibus demographic heterogeneity is robust across 7/7 specifications, whereas
7/7 robustness of each juvenile-adult pairwise contrast has not been directly
tested. The adult comparison is a failure to detect a difference, not evidence
that the two adult response functions are exactly equivalent.

## Manuscript-level claims

1. The juvenile wolf x jackal temporal response differed from the adult-female response in the primary joint model (Holm P = 0.0121): supported.
2. The juvenile wolf x jackal temporal response differed from the adult-male response in the primary joint model (Holm P = 0.0111): supported.
3. Adult-female and adult-male wolf x jackal temporal responses did not show a supported difference (Holm P = 0.2767): supported only as a statement of no detected difference, not equivalence.
4. Wolf-associated temporal responses differed among classes overall: unsupported (P = 0.0778).
5. The female wolf-associated temporal response differed specifically from juveniles or males: unsupported by direct joint-model testing. Pairwise wolf follow-ups were not warranted because the omnibus test was unsupported.
6. Jackal-associated temporal responses differed among classes: unsupported (P = 0.7066).

The current sentence, "the juvenile wolf x jackal temporal response differed
from both adult females and adult males, whereas the two adult classes did not
differ from one another," agrees with the primary direct joint-model tests. It
should not be described as pairwise-robust across all seven sensitivities,
because only the all-class omnibus interaction was audited in all seven.

## Reconciliation with separate models

The joint and separate models answer different questions. The separate models
support a robust juvenile wolf x jackal x time interaction, a robust
adult-female wolf x time effect, no adult-male predator-associated temporal
effect robust across all sensitivities, and a robust pooled wolf x time effect.
The joint model directly tests whether response functions differ among retained
demographic classes. It supports demographic heterogeneity only for the combined
wolf x jackal temporal component and locates that heterogeneity between
juveniles and each adult class. The supported female wolf effect in the separate
model does not establish a female-versus-other-class wolf difference; the joint
wolf heterogeneity test was unsupported. The pooled model is not a demographic
contrast and also includes unclassified and subadult gazelles.

## Manuscript-ready Results paragraph

Predator-associated diel responses differed among demographic classes for the
wolf x jackal temporal interaction (penalty-aware statistic = 6.110, EDF =
1.575, reference rank = 30, P = 0.006), with omnibus support in all seven
prespecified specifications. Planned follow-ups showed that the juvenile
response differed from those of adult females (statistic = 5.127, EDF = 1.504,
reference rank = 20, raw P = 0.006, Holm-adjusted P = 0.012) and adult males
(statistic = 6.724, EDF = 1.630, reference rank = 20, raw P = 0.004,
Holm-adjusted P = 0.011), whereas the two adult classes did not show a supported
difference (statistic < 0.001, EDF < 0.001, reference rank = 20, raw and
Holm-adjusted P = 0.277). Demographic heterogeneity was not supported for the
wolf-only (statistic = 2.081, EDF = 1.007, reference rank = 30, P = 0.078) or
jackal-only temporal response (statistic < 0.001, EDF < 0.001, reference rank =
30, P = 0.707).

## Manuscript-ready Discussion sentence

Direct joint-model contrasts indicate that juveniles, but not the two adult
classes, differed in how their diel schedules varied along the combined
wolf-jackal gradient, with no direct evidence of demographic heterogeneity in
responses to either predator alone.

