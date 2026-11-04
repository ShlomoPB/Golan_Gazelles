# Timestamp interpretation audit

Camera clocks were confirmed to have been set to Israel local time at installation. The primary unified spatiotemporal outputs therefore use:

```powershell
Rscript run_unified_spatiotemporal_analysis.R --timestamp_mode=local --bootstrap_reps=500 --prefix=publication_unified_spatiotemporal
```

## Robustness to the earlier UTC interpretation

The earlier UTC run is useful only as a sensitivity check. Spatial event-rate, hurdle, and fawn:female composition results are unchanged because they do not use time of day. Pairwise temporal-overlap conclusions were also qualitatively stable: adults overlapped most strongly with each other, fawns were more temporally distinct from adults, and fawn-female overlap increased later in the survey.

| Period | Pair | Earlier UTC overlap [95% CI] | Primary local-time overlap [95% CI] |
|---|---|---:|---:|
| All | Fawn vs Female | 0.777 [0.691, 0.859] | 0.789 [0.707, 0.860] |
| All | Fawn vs Male | 0.727 [0.653, 0.772] | 0.743 [0.661, 0.783] |
| All | Female vs Male | 0.944 [0.852, 0.954] | 0.948 [0.856, 0.958] |
| Early | Fawn vs Female | 0.713 [0.604, 0.807] | 0.712 [0.610, 0.810] |
| Late | Fawn vs Female | 0.848 [0.741, 0.883] | 0.851 [0.751, 0.889] |

## Daylight-use consequence

The main interpretation affected by the clock correction is daylight use. Under the confirmed local-clock interpretation, females showed greater daylight use in the late survey period, males showed a weaker positive shift, and fawns remained uncertain.

| Class | Earlier UTC late-vs-early daylight OR [95% CI] | Primary local-time late-vs-early daylight OR [95% CI] |
|---|---:|---:|
| Fawn | 0.60 [0.27, 1.35] | 0.52 [0.19, 1.47] |
| Female | 1.04 [0.79, 1.36] | 1.60 [1.16, 2.20] |
| Male | 0.94 [0.76, 1.17] | 1.27 [0.99, 1.64] |

The local-time results should be used in the manuscript. The early/late contrast remains a survey-period contrast and should not be described as a direct temperature effect unless temperature data are added.
