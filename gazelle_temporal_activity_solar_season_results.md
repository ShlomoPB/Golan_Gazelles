# Gazelle Temporal Activity Analysis in Apparent Solar Time

## Overview

- Event file: `Gazella gazella_independent_events_30min.csv`
- Time transformation: UTC timestamps were converted to apparent solar time using event date and event longitude.
- Seasonal split date (median event timestamp): `2020-10-18 06:31:00 UTC`
- Seasonal labels:
  - `Early (hotter)`: detections on or before the split date
  - `Late (mild-cold)`: detections after the split date
- Minimum event threshold per group: `20`

Monthly detection counts:

2020-01: 1 | 2020-02: 2 | 2020-03: 1 | 2020-04: 1 | 2020-05: 7 | 2020-06: 4 | 2020-07: 10 | 2020-08: 136 | 2020-09: 1117 | 2020-10: 1277 | 2020-11: 1018 | 2020-12: 358

## Overall class comparisons

| Period | Comparison | n1 | n2 | Overlap [95% CI] | Activity1 | Activity2 | Kernel comparison p |
|---|---|---:|---:|---:|---:|---:|---:|
| `All detections` | `Female vs Fawn` | `1479` | `180` | `0.777 [0.729, 0.824]` | `0.543` | `0.364` | `0.000` |
| `All detections` | `Female vs Male` | `1479` | `2273` | `0.945 [0.922, 0.968]` | `0.543` | `0.606` | `0.067` |
| `All detections` | `Fawn vs Male` | `180` | `2273` | `0.728 [0.678, 0.777]` | `0.364` | `0.606` | `0.000` |
| `Early (hotter)` | `Female vs Fawn` | `742` | `66` | `0.712 [0.630, 0.795]` | `0.642` | `0.439` | `0.000` |
| `Early (hotter)` | `Female vs Male` | `742` | `1158` | `0.956 [0.930, 0.982]` | `0.642` | `0.679` | `0.651` |
| `Early (hotter)` | `Fawn vs Male` | `66` | `1158` | `0.720 [0.648, 0.791]` | `0.439` | `0.679` | `0.000` |
| `Late (mild-cold)` | `Female vs Fawn` | `737` | `114` | `0.858 [0.798, 0.917]` | `0.450` | `0.333` | `0.208` |
| `Late (mild-cold)` | `Female vs Male` | `737` | `1115` | `0.892 [0.854, 0.929]` | `0.450` | `0.509` | `0.000` |
| `Late (mild-cold)` | `Fawn vs Male` | `114` | `1115` | `0.758 [0.695, 0.821]` | `0.333` | `0.509` | `0.000` |

## Significant high-vs-low quartile temporal shifts

| Period | Class | Covariate | Low events | High events | Overlap [95% CI] | Kernel comparison p |
|---|---|---|---:|---:|---:|---:|
| `Early (hotter)` | `Female` | `Orchards` | `201` | `50` | `0.613 [0.477, 0.748]` | `0.000` |
| `Early (hotter)` | `Female` | `Cropland` | `367` | `133` | `0.664 [0.588, 0.740]` | `0.000` |
| `Early (hotter)` | `Female` | `Firing zones` | `109` | `240` | `0.682 [0.599, 0.765]` | `0.000` |
| `Early (hotter)` | `Female` | `Human abundance` | `207` | `52` | `0.702 [0.600, 0.803]` | `0.000` |
| `Early (hotter)` | `Female` | `Built and roads` | `162` | `49` | `0.708 [0.579, 0.838]` | `0.000` |
| `Early (hotter)` | `Female` | `Nature reserves` | `111` | `290` | `0.803 [0.725, 0.880]` | `0.007` |
| `Early (hotter)` | `Female` | `Boar abundance` | `88` | `281` | `0.814 [0.711, 0.917]` | `0.007` |
| `Early (hotter)` | `Male` | `Built and roads` | `227` | `170` | `0.748 [0.665, 0.831]` | `0.000` |
| `Early (hotter)` | `Male` | `Jackal abundance` | `397` | `101` | `0.806 [0.718, 0.895]` | `0.000` |
| `Early (hotter)` | `Male` | `Boar abundance` | `232` | `337` | `0.820 [0.759, 0.882]` | `0.000` |
| `Early (hotter)` | `Male` | `Orchards` | `397` | `107` | `0.832 [0.757, 0.906]` | `0.007` |
| `Early (hotter)` | `Male` | `Wolf abundance` | `166` | `463` | `0.837 [0.769, 0.904]` | `0.000` |
| `Early (hotter)` | `Male` | `Human abundance` | `371` | `158` | `0.839 [0.772, 0.907]` | `0.000` |
| `Early (hotter)` | `Male` | `Firing zones` | `288` | `399` | `0.840 [0.782, 0.898]` | `0.000` |
| `Early (hotter)` | `Male` | `Minefields` | `346` | `207` | `0.844 [0.780, 0.908]` | `0.000` |
| `Early (hotter)` | `Male` | `Nature reserves` | `213` | `325` | `0.847 [0.785, 0.908]` | `0.000` |
| `Late (mild-cold)` | `Female` | `Built and roads` | `150` | `75` | `0.678 [0.567, 0.789]` | `0.000` |
| `Late (mild-cold)` | `Female` | `Orchards` | `148` | `115` | `0.761 [0.669, 0.853]` | `0.000` |
| `Late (mild-cold)` | `Female` | `Boar abundance` | `125` | `152` | `0.809 [0.728, 0.890]` | `0.027` |
| `Late (mild-cold)` | `Female` | `Wolf abundance` | `133` | `257` | `0.817 [0.740, 0.895]` | `0.000` |
| `Late (mild-cold)` | `Female` | `Minefields` | `266` | `123` | `0.818 [0.754, 0.882]` | `0.013` |
| `Late (mild-cold)` | `Male` | `Built and roads` | `161` | `148` | `0.769 [0.682, 0.857]` | `0.000` |
| `Late (mild-cold)` | `Male` | `Jackal abundance` | `458` | `125` | `0.774 [0.692, 0.855]` | `0.000` |
| `Late (mild-cold)` | `Male` | `Boar abundance` | `184` | `283` | `0.817 [0.744, 0.891]` | `0.000` |
| `Late (mild-cold)` | `Male` | `Human abundance` | `271` | `181` | `0.819 [0.750, 0.888]` | `0.000` |
| `Late (mild-cold)` | `Male` | `Wolf abundance` | `206` | `390` | `0.868 [0.810, 0.926]` | `0.020` |
