# Spatial predator-context count model

Camera-level class counts were modeled with a negative-binomial GLMM using an offset for the number of active rounds.

Model formula:

`count ~ class * predator_context_label + offset(log(active_rounds)) + (1 | Polygon.Name / camera_name)`

## Omnibus interaction test

Likelihood-ratio test comparing additive and interaction models: `Chisq = 7.523`, `df = 8`, `p = 0.481`.

## Predicted count per active round

| Class | Predator context | Predicted count [95% CI] |
|---|---|---:|
| `Female` | `Intermediate / mixed` | `1.97 [0.95, 4.08]` |
| `Male` | `Intermediate / mixed` | `3.00 [1.45, 6.19]` |
| `Fawn` | `Intermediate / mixed` | `0.19 [0.09, 0.42]` |
| `Female` | `Low both` | `11.62 [3.19, 42.38]` |
| `Male` | `Low both` | `12.28 [3.37, 44.75]` |
| `Fawn` | `Low both` | `0.98 [0.25, 3.77]` |
| `Female` | `High both` | `2.81 [0.78, 10.10]` |
| `Male` | `High both` | `3.83 [1.07, 13.73]` |
| `Fawn` | `High both` | `0.22 [0.05, 0.93]` |
| `Female` | `Low jackal / High wolf` | `4.51 [1.29, 15.82]` |
| `Male` | `Low jackal / High wolf` | `6.66 [1.92, 23.04]` |
| `Fawn` | `Low jackal / High wolf` | `0.45 [0.12, 1.63]` |
| `Female` | `High jackal / Low wolf` | `0.21 [0.06, 0.82]` |
| `Male` | `High jackal / Low wolf` | `0.94 [0.26, 3.33]` |
| `Fawn` | `High jackal / Low wolf` | `0.02 [0.00, 0.14]` |

## Context contrasts

| Class | Contrast | Rate ratio [95% CI] | p |
|---|---|---:|---:|
| `Female` | `High jackal / Low wolf vs Intermediate / mixed` | `0.11 [0.03, 0.41]` | `0.001` |
| `Male` | `High jackal / Low wolf vs Intermediate / mixed` | `0.31 [0.09, 1.10]` | `0.070` |
| `Fawn` | `High jackal / Low wolf vs Intermediate / mixed` | `0.11 [0.02, 0.72]` | `0.022` |
| `Female` | `Low jackal / High wolf vs Intermediate / mixed` | `2.30 [0.65, 8.11]` | `0.197` |
| `Male` | `Low jackal / High wolf vs Intermediate / mixed` | `2.22 [0.64, 7.70]` | `0.209` |
| `Fawn` | `Low jackal / High wolf vs Intermediate / mixed` | `2.33 [0.62, 8.76]` | `0.211` |
| `Female` | `High both vs Intermediate / mixed` | `1.43 [0.40, 5.09]` | `0.583` |
| `Male` | `High both vs Intermediate / mixed` | `1.28 [0.36, 4.54]` | `0.703` |
| `Fawn` | `High both vs Intermediate / mixed` | `1.16 [0.28, 4.89]` | `0.839` |
| `Female` | `Low both vs Intermediate / mixed` | `5.91 [1.63, 21.36]` | `0.007` |
| `Male` | `Low both vs Intermediate / mixed` | `4.10 [1.14, 14.77]` | `0.031` |
| `Fawn` | `Low both vs Intermediate / mixed` | `5.07 [1.29, 19.97]` | `0.020` |
| `Female` | `High jackal / Low wolf vs Low jackal / High wolf` | `0.05 [0.01, 0.29]` | `<0.001` |
| `Male` | `High jackal / Low wolf vs Low jackal / High wolf` | `0.14 [0.02, 0.80]` | `0.027` |
| `Fawn` | `High jackal / Low wolf vs Low jackal / High wolf` | `0.05 [0.00, 0.44]` | `0.008` |
