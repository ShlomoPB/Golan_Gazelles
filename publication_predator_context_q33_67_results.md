# Predator-context temporal analysis

## Camera classification

Cameras were classified using quartiles of `null.jackal.abundance` and `Null.Wolf.Abundance`: jackal Q1 = `0.416`, jackal Q3 = `0.601`, wolf Q1 = `0.490`, and wolf Q3 = `0.611`.

| Category | Cameras |
|---|---:|
| `High jackal / Low wolf` | `9` |
| `Low jackal / High wolf` | `8` |
| `High both` | `6` |
| `Low both` | `6` |
| `Intermediate / mixed` | `31` |

## Event totals by class and predator context

| Category | Class | Period | Events |
|---|---|---|---:|
| `High jackal / Low wolf` | `Fawn` | `Early (hotter)` | `1` |
| `High jackal / Low wolf` | `Fawn` | `Late (mild-cold)` | `1` |
| `High jackal / Low wolf` | `Female` | `Early (hotter)` | `7` |
| `High jackal / Low wolf` | `Female` | `Late (mild-cold)` | `10` |
| `High jackal / Low wolf` | `Male` | `Early (hotter)` | `42` |
| `High jackal / Low wolf` | `Male` | `Late (mild-cold)` | `45` |
| `Low jackal / High wolf` | `Fawn` | `Early (hotter)` | `31` |
| `Low jackal / High wolf` | `Fawn` | `Late (mild-cold)` | `27` |
| `Low jackal / High wolf` | `Female` | `Early (hotter)` | `302` |
| `Low jackal / High wolf` | `Female` | `Late (mild-cold)` | `169` |
| `Low jackal / High wolf` | `Male` | `Early (hotter)` | `339` |
| `Low jackal / High wolf` | `Male` | `Late (mild-cold)` | `243` |
| `High both` | `Fawn` | `Early (hotter)` | `4` |
| `High both` | `Fawn` | `Late (mild-cold)` | `8` |
| `High both` | `Female` | `Early (hotter)` | `42` |
| `High both` | `Female` | `Late (mild-cold)` | `66` |
| `High both` | `Male` | `Early (hotter)` | `71` |
| `High both` | `Male` | `Late (mild-cold)` | `109` |
| `Low both` | `Fawn` | `Early (hotter)` | `3` |
| `Low both` | `Fawn` | `Late (mild-cold)` | `27` |
| `Low both` | `Female` | `Early (hotter)` | `91` |
| `Low both` | `Female` | `Late (mild-cold)` | `190` |
| `Low both` | `Male` | `Early (hotter)` | `108` |
| `Low both` | `Male` | `Late (mild-cold)` | `227` |
| `Intermediate / mixed` | `Fawn` | `Early (hotter)` | `27` |
| `Intermediate / mixed` | `Fawn` | `Late (mild-cold)` | `51` |
| `Intermediate / mixed` | `Female` | `Early (hotter)` | `300` |
| `Intermediate / mixed` | `Female` | `Late (mild-cold)` | `302` |
| `Intermediate / mixed` | `Male` | `Early (hotter)` | `598` |
| `Intermediate / mixed` | `Male` | `Late (mild-cold)` | `491` |

## Temporal comparisons

Only comparisons with at least the minimum event threshold in both groups were tested formally.

| Period | Class | Contrast | n1 | n2 | Overlap [95% CI] | Kernel p |
|---|---|---|---:|---:|---:|---:|
| `All detections` | `Female` | `High jackal / Low wolf vs Low jackal / High wolf` | `17` | `471` | `NA` | `NA` |
| `All detections` | `Female` | `High both vs Low both` | `108` | `281` | `0.873 [0.818, 0.927]` | `0.273` |
| `All detections` | `Female` | `Extreme predator contexts vs Intermediate / mixed` | `877` | `602` | `0.845 [0.801, 0.888]` | `<0.001` |
| `All detections` | `Fawn` | `High jackal / Low wolf vs Low jackal / High wolf` | `2` | `58` | `NA` | `NA` |
| `All detections` | `Fawn` | `High both vs Low both` | `12` | `30` | `NA` | `NA` |
| `All detections` | `Fawn` | `Extreme predator contexts vs Intermediate / mixed` | `102` | `78` | `0.844 [0.739, 0.949]` | `0.263` |
| `All detections` | `Male` | `High jackal / Low wolf vs Low jackal / High wolf` | `87` | `582` | `0.801 [0.733, 0.869]` | `<0.001` |
| `All detections` | `Male` | `High both vs Low both` | `180` | `335` | `0.861 [0.794, 0.928]` | `0.051` |
| `All detections` | `Male` | `Extreme predator contexts vs Intermediate / mixed` | `1184` | `1089` | `0.862 [0.824, 0.900]` | `<0.001` |
| `Early (hotter)` | `Female` | `High jackal / Low wolf vs Low jackal / High wolf` | `7` | `302` | `NA` | `NA` |
| `Early (hotter)` | `Female` | `High both vs Low both` | `42` | `91` | `0.812 [0.694, 0.929]` | `0.202` |
| `Early (hotter)` | `Female` | `Extreme predator contexts vs Intermediate / mixed` | `442` | `300` | `0.761 [0.705, 0.816]` | `<0.001` |
| `Early (hotter)` | `Fawn` | `High jackal / Low wolf vs Low jackal / High wolf` | `1` | `31` | `NA` | `NA` |
| `Early (hotter)` | `Fawn` | `High both vs Low both` | `4` | `3` | `NA` | `NA` |
| `Early (hotter)` | `Fawn` | `Extreme predator contexts vs Intermediate / mixed` | `39` | `27` | `0.800 [0.643, 0.957]` | `0.354` |
| `Early (hotter)` | `Male` | `High jackal / Low wolf vs Low jackal / High wolf` | `42` | `339` | `0.764 [0.663, 0.865]` | `<0.001` |
| `Early (hotter)` | `Male` | `High both vs Low both` | `71` | `108` | `0.859 [0.760, 0.958]` | `0.202` |
| `Early (hotter)` | `Male` | `Extreme predator contexts vs Intermediate / mixed` | `560` | `598` | `0.838 [0.792, 0.884]` | `<0.001` |
| `Late (mild-cold)` | `Female` | `High jackal / Low wolf vs Low jackal / High wolf` | `10` | `169` | `NA` | `NA` |
| `Late (mild-cold)` | `Female` | `High both vs Low both` | `66` | `190` | `0.815 [0.733, 0.897]` | `0.091` |
| `Late (mild-cold)` | `Female` | `Extreme predator contexts vs Intermediate / mixed` | `435` | `302` | `0.904 [0.855, 0.953]` | `0.293` |
| `Late (mild-cold)` | `Fawn` | `High jackal / Low wolf vs Low jackal / High wolf` | `1` | `27` | `NA` | `NA` |
| `Late (mild-cold)` | `Fawn` | `High both vs Low both` | `8` | `27` | `NA` | `NA` |
| `Late (mild-cold)` | `Fawn` | `Extreme predator contexts vs Intermediate / mixed` | `63` | `51` | `0.855 [0.741, 0.969]` | `0.667` |
| `Late (mild-cold)` | `Male` | `High jackal / Low wolf vs Low jackal / High wolf` | `45` | `243` | `0.750 [0.642, 0.858]` | `0.030` |
| `Late (mild-cold)` | `Male` | `High both vs Low both` | `109` | `227` | `0.805 [0.730, 0.881]` | `<0.001` |
| `Late (mild-cold)` | `Male` | `Extreme predator contexts vs Intermediate / mixed` | `624` | `491` | `0.875 [0.824, 0.926]` | `<0.001` |

## Supported contrasts

- In `Early (hotter)`, `Female` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.761 [0.705, 0.816]`, `p = <0.001`).
- In `Early (hotter)`, `Male` detections differed between `High jackal / Low wolf` and `Low jackal / High wolf` cameras (overlap `0.764 [0.663, 0.865]`, `p = <0.001`).
- In `All detections`, `Male` detections differed between `High jackal / Low wolf` and `Low jackal / High wolf` cameras (overlap `0.801 [0.733, 0.869]`, `p = <0.001`).
- In `Late (mild-cold)`, `Male` detections differed between `High both` and `Low both` cameras (overlap `0.805 [0.730, 0.881]`, `p = <0.001`).
- In `Early (hotter)`, `Male` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.838 [0.792, 0.884]`, `p = <0.001`).
- In `All detections`, `Female` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.845 [0.801, 0.888]`, `p = <0.001`).
- In `All detections`, `Male` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.862 [0.824, 0.900]`, `p = <0.001`).
- In `Late (mild-cold)`, `Male` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.875 [0.824, 0.926]`, `p = <0.001`).
- In `Late (mild-cold)`, `Male` detections differed between `High jackal / Low wolf` and `Low jackal / High wolf` cameras (overlap `0.750 [0.642, 0.858]`, `p = 0.030`).
