# Predator-context temporal analysis

## Camera classification

Cameras were classified using quartiles of `null.jackal.abundance` and `Null.Wolf.Abundance`: jackal Q1 = `0.331`, jackal Q3 = `0.657`, wolf Q1 = `0.454`, and wolf Q3 = `0.667`.

| Category | Cameras |
|---|---:|
| `High jackal / Low wolf` | `6` |
| `Low jackal / High wolf` | `6` |
| `High both` | `3` |
| `Low both` | `2` |
| `Intermediate / mixed` | `43` |

## Event totals by class and predator context

| Category | Class | Period | Events |
|---|---|---|---:|
| `High jackal / Low wolf` | `Fawn` | `Early (hotter)` | `1` |
| `High jackal / Low wolf` | `Fawn` | `Late (mild-cold)` | `1` |
| `High jackal / Low wolf` | `Female` | `Early (hotter)` | `3` |
| `High jackal / Low wolf` | `Female` | `Late (mild-cold)` | `8` |
| `High jackal / Low wolf` | `Male` | `Early (hotter)` | `26` |
| `High jackal / Low wolf` | `Male` | `Late (mild-cold)` | `14` |
| `Low jackal / High wolf` | `Fawn` | `Early (hotter)` | `20` |
| `Low jackal / High wolf` | `Fawn` | `Late (mild-cold)` | `20` |
| `Low jackal / High wolf` | `Female` | `Early (hotter)` | `100` |
| `Low jackal / High wolf` | `Female` | `Late (mild-cold)` | `111` |
| `Low jackal / High wolf` | `Male` | `Early (hotter)` | `295` |
| `Low jackal / High wolf` | `Male` | `Late (mild-cold)` | `211` |
| `High both` | `Fawn` | `Early (hotter)` | `4` |
| `High both` | `Fawn` | `Late (mild-cold)` | `8` |
| `High both` | `Female` | `Early (hotter)` | `33` |
| `High both` | `Female` | `Late (mild-cold)` | `48` |
| `High both` | `Male` | `Early (hotter)` | `60` |
| `High both` | `Male` | `Late (mild-cold)` | `76` |
| `Low both` | `Fawn` | `Late (mild-cold)` | `5` |
| `Low both` | `Female` | `Early (hotter)` | `40` |
| `Low both` | `Female` | `Late (mild-cold)` | `87` |
| `Low both` | `Male` | `Early (hotter)` | `36` |
| `Low both` | `Male` | `Late (mild-cold)` | `114` |
| `Intermediate / mixed` | `Fawn` | `Early (hotter)` | `41` |
| `Intermediate / mixed` | `Fawn` | `Late (mild-cold)` | `80` |
| `Intermediate / mixed` | `Female` | `Early (hotter)` | `566` |
| `Intermediate / mixed` | `Female` | `Late (mild-cold)` | `483` |
| `Intermediate / mixed` | `Male` | `Early (hotter)` | `741` |
| `Intermediate / mixed` | `Male` | `Late (mild-cold)` | `700` |

## Temporal comparisons

Only comparisons with at least the minimum event threshold in both groups were tested formally.

| Period | Class | Contrast | n1 | n2 | Overlap [95% CI] | Kernel p |
|---|---|---|---:|---:|---:|---:|
| `All detections` | `Female` | `High jackal / Low wolf vs Low jackal / High wolf` | `11` | `211` | `NA` | `NA` |
| `All detections` | `Female` | `High both vs Low both` | `81` | `127` | `0.724 [0.623, 0.825]` | `0.010` |
| `All detections` | `Female` | `Extreme predator contexts vs Intermediate / mixed` | `430` | `1049` | `0.888 [0.840, 0.935]` | `<0.001` |
| `All detections` | `Fawn` | `High jackal / Low wolf vs Low jackal / High wolf` | `2` | `40` | `NA` | `NA` |
| `All detections` | `Fawn` | `High both vs Low both` | `12` | `5` | `NA` | `NA` |
| `All detections` | `Fawn` | `Extreme predator contexts vs Intermediate / mixed` | `59` | `121` | `0.873 [0.782, 0.963]` | `0.646` |
| `All detections` | `Male` | `High jackal / Low wolf vs Low jackal / High wolf` | `40` | `506` | `0.694 [0.610, 0.778]` | `<0.001` |
| `All detections` | `Male` | `High both vs Low both` | `136` | `150` | `0.846 [0.765, 0.928]` | `0.081` |
| `All detections` | `Male` | `Extreme predator contexts vs Intermediate / mixed` | `832` | `1441` | `0.904 [0.874, 0.934]` | `<0.001` |
| `Early (hotter)` | `Female` | `High jackal / Low wolf vs Low jackal / High wolf` | `3` | `100` | `NA` | `NA` |
| `Early (hotter)` | `Female` | `High both vs Low both` | `33` | `40` | `0.608 [0.450, 0.766]` | `<0.001` |
| `Early (hotter)` | `Female` | `Extreme predator contexts vs Intermediate / mixed` | `176` | `566` | `0.868 [0.809, 0.927]` | `0.040` |
| `Early (hotter)` | `Fawn` | `High jackal / Low wolf vs Low jackal / High wolf` | `1` | `20` | `NA` | `NA` |
| `Early (hotter)` | `Fawn` | `High both vs Low both` | `4` | `0` | `NA` | `NA` |
| `Early (hotter)` | `Fawn` | `Extreme predator contexts vs Intermediate / mixed` | `25` | `41` | `0.810 [0.664, 0.955]` | `0.556` |
| `Early (hotter)` | `Male` | `High jackal / Low wolf vs Low jackal / High wolf` | `26` | `295` | `0.619 [0.488, 0.750]` | `0.010` |
| `Early (hotter)` | `Male` | `High both vs Low both` | `60` | `36` | `0.939 [0.828, 1.050]` | `0.980` |
| `Early (hotter)` | `Male` | `Extreme predator contexts vs Intermediate / mixed` | `417` | `741` | `0.898 [0.851, 0.945]` | `<0.001` |
| `Late (mild-cold)` | `Female` | `High jackal / Low wolf vs Low jackal / High wolf` | `8` | `111` | `NA` | `NA` |
| `Late (mild-cold)` | `Female` | `High both vs Low both` | `48` | `87` | `0.651 [0.531, 0.771]` | `<0.001` |
| `Late (mild-cold)` | `Female` | `Extreme predator contexts vs Intermediate / mixed` | `254` | `483` | `0.863 [0.808, 0.917]` | `<0.001` |
| `Late (mild-cold)` | `Fawn` | `High jackal / Low wolf vs Low jackal / High wolf` | `1` | `20` | `NA` | `NA` |
| `Late (mild-cold)` | `Fawn` | `High both vs Low both` | `8` | `5` | `NA` | `NA` |
| `Late (mild-cold)` | `Fawn` | `Extreme predator contexts vs Intermediate / mixed` | `34` | `80` | `0.751 [0.620, 0.882]` | `0.071` |
| `Late (mild-cold)` | `Male` | `High jackal / Low wolf vs Low jackal / High wolf` | `14` | `211` | `NA` | `NA` |
| `Late (mild-cold)` | `Male` | `High both vs Low both` | `76` | `114` | `0.751 [0.645, 0.856]` | `0.010` |
| `Late (mild-cold)` | `Male` | `Extreme predator contexts vs Intermediate / mixed` | `415` | `700` | `0.905 [0.858, 0.951]` | `0.030` |

## Supported contrasts

- In `Early (hotter)`, `Female` detections differed between `High both` and `Low both` cameras (overlap `0.608 [0.450, 0.766]`, `p = <0.001`).
- In `Late (mild-cold)`, `Female` detections differed between `High both` and `Low both` cameras (overlap `0.651 [0.531, 0.771]`, `p = <0.001`).
- In `All detections`, `Male` detections differed between `High jackal / Low wolf` and `Low jackal / High wolf` cameras (overlap `0.694 [0.610, 0.778]`, `p = <0.001`).
- In `Late (mild-cold)`, `Female` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.863 [0.808, 0.917]`, `p = <0.001`).
- In `All detections`, `Female` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.888 [0.840, 0.935]`, `p = <0.001`).
- In `Early (hotter)`, `Male` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.898 [0.851, 0.945]`, `p = <0.001`).
- In `All detections`, `Male` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.904 [0.874, 0.934]`, `p = <0.001`).
- In `Early (hotter)`, `Male` detections differed between `High jackal / Low wolf` and `Low jackal / High wolf` cameras (overlap `0.619 [0.488, 0.750]`, `p = 0.010`).
- In `All detections`, `Female` detections differed between `High both` and `Low both` cameras (overlap `0.724 [0.623, 0.825]`, `p = 0.010`).
- In `Late (mild-cold)`, `Male` detections differed between `High both` and `Low both` cameras (overlap `0.751 [0.645, 0.856]`, `p = 0.010`).
- In `Late (mild-cold)`, `Male` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.905 [0.858, 0.951]`, `p = 0.030`).
- In `Early (hotter)`, `Female` detections differed between `Extreme predator contexts` and `Intermediate / mixed` cameras (overlap `0.868 [0.809, 0.927]`, `p = 0.040`).
