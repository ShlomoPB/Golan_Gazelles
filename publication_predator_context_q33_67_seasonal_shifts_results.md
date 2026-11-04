# Seasonal shifts within predator-context categories

Predator-context categories used the 33rd and 67th percentiles of jackal and wolf abundance. Season was split at the median gazelle detection timestamp (`2020-10-18 06:31:00` UTC).

## Event counts

| Predator context | Class | Early events | Late events |
|---|---|---:|---:|
| `High jackal / Low wolf` | `Fawn` | `1` | `1` |
| `High jackal / Low wolf` | `Female` | `7` | `10` |
| `High jackal / Low wolf` | `Male` | `42` | `45` |
| `Low jackal / High wolf` | `Fawn` | `31` | `27` |
| `Low jackal / High wolf` | `Female` | `302` | `169` |
| `Low jackal / High wolf` | `Male` | `339` | `243` |
| `High both` | `Fawn` | `4` | `8` |
| `High both` | `Female` | `42` | `66` |
| `High both` | `Male` | `71` | `109` |
| `Low both` | `Fawn` | `3` | `27` |
| `Low both` | `Female` | `91` | `190` |
| `Low both` | `Male` | `108` | `227` |
| `Intermediate / mixed` | `Fawn` | `27` | `51` |
| `Intermediate / mixed` | `Female` | `300` | `302` |
| `Intermediate / mixed` | `Male` | `598` | `491` |

## Early vs late activity comparisons

| Predator context | Class | Early n | Late n | Overlap [95% CI] | Kernel p |
|---|---|---:|---:|---:|---:|
| `High jackal / Low wolf` | `Fawn` | `1` | `1` | `NA` | `NA` |
| `High jackal / Low wolf` | `Female` | `7` | `10` | `NA` | `NA` |
| `High jackal / Low wolf` | `Male` | `42` | `45` | `0.772 [0.629, 0.915]` | `0.182` |
| `Low jackal / High wolf` | `Fawn` | `31` | `27` | `0.801 [0.642, 0.960]` | `0.354` |
| `Low jackal / High wolf` | `Female` | `302` | `169` | `0.651 [0.577, 0.726]` | `<0.001` |
| `Low jackal / High wolf` | `Male` | `339` | `243` | `0.849 [0.785, 0.912]` | `<0.001` |
| `High both` | `Fawn` | `4` | `8` | `NA` | `NA` |
| `High both` | `Female` | `42` | `66` | `0.744 [0.625, 0.863]` | `0.030` |
| `High both` | `Male` | `71` | `109` | `0.672 [0.559, 0.785]` | `<0.001` |
| `Low both` | `Fawn` | `3` | `27` | `NA` | `NA` |
| `Low both` | `Female` | `91` | `190` | `0.758 [0.675, 0.842]` | `<0.001` |
| `Low both` | `Male` | `108` | `227` | `0.753 [0.661, 0.844]` | `<0.001` |
| `Intermediate / mixed` | `Fawn` | `27` | `51` | `0.804 [0.669, 0.940]` | `0.384` |
| `Intermediate / mixed` | `Female` | `300` | `302` | `0.817 [0.760, 0.874]` | `<0.001` |
| `Intermediate / mixed` | `Male` | `598` | `491` | `0.849 [0.808, 0.891]` | `<0.001` |

## Supported seasonal shifts

- `Female` in `Low jackal / High wolf` shifted between early and late periods (overlap `0.651 [0.577, 0.726]`, `p = <0.001`).
- `Female` in `Low both` shifted between early and late periods (overlap `0.758 [0.675, 0.842]`, `p = <0.001`).
- `Female` in `Intermediate / mixed` shifted between early and late periods (overlap `0.817 [0.760, 0.874]`, `p = <0.001`).
- `Female` in `High both` shifted between early and late periods (overlap `0.744 [0.625, 0.863]`, `p = 0.030`).
- `Male` in `Low jackal / High wolf` shifted between early and late periods (overlap `0.849 [0.785, 0.912]`, `p = <0.001`).
- `Male` in `High both` shifted between early and late periods (overlap `0.672 [0.559, 0.785]`, `p = <0.001`).
- `Male` in `Low both` shifted between early and late periods (overlap `0.753 [0.661, 0.844]`, `p = <0.001`).
- `Male` in `Intermediate / mixed` shifted between early and late periods (overlap `0.849 [0.808, 0.891]`, `p = <0.001`).
