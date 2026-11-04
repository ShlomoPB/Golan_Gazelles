# Pairwise predator-context temporal overlap within gazelle classes

## Methods

All six context pairs were evaluated within each class using the all-camera median classification. `Dhat1` was used when the smaller sample contained fewer than 75 events; otherwise `Dhat4` was used. Intervals were obtained from 500 camera-cluster bootstrap replicates. Distributional differences were tested with `activity::compareCkern` using 999 null replicates and a +1 Monte Carlo correction. We report unadjusted P-values.

Because `compareCkern` resamples events rather than cameras, its P-values do not fully account for camera-level clustering. Comparisons with fewer than 20 events in either context were calculated descriptively but excluded from formal inference.

## Results

| Class | Context pair | Events | Estimator | Dhat [camera-bootstrap 95% CI] | compareCkern P |
|---|---|---:|---|---:|---:|
| Fawn | High jackal / low wolf vs High both | `3 / 38` | `Dhat1` | `0.672 [NA]` | `NA` |
| Fawn | High jackal / low wolf vs Low jackal / high wolf | `3 / 72` | `Dhat1` | `0.716 [NA]` | `NA` |
| Fawn | Low both vs High jackal / low wolf | `67 / 3` | `Dhat1` | `0.743 [NA]` | `NA` |
| Fawn | Low jackal / high wolf vs High both | `72 / 38` | `Dhat1` | `0.796 [0.477, 0.854]` | `0.179` |
| Fawn | Low both vs Low jackal / high wolf | `67 / 72` | `Dhat1` | `0.814 [0.516, 0.869]` | `0.167` |
| Fawn | Low both vs High both | `67 / 38` | `Dhat1` | `0.864 [0.591, 0.866]` | `0.762` |
| Female | Low jackal / high wolf vs High both | `601 / 276` | `Dhat4` | `0.788 [0.652, 0.883]` | `0.001` |
| Female | High jackal / low wolf vs Low jackal / high wolf | `69 / 601` | `Dhat1` | `0.794 [0.585, 0.859]` | `0.007` |
| Female | Low both vs Low jackal / high wolf | `531 / 601` | `Dhat4` | `0.799 [0.682, 0.874]` | `0.001` |
| Female | Low both vs High both | `531 / 276` | `Dhat4` | `0.814 [0.632, 0.896]` | `0.001` |
| Female | Low both vs High jackal / low wolf | `531 / 69` | `Dhat1` | `0.835 [0.599, 0.891]` | `0.179` |
| Female | High jackal / low wolf vs High both | `69 / 276` | `Dhat1` | `0.843 [0.590, 0.867]` | `0.112` |
| Male | High jackal / low wolf vs High both | `163 / 544` | `Dhat4` | `0.838 [0.652, 0.894]` | `0.002` |
| Male | Low both vs High jackal / low wolf | `791 / 163` | `Dhat4` | `0.877 [0.697, 0.898]` | `0.124` |
| Male | High jackal / low wolf vs Low jackal / high wolf | `163 / 770` | `Dhat4` | `0.900 [0.734, 0.922]` | `0.205` |
| Male | Low both vs High both | `791 / 544` | `Dhat4` | `0.902 [0.757, 0.922]` | `0.009` |
| Male | Low jackal / high wolf vs High both | `770 / 544` | `Dhat4` | `0.906 [0.781, 0.928]` | `0.018` |
| Male | Low both vs Low jackal / high wolf | `791 / 770` | `Dhat4` | `0.912 [0.802, 0.917]` | `0.012` |

## Lowest eligible overlap by class

- **Fawn:** Low jackal / high wolf vs High both, `Dhat = 0.796` (95% camera-bootstrap CI `0.477-0.854`); unadjusted `compareCkern P = 0.179`.
- **Female:** Low jackal / high wolf vs High both, `Dhat = 0.788` (95% camera-bootstrap CI `0.652-0.883`); unadjusted `compareCkern P = 0.001`.
- **Male:** High jackal / low wolf vs High both, `Dhat = 0.838` (95% camera-bootstrap CI `0.652-0.894`); unadjusted `compareCkern P = 0.002`.
