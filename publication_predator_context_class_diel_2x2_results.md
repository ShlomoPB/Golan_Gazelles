# Four-context gazelle diel activity figure

## Design

The figure uses independent 30-min events within documented active camera rounds. Camera timestamps were interpreted as Israel local civil time and converted to apparent solar time. Low and high predator contexts were defined by the 33rd and 67th percentiles: jackal cutoffs `0.416` and `0.601`; wolf cutoffs `0.490` and `0.611`. Intermediate cameras were excluded.

Kernel curves are normalized within each class and predator context. They describe temporal distribution, not abundance or event rate.

## Panel support

| Predator context | Cameras | Fawn events | Female events | Male events |
|---|---:|---:|---:|---:|
| Low both | `6` | `30` | `281` | `335` |
| High jackal / Low wolf | `9` | `2` | `17` | `87` |
| Low jackal / High wolf | `8` | `58` | `471` | `582` |
| High both | `6` | `12` | `108` | `180` |

The high-jackal/low-wolf category contained only two independent fawn events. No fawn density curve or fawn-adult overlap estimate was produced for that panel; the two event times are displayed as rug ticks. This absence is ecologically notable but does not provide enough information to estimate a diel distribution.

## Pairwise temporal overlap

Intervals use a camera-cluster bootstrap. Estimates requiring fewer than 20 events in either class were not calculated.

| Predator context | Pair | Events | Cameras with detections | Overlap [95% CI] | Status |
|---|---|---:|---:|---:|---|
| Low both | Fawn vs Female | `30 / 281` | `6 / 6` | `0.696 [0.414, 0.807]` | Estimated |
| Low both | Fawn vs Male | `30 / 335` | `6 / 6` | `0.575 [0.295, 0.766]` | Estimated |
| Low both | Female vs Male | `281 / 335` | `6 / 6` | `0.839 [0.637, 0.874]` | Estimated |
| High jackal / Low wolf | Fawn vs Female | `2 / 17` | `1 / 5` | `NA` | Fewer than 20 events in at least one class |
| High jackal / Low wolf | Fawn vs Male | `2 / 87` | `1 / 5` | `NA` | Fewer than 20 events in at least one class |
| High jackal / Low wolf | Female vs Male | `17 / 87` | `5 / 5` | `NA` | Fewer than 20 events in at least one class |
| Low jackal / High wolf | Fawn vs Female | `58 / 471` | `6 / 8` | `0.714 [0.622, 0.843]` | Estimated |
| Low jackal / High wolf | Fawn vs Male | `58 / 582` | `6 / 8` | `0.726 [0.630, 0.759]` | Estimated |
| Low jackal / High wolf | Female vs Male | `471 / 582` | `8 / 8` | `0.925 [0.781, 0.928]` | Estimated |
| High both | Fawn vs Female | `12 / 108` | `3 / 6` | `NA` | Fewer than 20 events in at least one class |
| High both | Fawn vs Male | `12 / 180` | `3 / 6` | `NA` | Fewer than 20 events in at least one class |
| High both | Female vs Male | `108 / 180` | `6 / 6` | `0.887 [0.742, 0.896]` | Estimated |
