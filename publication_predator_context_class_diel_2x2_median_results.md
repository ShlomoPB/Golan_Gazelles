# Four-context gazelle diel activity figure

## Design

The figure uses independent 30-min events within documented active camera rounds. Camera timestamps were interpreted as Israel local civil time and converted to apparent solar time. Low and high predator contexts were defined using a median split: jackal cutoffs `0.495` and `0.495`; wolf cutoffs `0.556` and `0.556`. All cameras were assigned to one of the four predator contexts.

Kernel curves are normalized within each class and predator context. They describe temporal distribution, not abundance or event rate.

## Panel support

| Predator context | Cameras | Fawn events | Female events | Male events |
|---|---:|---:|---:|---:|
| Low both | `14` | `67` | `531` | `791` |
| High jackal / Low wolf | `16` | `3` | `69` | `163` |
| Low jackal / High wolf | `16` | `72` | `601` | `770` |
| High both | `14` | `38` | `276` | `544` |

Class-context groups with fewer than 5 independent events were displayed as event-time ticks rather than assigned a kernel-density curve. Pairwise overlap was not estimated when either class had fewer than 20 events.

## Pairwise temporal overlap

Intervals use a camera-cluster bootstrap. Estimates requiring fewer than 20 events in either class were not calculated.

| Predator context | Pair | Events | Cameras with detections | Overlap [95% CI] | Status |
|---|---|---:|---:|---:|---|
| Low both | Fawn vs Female | `67 / 531` | `12 / 13` | `0.827 [0.628, 0.896]` | Estimated |
| Low both | Fawn vs Male | `67 / 791` | `12 / 13` | `0.739 [0.488, 0.841]` | Estimated |
| Low both | Female vs Male | `531 / 791` | `13 / 13` | `0.891 [0.768, 0.913]` | Estimated |
| High jackal / Low wolf | Fawn vs Female | `3 / 69` | `2 / 10` | `NA` | Fewer than 20 events in at least one class |
| High jackal / Low wolf | Fawn vs Male | `3 / 163` | `2 / 11` | `NA` | Fewer than 20 events in at least one class |
| High jackal / Low wolf | Female vs Male | `69 / 163` | `10 / 11` | `0.760 [0.611, 0.828]` | Estimated |
| Low jackal / High wolf | Fawn vs Female | `72 / 601` | `9 / 14` | `0.715 [0.623, 0.827]` | Estimated |
| Low jackal / High wolf | Fawn vs Male | `72 / 770` | `9 / 15` | `0.734 [0.634, 0.774]` | Estimated |
| Low jackal / High wolf | Female vs Male | `601 / 770` | `14 / 15` | `0.936 [0.829, 0.938]` | Estimated |
| High both | Fawn vs Female | `38 / 276` | `8 / 13` | `0.761 [0.577, 0.815]` | Estimated |
| High both | Fawn vs Male | `38 / 544` | `8 / 13` | `0.719 [0.539, 0.751]` | Estimated |
| High both | Female vs Male | `276 / 544` | `13 / 13` | `0.887 [0.781, 0.911]` | Estimated |
