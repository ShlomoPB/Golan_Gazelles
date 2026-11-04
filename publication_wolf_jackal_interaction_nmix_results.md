# Wolf x Jackal interaction N-mixture comparison

Models were fit separately for fawns, females, and males. The `unmarked` models used the camera- and habitat-based detection model:

`detection ~ Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera`

The `spAbundance` models used the simpler effort-based detection model:

`detection ~ log_days_operational.z + cattle_sightings.z`

Both frameworks used the predator-only abundance model:

`abundance ~ wolf * jackal`

Unmarked models used `pcount(..., mixture = "NB")`. spAbundance models used `NMix(..., family = "NB")`.

## Coefficient summary

| Species | Framework | Term | Estimate or posterior mean | Interval / p |
|---|---|---|---:|---:|
| `Fawn` | `unmarked` | `(Intercept)` | `3.173` | `p = <0.001` |
| `Fawn` | `unmarked` | `wolf_z` | `0.265` | `p = 0.237` |
| `Fawn` | `unmarked` | `jackal_z` | `-0.742` | `p = 0.011` |
| `Fawn` | `unmarked` | `wolf_z:jackal_z` | `0.485` | `p = 0.092` |
| `Fawn` | `spAbundance` | `Null.Wolf.Abundance.z` | `0.265` | `0.27 [-0.27, 0.85]` |
| `Fawn` | `spAbundance` | `null.jackal.abundance.z` | `-0.803` | `-0.80 [-1.53, -0.13]` |
| `Fawn` | `spAbundance` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.497` | `0.50 [-0.21, 1.30]` |
| `Female` | `unmarked` | `(Intercept)` | `4.298` | `p = <0.001` |
| `Female` | `unmarked` | `wolf_z` | `0.355` | `p = 0.027` |
| `Female` | `unmarked` | `jackal_z` | `-0.547` | `p = 0.001` |
| `Female` | `unmarked` | `wolf_z:jackal_z` | `0.444` | `p = 0.017` |
| `Female` | `spAbundance` | `Null.Wolf.Abundance.z` | `0.540` | `0.54 [0.13, 0.97]` |
| `Female` | `spAbundance` | `null.jackal.abundance.z` | `-0.721` | `-0.72 [-1.19, -0.27]` |
| `Female` | `spAbundance` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.521` | `0.52 [0.04, 1.02]` |
| `Male` | `unmarked` | `(Intercept)` | `4.007` | `p = <0.001` |
| `Male` | `unmarked` | `wolf_z` | `0.318` | `p = 0.021` |
| `Male` | `unmarked` | `jackal_z` | `-0.542` | `p = <0.001` |
| `Male` | `unmarked` | `wolf_z:jackal_z` | `0.423` | `p = 0.007` |
| `Male` | `spAbundance` | `Null.Wolf.Abundance.z` | `0.352` | `0.35 [0.01, 0.71]` |
| `Male` | `spAbundance` | `null.jackal.abundance.z` | `-0.607` | `-0.61 [-0.96, -0.26]` |
| `Male` | `spAbundance` | `Null.Wolf.Abundance.z:null.jackal.abundance.z` | `0.395` | `0.39 [0.03, 0.77]` |

## Interpretation

Both frameworks recovered the same qualitative predator pattern: jackal effects were negative, wolf effects were positive especially for adults, and wolf x jackal interactions were positive.

The adult interaction terms were supported in both `unmarked` and `spAbundance`, which makes predator context a stronger result than if it appeared in only one model family. The ecological interpretation is that the negative jackal association weakened where wolf activity was high, or equivalently that the positive wolf association was strongest where jackal activity was also high.

This supports a two-predator-context hypothesis for gazelle spatial abundance. It should not be interpreted as direct evidence of wolf-mediated jackal suppression without temporally resolved predator data.
