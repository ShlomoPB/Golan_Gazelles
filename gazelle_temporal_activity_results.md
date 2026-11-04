# Gazelle Temporal Activity Analysis

## Overview

- Event file: `Gazella gazella_independent_events_30min.csv`
- Site covariates: `SiteCovs_060324.csv`
- Local timezone used for clock time conversion: `Asia/Jerusalem`
- Each independent event timestamp was treated as one detection; `event_count` was not used as a weight.
- High and low pressure groups were defined as the top and bottom quartiles of each site covariate across cameras.
- Minimum event threshold per group for formal comparison: `20`

## Overall class comparisons

| Comparison | n1 | n2 | Overlap estimator | Overlap [95% CI] | Activity1 | Activity2 | Kernel comparison p |
|---|---:|---:|---|---:|---:|---:|---:|
| `Female vs Fawn` | `1479` | `180` | `Dhat4` | `0.779 [0.731, 0.828]` | `0.500` | `0.350` | `0.000` |
| `Female vs Male` | `1479` | `2273` | `Dhat4` | `0.944 [0.920, 0.968]` | `0.500` | `0.558` | `0.045` |
| `Fawn vs Male` | `180` | `2273` | `Dhat4` | `0.728 [0.676, 0.781]` | `0.350` | `0.558` | `0.000` |

## High vs low pressure comparisons

| Class | Covariate | Low cameras | High cameras | Low events | High events | Overlap [95% CI] | Kernel comparison p |
|---|---|---:|---:|---:|---:|---:|---:|
| `Female` | `Human abundance` | `15` | `15` | `344` | `126` | `0.797 [0.726, 0.868]` | `0.000` |
| `Female` | `Wolf abundance` | `15` | `15` | `216` | `611` | `0.894 [0.843, 0.945]` | `0.060` |
| `Female` | `Jackal abundance` | `15` | `15` | `465` | `107` | `0.847 [0.768, 0.927]` | `0.095` |
| `Female` | `Boar abundance` | `15` | `15` | `213` | `433` | `0.795 [0.726, 0.864]` | `0.000` |
| `Female` | `Cropland` | `15` | `15` | `564` | `341` | `0.737 [0.683, 0.792]` | `0.000` |
| `Female` | `Built and roads` | `15` | `15` | `312` | `124` | `0.697 [0.613, 0.780]` | `0.000` |
| `Female` | `Firing zones` | `15` | `15` | `245` | `491` | `0.851 [0.793, 0.909]` | `0.000` |
| `Female` | `Minefields` | `15` | `15` | `659` | `223` | `0.877 [0.822, 0.932]` | `0.050` |
| `Female` | `Nature reserves` | `15` | `15` | `319` | `450` | `0.795 [0.738, 0.852]` | `0.000` |
| `Female` | `Open spaces` | `15` | `15` | `149` | `733` | `0.849 [0.779, 0.919]` | `0.000` |
| `Female` | `Orchards` | `15` | `15` | `349` | `165` | `0.704 [0.622, 0.785]` | `0.000` |
| `Fawn` | `Human abundance` | `15` | `15` | `41` | `13` | `NA` | `NA` |
| `Fawn` | `Wolf abundance` | `15` | `15` | `30` | `73` | `0.867 [0.751, 0.984]` | `0.824` |
| `Fawn` | `Jackal abundance` | `15` | `15` | `54` | `14` | `NA` | `NA` |
| `Fawn` | `Boar abundance` | `15` | `15` | `33` | `60` | `0.850 [0.731, 0.968]` | `0.709` |
| `Fawn` | `Cropland` | `15` | `15` | `59` | `30` | `0.785 [0.663, 0.907]` | `0.236` |
| `Fawn` | `Built and roads` | `15` | `15` | `30` | `41` | `0.777 [0.637, 0.918]` | `0.055` |
| `Fawn` | `Firing zones` | `15` | `15` | `39` | `45` | `0.828 [0.700, 0.955]` | `0.477` |
| `Fawn` | `Minefields` | `15` | `15` | `50` | `15` | `NA` | `NA` |
| `Fawn` | `Nature reserves` | `15` | `15` | `49` | `68` | `0.831 [0.705, 0.956]` | `0.201` |
| `Fawn` | `Open spaces` | `15` | `15` | `5` | `90` | `NA` | `NA` |
| `Fawn` | `Orchards` | `15` | `15` | `33` | `19` | `NA` | `NA` |
| `Male` | `Human abundance` | `15` | `15` | `642` | `339` | `0.847 [0.795, 0.899]` | `0.000` |
| `Male` | `Wolf abundance` | `15` | `15` | `372` | `853` | `0.899 [0.854, 0.944]` | `0.030` |
| `Male` | `Jackal abundance` | `15` | `15` | `855` | `226` | `0.866 [0.809, 0.924]` | `0.010` |
| `Male` | `Boar abundance` | `15` | `15` | `416` | `620` | `0.869 [0.824, 0.913]` | `0.000` |
| `Male` | `Cropland` | `15` | `15` | `721` | `560` | `0.906 [0.864, 0.948]` | `0.020` |
| `Male` | `Built and roads` | `15` | `15` | `388` | `318` | `0.766 [0.710, 0.822]` | `0.000` |
| `Male` | `Firing zones` | `15` | `15` | `512` | `771` | `0.900 [0.855, 0.944]` | `0.005` |
| `Male` | `Minefields` | `15` | `15` | `743` | `392` | `0.911 [0.871, 0.950]` | `0.136` |
| `Male` | `Nature reserves` | `15` | `15` | `474` | `626` | `0.910 [0.870, 0.949]` | `0.095` |
| `Male` | `Open spaces` | `15` | `15` | `354` | `912` | `0.891 [0.844, 0.937]` | `0.005` |
| `Male` | `Orchards` | `15` | `15` | `732` | `301` | `0.862 [0.810, 0.914]` | `0.000` |
