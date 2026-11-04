# Dynamic predator circular GAMM diagnostics

## Predictor construction

Predator timestamps from cameras 22 and 24 were corrected by +148 and +238 days, respectively, before active-round assignment and apparent-solar-time conversion. Records were collapsed to 30-min independent events within camera, Camelot session, and taxon.
One-hour predator event rates were decomposed into between-camera means, between-round deviations within cameras, and hourly deviations within camera-rounds. The cyclic predator smooths used the final hourly component.
A three-hour circular moving-window rate was fitted as a sensitivity to sparse one-hour predator counts.

## Model diagnostics

| Class | Window | Structure | Theta | Deviance explained | Dispersion | Lag-1 residual r | Min k-index | Min k P | Max concurvity |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| Juvenile | 1h | random intercepts | 2.932 | 45.6% | 0.417 | 0.725 | 0.854 | <0.001 | 0.425 |
| Adult female | 1h | camera-specific cyclic deviations | 3.984 | 53.1% | 0.646 | 0.312 | 0.806 | <0.001 | 0.419 |
| Adult male | 1h | camera-specific cyclic deviations | 8.871 | 45.3% | 0.764 | 0.168 | 0.850 | <0.001 | 0.419 |
| Juvenile | 3h | random intercepts | 3.299 | 45.8% | 0.411 | 0.723 | 0.856 | <0.001 | 0.577 |
| Adult female | 3h | camera-specific cyclic deviations | 3.982 | 53.0% | 0.647 | 0.312 | 0.806 | <0.001 | 0.577 |
| Adult male | 3h | camera-specific cyclic deviations | 9.318 | 45.4% | 0.767 | 0.170 | 0.850 | <0.001 | 0.577 |

## Dynamic predator smooths

| Class | Window | Effect | EDF | Approximate P |
|---|---|---|---:|---:|
| Juvenile | 1h | Dynamic wolf by time | 0.00 | 0.668 |
| Juvenile | 1h | Dynamic jackal by time | 0.12 | 0.305 |
| Juvenile | 1h | Dynamic wolf x jackal by time | 0.00 | 0.508 |
| Adult female | 1h | Dynamic wolf by time | 0.00 | 0.799 |
| Adult female | 1h | Dynamic jackal by time | 0.00 | 0.820 |
| Adult female | 1h | Dynamic wolf x jackal by time | 0.32 | 0.280 |
| Adult male | 1h | Dynamic wolf by time | 0.71 | 0.208 |
| Adult male | 1h | Dynamic jackal by time | 0.00 | 0.455 |
| Adult male | 1h | Dynamic wolf x jackal by time | 0.00 | 0.424 |
| Juvenile | 3h | Dynamic wolf by time | 1.32 | 0.085 |
| Juvenile | 3h | Dynamic jackal by time | 0.00 | 0.793 |
| Juvenile | 3h | Dynamic wolf x jackal by time | 0.00 | 0.454 |
| Adult female | 3h | Dynamic wolf by time | 0.00 | 0.780 |
| Adult female | 3h | Dynamic jackal by time | 0.00 | 0.572 |
| Adult female | 3h | Dynamic wolf x jackal by time | 0.00 | 0.617 |
| Adult male | 3h | Dynamic wolf by time | 0.00 | 0.615 |
| Adult male | 3h | Dynamic jackal by time | 0.00 | 0.561 |
| Adult male | 3h | Dynamic wolf x jackal by time | 1.35 | 0.081 |

Smooth-term P values are approximate. Interpretation requires agreement across temporal windows, acceptable concurvity and basis checks, and biologically coherent prediction curves.
