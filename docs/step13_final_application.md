# Step 13: final Madrid application (20 Oct 2026)

Script: `application/05_final_application.R` (about 13 min). Results `application/results/final_*.csv`; figures `final_summers.png`, `final_grid_sizes.png`; table `paper/tables/madrid_final_out_of_sample.tex`.
Final models: M1 exposure (log residential, Densidad, log area); M2 = M1 + dist_verd, AreaVer, YCons, dist_agua, dist_fuent, D_GatoOK, Desempleo; M3 = M2 + log rat reports.

## Per summer (1 km)
- Gate passes every summer and model (Moran z: M1 9.0–10.2, M2 7.2–9.4, M3 4.6–6.2).
- Hot cells per summer: M1 52–60, M2 54–66, M3 51–52; they hold 40–46% (M2) and 29–36% (M3) of reports.
- Hot again next summer: M2 61–77%, M3 52–59%.
- Persistence: M2 27 cells hot in all four summers, 70 in two or more; M3 12 and 60.

## Out-of-sample (mean of 2010→11, 11→12, 12→13, 2010–12→13)
| ranking | 25 cells: sightings / excess | 50 cells: sightings / excess |
|---|---|---|
| ESF M1 | 21% / 19% | 36% / 26% |
| ESF M2 | 17% / 16% | 32% / 30% |
| ESF M3 | 15% / 18% | 29% / 28% |
| Gi* (M2 residuals) | 17% / 12% | 30% / 24% |
| LISA (M2 residuals) | 21% / 29% | 35% / 41% |
| Top residuals (M2) | 21% / 33% | 38% / 48% |
| Top counts | 29% / 32% | 49% / 43% |
| Random cells | 5% / 0% | 10% / 0% |

## Grid size (2010–2012 pooled)
| grid | cells | cells without reports | Moran z M1 / M2 / M3 | hot M2 | share of city hot (M2) |
|---|---|---|---|---|---|
| 500 m | 1901 | 65% | 21.8 / 18.9 / 14.7 | 171 | 9% |
| 1 km | 493 | 53% | 10.8 / 9.1 / 5.6 | 70 | 14% |
| 2 km | 136 | 39% | 4.6 / 4.0 / 2.1 | 25 | 22% |

Overlap of hot areas on a common 250 m lattice (Jaccard), M2: 500 m–1 km 0.29, 1 km–2 km 0.36, 500 m–2 km 0.25 (M1 0.32 / 0.45 / 0.27; M3 0.28 / 0.19 / 0.14).
The main zones (south-central band, north-central area) appear at all three sizes, but their extent changes with the cell size; at 2 km the M3 signal is weak (z = 2.1, 16 hot cells).
