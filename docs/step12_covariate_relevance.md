# Step 12: models with the documented covariates, and covariate relevance (20 Oct 2026)

Script: `application/04_models_covariates.R` (2 min). Results: `application/results/madrid_cov_*.csv`; figures `madrid_cov_models.png`, `madrid_cov_relevance.png`.

Cockroaches April–September, 2010–2012 pooled (2013 out-of-sample); rats all months 2010–2012.
- M1 exposure: log residential points, Densidad, log cell area.
- M2: M1 + dist_verd, AreaVer, YCons, dist_agua, dist_fuent, D_Mercados, D_GatoOK, Euros, Desempleo.
- M2 variant with D_Gatos; M3 = M2 + log rat reports; M2 + field-campaign variables (sensitivity); rats with M2.

| Model | covariates explain | Moran z | hot | cold | reports in hot |
|---|---|---|---|---|---|
| M1 exposure | 80.0% | 10.8 | 66 | 64 | 45% |
| M2 + environment | 82.8% | 8.7 | 66 | 63 | 37% |
| M2, D_Gatos variant | 82.8% | 8.7 | 68 | 61 | 36% |
| M3 + rat reports | 87.5% | 5.5 | 58 | 63 | 33% |
| M2 + field variables | 83.4% | 8.4 | 59 | 61 | 36% |
| Rats, M2 | 86.0% | 8.6 | 57 | 60 | 32% |

Overlap of hot cells (Jaccard): M2 vs D_Gatos variant 0.91; M2 vs M2 + field 0.81; M1 vs M2 0.55; M2 vs M3 0.53; cockroach M2 vs rat M2 0.28.

## Relevance (M2, cockroaches; coefficients per SD of log reports, ESF-adjusted)

| Covariate | coef (SE) | p | same sign in all 4 summers | drop-one overlap | rats coef (p) | verdict |
|---|---|---|---|---|---|---|
| YCons | −0.11 (0.02) | <0.001 | yes | 0.80 | −0.15 (<0.001) | relevant: older buildings, more reports |
| AreaVer | +0.11 (0.02) | <0.001 | yes | 0.81 | +0.10 (<0.001) | relevant |
| dist_fuent | −0.12 (0.03) | <0.001 | yes | 0.81 | −0.02 (0.59) | relevant (cockroaches only) |
| Desempleo | +0.12 (0.04) | 0.001 | yes | 0.93 | −0.02 (0.60) | relevant (cockroaches only) |
| dist_verd | +0.07 (0.04) | 0.06 | yes | 0.83 | +0.08 (0.02) | borderline; relevant for rats |
| D_GatoOK | −0.06 (0.03) | 0.07 | yes | 0.91 | −0.03 (0.39) | borderline |
| dist_agua | −0.03 (0.03) | 0.32 | yes | 0.86 | −0.10 (<0.001) | not for cockroaches; relevant for rats |
| D_Mercados | +0.01 (0.04) | 0.69 | no | 0.97 | −0.07 (0.06) | not relevant |
| Euros | +0.01 (0.04) | 0.81 | no | 0.91 | −0.06 (0.15) | not relevant |

Exposure: log residential 0.63 (p ≈ 0), Densidad 0.14 (p < 0.001); log area not significant (kept as a design control).
Plain OLS p-values are much smaller (spatial autocorrelation inflates them): e.g. dist_agua p = 0.0005, Euros 0.04.

Field-campaign variables add 0.6 points of explained variance and change the hot map little (overlap 0.81).

Out-of-sample 2013, 25 / 50 cells: M1 19% / 33%, M2 13% / 28%, M3 11% / 25%, M2 + field 14% / 29%, top counts 29% / 50% of sightings.

## Final M2 (20 Oct 2026): D_Mercados and Euros removed

M2 = exposure (log residential, Densidad, log area) + dist_verd, AreaVer, YCons, dist_agua, dist_fuent, D_GatoOK, Desempleo.
- M2: covariates explain 82.7%, Moran z 9.1, 70 hot / 61 cold, hot cells hold 38% of reports.
- Coefficients barely change; dist_verd becomes significant (+0.08, p = 0.04); all seven keep the same sign in every summer.
- M3: 87.5%, z 5.6, 59 hot. Rats M2: 85.9%, z 8.8, 62 hot. D_Gatos variant: 69 hot.
- Out-of-sample 2013 (25 / 50 cells): M2 13% / 29% of sightings.
