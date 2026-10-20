# Step 11: Madrid, nested models with exposure and rat reports (20 Oct 2026)

Script: `application/03_models_exposure.R` (40 s). Results: `application/results/madrid_models_*.csv`; figure `application/figures/madrid_models.png`.

Cockroaches April–September, main map 2010–2012 pooled (5795 reports), 2013 for out-of-sample. Rats: all months, 2010–2012 (5030).

- M1 exposure: log(1 + residential points of TODO_MADRID119m in the cell), population density, log cell area.
- M2: M1 + provisional risk factors (D_AlcSuc, D_AlcDef, D_RsuDefi, D_Mercados, D_GatoOK, dist_agua, dist_verd, AreaVer, YCons, Euros, Desempleo).
- M3: M2 + log(1 + rat reports).
- Rats with M2, for comparison.

| | covariates explain | Moran z | hot | cold | reports in hot cells |
|---|---|---|---|---|---|
| M1 exposure | 80% | 10.8 | 66 | 64 | 45% |
| M2 + environment | 83% | 8.5 | 66 | 61 | 39% |
| M3 + rat reports | 87% | 5.1 | 56 | 56 | 32% |
| Rats, M2 | 86% | 8.9 | 55 | 62 | 30% |

- Exposure alone explains 80% of the variation in log cockroach reports (residential points: standardised coefficient 0.80). The original 19 covariates did not include it.
- Clear spatial structure remains in every model, so hot spots are not only exposure.
- Overlap of hot cells (Jaccard): M1–M2 0.63, M2–M3 0.58, M1–M3 0.47; cockroach M2 vs rat M2 only 0.22. Once exposure is removed, the two species' hot spots differ.
- Environmental risk factors add little (3 points); the largest coefficients are dist_verd, AreaVer, dist_agua.
- Out-of-sample (2013 sightings in 25 / 50 cells): M1 19% / 33%, M2 14% / 30%, M3 12% / 26%, top counts 2010–12 29% / 50%. Raw sightings follow exposure, so counts predict them best; hot spots measure excess, a different target.
