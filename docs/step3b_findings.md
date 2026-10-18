# Step 3b: signal strength × cell size × penalty exponent × response scale (17 Oct 2026)

Script: `simulation/02_signal_test.R`, run on Yandira's Mac (120 datasets, 51.7 min).
Raw results: `simulation/results/02_signal_test.csv`. Budget m = 5% of the cells.
Clusters have radius 100. k = cells per side (cell side 100 / 50 / 25 for k = 10 / 20 / 40), s = points per cluster.

## Precision in S2 (share of the top-m cells inside a true cluster)

| Method | k10 s50 | k10 s200 | k10 s500 | k20 s50 | k20 s200 | k20 s500 | k40 s50 | k40 s200 | k40 s500 | mean |
|---|---|---|---|---|---|---|---|---|---|---|
| B0 top counts | 0.00 | 0.04 | 0.54 | 0.00 | 0.08 | 0.52 | 0.01 | 0.17 | 0.60 | 0.22 |
| B1 top residuals (std) | 0.32 | 0.72 | 0.96 | 0.30 | 0.74 | 0.96 | 0.21 | 0.64 | 0.90 | 0.64 |
| **B1 top residuals (log)** | **0.76** | **1.00** | 0.98 | **0.50** | 0.98 | 1.00 | 0.30 | 0.79 | 0.99 | **0.81** |
| ESF log, a=2, q75, R2 | 0.30 | 1.00 | 0.98 | 0.24 | 0.98 | 1.00 | 0.09 | 0.80 | 0.99 | 0.71 |
| ESF log, a=2, no pruning, R2 | 0.20 | 1.00 | 0.98 | 0.24 | 0.98 | 1.00 | 0.09 | 0.78 | 0.99 | 0.70 |
| ESF log, a=1, no pruning, R2 | 0.00 | 0.48 | 0.90 | 0.00 | 0.96 | 1.00 | 0.01 | 0.89 | 0.98 | 0.58 |
| **ESF std, a=2, no pruning, R3 (zones)** | 0.34 | 0.68 | 0.84 | 0.46 | 0.77 | 0.85 | **0.62** | **0.89** | 0.79 | 0.69 |
| published (std, a=2, q75, R0) | 0.10 | 0.20 | 0.38 | 0.12 | 0.24 | 0.22 | 0.10 | 0.25 | 0.08 | 0.19 |

## S3 (nothing is hot): share of the top-m cells in the covariate bump

- B0 0.98, B1 std 0.56, **B1 log 0.01**.
- With `response = "log"`, **Step 0 correctly says "no residual spatial structure" in 26 of 30 datasets**. With standardised counts it never does (0 of 30).

## Findings

1. **The log scale is essential.** On the count scale the covariate effect is non-linear, and the leftover looks like a hot spot. On log(Y + 1) the covariate is absorbed and the method correctly finds nothing in S3. Recommendation: `response = "log"` becomes the default, and the application is rerun on log counts.
2. **The published pipeline (R0) is no better than top counts** (precision 0.19 vs 0.22).
3. **Iterative pruning adds nothing measurable.** q75 vs no pruning differ by at most 0.01 in mean precision, in every configuration.
4. **ESF beats the residual baseline only on fine grids with weak clusters** (k = 40, clusters spread over about 50 cells). There, zone ranking R3 reaches 0.62 / 0.89 vs 0.30 / 0.79 for B1 log. Elsewhere B1 log is as good or better, sometimes by a lot (k10 s50: 0.76 vs 0.30). This matches the hypothesis: spatial smoothing helps only when single-cell residuals are noisy *and* clusters cover many cells.
5. **Exponent a:** a = 0.5 selects no eigenvectors (useless). a = 1 helps on fine grids and hurts on coarse ones. a = 2 is the most robust. The number of eigenvectors selected grows with n (a = 2, log: 47 / 247 / 954 for n = 100 / 400 / 1600).

## Implications (to discuss with co-authors)

- The contribution has to be reframed. The value is (i) covariate adjustment on the log scale, (ii) the hot/cold labels, and (iii) ESF spatial smoothing with zone ranking, which wins where data are noisy and fine-grained. The graph reconstruction is the original core of the paper, but it shows no benefit in ranking. It could stay as a zone-delineation tool or be dropped.
- Before claiming (iii), Getis-Ord Gi* and local Moran (LISA) on the log residuals must be added as baselines, since they also smooth spatially.
