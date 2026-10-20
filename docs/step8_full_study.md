# Step 8: full simulation study (19–20 Oct 2026)

Script: `simulation/17_full_study.R` (2400 datasets, 157 min on 2 cloud cores). Summary: `simulation/18_full_study_summary.R` → `simulation/results/17_summary.csv`, `paper/tables/sim_ranking.tex`, `paper/tables/sim_labels.tex`, figures fig15_full_ranking and fig16_full_labels.

## Design

- Proposed method (ESF): log(Y + 1), covariate adjustment, Moran gate (5%), ESF with smooth eigenvectors (eigenvalue ≥ 0.25 × largest), ridge, penalty by 10-fold CV, permutation band (hot above the 95% quantile, cold below the 5%), ranking by largest excess.
- Benchmarks on the same residuals: Gi* (no gate), Gi* + gate, Gi* + gate + FDR, LISA + gate, top residuals, top counts.
- Scenarios (2000 random points + a covariate area of 3000 points in all, model receives the covariate except in U): B compact clusters; L large diffuse; E elongated; M mixed sizes; S many small; K clusters + a cold area (keeps 50% / 20% of points); U covariate not in the model; C covariate area only; N random points. Weak and strong versions.
- Grids 20×20 (100 replicates, 200 permutations) and 40×40 (50 replicates, 100 permutations). ESF run time: median 1.2 s (20×20), 16 s (40×40).
- Simulator: `sim_thomas()` gained cold areas (`n_cold`, `cold_radius`, `cold_keep`); old seeds still give the same data.

## Ranking (precision of the top 5% of cells)

- Strong clusters (B, L, M, K, 40×40 E): ESF, Gi* and LISA all 0.96–1.00.
- Weak compact clusters (B weak): Gi* is best (0.81 / 0.59 on 20×20 / 40×40), then LISA (0.75 / 0.52), then ESF (0.70 / 0.43).
- Mixed sizes (M): all ≥ 0.94; ESF finds the radius-50 cluster on 40×40 more often (0.96 strong, 0.88 weak) than Gi* (0.68 / 0.62) and LISA (0.76 / 0.86), and the radius-100 cluster on 20×20 (0.97–1.00 vs Gi* 0.79–0.90).
- Small clusters (S) and elongated clusters on 20×20: no smoothing (top residuals) or LISA does best; ESF and Gi* are lower.
- Cold area (K): ranking unaffected (0.99–1.00).
- Covariate missing from the model (U): every method fails (top 5% precision ≤ 0.16); the covariate area dominates the ranking.
- Top counts fails everywhere (≤ 0.32): covariate adjustment is essential.

## Labels

- Hot cells inside a cluster: ESF 0.77 vs Gi* 0.67 for strong compact clusters on 20×20; otherwise ESF, Gi* and LISA are within about 0.05, with ESF slightly lower on weak clusters. Gi* + gate + FDR has the most precise labels (0.88–1.00) but finds few cells with weak clusters (0.01–0.15).
- Cells labelled hot with no clusters (% of the map): ESF 0.2–2.2, Gi* + gate 0.2–1.9, LISA + gate 0.1–1.3, Gi* without gate 4.0–5.6, FDR 0.0.
- Cold area (K strong): ESF labels 97% / 88% of the cold area cold (Gi* 100% / 96%); 79% / 58% of ESF's cold cells are inside the cold area (Gi* 70% / 58%, LISA 81% / 74%).
- U: about half of every method's hot cells are in the covariate area.

## Moran gate: how often it passes

- Scenarios with clusters: 86–100%, except S weak (64% / 76%) and B weak (89% / 86%).
- N (random points): 3% on 20×20, 16% on 40×40 (8 of 50). The 40×40 rate is above the nominal 5%; with about 2 points per cell, log(Y + 1) is far from normal, and the test uses moments derived under normality.
- C (covariate area only): 19% / 32%. The covariate enters linearly, but its effect on log counts is not linear, so the residuals keep some structure around the covariate area.

## Conclusions for the paper

1. ESF performs as well as Gi* and LISA when clusters are strong, and controls false hot cells as well as the gated benchmarks.
2. ESF's distinctive strength is clusters of different sizes: it finds the smaller clusters more often.
3. Gi* is stronger with weak compact clusters; no smoothing is best for very small clusters.
4. Covariate adjustment is the single most important step (top counts and U fail).
5. Limitations to report: the gate's false-pass rate is above 5% with sparse counts (N 40×40) and when the covariate effect is not linear (C). Possible fixes: a permutation version of the Moran test; a spline of each covariate instead of a linear term.
