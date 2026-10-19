# Step 6: Moran gate for Gi*, spatial CV and ridge for ESF (19 Oct 2026)

Code: `R/esf.R` (`alpha = 0` ridge option), `R/penalty.R` (`block_folds()`, `foldid`, `mix`), `hot_cold_map(fit_type, foldid)`.
Experiment: `simulation/13_gate_and_cv.R` (same 150 datasets, 19 min in the cloud). Figure: `simulation/14_gate_and_cv_figure.R` → fig13_gate_and_cv.

## Test 1: Gi* on residuals with the Moran gate and/or FDR

| | Gi* | Gi* + gate | Gi* + FDR | Gi* + gate + FDR |
|---|---|---|---|---|
| Top 5% precision, 20×20 B50 / 40×40 B50 | 0.81 / 0.58 | same | same | same |
| Hot cells inside a cluster, 20×20 A / B500 / B200 | 0.61 / 0.63 / 0.68 | same | 0.76 / 0.80 / 0.88 | same as FDR |
| Cluster cells labelled hot, 20×20 B50 / 40×40 B50 | 0.62 / 0.40 | 0.55 / 0.40 | 0.14 / 0.00 | 0.13 / 0.00 |
| Hot cells with no clusters, 20×20 C / N, 40×40 C | 20 / 17 / 81 | **0 / 0 / 18** | 0 / 0 / 0.2 | 0 / 0 / 0.2 |

The gate passed in 0 of 20 datasets for 20×20 C and N, 2 of 10 for 40×40 C, and 17 of 20 for 20×20 B50.

**The Moran gate alone gives Gi* the same false-positive control as the ESF procedure** (ESF with CV: 0 / 0 / 20–46). The gate does not change the Gi* ranking, so Gi* + gate keeps its advantage with weak clusters. FDR removes almost all false hot cells but also almost all weak-cluster cells.

## Test 2: ESF with smooth eigenvectors, penalty by CV

| | Lasso, random folds | Lasso, spatial folds | Ridge, random folds | Ridge, spatial folds |
|---|---|---|---|---|
| Top 5% precision, strong clusters | 0.98–1.00 | 0.93–1.00 | 0.99–1.00 | 0.99–1.00 |
| Top 5% precision, 20×20 B50 / 40×40 B50 | 0.52 / 0.36 | 0.36 / 0.17 | **0.67 / 0.47** | 0.66 / 0.47 |
| Hot cells inside a cluster, 20×20 strong | 0.72–0.77 | 0.65–0.77 | 0.76–0.79 | 0.77–0.78 |
| Hot cells with no clusters, 40×40 C | 46 | 20 | 22 | 22 |

- Spatial folds choose a larger penalty: the Lasso then selects 0–1 eigenvectors with weak clusters and loses power.
- Ridge keeps all smooth eigenvectors and shrinks them. It is the best ESF version so far: weak clusters 0.67 / 0.47 (Lasso CV 0.52 / 0.36; current rule 0.19 / 0.11). Fold type makes no difference for ridge.
- Ridge still trails Gi* with weak clusters (0.81 / 0.58).

## What this means

With the same covariate adjustment and the same Moran gate, Gi* on residuals is at least as good as the best ESF version on every measure except one: with strong clusters, more of ESF's hot cells lie inside a true cluster (0.76–0.79 vs 0.61–0.68). Gi* with FDR beats that too (0.76–0.88), at the cost of sensitivity with weak clusters.

So, in these simulations, the gains come from the pipeline (covariate adjustment, log scale, Moran gate, a local statistic, a ranking), not from ESF specifically. This needs a decision with the co-authors on how the paper frames its contribution.
