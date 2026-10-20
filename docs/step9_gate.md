# Step 9: the Moran gate (20 Oct 2026)

Problem from the full study (step 8): with no clusters the gate passed 16% of datasets for random points on 40×40 (N) and 19% / 32% for the covariate area alone (C, 20×20 / 40×40).

Code: `moran_perm()` (Freedman–Lane permutation version of the Moran test), `covariate_basis()` (linear terms or a natural spline per covariate), and `hot_cold_map(gate = c("normal", "permutation"), spline_df = 0)`. Defaults are unchanged (normal test, linear covariates).
Experiments: `simulation/19_gate_fixes.R` (900 datasets of study 17, 63 min) and `simulation/20_gate_checks.R` (gate only).

## Findings

1. **Random points (N), 40×40: not a real problem.** With 500 replicates the normal gate passes 5.0% (95% CI 3.3–7.3%) and the permutation gate 5.4%. The 16% in the full study was sampling variation over 50 replicates. The permutation test brings no gain; keep the normal test.

2. **Covariate area alone (C): the cause is the linear covariate term.** Pass rate of the gate (normal test) by spline degrees of freedom:

| | linear | df 2 | df 3 | df 4 |
|---|---|---|---|---|
| C 20×20 (should be low) | 0.19 | 0.04 | 0.04 | 0.06 |
| C 40×40 | 0.32 | 0.04 | 0.04 | 0.04 |
| B weak 20×20 (should be high) | 0.89 | 0.82 | 0.80 | 0.76 |
| B weak 40×40 | 0.86 | 0.66 | 0.64 | 0.60 |
| S weak 20×20 | 0.64 | 0.30 | 0.28 | 0.30 |
| S weak 40×40 | 0.76 | 0.44 | 0.40 | 0.42 |

   Splines fix the false passes in C, but cost power with weak clusters: a flexible function of a smooth covariate can absorb part of the cluster signal (spatial confounding). Already with df = 2 the pass rate for weak small clusters halves.

3. **Effect on the full procedure (spline df 4 + permutation gate, experiment 19).**
   - False hot cells with C: ESF 5.8 → 1.0 cells (20×20), 35 → 4.5 (40×40).
   - ESF ranking with weak clusters drops: B weak 0.70 → 0.59 (20×20), 0.43 → 0.29 (40×40); S weak 0.17 → 0.08 and 0.25 → 0.16. Gi* and LISA lose less.
   - Strong clusters and the cold area: unchanged (0.95–1.00).

## Recommendation

Keep the defaults (normal test, linear covariates) for the main results and report the trade-off: when the covariate effect is not linear, the gate lets through 19–32% of datasets with no clusters. These false passes give about 1–2% of the map as false hot cells. A spline removes them, at the cost of power with weak clusters. In the application, a spline version can be run as a sensitivity analysis.
