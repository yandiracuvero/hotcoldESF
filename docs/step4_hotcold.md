# Step 4: hot/cold map with a permutation band vs Gi* and LISA (18 Oct 2026)

Code: `R/hotcold.R` (`hot_cold_map()`, `local_gi_star()`, `local_moran()`).
Experiment: `simulation/07_hotcold_benchmarks.R` (150 datasets, 21 min). Figures: `simulation/08_benchmark_figures.R` → fig6_labels, fig7_ranking.

## The procedure tested

1. y = standardised log(Y + 1); Z standardised.
2. Moran test of the residuals of y on Z (A, 5%). Not significant → no hot cells.
3. MI-Lasso ESF fit on all eigenvectors of A, λ = |z0|^(-2). The result is f̂, the spatial component.
4. Permutation band: 200 Freedman–Lane permutations (100 on 40×40), same λ and eigenvectors.
   - *pointwise*: hot if f̂_i > 95th percentile of its null, cold if below the 5th;
   - *max*: hot if above the 95th percentile of the null maximum over the map.
5. Zones = connected hot (cold) cells in A.

## Labels (hot vs not hot)

| | ESF band (pointwise) | Gi* residuals | LISA residuals | Gi* counts | LISA counts |
|---|---|---|---|---|---|
| Hot cells inside a cluster, 20×20 A / B500 / B200 / B50 | 0.98 / 0.98 / 0.96 / 0.34 | 0.61 / 0.63 / 0.68 / 0.64 | 0.68 / 0.73 / 0.75 / 0.71 | 0.66 / 0.46 / 0.25 / 0.00 | 0.76 / 0.49 / 0.17 / 0.00 |
| Cluster cells labelled hot, same | 0.61 / 0.59 / 0.58 / 0.17 | 1.00 / 1.00 / 0.98 / 0.62 | 1.00 / 1.00 / 0.96 / 0.47 | 1.00 / 0.99 / 0.42 / 0.00 | 1.00 / 0.91 / 0.23 / 0.00 |
| Hot cells inside a cluster, 40×40 B200 / B50 | 0.74 / 0.19 | 0.78 / 0.50 | 0.83 / 0.53 | 0.37 / 0.01 | 0.30 / 0.00 |
| Cluster cells labelled hot, 40×40 B200 / B50 | 0.48 / 0.09 | 0.93 / 0.40 | 0.85 / 0.31 | 0.63 / 0.01 | 0.40 / 0.00 |
| Cells labelled hot with no clusters, 20×20 C / N, 40×40 C | **0 / 0 / 0** | 20 / 17 / 81 | 12 / 9 / 56 | 51 / 21 / 179 | 46 / 8 / 146 |

The max band labelled almost no cells hot in any scenario (0–0.7 cells on average), so it was dropped.

## Ranking (precision of the top 5% of cells)

| Scenario | ESF largest excess | ESF hot zones first | Gi* residuals | LISA residuals | Top residuals | Top counts |
|---|---|---|---|---|---|---|
| 20×20 A | 1.00 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 |
| 20×20 B500 | 1.00 | 0.99 | 1.00 | 1.00 | 1.00 | 0.52 |
| 20×20 B200 | 0.98 | 0.96 | 0.99 | 1.00 | 0.98 | 0.04 |
| 20×20 B50 | 0.19 | 0.19 | **0.81** | 0.75 | 0.47 | 0.00 |
| 40×40 B200 | 0.82 | 0.87 | **1.00** | 1.00 | 0.83 | 0.16 |
| 40×40 B50 | 0.11 | 0.12 | **0.58** | 0.51 | 0.27 | 0.01 |

A quick check (20×20, B50, 10 replicates) of other ESF settings: log/counts × a = 1/2 × all/Griffith candidates. Precision was 0.00–0.42, against 0.80 for Gi* on residuals.

## Findings

1. **The ESF band controls false hot cells.** In the scenarios with no clusters, it labelled 0 cells hot, because the Moran test stops the procedure. Gi* and LISA on residuals labelled 9–81 cells hot there. With strong clusters, 96–98% of ESF hot cells are inside a cluster, against 61–75% for Gi*/LISA on residuals.
2. **The ESF band is conservative.** It labels about 60% of the cluster cells hot with strong clusters, and 9–17% with weak ones. Gi*/LISA on residuals label 93–100% with strong clusters.
3. **For ranking, Gi* on covariate-adjusted residuals equals or exceeds ESF in every scenario.** The gap is large with weak clusters: 0.81 vs 0.19 on 20×20, and 0.58 vs 0.11 on 40×40.
4. **Gi* and LISA on raw counts flag the covariate area as hot** (45–94% of their hot cells), and flag hot cells when no clusters exist. Covariate adjustment is the step that matters most.

## Open decision (for the co-authors)

ESF's distinctive contribution in these simulations is error control in the hot/cold labels: a Moran gate plus a permutation band. It is not ranking power. Options:

- (a) Keep ESF for the map and use Gi* on residuals for the priority list.
- (b) Gate Gi* on residuals with the same Moran test, and compare its labels to the ESF band.
- (c) Investigate why ESF loses power with weak clusters (with λ = |z0|^(-2), many eigenvectors are selected and f̂ also follows noise).
