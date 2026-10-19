# Step 7: cluster shapes and sizes, ESF vs Gi* (19 Oct 2026)

Decision (Yandira, 19 Oct): keep ESF as the main method and report Gi* on residuals as a benchmark.
Proposed ESF version: smooth eigenvectors (eigenvalue ≥ 0.25 × largest), ridge, penalty by 10-fold CV, Moran gate, permutation band.

Code: `sim_thomas()` now accepts `aspect` (elliptical clusters, random angle) and per-cluster `radius` / `sons`; old seeds give the same data.
Experiment: `simulation/15_cluster_shapes.R` (240 datasets, 31 min). Figure: `simulation/16_cluster_shapes_figure.R` → fig14_cluster_shapes.

## Scenarios (all with 2000 random points and a covariate area of 3000 points)

- Test 3 (expected to favour ESF): **L** 2 large diffuse clusters, radius 200; **E** 3 elongated clusters (ellipses 300 × 33, same area as radius 100); **M** mixed sizes, radii 50 / 100 / 200, same intensity.
- Test 4 (expected to favour Gi*): **S** 8 small clusters, radius 40.
- Weak and strong versions; 20×20 (20 replicates) and 40×40 (10 replicates).

## Results

| | ESF ridge | Gi* + gate | Top residuals |
|---|---|---|---|
| L, top 5% precision (20 / 40 weak) | 0.99 / 0.84 | 1.00 / 0.91 | 0.85 / 0.57 |
| E, top 5% precision (20 / 40 weak) | 0.70 / 0.72 | 0.70 / 0.76 | 0.71 / 0.46 |
| M, clusters with a cell in the top 5% (20 weak / 40 weak) | 0.70 / 0.90 | 0.60 / 0.83 | 0.85 / 1.00 |
| M, radius-50 cluster found (20 strong / 40 strong / 40 weak) | 0.20 / 1.00 / 0.70 | 0.20 / 0.70 / 0.50 | 0.80 / 1.00 / 1.00 |
| M, radius-100 cluster found (20 strong / 20 weak) | 1.00 / 0.95 | 0.85 / 0.75 | 1.00 / 0.95 |
| S, top 5% precision (20 strong / 40 weak) | 0.49 / 0.26 | 0.44 / 0.32 | 0.67 / 0.19 |
| S, clusters with a cell in the top 5% (40 weak) | 0.72 | 0.84 | 0.88 |

Labels (hot cells inside a cluster, cluster cells labelled hot) are within a few points of each other for ESF ridge and Gi* + gate in every scenario. Gi* + gate + FDR has the most precise labels but finds few cells with weak clusters.

## Findings

1. **ESF and Gi* are close in most settings.** Large, elongated and small clusters give differences of 0.00–0.07 in most cells of the table.
2. **ESF's advantage: clusters of different sizes.** When small and large clusters coexist, ESF finds the smaller ones more often (radius 100 on 20×20: 0.95–1.00 vs 0.75–0.85; radius 50 on 40×40: 0.70–1.00 vs 0.50–0.70). Gi*'s fixed window is dominated by the large cluster's cells.
3. **Gi*'s advantage: weak signal on the fine grid** (L weak 40×40: 0.91 vs 0.84; S weak 40×40: 0.32 vs 0.26, clusters found 0.84 vs 0.72), and weak clusters in experiment 13.
4. **Very small clusters favour no smoothing at all.** Top residuals finds the most small clusters (S, and the radius-50 cluster on 20×20), because both ESF and Gi* average over a neighbourhood larger than the cluster.
5. **Elongated clusters favour neither method.**
