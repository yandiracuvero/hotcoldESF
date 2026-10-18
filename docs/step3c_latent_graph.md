# Step 3c: does the reconstruction recover a latent graph? (17 Oct 2026)

Script: `simulation/03_latent_graph.R` (180 datasets, 9.5 min). Results: `simulation/results/03_latent_graph.csv`.

**Data-generating process.** Queen grid 20 × 20 (A, 1482 edges). Y is a SAR process on a latent graph W* ⊂ A:
- *barriers*: a river with one bridge and a highway, about 85 edges cut;
- *random*: 20% of edges missing;
- *none*: W* = A.

ρ ∈ {0.5, 0.8, 0.95}, with a covariate gradient. Truth per edge: in W* (keep) or barrier (cut).

## Threshold-free signal: AUC for separating W* edges from barrier edges (0.5 = chance)

| ρ | b_ij from ESF (Δ) | b_ij from ESF (γ̂²) | trivial: −\|e_i − e_j\| (OLS residuals) |
|---|---|---|---|
| 0.50 | 0.55 | 0.55 | 0.55 |
| 0.80 | 0.55 | 0.53 | 0.65 |
| 0.95 | **0.45** | 0.61 | **0.96** |

The random-removal case is at 0.50 for every method, so it is unidentifiable, as expected.

## Full algorithm, barriers (J = cut rate + keep rate − 1; 0 = chance, 1 = perfect)

- The best variant (Griffith 0.25, Δ, mixture) reaches J = 0.05 / 0.11 / 0.14 for ρ = 0.5 / 0.8 / 0.95. The trivial baseline cutting the same number of edges gets 0.02 / 0.11 / **0.32**.
- The other variants are at J ≈ 0.
- **Control (no barriers): the algorithm still removes 0–99% of the edges**, depending on the variant (q75 with all eigenvectors keeps only 1–27%).

## Why

With MI-Lasso at a = 2, about 130–300 of the 400 eigenvectors are selected. B = Σ_selected λ_j e_j e_jᵀ is then A minus a handful of unselected global patterns. Its values on edges are unrelated to local differences in Y: cor(b_ij, |e_i − e_j|) = 0.05 in an example with ρ = 0.95. The eigenvectors describe the geometry of A, not the local discontinuities that a barrier creates.

## Conclusion

In its current form, the reconstruction does **not** recover latent spatial relationships, not even strong barriers that a trivial neighbour-difference rule detects with AUC 0.96. The paper cannot claim that the reconstructed graph reflects the connections pests use.
