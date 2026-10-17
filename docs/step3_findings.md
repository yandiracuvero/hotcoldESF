# Step 3: small test, findings (17 Oct 2026)

Script: `simulation/01_small_test.R` (20 replicates × 3 scenarios, 400 cells, top m = 20).
Tables: `Rscript simulation/01_summary.R`. Raw results: `simulation/results/01_small_test.csv`.

Scenarios:
- S1: 3 clusters + uniform noise, no covariates.
- S2: S1 + points driven by an environmental covariate; Z is given to the method.
- S3: covariate-driven points only. Nothing should be hot once Z is used.

## Main table (means over 20 replicates)

precision = share of the top 20 cells inside a true cluster; recall = share of parents found;
cov_fp = share of the top 20 cells in the covariate bump; cold = share of the top 20 cells labelled cold.

| Variant | S1 precision | S1 recall | S1 cold | S2 precision | S2 recall | S2 cov_fp | S2 cold | S3 cov_fp |
|---|---|---|---|---|---|---|---|---|
| B0 top counts | 1.00 | 1.00 | – | 0.57 | 0.98 | 0.44 | – | 1.00 |
| B1 top OLS residuals on Z | 1.00 | 1.00 | – | 0.96 | 1.00 | 0.04 | – | 0.54 |
| **published** (all eigvecs, Δ, q75, R0) | **0.28** | 0.60 | **0.57** | **0.24** | 0.52 | 0.10 | **0.50** | 0.43 |
| no pruning, R2 (hot cells by f̂) | 1.00 | 1.00 | 0 | 0.92 | 1.00 | 0.07 | 0 | 0.54 |
| no pruning, R3 | 0.76 | 0.58 | 0 | 0.76 | 0.63 | 0.06 | 0 | 0.00 |
| Δ, q75, R2 | 1.00 | 1.00 | 0 | 0.96 | 1.00 | 0.04 | 0 | 0.74 |
| Δ, q75, R3 | 0.84 | 0.90 | 0 | 0.85 | 0.97 | 0.03 | 0 | 0.69 |
| γ̂², q75, R3 | 0.90 | 0.80 | 0 | 0.94 | 0.97 | 0.04 | 0 | 0.76 |
| Griffith 0.25, Δ, q50, R3 | 0.78 | 0.58 | 0 | 0.75 | 0.60 | 0.00 | 0 | 0.83 |

## Findings

1. **The published pipeline fails.** About half of its top 20 cells are *cold*. This confirms the sign-blindness of B. Adding the hot/cold labels raises precision from about 0.25 to 0.85–1.00.
2. **The Lasso hardly selects.** z0 ≈ 25, so λ = |z0|^(-2) ≈ 0.0015, and 360 of 400 eigenvectors are selected. The model is almost saturated, so the Moran check almost always passes and the graph shrinks to 7–12% of its edges. The open question is whether glmnet's scaling (1/2n) matches Barde et al.'s theory. **Ask Rowan.**
3. **Griffith's candidate rule** (eigenvalue ≥ 0.25 × max) avoids the saturation: about 80 eigenvectors are selected. But then the Moran test rejects the first pruning in most runs, so the graph stays close to A and R3 gets worse.
4. **Pruning adds a little; ranking by f̂ adds the most.** Without any pruning, R2 already reaches 0.92 (S2). Pruning with q75 brings it to 0.96. **A trivial baseline (B1, top OLS residuals) also reaches 0.96.** In these easy scenarios (500 points per cluster), single-cell residuals are clean, so spatial smoothing has nothing to add.
5. **S3: no method is "conditional on Z" in practice.** The covariate acts as exp(3·bump), but the model is linear in Z. The leftover non-linearity looks like a hot spot: 50–85% of every method's top 20 falls in the bump, B1 included, and Step 0 never says "no structure". Options: model log(Y + 1), or add flexible terms in Z.

## Proposed next experiments (to discuss)

- **Signal strength:** offspring per cluster 50 / 100 / 200 / 500, and cluster radius vs cell size. With weak, noisy clusters, single-cell residuals (B1) should degrade, and the spatially smoothed f̂ is where ESF should show its value. This is the scenario the paper needs.
- **Response scale:** standardised Y (published) vs log(Y + 1), for S3.
- **λ:** the published |z0|^(-2) vs a scaling agreed with Rowan.
