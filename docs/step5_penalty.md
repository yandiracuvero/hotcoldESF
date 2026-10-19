# Step 5: how the Lasso penalty is chosen (18 Oct 2026)

Code: `R/penalty.R` (`select_lambda()`), option `penalty` in `hot_cold_map()` (default unchanged: |z0|^-2).
Experiment: `simulation/10_penalty_rules.R` (same 150 datasets as experiment 07, 135 min in the cloud). Figures: `simulation/11_penalty_figures.R` → fig10_penalty_ranking, fig11_penalty_labels.

## What was tested

8 penalty rules × 2 candidate sets = 16 variants, plus Gi* on residuals as the reference.

- Rules: |z0|^-a with a = 1, 2 (default), 3, 4; 10-fold cross-validation (lowest error, and largest penalty within 1 SE); BIC along the Lasso path; "Moran stop" (largest penalty whose post-Lasso residuals pass the Moran test).
- Candidates: all eigenvectors, or smooth ones only (eigenvalue ≥ 0.25 × the largest, Griffith's rule; about 21% of them).
- The permutation band reuses the chosen penalty in every permutation. "Moran stop" with all eigenvectors was skipped on 40×40 (too slow).

## Precision of the top 5% (largest excess)

| | A 20 | B500 20 | B200 20 | B50 20 | B200 40 | B50 40 |
|---|---|---|---|---|---|---|
| Gi* residuals | 0.99 | 0.99 | 0.99 | **0.81** | 1.00 | **0.58** |
| \|z\|^-2, all (default) | 0.99 | 0.99 | 0.98 | 0.19 | 0.80 | 0.11 |
| \|z\|^-2, smooth | 1.00 | 0.99 | 0.97 | 0.27 | 0.98 | 0.15 |
| \|z\|^-3, smooth | 1.00 | 0.79 | 0.91 | 0.55 | 0.99 | 0.40 |
| \|z\|^-4, smooth | 1.00 | 0.79 | 0.89 | 0.53 | 0.91 | 0.45 |
| **CV lowest error, smooth** | 1.00 | 0.99 | 1.00 | 0.53 | 0.98 | 0.36 |
| Moran stop, smooth | 0.99 | 1.00 | 0.96 | 0.41 | 0.94 | 0.32 |
| CV 1 SE / BIC (either set) | 0.97–1.00 | 0.96–0.98 | 0.82–0.97 | 0.13–0.39 | 0.90–0.93 | 0.02–0.07 |

## Labels (CV lowest error, smooth vs Gi* residuals vs default)

| | CV min, smooth | Gi* residuals | default |
|---|---|---|---|
| Hot cells inside a cluster, 20×20 A / B500 / B200 / B50 | 0.77 / 0.73 / 0.73 / 0.44 | 0.61 / 0.63 / 0.68 / 0.64 | 0.98 / 0.99 / 0.97 / 0.33 |
| Cluster cells labelled hot, same | 1.00 / 1.00 / 0.95 / 0.46 | 1.00 / 1.00 / 0.98 / 0.62 | 0.62 / 0.59 / 0.56 / 0.16 |
| 40×40 B200 / B50: inside a cluster; labelled hot | 0.69 / 0.32; 0.89 / 0.34 | 0.78 / 0.50; 0.93 / 0.40 | 0.74 / 0.18; 0.47 / 0.08 |
| Hot cells with no clusters, 20×20 C / N, 40×40 C | 0 / 0 / 30 | 20 / 17 / 81 | 0 / 0 / 0 |

## Findings

1. **Smooth candidates help almost everywhere.** Restricting to eigenvectors with large positive eigenvalues removes the patterns that follow cell-level noise.
2. **Cross-validation (lowest error) with smooth candidates is the most stable rule.** It matches Gi* with strong clusters (0.98–1.00), and with weak clusters it is about 3 times the default (0.53 vs 0.19; 0.36 vs 0.11). It still trails Gi* there (0.81 and 0.58).
3. **Rules that penalise more (|z|^-1, CV 1 SE, BIC) fail with weak clusters**: they select 0 eigenvectors.
4. **Small penalties with smooth candidates (|z|^-3, |z|^-4) absorb the covariate.** They select almost all 84 smooth eigenvectors; the covariate bump lies in their span, so f̂ takes part of its effect (B500: 0.79, sensitivity 0.56).
5. **Labels.** With CV the band finds almost all cluster cells (sensitivity 0.95–1.00 with strong clusters) and keeps more of its hot cells inside clusters than Gi* on 20×20 (0.73–0.77 vs 0.61–0.68). It is lower than Gi* on 40×40.
6. **No-cluster scenarios.** In 20×20 C and N no variant labels any cell hot (Moran gate). On 40×40 C the gate passed in some datasets, and the variants with many eigenvectors then label about 30 cells hot (Gi*: 81).

## Questions for the co-authors

- Use cross-validation with smooth candidates as the penalty rule in the paper, instead of |z0|^-2?
- With weak clusters Gi* on residuals still ranks better. Is the paper's claim the ESF map with error control (Moran gate + band), with Gi* as a complementary ranking?
