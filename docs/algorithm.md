# Algorithm: merged version (draft 2.1, 17 Oct 2026). Implemented in R/

This is the single algorithm that the package, the simulation, the application and the paper will all use.
- **[TEST]** marks a choice that the small test (plan step 3) or the simulation will settle.
- **[DISCUSS]** marks a choice that is still open in conversation.

Sources: nested pruning comes from `Paquete/old/R/matrix2sf.R`. The fixed λ and the Moran test with the original A come from `Codigo/50Location.Rmd`. The rounding rule is a pluggable argument.

**Decisions taken (17 Oct):** Q2 uses current edges only. Q3 uses queen everywhere. Q4 continues with a warning flag. Nested pruning confirmed. Hot/cold labelling added (Step 4).

---

## Inputs

| Symbol | Meaning |
|---|---|
| `Y` (n) | counts per cell, standardised |
| `Z` (n × k) | covariates, standardised (may be empty: intercept only) |
| `A` (n × n) | binary, symmetric queen contiguity matrix, zero diagonal |
| `α` | significance level for Moran's I (default 0.05) |
| `weight` | how B is built: `"lambda"` (Δ, published) or `"gamma2"` (γ̂²) **[TEST]** |
| `R(·)` | rounding rule: takes a vector of values and returns a threshold τ **[TEST]** |
| `m` | budget: number of cells to return |
| `K_max` | safety cap on iterations (default: number of edges in A) |

## Output

- the final graph `D*`;
- the fitted spatial component `f̂` (n) and the label of every cell: **hot / cold / neutral**;
- **hot zones** and **cold zones**: connected components of the labelled subgraphs;
- a priority list of m cells (Step 5, **[DISCUSS]**);
- a trace of every iteration: number of edges, |S|, Moran p-value, τ, and the reason for stopping;
- the flag `esf_inadequate` (Q4).

---

## Step 0: Is there anything to explain?

1. `u = M_Z Y`, the OLS residuals of Y on (1, Z).
2. Compute the standardised Moran's I of `u` with **A**. This gives `z₀` and `p₀`.
3. If `p₀ > α`, **stop**. There is no residual spatial structure, so there are no hot or cold spots to report.
4. Fix `λ = |z₀|^(-2)` **once**. It is never recomputed.

## Step 1: ESF fit for a graph D (the subroutine `ESF(D)`)

1. Compute the eigendecomposition `D = E Δ Eᵀ`.
2. Run the Lasso (glmnet) of Y on (Z unpenalised, E penalised) with the fixed λ. The selected set is `S = {j : γ̂ⱼ ≠ 0}`.
3. Run post-Lasso OLS of Y on (1, Z, E_S). This gives `γ̂_S` (post-Lasso) and residuals `e`.
4. Compute Moran's I of `e` **with the original A**, giving `p`.
5. Return `(E_S, Δ_S, γ̂_S, p)`.

## Step 2: Initialise

`D⁽⁰⁾ = A`, `fit⁽⁰⁾ = ESF(A)`.
- If `S` is empty, stop with "no eigenvectors selected for A".
- If `p⁽⁰⁾ ≤ α`, set `esf_inadequate = TRUE` and continue (Q4). The simulation counts how often this happens.

## Step 3: Iterate (k = 0, 1, 2, …)

1. **Reconstruct**
   - `weight = "lambda"`: `B = E_S Δ_S E_Sᵀ` (published version);
   - `weight = "gamma2"`: `B = E_S diag(γ̂_S²) E_Sᵀ`, which weights each pattern by its importance for Y.

   Both are invariant to flipping the sign of an eigenvector. `diag(γ̂)` without the square is **not** used, because it would depend on the arbitrary sign that `eigen()` returns.
2. **Values to round**: `v = { b_ij : i < j, d_ij⁽ᵏ⁾ = 1, b_ij > 0 }`. These are only the edges still present, and only positive values.
3. **Threshold**: `τ⁽ᵏ⁾ = R(v)`.
4. **Nested pruning**: `D_c = D⁽ᵏ⁾ ∘ 1(B ≥ τ⁽ᵏ⁾)`. The *previous* graph is multiplied by the new 0/1 matrix, so an edge can never come back.
   (Bug in `50Location.Rmd`: it multiplies by `morg`, the original A, instead of the previous graph.)
5. **Stop checks**: if `D_c = D⁽ᵏ⁾`, stop with "converged". If `D_c` has no edges, stop with "empty graph". In both cases return `D⁽ᵏ⁾`.
6. **Adequacy**: `fit_c = ESF(D_c)`. If `S` is empty, stop with "no eigenvectors". If `p_c ≤ α`, stop with "Moran rejects". In both cases return `D⁽ᵏ⁾` and `fit⁽ᵏ⁾`.
7. **Accept**: `D⁽ᵏ⁺¹⁾ = D_c`, `fit⁽ᵏ⁺¹⁾ = fit_c`, and go to 3.1.

**Proposition (termination).** Every accepted step removes at least one edge, so the loop stops in at most |E(A)| iterations, for any rounding rule.

**Remark.** With a quantile rule, τ is above the minimum, so every step removes edges (unless many b_ij are tied, as happens in small symmetric components). In practice, the Moran test is what stops the loop. Mixture and permutation rules can keep everything, so real convergence is possible with them.

## Step 4: Hot / cold labelling (new)

**Why it is needed.** B is blind to sign. `b_ik > 0` means that i and k are on the *same side* of the selected patterns: both high or both low. So pruning keeps coherent zones, **hot and cold alike**, and cuts the edges on the borders between them. The sign lives in the fitted spatial component, which is well defined: if eⱼ flips, so does γ̂ⱼ.

```
f̂ = E_S γ̂_S − mean(E_S γ̂_S)   (from the last accepted fit)
```

f̂ is centred because the leading eigenvector of a connected graph is nearly constant and therefore collinear with the intercept. Only deviations from the mean are meaningful. (Found while implementing.)

`f̂ᵢ` is the number of sightings in cell i that the covariates do not explain and that has spatial structure.

1. **Label** each cell:
   - hot if `f̂ᵢ > c`;
   - cold if `f̂ᵢ < −c`;
   - neutral otherwise.

   Default `c = 0`. **[TEST]** Alternatives are c = sd(f̂)/2 or a permutation band.
2. **Hot subgraph** `D_hot = D*` restricted to the hot cells. Its connected components (size ≥ 2) are the **hot zones**. Cold zones are built the same way.
3. **Check (simulation)**: with no labelling (the published version), how many cold cells enter the top m? This measures how much cold zones contaminated the old ranking.

## Step 5: Priority list for a budget of m cells **[DISCUSS]**

These candidates are compared in the simulation with the same top-m metric:

| Variant | Rule | Meaning |
|---|---|---|
| R0 published | harmonic centrality on `D*`, no labels | position in the whole retained graph |
| R1 | harmonic centrality on `D_hot` | position inside hot zones only |
| R2 | rank hot cells by `f̂` | largest unexplained excess first |
| R3 zones first | rank hot zones by total excess `Σ f̂`; within each zone, order cells by harmonic centrality | allocate to the most important zone first; inside a zone, start from its most central cell |
| B0 baseline | rank cells by Y | naive strategy |

---

## Rounding rules (the argument `R`) **[TEST]**

| Name | Rule | Note |
|---|---|---|
| `quantile(q)` | τ = q-quantile of v | published version; q is a free choice |
| `mixture` | τ = antimode between the two components of a 2-component Gaussian mixture on v; if one component is preferred by BIC, keep everything | stable, formal version of the valley rule; replaces `locmodes(mod0 = 5)` |
| `permutation(B, 0.95)` | Freedman–Lane: permute the residuals of Y on Z, recompute B, and take the pooled 95th percentile as τ | τ means "stronger than chance" |
| `stability(B, π)` | run the algorithm on B bootstrap versions and keep edges present in ≥ π of the final graphs | edge-inclusion probabilities; run once at the end |

## Paper fixes that follow from this document

- The pseudocode loop condition is reversed (`or p.val ≤ α/2`). Replace it with Step 3.
- The Moran test is done with **A**, not D⁽ᵏ⁾.
- Write the Lasso objective as glmnet and Barde et al.'s reference code define it: `(1/2n)‖Y − Zβ − Eγ‖² + λ‖γ‖₁`, with λ = |z₀|^(-2).
- B entries are not restricted to (−1, 1).
- Add the termination proposition.
- Explain why B is sign-blind and why the labelling step is needed. This answers the editor's point 2 (the rationale for each step).
