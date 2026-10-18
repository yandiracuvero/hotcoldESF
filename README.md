# hotspots

The clean project for the paper *Reconstructing latent spatial relationships in urban pests to control for potential public health risks* (resubmission to SERRA).

The repository root is also an R package (`hotspots`). The simulation, the application and the paper all use **one** algorithm, described in `docs/algorithm.md`.

```
Hotspot/
├── README.md          this file
├── DESCRIPTION        R package metadata (dependencies: glmnet, igraph)
├── NAMESPACE
├── R/                 the package code: the single implementation of the algorithm
│   ├── moran.R          Moran's I for regression residuals (exact moments)
│   ├── esf.R            ESF(D): MI-Lasso selection, post-Lasso fit, reconstruct_B()
│   ├── rounding.R       rounding rules: none, quantile, mixture, permutation
│   ├── find_hotspots.R  main algorithm (Steps 0-4) + hot/cold labels + zones
│   ├── rank.R           priority lists R0-R3 and baseline B0 (Step 5)
│   └── simulate.R       queen grid + modified Thomas process
├── tests/testthat/    unit tests (nesting, termination, sign invariance, ...)
├── docs/              algorithm spec, decisions, notes
├── simulation/        simulation study scripts and results
├── application/       Madrid cockroach application
├── paper/             LaTeX source of the resubmission + response letter
└── archive/           pointers to old code (originals stay in 02SepTest)
```

## How to use it (in RStudio, from this folder)

```r
devtools::load_all()          # load the package without installing
devtools::test()              # run the unit tests

s   <- sim_thomas(k = 20, n_parents = 3, radius = 100, sons = 1000)
res <- find_hotspots(s$Y, Z = NULL, A = s$A,
                     weight = "lambda",              # or "gamma2"
                     rounding = round_quantile(0.75)) # or round_mixture(), round_permutation()
res                             # summary
rank_cells(res, "R3")[1:20]     # top-20 cells (zones first, then centrality)
```

## Plan and status

- [x] **Step 1:** merged algorithm written down (`docs/algorithm.md`, draft 2)
- [x] **Step 2:** core function in `R/` + unit tests (all passing)
- [x] **Step 3:** small test (`simulation/01_small_test.R`). Findings and open questions are in `docs/step3_findings.md`
- [x] **Step 3b:** signal strength × cell size × exponent a × log scale (`simulation/02_signal_test.R`, run on the Mac). Findings: `docs/step3b_findings.md`
- [ ] **Decision:** how to frame the contribution (see step3b findings) + add Gi*/LISA baselines
- [ ] **Step 4:** full simulation study: top-m metric, ranking variants R0–R3, baselines (top counts B0, Getis-Ord Gi*, LISA), CSR null scenario, scenario **with covariates**
- [ ] **Step 5:** Madrid application, rerun with the package:
  - hot/cold map per summer (2010–2013)
  - **persistence map**: zones that are hot in several summers
  - **out-of-sample table**: choose m cells with year t and count the sightings of year t+1 (2010→11, 11→12, 12→13), for each ranking variant and baseline. Caveat for the paper: possible interventions between years are unknown and affect all methods equally, so the comparison between methods stays fair, but absolute capture rates may be underestimated.
- [ ] **Step 6:** revise the paper + point-by-point response (manuscript number, tracked changes)

## Where the old code is

The originals are untouched in `~/Documents/UJI/02SepTest`:
`Codigo/50Location.Rmd` (application), `Paquete/old/R` (first package + simulation),
`Paquete/SpatialStructure` (new package skeleton), and `CodigoV1`–`V5` (older versions).
