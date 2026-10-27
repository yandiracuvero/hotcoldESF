# hotspots

The clean project for the paper *Reconstructing latent spatial relationships in urban pests to control for potential public health risks* (resubmission to SERRA).

The repository root is also an R package (`hotcoldESF`). The simulation, the application and the paper all use **one** algorithm, described in `docs/algorithm.md`.

```
Hotspot/
├── README.md          this file
├── DESCRIPTION        R package metadata (hotcoldESF; dependencies: glmnet, igraph)
├── NAMESPACE
├── R/                 the package code: the single implementation of the algorithm
│   ├── moran.R          Moran's I for regression residuals (exact moments)
│   ├── esf.R            ESF(D): MI-Lasso selection, post-Lasso fit, reconstruct_B()
│   ├── rounding.R       rounding rules: none, quantile, mixture, permutation
│   ├── find_hotspots.R  main algorithm (Steps 0-4) + hot/cold labels + zones
│   ├── rank.R           priority lists R0-R3 and baseline B0 (Step 5)
│   └── simulate.R       queen grid, modified Thomas process, SAR on a latent graph
├── tests/testthat/    unit tests (nesting, termination, sign invariance, ...)
├── docs/              algorithm spec, decisions, notes
├── simulation/        simulation study scripts and results
├── application/       Madrid cockroach application
├── paper/             LaTeX source of the resubmission + response letter
└── archive/           pointers to old code (originals stay in 02SepTest)
```

## How to use it (in RStudio, from this folder)

The package is called `hotcoldESF` (the name `hotspots` is already taken on CRAN).

```r
devtools::load_all()          # load the package without installing
devtools::test()              # run the unit tests
devtools::check()             # R CMD check

s <- sim_thomas(k = 20, n_parents = 3, radius = 100, sons = 200,
                n_background = 2000, n_cov = 3000)
h <- hot_cold_map(s$Y, s$Z, s$A)   # defaults = the procedure of the paper
h                                   # summary
rank_cells(h)[1:20]                 # top-20 cells of the priority list
```

The functions of the earlier graph reconstruction (`find_hotspots()`, the `round_*()` rules,
`reconstruct_B()`) are no longer exported; they stay in `R/` so that the appendix scripts can
be rerun with `devtools::load_all()`.

## Plan and status

- [x] Steps 1-3c: algorithm written down, core functions + tests, small tests, signal test, latent-graph test (the reconstruction does not recover W*: moved to an appendix)
- [x] Step 4: one-pass hot/cold map with a permutation band; Gi*/LISA benchmarks (`docs/step4_hotcold.md`)
- [x] Steps 5-7: penalty choice (ridge, smooth eigenvectors, CV), Moran gate, cluster shapes (`docs/step5*`-`step7*`)
- [x] Step 8: full simulation study, 2400 datasets (`simulation/17_full_study.R`, `docs/step8_full_study.md`)
- [x] Step 9: gate checks (`docs/step9_gate.md`)
- [x] Steps 10-12: Madrid application; nested models (exposure, environment, rat reports); covariate relevance (`application/`, `docs/step10*`-`step12*`)
- [ ] Open data questions: D_GatoOK vs D_Gatos; units of Densidad; deduplication of repeated reports
- [ ] Grid-size sensitivity (500 m, 1 km, 2 km) for the Madrid models
- [ ] Revise the paper + point-by-point response (manuscript number, tracked changes)

## Where the old code is

The originals are untouched in `~/Documents/UJI/02SepTest`:
`Codigo/50Location.Rmd` (application), `Paquete/old/R` (first package + simulation),
`Paquete/SpatialStructure` (new package skeleton), and `CodigoV1`–`V5` (older versions).
