# hotcoldESF 0.1.0

* Initial CRAN submission.
* `hot_cold_map()` maps covariate-adjusted hot and cold spots: Moran test of the
  regression residuals, eigenvector spatial filter with a cross-validated ridge
  penalty, and a pointwise permutation band.
* `rank_cells()` turns the map into a priority list.
* Benchmarks `local_gi_star()` and `local_moran()`, and the simulator `sim_thomas()`.
* Datasets `madrid`, `madrid_A` and `madrid_boundary`: cockroach and rat reports in
  Madrid, 2010-2013, on a 1 km grid, with covariates.
* Vignettes: a worked example on a simulated dataset, and the Madrid application.
