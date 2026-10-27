#' hotcoldESF: covariate-adjusted hot and cold spots
#'
#' The main function is [hot_cold_map()], which tests the residuals of a
#' regression of area-level counts on covariates for spatial autocorrelation,
#' estimates their spatial component by eigenvector spatial filtering, and
#' labels each area hot, cold or neutral with a permutation band.
#' [rank_cells()] turns the result into a priority list. [local_gi_star()] and
#' [local_moran()] compute the usual local statistics, and [sim_thomas()]
#' simulates clustered data on a grid for testing.
#'
#' @keywords internal
"_PACKAGE"
