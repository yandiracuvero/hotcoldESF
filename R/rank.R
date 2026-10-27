# Step 5: priority list for a budget of m cells ----------------------------------
#
# Every method returns a full ordering of the n cells (best first), so any
# budget m can be evaluated. Ties are broken by Y.

#' @noRd
harmonic <- function(D, cells = seq_len(nrow(D))) {
  Ds <- D[cells, cells, drop = FALSE]
  if (sum(Ds) == 0) return(rep(0, length(cells)))
  g <- igraph::graph_from_adjacency_matrix(Ds, mode = "undirected")
  as.numeric(igraph::harmonic_centrality(g))
}

#' Priority list of cells
#'
#' Orders all cells from the highest to the lowest priority, so that any budget
#' of m cells can be served by the first m. Ties are broken by the counts.
#' @param obj an object of class `"hotcold"` from [hot_cold_map()].
#' @param method `"R2"` (default): hot cells first, ordered by their estimated
#'   excess \eqn{\hat f_i}, then all other cells by \eqn{\hat f_i};
#'   `"R3"`: hot zones by their total excess, then cells by harmonic
#'   centrality within the zone; `"R1"`: hot cells by harmonic centrality within
#'   the hot cells; `"R0"`: harmonic centrality in the contiguity graph, without
#'   labels; `"B0"`: the counts (a naive baseline).
#' @return Integer vector of cell indices, from the highest to the lowest priority.
#' @examples
#' set.seed(1)
#' s <- sim_thomas(k = 15, n_parents = 2, radius = 120, sons = 300,
#'                 n_background = 1500, n_cov = 1500)
#' h <- hot_cold_map(s$Y, s$Z, s$A, nperm = 199)
#' top <- rank_cells(h)[1:10]
#' h$labels[top]
#' @export
rank_cells <- function(obj, method = c("R2", "R3", "R1", "R0", "B0")) {
  method <- match.arg(method)
  n <- obj$n; Y <- obj$Y; f <- obj$fhat; hot <- obj$labels == "hot"

  # order by (primary, secondary, ...) descending; Y is always the last tie-breaker
  ord <- function(...) do.call(order, c(lapply(list(...), function(x) -x),
                                         list(-Y)))
  switch(method,
    B0 = ord(Y),
    R0 = ord(harmonic(obj$D)),
    R1 = {
      h <- rep(0, n); h[hot] <- harmonic(obj$D, which(hot))
      ord(hot, h, f)                                 # hot cells first
    },
    R2 = ord(hot, f),
    R3 = {
      z  <- obj$zones$hot
      zs <- tapply(f, z, sum)                        # total excess per zone
      zone_score <- ifelse(is.na(z), -Inf, zs[as.character(z)])
      h <- rep(0, n)
      for (id in unique(stats::na.omit(z))) {
        cells <- which(z == id)
        h[cells] <- harmonic(obj$D, cells)
      }
      ord(hot, zone_score, h, f)
    })
}
