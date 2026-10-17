# Step 5: priority list for a budget of m cells ----------------------------------
#
# Spec: docs/algorithm.md, Step 5. Every method returns a full ordering of the
# n cells (best first), so any budget m can be evaluated. Ties are broken by Y.

#' @noRd
harmonic <- function(D, cells = seq_len(nrow(D))) {
  Ds <- D[cells, cells, drop = FALSE]
  if (sum(Ds) == 0) return(rep(0, length(cells)))
  g <- igraph::graph_from_adjacency_matrix(Ds, mode = "undirected")
  as.numeric(igraph::harmonic_centrality(g))
}

#' Rank cells
#'
#' @param obj output of [find_hotspots()].
#' @param method one of
#'   "R0" harmonic centrality on the final graph, no labels (published);
#'   "R1" harmonic centrality on the hot subgraph;
#'   "R2" hot cells by fhat;
#'   "R3" hot zones by total excess, then cells by centrality within the zone;
#'   "B0" raw counts Y (naive baseline).
#' @return integer vector: cell indices ordered from highest to lowest priority.
#' @export
rank_cells <- function(obj, method = c("R0", "R1", "R2", "R3", "B0")) {
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
