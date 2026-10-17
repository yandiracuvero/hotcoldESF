# Simulation helpers ---------------------------------------------------------------
#
# Ported from Paquete/SpatialStructure (sim_thomas_modified, processSim), without
# the sf/spatstat dependency: the window is a square and the grid is regular.
# Cells are numbered row by row from the bottom-left, like the application.

#' Queen contiguity matrix of a k x k regular grid
#' @param k number of cells per side.
#' @param queen TRUE for queen (8 neighbours), FALSE for rook (4).
#' @export
grid_adjacency <- function(k, queen = TRUE) {
  n  <- k * k
  rc <- cbind(row = (seq_len(n) - 1) %/% k, col = (seq_len(n) - 1) %% k)
  dr <- abs(outer(rc[, "row"], rc[, "row"], "-"))
  dc <- abs(outer(rc[, "col"], rc[, "col"], "-"))
  A  <- if (queen) (pmax(dr, dc) == 1) else (dr + dc == 1)
  A  <- A * 1
  dimnames(A) <- NULL
  A
}

#' Simulate a modified Thomas process and aggregate it on a grid
#'
#' Parents are uniform in the window, at least `radius` from the border.
#' Each parent has `sons` offspring uniform in a disc of radius `radius`
#' (as in nclust / sim_thomas_modified).
#'
#' @param k cells per side (the grid has k^2 cells). NOTE: k is the number of
#'   cells, not the cell side. Cell side = size / k.
#' @param n_parents number of clusters (0 gives complete spatial randomness
#'   with `n_background` points).
#' @param radius cluster radius (same units as `size`).
#' @param sons offspring per parent.
#' @param n_background extra uniform points (noise).
#' @param size side of the square window.
#' @param queen contiguity type.
#' @return list: `Y` counts per cell, `A`, `k`, `cell_side`, `parents`
#'   (coordinates), `parent_cells`, `points`.
#' @export
sim_thomas <- function(k = 20, n_parents = 3, radius = 100, sons = 1000,
                       n_background = 0, size = 1000, queen = TRUE) {
  px <- stats::runif(n_parents, radius, size - radius)
  py <- stats::runif(n_parents, radius, size - radius)
  pts <- do.call(rbind, lapply(seq_len(n_parents), function(i) {
    r  <- radius * sqrt(stats::runif(sons))
    th <- stats::runif(sons, 0, 2 * pi)
    cbind(px[i] + r * cos(th), py[i] + r * sin(th))
  }))
  if (n_background > 0)
    pts <- rbind(pts, cbind(stats::runif(n_background, 0, size),
                            stats::runif(n_background, 0, size)))
  if (is.null(pts)) pts <- matrix(numeric(0), ncol = 2)
  cell_of <- function(x, y) {
    cx <- pmin(floor(x / (size / k)), k - 1)
    cy <- pmin(floor(y / (size / k)), k - 1)
    cy * k + cx + 1
  }
  Y <- tabulate(cell_of(pts[, 1], pts[, 2]), nbins = k * k)
  list(Y = Y, A = grid_adjacency(k, queen), k = k, cell_side = size / k,
       parents = cbind(x = px, y = py),
       parent_cells = if (n_parents > 0) cell_of(px, py) else integer(0),
       points = pts)
}
