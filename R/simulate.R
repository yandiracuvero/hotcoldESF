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
#' @param n_cov points whose intensity follows an environmental covariate
#'   (a Gaussian bump placed at least `min_dist` from every parent). These
#'   points are explained by Z, so they are NOT hot spots once Z is used.
#' @param cov_strength,cov_sd strength and width of the bump.
#' @return list: `Y` counts per cell, `Z` (covariates: `env` bump + `noise`;
#'   NULL if n_cov = 0), `A`, `k`, `cell_side`, `centres`, `footprint` (cells
#'   inside a cluster disc), `cov_cells` (cells where the bump > 0.5),
#'   `parents`, `parent_cells`, `points`.
#' @export
sim_thomas <- function(k = 20, n_parents = 3, radius = 100, sons = 1000,
                       n_background = 0, size = 1000, queen = TRUE,
                       n_cov = 0, cov_strength = 3, cov_sd = 150,
                       min_dist = 300) {
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

  # covariate-driven points: an environmental "bump" away from the parents
  cov_center <- NULL
  bump <- function(x, y) exp(-((x - cov_center[1])^2 + (y - cov_center[2])^2) /
                               (2 * cov_sd^2))
  if (n_cov > 0) {
    repeat {
      cov_center <- stats::runif(2, cov_sd, size - cov_sd)
      if (n_parents == 0 ||
          all(sqrt((px - cov_center[1])^2 + (py - cov_center[2])^2) > min_dist)) break
    }
    cp <- matrix(numeric(0), ncol = 2)
    while (nrow(cp) < n_cov) {                      # rejection sampling
      x <- stats::runif(4 * n_cov, 0, size); y <- stats::runif(4 * n_cov, 0, size)
      keep <- stats::runif(4 * n_cov) < exp(cov_strength * (bump(x, y) - 1))
      cp <- rbind(cp, cbind(x, y)[keep, , drop = FALSE])
    }
    pts <- rbind(pts, cp[seq_len(n_cov), , drop = FALSE])
  }

  cell_of <- function(x, y) {
    cx <- pmin(floor(x / (size / k)), k - 1)
    cy <- pmin(floor(y / (size / k)), k - 1)
    cy * k + cx + 1
  }
  Y <- tabulate(cell_of(pts[, 1], pts[, 2]), nbins = k * k)

  # cell centres; observed covariates = bump value + an irrelevant noise covariate
  h  <- size / k
  cc <- cbind(x = ((seq_len(k * k) - 1) %% k + 0.5) * h,
              y = ((seq_len(k * k) - 1) %/% k + 0.5) * h)
  Z  <- NULL; cov_cells <- integer(0)
  if (n_cov > 0) {
    b  <- bump(cc[, "x"], cc[, "y"])
    Z  <- cbind(env = b, noise = stats::rnorm(k * k))
    cov_cells <- which(b > 0.5)
  }
  # cells whose centre lies inside some cluster disc (the true hot footprint)
  foot <- if (n_parents > 0)
    which(apply(cc, 1, function(p) any((p[1] - px)^2 + (p[2] - py)^2 <= radius^2)))
  else integer(0)

  list(Y = Y, Z = Z, A = grid_adjacency(k, queen), k = k, cell_side = h,
       centres = cc, footprint = foot, cov_center = cov_center, cov_cells = cov_cells,
       parents = cbind(x = px, y = py),
       parent_cells = if (n_parents > 0) cell_of(px, py) else integer(0),
       points = pts)
}
