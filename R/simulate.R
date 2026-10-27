# Simulation helpers ---------------------------------------------------------------
#
# Ported from Paquete/SpatialStructure (sim_thomas_modified, processSim), without
# the sf/spatstat dependency: the window is a square and the grid is regular.
# Cells are numbered row by row from the bottom-left, like the application.

#' Contiguity matrix of a regular grid
#'
#' Binary contiguity matrix of a k x k grid, with cells numbered row by row from
#' the bottom-left corner.
#' @param k number of cells per side.
#' @param queen `TRUE` (default) for queen contiguity (cells sharing a border or a
#'   corner), `FALSE` for rook contiguity (sharing a border).
#' @return A k^2 x k^2 binary symmetric matrix with zero diagonal.
#' @examples
#' A <- grid_adjacency(4)
#' table(rowSums(A))    # corners 3, borders 5, interior 8 neighbours
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

#' Simulate clustered points and count them on a grid
#'
#' Simulates a modified Thomas cluster process in a square window and counts the
#' points in the cells of a k x k grid. Parents are uniform in the window, at
#' least `radius` from the border; each parent has `sons` offspring uniform in a
#' disc of radius `radius`. Optionally adds uniform background points, points
#' that follow an environmental covariate, and cold areas where points are
#' thinned. Used in the simulation study of the package's paper.
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
#' @param aspect elongation of the clusters (1 = discs). Each cluster becomes an
#'   ellipse with the same area, semi-axes radius*sqrt(aspect) and
#'   radius/sqrt(aspect), at a random angle. `radius` and `sons` may be vectors
#'   (one value per cluster).
#' @param n_cold,cold_radius,cold_keep cold areas: `n_cold` discs of radius
#'   `cold_radius` where only a share `cold_keep` of the points is kept (placed
#'   away from the clusters and the covariate area). Returned as `cold_cells`.
#' @param n_cov points whose intensity follows an environmental covariate
#'   (a Gaussian bump placed at least `min_dist` from every parent). These
#'   points are explained by Z, so they are NOT hot spots once Z is used.
#' @param cov_strength,cov_sd strength and width of the bump.
#' @param min_dist minimum distance between parents, and between the parents and
#'   the centre of the covariate bump.
#' @return list: `Y` counts per cell, `Z` (covariates: `env` bump + `noise`;
#'   NULL if n_cov = 0), `A`, `k`, `cell_side`, `centres`, `footprint` (cells
#'   inside a cluster disc), `cov_cells` (cells where the bump > 0.5),
#'   `parents`, `parent_cells`, `points`.
#' @examples
#' set.seed(7)
#' s <- sim_thomas(k = 20, n_parents = 3, radius = 100, sons = 200,
#'                 n_background = 2000, n_cov = 3000)
#' sum(s$Y)                      # 2000 + 3000 + 3 * 200 points
#' image(matrix(log1p(s$Y), 20, 20), asp = 1, axes = FALSE)
#' @export
sim_thomas <- function(k = 20, n_parents = 3, radius = 100, sons = 1000,
                       n_background = 0, size = 1000, queen = TRUE,
                       n_cov = 0, cov_strength = 3, cov_sd = 150,
                       min_dist = 300, aspect = 1,
                       n_cold = 0, cold_radius = 150, cold_keep = 0.2) {
  radius <- rep_len(radius, max(n_parents, 1)); sons <- rep_len(sons, max(n_parents, 1))
  ax <- radius * sqrt(aspect); ay <- radius / sqrt(aspect)   # semi-axes, same area
  px <- stats::runif(n_parents, ax, size - ax)
  py <- stats::runif(n_parents, ax, size - ax)
  ang <- if (aspect != 1) stats::runif(n_parents, 0, pi) else rep(0, n_parents)
  pts <- do.call(rbind, lapply(seq_len(n_parents), function(i) {
    r  <- sqrt(stats::runif(sons[i]))
    th <- stats::runif(sons[i], 0, 2 * pi)
    u <- ax[i] * r * cos(th); v <- ay[i] * r * sin(th)
    cbind(px[i] + u * cos(ang[i]) - v * sin(ang[i]),
          py[i] + u * sin(ang[i]) + v * cos(ang[i]))
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
    tries <- 0
    repeat {
      tries <- tries + 1
      if (tries > 5000) stop("no room for the covariate area: lower min_dist")
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

  # cold areas: discs where only a share `cold_keep` of the points is kept, placed
  # away from the clusters and from the covariate area
  cold_centers <- matrix(numeric(0), ncol = 2)
  if (n_cold > 0) {
    for (c0 in seq_len(n_cold)) {
      tries <- 0
      repeat {
        tries <- tries + 1
        if (tries > 5000) stop("no room for the cold area")
        cc0 <- stats::runif(2, cold_radius, size - cold_radius)
        far_p <- n_parents == 0 || all(sqrt((px - cc0[1])^2 + (py - cc0[2])^2) > ax + cold_radius + 50)
        far_c <- is.null(cov_center) || sqrt(sum((cov_center - cc0)^2)) > cold_radius + 2 * cov_sd
        far_o <- nrow(cold_centers) == 0 ||
          all(sqrt((cold_centers[, 1] - cc0[1])^2 + (cold_centers[, 2] - cc0[2])^2) > 2 * cold_radius)
        if (far_p && far_c && far_o) break
      }
      cold_centers <- rbind(cold_centers, cc0)
      inside <- (pts[, 1] - cc0[1])^2 + (pts[, 2] - cc0[2])^2 <= cold_radius^2
      drop <- inside & stats::runif(nrow(pts)) > cold_keep
      pts <- pts[!drop, , drop = FALSE]
    }
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
  inside <- function(i) {                          # cells whose centre is in cluster i
    dx <- cc[, 1] - px[i]; dy <- cc[, 2] - py[i]
    u <- dx * cos(ang[i]) + dy * sin(ang[i]); v <- -dx * sin(ang[i]) + dy * cos(ang[i])
    which((u / ax[i])^2 + (v / ay[i])^2 <= 1)
  }
  clusters <- lapply(seq_len(n_parents), inside)
  foot <- sort(unique(unlist(clusters)))

  list(Y = Y, Z = Z, A = grid_adjacency(k, queen), k = k, cell_side = h,
       centres = cc, footprint = foot, cov_center = cov_center, cov_cells = cov_cells,
       parents = cbind(x = px, y = py), radius = radius, aspect = aspect, angle = ang,
       clusters = clusters,
       cold_cells = if (n_cold > 0) which(apply(cc, 1, function(p)
         any((p[1] - cold_centers[, 1])^2 + (p[2] - cold_centers[, 2])^2 <= cold_radius^2))) else integer(0),
       parent_cells = if (n_parents > 0) cell_of(px, py) else integer(0),
       points = pts)
}

#' Simulate a spatial autoregressive process on a latent graph with barriers
#'
#' The observed contiguity A is a queen grid. The process only flows through a
#' latent graph W* (a subgraph of A): A minus "barrier" edges. Y follows a SAR
#' model on W*: Y = rho W Y + Z beta + e, with W the row-standardised W*.
#'
#' @param k cells per side.
#' @param rho spatial autoregressive parameter.
#' @param type "barriers": a vertical river with one bridge plus a horizontal
#'   highway on the left half; "random": a random share `p_remove` of edges is
#'   missing; "none": W* = A (control).
#' @param p_remove share of edges removed when type = "random".
#' @param beta coefficients of Z = (smooth gradient, noise).
#' @return list: `Y`, `Z`, `A`, `Wstar` (binary latent graph), `k`.
#' @noRd
sim_sar_barrier <- function(k = 20, rho = 0.8, type = c("barriers", "random", "none"),
                            p_remove = 0.2, beta = c(1, 0)) {
  type <- match.arg(type)
  A  <- grid_adjacency(k)
  n  <- k * k
  rc <- cbind(row = (seq_len(n) - 1) %/% k, col = (seq_len(n) - 1) %% k)
  W  <- A
  if (type == "barriers") {
    river <- floor(k / 2) - 1                      # cut between col river and river+1
    bridge <- floor(k / 4)                         # one row keeps its crossing
    hwy   <- floor(0.6 * k)                        # cut between row hwy and hwy+1
    for (i in seq_len(n)) for (j in which(A[i, ] == 1)) {
      ci <- rc[i, "col"]; cj <- rc[j, "col"]; ri <- rc[i, "row"]; rj <- rc[j, "row"]
      cross_river <- (min(ci, cj) == river && max(ci, cj) == river + 1) &&
                     !(ri == bridge && rj == bridge)
      cross_hwy   <- (min(ri, rj) == hwy && max(ri, rj) == hwy + 1) &&
                     max(ci, cj) <= river
      if (cross_river || cross_hwy) W[i, j] <- 0
    }
  } else if (type == "random") {
    up <- which(upper.tri(A) & A == 1)
    drop <- sample(up, round(p_remove * length(up)))
    W[drop] <- 0
    W <- W * t(W)                                  # keep symmetry
  }
  Wr <- W / pmax(rowSums(W), 1)
  Z  <- cbind(grad = as.numeric(scale(rc[, "col"] + rc[, "row"])),
              noise = stats::rnorm(n))
  Y  <- solve(diag(n) - rho * Wr, Z %*% beta + stats::rnorm(n))
  list(Y = as.numeric(Y), Z = Z, A = A, Wstar = W, k = k)
}
