# Main algorithm -----------------------------------------------------------------
#
# Spec: docs/algorithm.md, Steps 0-4. Ranking (Step 5) is in rank.R.

#' Find hot and cold spots by iterative ESF graph reconstruction
#'
#' @param Y response per cell (counts). It is standardised internally.
#' @param Z covariates per cell (n x k) or NULL. Standardised internally.
#' @param A binary symmetric contiguity matrix (queen), zero diagonal.
#' @param alpha significance level of the Moran's I tests.
#' @param weight "lambda" (published) or "gamma2"; see [reconstruct_B()].
#' @param rounding a rounding rule, e.g. [round_quantile()], [round_mixture()],
#'   [round_permutation()].
#' @param a exponent of the Lasso penalty, lambda = |z0|^(-a). 2 is the
#'   MI-Lasso default (Barde et al.); any positive value is allowed.
#' @param response "std": standardised counts (published); "log": standardised
#'   log(Y + 1).
#' @param label_c cells with fhat > label_c are hot, fhat < -label_c are cold.
#' @param max_iter cap on iterations (default: number of edges of A).
#' @param standardize passed to glmnet.
#' @param mc_min candidate eigenvectors; see [esf_fit()]. NULL = all (published).
#' @param verbose print one line per iteration.
#' @param keep_graphs store every accepted graph (and its fhat) in `out$graphs`,
#'   to plot the iterations.
#' @return an object of class "hotspots"; see the README for its fields.
#' @export
find_hotspots <- function(Y, Z = NULL, A, alpha = 0.05,
                          weight = c("lambda", "gamma2"),
                          a = 2, response = c("std", "log"),
                          rounding = round_quantile(0.75),
                          label_c = 0, max_iter = NULL,
                          standardize = TRUE, mc_min = 0.25,
                          verbose = FALSE, keep_graphs = FALSE) {
  weight <- match.arg(weight)
  stopifnot(isSymmetric(unname(A)), all(A %in% c(0, 1)), all(diag(A) == 0))
  A  <- unname(A)
  n  <- length(Y)
  response <- match.arg(response)
  Ys <- as.numeric(scale(if (response == "log") log1p(Y) else Y))
  Zs <- if (is.null(Z)) NULL else scale(as.matrix(Z))
  if (!is.null(Zs)) Zs <- Zs[, apply(Zs, 2, function(z) all(is.finite(z))), drop = FALSE]
  if (!is.null(Zs) && ncol(Zs) == 0) Zs <- NULL
  n_edges <- function(D) sum(D[upper.tri(D)])
  if (is.null(max_iter)) max_iter <- n_edges(A)

  # Step 0: is there residual spatial structure? fix lambda once ------------
  X0  <- design0(n, Zs)
  m0  <- moran_resid(stats::lm.fit(X0, Ys)$residuals, A, X0)
  out <- list(call = match.call(), a = a, response = response, n = n, A = A, D = A, Y = Y, z0 = m0$z,
              p0 = m0$p, lambda = NA, weight = weight,
              rounding = rule_label(rounding), alpha = alpha, mc_min = mc_min,
              esf_inadequate = NA, fit = NULL, trace = NULL,
              stop_reason = NA_character_)
  class(out) <- "hotspots"
  if (m0$p > alpha) {
    out$stop_reason <- "no residual spatial structure (step 0)"
    out$fhat   <- rep(0, n)
    out$labels <- label_cells(out$fhat, label_c)
    out$zones  <- find_zones(A, out$labels)
    return(out)
  }
  lambda     <- abs(m0$z)^(-a)
  out$lambda <- lambda

  # Step 2: initialise -------------------------------------------------------
  D   <- A
  graphs <- list()
  fit <- esf_fit(Ys, Zs, D, A, lambda, standardize, mc_min)
  out$esf_inadequate <- !isTRUE(fit$moran$p > alpha)
  trace <- data.frame(iter = 0, edges = n_edges(D), n_eig = length(fit$S),
                      moran_p = fit$moran$p, tau = NA_real_, accepted = TRUE)
  ctx <- list(Y = Ys, Z = Zs, A = A, lambda = lambda, weight = weight,
              standardize = standardize, mc_min = mc_min)
  stop_reason <- if (length(fit$S) == 0) "no eigenvectors selected for A" else NA

  # Step 3: iterate ------------------------------------------------------------
  k <- 0
  while (is.na(stop_reason)) {
    if (k >= max_iter) { stop_reason <- "max_iter"; break }
    B   <- reconstruct_B(fit, weight)
    up  <- upper.tri(D) & D == 1
    v   <- B[up]
    v   <- v[v > 0]
    if (length(v) == 0) { stop_reason <- "no positive values on edges"; break }
    ctx$D <- D
    tau <- rounding(v, ctx)
    Dc  <- D * (B >= tau)                       # nested pruning
    Dc  <- pmax(Dc, t(Dc)) * D                  # keep exact symmetry
    if (identical(Dc, D))  { stop_reason <- "converged"; break }
    if (n_edges(Dc) == 0)  { stop_reason <- "empty graph"; break }

    fc <- esf_fit(Ys, Zs, Dc, A, lambda, standardize, mc_min)
    ok <- length(fc$S) > 0 && isTRUE(fc$moran$p > alpha)
    trace <- rbind(trace, data.frame(iter = k + 1, edges = n_edges(Dc),
                                     n_eig = length(fc$S), moran_p = fc$moran$p,
                                     tau = tau, accepted = ok))
    if (verbose) message(sprintf("iter %d: edges %d -> %d, |S| = %d, p = %.3f%s",
                                 k + 1, n_edges(D), n_edges(Dc), length(fc$S),
                                 fc$moran$p, if (ok) "" else "  (rejected)"))
    if (length(fc$S) == 0) { stop_reason <- "no eigenvectors"; break }
    if (is.na(fc$moran$p)) { stop_reason <- "Moran undefined (saturated)"; break }
    if (fc$moran$p <= alpha) { stop_reason <- "Moran rejects"; break }
    if (keep_graphs) graphs[[length(graphs) + 1]] <- list(D = D, fhat = fit$fhat)
    D <- Dc; fit <- fc; k <- k + 1
  }

  # Step 4: hot / cold labelling ----------------------------------------------
  out$D <- D
  out$fit <- fit
  out$trace <- trace
  if (keep_graphs) out$graphs <- c(graphs, list(list(D = D, fhat = fit$fhat)))
  out$stop_reason <- stop_reason
  out$fhat   <- fit$fhat
  out$labels <- label_cells(fit$fhat, label_c)
  out$zones  <- find_zones(D, out$labels)
  out
}

#' @noRd
label_cells <- function(fhat, c = 0) {
  factor(ifelse(fhat > c, "hot", ifelse(fhat < -c, "cold", "neutral")),
         levels = c("hot", "neutral", "cold"))
}

#' Connected components (size >= 2) of the hot and of the cold subgraphs.
#' @return integer vectors `hot` and `cold` (zone id per cell, NA if none).
#' @noRd
find_zones <- function(D, labels) {
  zone_of <- function(type) {
    keep <- labels == type
    Dt <- D * outer(keep, keep)
    g  <- igraph::graph_from_adjacency_matrix(Dt, mode = "undirected")
    cm <- igraph::components(g)
    id <- cm$membership
    id[!keep | cm$csize[id] < 2] <- NA
    as.integer(factor(id))                         # renumber 1, 2, ...
  }
  list(hot = zone_of("hot"), cold = zone_of("cold"))
}

#' @export
print.hotspots <- function(x, ...) {
  e0 <- sum(x$A[upper.tri(x$A)]); e1 <- sum(x$D[upper.tri(x$D)])
  cat("hotspots: n =", x$n, " weight =", x$weight, " rounding =", x$rounding, "\n")
  cat("  step 0: z0 =", round(x$z0, 2), " p0 =", signif(x$p0, 3),
      " lambda =", signif(x$lambda, 3), "\n")
  cat("  edges:", e0, "->", e1, sprintf("(%.1f%% kept)", 100 * e1 / e0), "\n")
  cat("  stop:", x$stop_reason, if (isTRUE(x$esf_inadequate)) " [ESF with A inadequate]", "\n")
  tab <- table(x$labels)
  cat("  cells: hot", tab["hot"], " neutral", tab["neutral"], " cold", tab["cold"], "\n")
  cat("  zones: hot", max(c(0, x$zones$hot), na.rm = TRUE),
      " cold", max(c(0, x$zones$cold), na.rm = TRUE), "\n")
  invisible(x)
}
