# One-pass hot/cold map ----------------------------------------------------------------
#
# The proposed procedure (no graph reconstruction):
#   1. y = standardised log(Y + 1) (or counts); Z standardised
#   2. Moran's I of the residuals of y on Z (with A); stop if not significant
#   3. MI-Lasso ESF fit on the eigenvectors of A, lambda = |z0|^(-a):
#      fhat = E_S gamma_S, the spatial component not explained by Z
#   4. permutation band: Freedman-Lane permutations of the residuals of y on Z,
#      refit (same lambda, same eigenvectors) -> fhat under no spatial structure;
#      hot  = fhat above the upper quantile of that null distribution,
#      cold = below the lower quantile, neutral otherwise
#   5. hot (cold) zones = connected components of hot (cold) cells in A

#' Hot/cold map with a permutation band
#'
#' @param Y counts per cell.
#' @param Z covariates per cell (n x k) or NULL.
#' @param A binary symmetric contiguity matrix.
#' @param response "log" (log(Y + 1), default) or "std" (counts).
#' @param a exponent of the penalty, lambda = |z0|^(-a).
#' @param mc_min candidate eigenvectors (NULL = all; see [esf_fit()]).
#' @param nperm number of permutations for the band.
#' @param level the band is [1 - level, level] quantiles (0.95 -> 5% / 95%).
#' @param band "pointwise": a quantile per cell; "max": the quantile of the
#'   maximum (minimum) of fhat over all cells, which controls the chance of
#'   any false hot (cold) cell in the map.
#' @param alpha level of the Moran test in step 2.
#' @param penalty how the Lasso penalty is chosen (see [select_lambda()]).
#'   The permutation band reuses the chosen penalty in every permutation.
#' @param fit_type "lasso" (select eigenvectors, then refit by OLS) or "ridge"
#'   (shrink all candidate eigenvectors, no selection).
#' @param foldid folds for the CV penalty rules (NULL = random; see [block_folds()]).
#' @param gate Moran test of step 2: "normal" (exact moments, normal
#'   approximation) or "permutation" (Freedman-Lane, see [moran_perm()]).
#' @param gate_nperm permutations for the permutation gate.
#' @param spline_df 0: covariates enter linearly; > 0: each covariate enters
#'   as a natural spline with this many degrees of freedom (see [covariate_basis()]).
#' @return an object of class "hotspots" (works with [rank_cells()]) with extra
#'   fields `upper`, `lower` (band per cell), `null_hot` (expected number of hot
#'   cells under the null), `band`.
#' @export
hot_cold_map <- function(Y, Z = NULL, A, response = c("log", "std"), a = 2,
                         mc_min = NULL, nperm = 999, level = 0.95,
                         band = c("pointwise", "max"), alpha = 0.05,
                         penalty = c("z", "cv_min", "cv_1se", "bic", "moran"),
                         fit_type = c("lasso", "ridge"), foldid = NULL,
                         gate = c("normal", "permutation"), gate_nperm = 999,
                         spline_df = 0) {
  gate <- match.arg(gate)
  fit_type <- match.arg(fit_type); mix <- if (fit_type == "ridge") 0 else 1
  penalty <- match.arg(penalty)
  response <- match.arg(response); band <- match.arg(band)
  A  <- unname(A); n <- length(Y)
  y  <- as.numeric(scale(if (response == "log") log1p(Y) else Y))
  Zs <- covariate_basis(Z, spline_df)
  X0  <- design0(n, Zs)
  ols <- stats::lm.fit(X0, y)
  m0  <- if (gate == "normal") moran_resid(ols$residuals, A, X0) else
    moran_perm(ols$residuals, A, X0, nperm = gate_nperm)

  out <- structure(list(n = n, A = A, D = A, Y = Y, z0 = m0$z, p0 = m0$p,
                        response = response, a = a, band = band, level = level,
                        fhat = rep(0, n), upper = rep(Inf, n), lower = rep(-Inf, n),
                        stop_reason = NA_character_), class = "hotspots")
  if (!isTRUE(m0$p <= alpha)) {
    out$stop_reason <- "no residual spatial structure (step 0)"
    out$labels <- factor(rep("neutral", n), levels = c("hot", "neutral", "cold"))
    out$zones  <- find_zones(A, out$labels)
    return(out)
  }
  eig    <- eigen(A, symmetric = TRUE)
  lambda <- select_lambda(y, Zs, A, eig, mc_min = mc_min, rule = penalty,
                          z0 = m0$z, a = a, alpha = alpha, mix = mix, foldid = foldid)
  fit    <- esf_fit(y, Zs, A, A, lambda, mc_min = mc_min, eig = eig, alpha = mix)

  # permutation band: what fhat looks like when the residuals carry no structure
  yhat <- y - ols$residuals
  if (fit_type == "ridge") {
    # ridge is a linear smoother, fhat = S y with S X = 0, so each permutation is
    # fhat^(b) = S Pi_b e: build S once and permute the residuals (no refits)
    Sm   <- ridge_smoother(y, Zs, eig, lambda, mc_min = mc_min)
    fit$fhat <- as.numeric(Sm %*% y)
    null <- Sm %*% vapply(seq_len(nperm), function(b) sample(ols$residuals), numeric(n))
  } else {
    null <- vapply(seq_len(nperm), function(b) {
      yb <- yhat + sample(ols$residuals)
      esf_fit(yb, Zs, A, A, lambda, mc_min = mc_min, eig = eig, test = FALSE, alpha = mix)$fhat
    }, numeric(n))                                      # n x nperm
  }
  up_pw <- apply(null, 1, stats::quantile, probs = level, names = FALSE)
  lo_pw <- apply(null, 1, stats::quantile, probs = 1 - level, names = FALSE)
  up_mx <- rep(stats::quantile(apply(null, 2, max), level, names = FALSE), n)
  lo_mx <- rep(stats::quantile(apply(null, 2, min), 1 - level, names = FALSE), n)
  upper <- if (band == "pointwise") up_pw else up_mx
  lower <- if (band == "pointwise") lo_pw else lo_mx
  f <- fit$fhat
  labels <- factor(ifelse(f > upper, "hot", ifelse(f < lower, "cold", "neutral")),
                   levels = c("hot", "neutral", "cold"))
  out$lambda <- lambda; out$fit <- fit; out$fhat <- f
  out$upper <- upper; out$lower <- lower
  out$null_hot <- mean(colSums(null > upper))
  lab <- function(u, l) factor(ifelse(f > u, "hot", ifelse(f < l, "cold", "neutral")),
                               levels = c("hot", "neutral", "cold"))
  out$labels_pointwise <- lab(up_pw, lo_pw)
  out$labels_max       <- lab(up_mx, lo_mx)
  out$labels <- labels
  out$zones  <- find_zones(A, labels)
  out$stop_reason <- "ok"
  out$esf_moran_p <- fit$moran$p
  out
}

# Benchmarks: Getis-Ord Gi* and local Moran (LISA) -------------------------------------

#' Getis-Ord Gi* z-scores (binary weights including the cell itself)
#' @param x values per cell (counts, rates or residuals).
#' @param A binary contiguity matrix.
#' @return z-scores; large positive = high values surrounded by high values.
#' @export
local_gi_star <- function(x, A) {
  n  <- length(x)
  W  <- A + diag(n)                       # Gi*: the cell is its own neighbour
  wi <- rowSums(W)
  xb <- mean(x); s <- sqrt(sum(x^2) / n - xb^2)
  num <- as.numeric(W %*% x) - xb * wi
  den <- s * sqrt((n * rowSums(W^2) - wi^2) / (n - 1))
  num / den
}

#' Local Moran's I (LISA) with conditional permutation p-values
#' @param x values per cell. @param A binary contiguity matrix.
#' @param nperm permutations (neighbour values drawn with replacement from the
#'   other cells). @return data.frame: `I`, `p` (one-sided, positive
#'   association), `quad` ("HH", "LL", "HL", "LH").
#' @export
local_moran <- function(x, A, nperm = 999) {
  n  <- length(x); z <- x - mean(x); m2 <- sum(z^2) / n
  Wr <- A / pmax(rowSums(A), 1)                 # row-standardised
  lag <- as.numeric(Wr %*% z)
  I  <- z / m2 * lag
  k  <- rowSums(A)
  p  <- vapply(seq_len(n), function(i) {
    if (k[i] == 0) return(NA_real_)
    zo <- z[-i]
    lagp <- rowMeans(matrix(zo[sample.int(n - 1, nperm * k[i], replace = TRUE)], nperm, k[i]))
    Ip <- z[i] / m2 * lagp
    (sum(Ip >= I[i]) + 1) / (nperm + 1)
  }, numeric(1))
  quad <- ifelse(z > 0, ifelse(lag > 0, "HH", "HL"), ifelse(lag > 0, "LH", "LL"))
  data.frame(I = I, p = p, quad = quad)
}

#' Covariate design used by the procedure
#'
#' @param Z covariates (n x k) or NULL.
#' @param spline_df 0: standardised covariates (linear terms); > 0: a natural
#'   spline basis with `spline_df` degrees of freedom for each covariate,
#'   standardised. Columns that are constant or not finite are dropped.
#' @return matrix or NULL.
#' @export
covariate_basis <- function(Z, spline_df = 0) {
  if (is.null(Z)) return(NULL)
  Z <- as.matrix(Z)
  B <- if (spline_df > 0)
    do.call(cbind, lapply(seq_len(ncol(Z)), function(j) splines::ns(Z[, j], df = spline_df))) else Z
  B <- scale(B)
  B <- B[, apply(B, 2, function(z) all(is.finite(z))), drop = FALSE]
  if (ncol(B) == 0) NULL else B
}
