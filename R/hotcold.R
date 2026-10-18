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
#' @return an object of class "hotspots" (works with [rank_cells()]) with extra
#'   fields `upper`, `lower` (band per cell), `null_hot` (expected number of hot
#'   cells under the null), `band`.
#' @export
hot_cold_map <- function(Y, Z = NULL, A, response = c("log", "std"), a = 2,
                         mc_min = NULL, nperm = 200, level = 0.95,
                         band = c("pointwise", "max"), alpha = 0.05) {
  response <- match.arg(response); band <- match.arg(band)
  A  <- unname(A); n <- length(Y)
  y  <- as.numeric(scale(if (response == "log") log1p(Y) else Y))
  Zs <- if (is.null(Z)) NULL else scale(as.matrix(Z))
  if (!is.null(Zs)) Zs <- Zs[, apply(Zs, 2, function(z) all(is.finite(z))), drop = FALSE]
  if (!is.null(Zs) && ncol(Zs) == 0) Zs <- NULL
  X0  <- design0(n, Zs)
  ols <- stats::lm.fit(X0, y)
  m0  <- moran_resid(ols$residuals, A, X0)

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
  lambda <- abs(m0$z)^(-a)
  eig    <- eigen(A, symmetric = TRUE)
  fit    <- esf_fit(y, Zs, A, A, lambda, mc_min = mc_min, eig = eig)

  # permutation band: what fhat looks like when the residuals carry no structure
  yhat <- y - ols$residuals
  null <- vapply(seq_len(nperm), function(b) {
    yb <- yhat + sample(ols$residuals)
    esf_fit(yb, Zs, A, A, lambda, mc_min = mc_min, eig = eig, test = FALSE)$fhat
  }, numeric(n))                                        # n x nperm
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
