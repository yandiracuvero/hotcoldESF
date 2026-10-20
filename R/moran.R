# Moran's I for regression residuals ------------------------------------------
#
# Exact moments under the null of no spatial correlation, for residuals
# e = M_X y with M_X = I - X (X'X)^{-1} X' (Tiefelsdorf & Griffith 2007,
# eqs. 8-9). This is the same statistic used in Barde et al. (2025) and in
# the reference MI-Lasso code.

#' Projection matrix that removes the column space of X.
#' Uses a rank-revealing QR, so collinear columns are allowed (e.g. when the
#' graph splits into components and some eigenvectors add up to the intercept).
#' @noRd
proj_out <- function(X) {
  q <- qr(as.matrix(X), tol = 1e-9)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  diag(nrow(Q)) - tcrossprod(Q)
}

#' Moran's I test for regression residuals
#'
#' @param e residual vector (length n).
#' @param W symmetric spatial weights matrix (n x n), e.g. binary contiguity.
#' @param X design matrix that produced the residuals (including the intercept).
#' @return list with `I`, `z` (standardised), `p` (one-sided, H1: positive
#'   spatial correlation). `z` and `p` are NA when the model is (nearly)
#'   saturated and the test is undefined.
#' @export
moran_resid <- function(e, W, X) {
  M   <- proj_out(X)
  MWM <- M %*% W %*% M
  df  <- nrow(X) - qr(as.matrix(X), tol = 1e-9)$rank
  t1  <- sum(diag(MWM))
  t2  <- sum(MWM * t(MWM))                     # trace(MWM %*% MWM), cheaper
  EI  <- t1 / df
  VI  <- 2 * (df * t2 - t1^2) / (df^2 * (df + 2))
  I   <- as.numeric(crossprod(e, W %*% e) / crossprod(e))
  if (df < 3 || !is.finite(VI) || VI <= 0)        # (nearly) saturated model
    return(list(I = I, z = NA_real_, p = NA_real_))
  z   <- (I - EI) / sqrt(VI)
  list(I = I, z = z, p = stats::pnorm(z, lower.tail = FALSE))
}

#' Intercept + covariates design matrix (works when Z is NULL).
#' @noRd
design0 <- function(n, Z) if (is.null(Z)) matrix(1, n, 1) else cbind(1, Z)

#' Moran's I test for regression residuals, permutation version
#'
#' Freedman-Lane: the OLS residuals are permuted, added back to the fitted
#' values and the model is refitted; since the fitted part is removed again,
#' this amounts to projecting the permuted residuals. Exact under
#' exchangeable errors, whatever their distribution (e.g. sparse counts).
#' @param e residuals of y on X. @param W weights. @param X design matrix.
#' @param nperm number of permutations.
#' @return list with `I`, `z` (standardised by the permutation mean and sd),
#'   `p` (one-sided, H1: positive spatial correlation).
#' @export
moran_perm <- function(e, W, X, nperm = 999) {
  q <- qr(as.matrix(X), tol = 1e-9)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  I <- as.numeric(crossprod(e, W %*% e) / crossprod(e))
  P <- vapply(seq_len(nperm), function(b) sample(e), numeric(length(e)))
  R <- P - Q %*% crossprod(Q, P)
  Ib <- colSums(R * (W %*% R)) / colSums(R^2)
  list(I = I, z = (I - mean(Ib)) / stats::sd(Ib), p = (1 + sum(Ib >= I)) / (nperm + 1))
}
