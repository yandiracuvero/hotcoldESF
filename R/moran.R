# Moran's I for regression residuals ------------------------------------------
#
# Exact moments under the null of no spatial correlation, for residuals
# e = M_X y with M_X = I - X (X'X)^{-1} X' (Tiefelsdorf & Griffith 2007,
# eqs. 8-9). This is the same statistic used in Barde et al. (2025).

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
#' Moran's I of OLS residuals, standardised with its exact mean and variance
#' under the null hypothesis of no spatial correlation (Tiefelsdorf and
#' Griffith, 2007), which depend only on the design matrix and the weights.
#' @param e residual vector (length n) of the OLS regression on `X`.
#' @param W symmetric spatial weights matrix (n x n), e.g. binary contiguity.
#' @param X design matrix that produced the residuals, including the intercept.
#' @return List with `I`, `z` (standardised statistic) and `p` (one-sided
#'   p-value; the alternative is positive spatial correlation). `z` and `p` are
#'   `NA` when the model is (nearly) saturated and the test is undefined.
#' @references Tiefelsdorf, M. and Griffith, D. A. (2007). Semiparametric
#'   filtering of spatial autocorrelation: the eigenvector approach.
#'   *Environment and Planning A*, 39(5), 1193--1221. \doi{10.1068/a37378}
#' @examples
#' set.seed(5)
#' A <- grid_adjacency(10)
#' X <- cbind(1, rnorm(100))
#' e <- lm.fit(X, rnorm(100))$residuals
#' moran_resid(e, A, X)
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
#' Permutation test in the manner of Freedman and Lane (1983): the OLS residuals
#' are permuted, added back to the fitted values and the model is refitted.
#' Since the fitted part is removed again, this amounts to projecting the
#' permuted residuals. It does not rely on normality, which helps with sparse
#' counts.
#' @param e residual vector (length n) of the OLS regression on `X`.
#' @param W symmetric spatial weights matrix (n x n).
#' @param X design matrix that produced the residuals, including the intercept.
#' @param nperm number of permutations (default 999).
#' @return List with `I`, `z` (standardised by the permutation mean and standard
#'   deviation) and `p` (one-sided permutation p-value; the alternative is
#'   positive spatial correlation).
#' @examples
#' set.seed(6)
#' A <- grid_adjacency(10)
#' X <- cbind(1, rnorm(100))
#' e <- lm.fit(X, rnorm(100))$residuals
#' moran_perm(e, A, X, nperm = 199)$p
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
