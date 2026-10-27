# ESF(D): MI-Lasso eigenvector selection for a given graph D -------------------
#
# Spec: docs/algorithm.md, Step 1.

#' Fit the ESF model for one graph
#'
#' @param Y standardised response (length n).
#' @param Z standardised covariates (n x k matrix) or NULL.
#' @param D candidate graph (binary symmetric n x n). Its eigenvectors are used.
#' @param A original graph. Moran's I of the residuals is ALWAYS computed with A.
#' @param lambda fixed Lasso penalty, |z0|^(-2) (glmnet scaling).
#' @param standardize passed to glmnet (TRUE matches the reference code).
#' @param mc_min candidate set: only eigenvectors with eigenvalue >=
#'   mc_min * (largest eigenvalue) enter the Lasso (Griffith's rule, 0.25 is
#'   the usual choice). NULL uses all n eigenvectors (published version).
#' @param test compute Moran's I of the residuals (skipped inside permutations).
#' @param alpha 1 (default): Lasso, then post-Lasso OLS on the selected
#'   eigenvectors. 0: ridge on all candidate eigenvectors, no selection and no
#'   refit; fhat uses the shrunk coefficients.
#' @param eig optional precomputed `eigen(D, symmetric = TRUE)` (reused by the
#'   permutation rule, where D is the same for every permutation).
#' @return list: `E_S`, `values_S` (eigenvalues), `gamma_S` (post-Lasso
#'   coefficients), `S` (indices), `fhat` (= E_S gamma_S), `resid`, `moran`.
#' @noRd
esf_fit <- function(Y, Z, D, A, lambda, standardize = TRUE, mc_min = 0.25,
                    eig = NULL, test = TRUE, alpha = 1) {
  n   <- length(Y)
  if (is.null(eig)) eig <- eigen(D, symmetric = TRUE)
  cand <- if (is.null(mc_min)) seq_len(n) else
    which(eig$values >= mc_min * eig$values[1] & eig$values > 1e-8)
  eig$vectors <- eig$vectors[, cand, drop = FALSE]
  eig$values  <- eig$values[cand]
  E   <- eig$vectors
  nE  <- ncol(E)

  if (is.null(Z)) {
    XE  <- E
    pen <- rep(1, nE)
    k   <- 0
  } else {
    Z   <- as.matrix(Z)
    k   <- ncol(Z)
    XE  <- cbind(Z, E)
    pen <- c(rep(0, k), rep(1, nE))
  }

  if (alpha == 0 && nE > 0) {                       # ridge: shrink, do not select
    if (ncol(XE) < 2) { XE <- cbind(XE, 0); pen <- c(pen, 1) }
    fit <- glmnet::glmnet(XE, Y, alpha = 0, lambda = lambda, penalty.factor = pen,
                          standardize = standardize)
    g   <- as.numeric(fit$beta[(k + 1):(k + nE), 1])
    res <- Y - as.numeric(stats::predict(fit, XE, s = lambda))
    return(list(E_S = E, values_S = eig$values, gamma_S = g, S = seq_len(nE),
                fhat = as.numeric(scale(E %*% g, scale = FALSE)), resid = res,
                moran = NULL))
  }
  if (nE == 0) {
    S <- integer(0)
  } else {
    if (ncol(XE) < 2) {                  # glmnet needs >= 2 columns
      XE  <- cbind(XE, 0)
      pen <- c(pen, 1)
    }
    fit  <- glmnet::glmnet(XE, Y, lambda = lambda, penalty.factor = pen,
                           standardize = standardize)
    beta <- as.numeric(fit$beta[, 1])
    S    <- which(beta[(k + 1):(k + nE)] != 0)
  }

  X0 <- design0(n, Z)
  if (length(S) == 0) {
    ols <- stats::lm.fit(X0, Y)
    return(list(E_S = NULL, values_S = numeric(0), gamma_S = numeric(0),
                S = S, fhat = rep(0, n), resid = ols$residuals,
                moran = if (test) moran_resid(ols$residuals, A, X0) else NULL))
  }

  E_S  <- E[, S, drop = FALSE]
  Xp   <- cbind(X0, E_S)
  ols  <- stats::lm.fit(Xp, Y)                       # post-Lasso OLS
  coef <- ols$coefficients
  coef[is.na(coef)] <- 0                             # guard against collinearity
  gamma_S <- coef[(ncol(X0) + 1):ncol(Xp)]

  list(E_S      = E_S,
       values_S = eig$values[S],
       gamma_S  = unname(gamma_S),
       S        = S,
       # centred: the level is absorbed by the intercept (the leading
       # eigenvector of a connected graph is nearly constant and collinear
       # with it), so only deviations from the mean are meaningful.
       fhat     = as.numeric(scale(E_S %*% gamma_S, scale = FALSE)),
       resid    = ols$residuals,
       moran    = if (test) moran_resid(ols$residuals, A, Xp) else NULL)
}

#' Reconstruct the matrix B from a fit
#'
#' Both versions are invariant to the sign of each eigenvector.
#' @param fit output of [esf_fit()].
#' @param weight "lambda" (eigenvalues; published version) or "gamma2"
#'   (squared post-Lasso coefficients: importance of each pattern for Y).
#' @noRd
reconstruct_B <- function(fit, weight = c("lambda", "gamma2")) {
  weight <- match.arg(weight)
  w <- switch(weight,
              lambda = fit$values_S,
              gamma2 = fit$gamma_S^2)
  fit$E_S %*% (w * t(fit$E_S))
}

#' Ridge eigenvector spatial filter as a linear smoother
#'
#' With a fixed set of candidate eigenvectors and a fixed penalty, the ridge
#' estimate of the spatial component is linear in the response,
#' \eqn{\hat f = S y}. Because \eqn{S X = 0} for the covariate design \eqn{X},
#' the fit to permuted residuals is \eqn{\hat f^{(b)} = S \Pi_b u}, so the
#' permutation band of [hot_cold_map()] needs one matrix product instead of one
#' refit per permutation.
#'
#' `S` reproduces the ridge fit of [glmnet::glmnet()] with `alpha = 0`: an
#' unpenalised intercept and covariates, standardised columns, penalty factors
#' rescaled to sum to the number of columns, and the penalty divided by the
#' standard deviation of the response (as `glmnet` does internally). The scale
#' of the response is taken from `y` and then held fixed.
#' @param y standardised response (sets the scale of the penalty).
#' @param Z standardised covariates (n x k matrix) or `NULL`.
#' @param eig eigendecomposition of the contiguity matrix,
#'   `eigen(A, symmetric = TRUE)`.
#' @param lambda penalty on the `glmnet` scale, e.g. from [select_lambda()].
#' @param mc_min candidate eigenvectors: eigenvalue at least `mc_min` times the
#'   largest; `NULL` uses all of them.
#' @return The n x n smoothing matrix `S`.
#' @examples
#' set.seed(3)
#' A <- grid_adjacency(8)
#' y <- as.numeric(scale(rnorm(64)))
#' S <- ridge_smoother(y, NULL, eigen(A, symmetric = TRUE), lambda = 0.5)
#' fhat <- S %*% y
#' max(abs(S %*% rep(1, 64)))   # constants are removed: S 1 = 0
#' @export
ridge_smoother <- function(y, Z, eig, lambda, mc_min = 0.25) {
  n <- length(y)
  cand <- if (is.null(mc_min)) seq_len(n) else
    which(eig$values >= mc_min * eig$values[1] & eig$values > 1e-8)
  E  <- eig$vectors[, cand, drop = FALSE]
  k  <- if (is.null(Z)) 0 else ncol(Z)
  X  <- if (k) cbind(as.matrix(Z), E) else E
  p  <- ncol(X)
  pf <- c(rep(0, k), rep(1, ncol(E))) * p / ncol(E)   # glmnet rescales penalty.factor
  Xc <- sweep(X, 2, colMeans(X))
  sj2 <- colMeans(Xc^2)                                # glmnet standardize = TRUE
  sy  <- sqrt(mean((y - mean(y))^2))                   # glmnet scales lambda by sd(y)
  lam <- lambda / sy
  G  <- solve(crossprod(Xc) + n * lam * diag(pf * sj2, p), t(Xc))  # coefficients = G (y - ybar)
  Ec <- sweep(E, 2, colMeans(E))
  S  <- Ec %*% G[(k + 1):p, , drop = FALSE]
  sweep(S, 2, colMeans(S))                            # right-multiply by the centring matrix
}
