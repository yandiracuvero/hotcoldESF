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
#' @return list: `E_S`, `values_S` (eigenvalues), `gamma_S` (post-Lasso
#'   coefficients), `S` (indices), `fhat` (= E_S gamma_S), `resid`, `moran`.
#' @export
esf_fit <- function(Y, Z, D, A, lambda, standardize = TRUE) {
  n   <- length(Y)
  eig <- eigen(D, symmetric = TRUE)
  E   <- eig$vectors

  if (is.null(Z)) {
    XE  <- E
    pen <- rep(1, n)
    k   <- 0
  } else {
    Z   <- as.matrix(Z)
    k   <- ncol(Z)
    XE  <- cbind(Z, E)
    pen <- c(rep(0, k), rep(1, n))
  }

  fit <- glmnet::glmnet(XE, Y, lambda = lambda, penalty.factor = pen,
                        standardize = standardize)
  beta <- as.numeric(fit$beta[, 1])
  S    <- which(beta[(k + 1):(k + n)] != 0)

  X0 <- design0(n, Z)
  if (length(S) == 0) {
    ols <- stats::lm.fit(X0, Y)
    return(list(E_S = NULL, values_S = numeric(0), gamma_S = numeric(0),
                S = S, fhat = rep(0, n), resid = ols$residuals,
                moran = moran_resid(ols$residuals, A, X0)))
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
       moran    = moran_resid(ols$residuals, A, Xp))
}

#' Reconstruct the matrix B from a fit
#'
#' Both versions are invariant to the sign of each eigenvector.
#' @param fit output of [esf_fit()].
#' @param weight "lambda" (eigenvalues; published version) or "gamma2"
#'   (squared post-Lasso coefficients: importance of each pattern for Y).
#' @export
reconstruct_B <- function(fit, weight = c("lambda", "gamma2")) {
  weight <- match.arg(weight)
  w <- switch(weight,
              lambda = fit$values_S,
              gamma2 = fit$gamma_S^2)
  fit$E_S %*% (w * t(fit$E_S))
}
