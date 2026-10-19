# Choice of the Lasso penalty for the ESF fit -------------------------------------------
#
# The default rule, lambda = |z0|^(-a), ties the amount of smoothing to the Moran z of
# the OLS residuals: weak signal -> large lambda -> few or no eigenvectors; strong
# signal -> tiny lambda -> hundreds of eigenvectors. The alternatives below choose
# lambda from the data along the glmnet path.
#
#   "z"      lambda = |z0|^(-a)                       (Barde et al. default, a = 2)
#   "cv_min" 10-fold cross-validation, lambda with the lowest CV error
#   "cv_1se" 10-fold cross-validation, largest lambda within 1 SE of the minimum
#   "bic"    BIC along the path: n log(RSS / n) + log(n) df
#   "moran"  largest lambda whose post-Lasso residuals show no spatial
#            autocorrelation (Moran test with A, p > alpha): the classical
#            ESF stopping rule, applied along the Lasso path

#' Choose the Lasso penalty for the ESF fit
#'
#' @param y standardised response. @param Zs standardised covariates or NULL.
#' @param A contiguity matrix (used by the "moran" rule).
#' @param eig `eigen(A)`; @param mc_min candidate eigenvectors (NULL = all).
#' @param rule one of "z", "cv_min", "cv_1se", "bic", "moran".
#' @param z0 Moran z of the OLS residuals (rule "z"). @param a exponent (rule "z").
#' @param alpha level of the Moran test (rule "moran").
#' @param mix 1 = Lasso (default), 0 = ridge (glmnet's `alpha`).
#' @param foldid fold of each cell for the CV rules (NULL = 10 random folds);
#'   see [block_folds()] for spatial folds.
#' @return the penalty, in glmnet scaling.
#' @export
select_lambda <- function(y, Zs, A, eig, mc_min = NULL,
                          rule = c("z", "cv_min", "cv_1se", "bic", "moran"),
                          z0 = NULL, a = 2, alpha = 0.05, mix = 1, foldid = NULL) {
  rule <- match.arg(rule)
  if (rule == "z") return(abs(z0)^(-a))
  n <- length(y)
  cand <- if (is.null(mc_min)) seq_len(n) else
    which(eig$values >= mc_min * eig$values[1] & eig$values > 1e-8)
  E <- eig$vectors[, cand, drop = FALSE]
  k <- if (is.null(Zs)) 0 else ncol(Zs)
  X <- if (k) cbind(Zs, E) else E
  pen <- c(rep(0, k), rep(1, ncol(E)))
  if (rule %in% c("cv_min", "cv_1se")) {
    cv <- if (is.null(foldid))
      glmnet::cv.glmnet(X, y, penalty.factor = pen, nfolds = 10, alpha = mix) else
      glmnet::cv.glmnet(X, y, penalty.factor = pen, foldid = foldid, alpha = mix)
    return(if (rule == "cv_min") cv$lambda.min else cv$lambda.1se)
  }
  path <- glmnet::glmnet(X, y, penalty.factor = pen)
  if (rule == "bic") {
    rss <- (1 - path$dev.ratio) * path$nulldev
    bic <- n * log(rss / n) + log(n) * (path$df + 1)
    return(path$lambda[which.min(bic)])
  }
  # rule "moran": walk down the path, stop at the first support whose
  # post-Lasso residuals pass the Moran test
  X0 <- design0(n, Zs); last <- NULL; best <- path$lambda[1]; bestp <- -1
  for (j in seq_along(path$lambda)) {
    S <- which(path$beta[(k + 1):(k + ncol(E)), j] != 0)
    key <- paste(S, collapse = ",")
    if (identical(key, last)) next
    last <- key
    Xp <- cbind(X0, E[, S, drop = FALSE])
    if (ncol(Xp) > n - 5) break
    m <- moran_resid(stats::lm.fit(Xp, y)$residuals, A, Xp)
    if (is.na(m$p)) break
    if (m$p > bestp) { bestp <- m$p; best <- path$lambda[j] }
    if (m$p > alpha) return(path$lambda[j])
  }
  best
}

#' Spatial folds for cross-validation: square blocks of cells
#'
#' The grid is cut into square blocks of `block` x `block` cells; blocks are
#' assigned at random to `nfolds` folds, so neighbouring cells are mostly left
#' out together.
#' @param k cells per side (cells numbered row by row, as in [grid_adjacency()]).
#' @return integer vector of fold ids, length k^2.
#' @export
block_folds <- function(k, block = 4, nfolds = 10) {
  rc <- cbind((seq_len(k^2) - 1) %/% k, (seq_len(k^2) - 1) %% k)
  b  <- (rc[, 1] %/% block) * ceiling(k / block) + rc[, 2] %/% block + 1
  nb <- max(b)
  fold_of_block <- sample(rep_len(seq_len(nfolds), nb))
  fold_of_block[b]
}
