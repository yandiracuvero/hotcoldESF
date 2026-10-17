# Rounding rules ---------------------------------------------------------------
#
# Spec: docs/algorithm.md, "Rounding rules".
# A rule is a function(v, ctx) returning a threshold tau. An edge is kept when
# b_ij >= tau. `v` are the positive b_ij on the edges still present; `ctx` holds
# what some rules need (Y, Z, D, A, lambda, weight, ...).
# Returning -Inf means "keep every edge".

#' Quantile rule (published version)
#' @param q quantile in (0, 1).
#' @export
round_quantile <- function(q = 0.75) {
  force(q)
  structure(function(v, ctx) as.numeric(stats::quantile(v, q, names = FALSE)),
            label = paste0("quantile(", q, ")"))
}

#' Two-component Gaussian mixture rule
#'
#' Fits 1 and 2 component Gaussian mixtures (unequal variances, EM) to v. If
#' BIC prefers one component, there is no separation and every edge is kept.
#' Otherwise tau is the point between the two means where the weighted
#' component densities are equal (the classification boundary). This is a
#' formal version of the "valley" rule.
#' @export
round_mixture <- function() {
  structure(function(v, ctx) {
    if (length(unique(v)) < 5) return(-Inf)
    m <- em2(v)
    if (m$bic2 <= m$bic1) return(-Inf)
    mu <- m$mu; sd <- m$sd; pr <- m$pr
    g  <- function(x) log(pr[1]) + stats::dnorm(x, mu[1], sd[1], log = TRUE) -
                     (log(pr[2]) + stats::dnorm(x, mu[2], sd[2], log = TRUE))
    if (sign(g(mu[1])) == sign(g(mu[2]))) return(mean(mu))
    stats::uniroot(g, mu)$root
  }, label = "mixture")
}

#' EM for a two-component univariate Gaussian mixture, plus BIC of 1 vs 2.
#' @noRd
em2 <- function(x, iter = 500, tol = 1e-8) {
  n  <- length(x)
  ll1 <- sum(stats::dnorm(x, mean(x), stats::sd(x), log = TRUE))
  q  <- stats::quantile(x, c(.25, .75), names = FALSE)
  mu <- q; sd <- rep(stats::sd(x) / 2, 2); pr <- c(.5, .5); ll_old <- -Inf
  for (i in seq_len(iter)) {
    d1 <- pr[1] * stats::dnorm(x, mu[1], sd[1])
    d2 <- pr[2] * stats::dnorm(x, mu[2], sd[2])
    tot <- d1 + d2 + 1e-300
    w  <- d1 / tot
    ll <- sum(log(tot))
    pr <- c(mean(w), 1 - mean(w))
    mu <- c(sum(w * x) / sum(w), sum((1 - w) * x) / sum(1 - w))
    sd <- sqrt(c(sum(w * (x - mu[1])^2) / sum(w),
                 sum((1 - w) * (x - mu[2])^2) / sum(1 - w)))
    sd <- pmax(sd, 1e-6 * stats::sd(x))
    if (abs(ll - ll_old) < tol * abs(ll)) break
    ll_old <- ll
  }
  o <- order(mu)
  list(mu = mu[o], sd = sd[o], pr = pr[o],
       bic1 = 2 * ll1 - 2 * log(n), bic2 = 2 * ll - 5 * log(n))
}

#' Permutation rule (Freedman-Lane)
#'
#' Permutes the residuals of Y on Z, refits ESF on the current graph, rebuilds
#' B, and pools the positive values on the current edges. tau is the `level`
#' quantile of that null distribution: kept edges are stronger than chance.
#' @param nperm number of permutations.
#' @param level quantile of the null distribution.
#' @export
round_permutation <- function(nperm = 50, level = 0.95) {
  force(nperm); force(level)
  structure(function(v, ctx) {
    X0   <- design0(length(ctx$Y), ctx$Z)
    ols  <- stats::lm.fit(X0, ctx$Y)
    yhat <- ctx$Y - ols$residuals
    up   <- upper.tri(ctx$D) & ctx$D == 1
    null <- unlist(lapply(seq_len(nperm), function(b) {
      Yb <- yhat + sample(ols$residuals)
      fb <- esf_fit(Yb, ctx$Z, ctx$D, ctx$A, ctx$lambda, ctx$standardize)
      if (length(fb$S) == 0) return(numeric(0))
      Bb <- reconstruct_B(fb, ctx$weight)[up]
      Bb[Bb > 0]
    }))
    if (length(null) == 0) return(-Inf)
    as.numeric(stats::quantile(null, level, names = FALSE))
  }, label = paste0("permutation(", nperm, ",", level, ")"))
}

rule_label <- function(rule) {
  lab <- attr(rule, "label")
  if (is.null(lab)) "custom" else lab
}
