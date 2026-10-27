# Covariate-adjusted hot/cold map (the procedure of the paper) ------------------------

#' Covariate-adjusted hot and cold spots
#'
#' Maps hot and cold spots of area-level counts after adjusting for
#' covariates. A hot (cold) spot is a cell whose spatially structured residual
#' component \eqn{f_i} is positive (negative) and larger in absolute value than
#' what spatially unstructured noise would produce. The procedure has four steps:
#' \enumerate{
#'   \item regress \eqn{y}, the standardised \eqn{\log(1 + Y)}, on the covariates and
#'     test the residuals for positive spatial autocorrelation with Moran's I and
#'     its exact moments (Tiefelsdorf and Griffith, 2007). If the test does not
#'     reject, every cell is labelled neutral and the procedure stops;
#'   \item estimate \eqn{f} by eigenvector spatial filtering: the eigenvectors of
#'     `A` whose eigenvalue is at least `mc_min` times the largest enter a ridge
#'     regression together with the (unpenalised) covariates, with the penalty
#'     chosen by ten-fold cross-validation;
#'   \item compare each \eqn{\hat f_i} with its distribution when the OLS residuals
#'     are permuted among the cells (Freedman and Lane, 1983): a cell is hot if
#'     \eqn{\hat f_i} is above the `level` quantile, cold if it is below the
#'     `1 - level` quantile, and neutral otherwise. For ridge the fit is linear in
#'     \eqn{y}, so the permuted fits are computed without refitting (see
#'     [ridge_smoother()]);
#'   \item hot (cold) zones are the connected components of the hot (cold) cells.
#' }
#' The defaults reproduce the procedure of Cuvero et al. Use [rank_cells()] to
#' turn the result into a priority list.
#'
#' @param Y numeric vector of counts per cell (length n).
#' @param Z covariates per cell: a numeric matrix or data frame with n rows, or
#'   `NULL` for no covariates. Covariates are standardised internally.
#' @param A binary, symmetric n x n contiguity matrix with zero diagonal (for
#'   example from [grid_adjacency()] or from `spdep::nb2mat(style = "B")`).
#' @param response `"log"` (default) models the standardised \eqn{\log(1 + Y)};
#'   `"std"` models the standardised counts.
#' @param a exponent of the Moran's I Lasso penalty \eqn{|z_0|^{-a}}, used only
#'   when `penalty = "z"`.
#' @param mc_min candidate eigenvectors: those with eigenvalue at least `mc_min`
#'   times the largest (Griffith's rule; default 0.25). `NULL` uses all n
#'   eigenvectors.
#' @param nperm number of permutations for the band (default 999).
#' @param level the band is given by the `1 - level` and `level` quantiles of the
#'   permuted fits (default 0.95, i.e. the 5\% and 95\% quantiles).
#' @param band `"pointwise"` (default): one band per cell; `"max"`: the quantile of
#'   the maximum (minimum) of \eqn{\hat f} over the map, which controls the
#'   chance of any false hot (cold) cell but is very conservative.
#' @param alpha level of the Moran test of step 1 (default 0.05).
#' @param penalty how the penalty is chosen: `"cv_min"` (default; ten-fold
#'   cross-validation, smallest error), `"cv_1se"` (largest penalty within one
#'   standard error), `"bic"`, `"moran"` (largest penalty whose residuals pass
#'   the Moran test), or `"z"` (\eqn{|z_0|^{-a}}, Barde et al., 2025). See
#'   [select_lambda()]. The permutation band reuses the chosen penalty.
#' @param fit_type `"ridge"` (default; shrinks all candidate eigenvectors) or
#'   `"lasso"` (selects eigenvectors and refits them by OLS).
#' @param foldid optional fold of each cell for the cross-validation (default:
#'   ten random folds); see [block_folds()] for spatial folds.
#' @param gate Moran test of step 1: `"normal"` (default; exact moments and a
#'   normal approximation, [moran_resid()]) or `"permutation"` ([moran_perm()]).
#' @param gate_nperm permutations for the permutation version of the test.
#' @param spline_df 0 (default): covariates enter linearly; a positive value
#'   enters each covariate as a natural spline with that many degrees of freedom
#'   (see [covariate_basis()]).
#' @param keep_null if `TRUE`, the fits to the permuted residuals are returned
#'   in element `null` (an n x `nperm` matrix), for example to plot the band.
#' @return An object of class `"hotcold"`, a list with elements
#'   \describe{
#'     \item{`labels`}{factor with levels `"hot"`, `"neutral"`, `"cold"`, one value per cell;}
#'     \item{`fhat`}{estimated spatial component \eqn{\hat f} (centred), in standard
#'       deviations of the modelled response;}
#'     \item{`lower`, `upper`}{the band for each cell;}
#'     \item{`zones`}{list with integer vectors `hot` and `cold`: zone number of each
#'       cell (`NA` outside a zone of at least two cells);}
#'     \item{`z0`, `p0`}{standardised Moran's I of the OLS residuals and its p-value;}
#'     \item{`lambda`}{the penalty (glmnet scale);}
#'     \item{`stop_reason`}{`"ok"`, or the reason why no band was computed;}
#'     \item{`null`}{the permuted fits, only if `keep_null = TRUE`;}
#'   }
#'   and further elements used by [rank_cells()] and [print.hotcold()].
#' @references
#' Barde, S., Cherodian, R. and Tchuente, G. (2025). Moran's I Lasso for models
#' with spatially correlated data. *The Econometrics Journal*.
#' \doi{10.1093/ectj/utaf008}
#'
#' Freedman, D. and Lane, D. (1983). A nonstochastic interpretation of reported
#' significance levels. *Journal of Business & Economic Statistics*, 1(4), 292--298.
#'
#' Griffith, D. A. (2003). *Spatial Autocorrelation and Spatial Filtering*.
#' Springer, Berlin.
#'
#' Tiefelsdorf, M. and Griffith, D. A. (2007). Semiparametric filtering of
#' spatial autocorrelation: the eigenvector approach. *Environment and
#' Planning A*, 39(5), 1193--1221. \doi{10.1068/a37378}
#' @seealso [rank_cells()] for the priority list; [local_gi_star()] and
#'   [local_moran()] for the usual local statistics.
#' @examples
#' set.seed(1)
#' s <- sim_thomas(k = 15, n_parents = 2, radius = 120, sons = 300,
#'                 n_background = 1500, n_cov = 1500)
#' h <- hot_cold_map(s$Y, s$Z, s$A, nperm = 199)
#' h
#' table(h$labels)
#' head(rank_cells(h), 10)          # the ten cells with the highest priority
#' @export
hot_cold_map <- function(Y, Z = NULL, A, response = c("log", "std"), a = 2,
                         mc_min = 0.25, nperm = 999, level = 0.95,
                         band = c("pointwise", "max"), alpha = 0.05,
                         penalty = c("cv_min", "cv_1se", "bic", "moran", "z"),
                         fit_type = c("ridge", "lasso"), foldid = NULL,
                         gate = c("normal", "permutation"), gate_nperm = 999,
                         spline_df = 0, keep_null = FALSE) {
  check_inputs(Y, Z, A)
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
                        stop_reason = NA_character_), class = "hotcold")
  if (!isTRUE(m0$p <= alpha)) {
    out$stop_reason <- "no residual spatial structure (step 1)"
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
  if (keep_null) out$null <- null
  out
}

#' @noRd
check_inputs <- function(Y, Z, A) {
  if (!is.numeric(Y) || anyNA(Y)) stop("`Y` must be a numeric vector without missing values.")
  n <- length(Y)
  A <- as.matrix(A)
  if (!all(dim(A) == n)) stop("`A` must be an n x n matrix, with n = length(Y).")
  if (!isSymmetric(unname(A)) || !all(A %in% c(0, 1)) || any(diag(A) != 0))
    stop("`A` must be a binary symmetric matrix with zero diagonal.")
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z) || nrow(Z) != n || anyNA(Z))
      stop("`Z` must be numeric, with one row per cell and no missing values.")
  }
  invisible(TRUE)
}

#' Print a hot/cold map
#'
#' @param x an object of class `"hotcold"`, from [hot_cold_map()].
#' @param ... not used.
#' @return `x`, invisibly.
#' @export
print.hotcold <- function(x, ...) {
  tab <- table(x$labels)
  cat("Covariate-adjusted hot and cold spots:", x$n, "cells\n")
  cat(sprintf("  Moran test of the residuals: z0 = %.2f, p = %s\n", x$z0, format.pval(x$p0, digits = 3)))
  if (!identical(x$stop_reason, "ok")) {
    cat("  ", x$stop_reason, ": all cells neutral\n", sep = "")
    return(invisible(x))
  }
  cat(sprintf("  band: %s, %g%% quantiles\n", x$band, 100 * x$level))
  cat("  cells: ", tab["hot"], " hot, ", tab["neutral"], " neutral, ", tab["cold"], " cold\n", sep = "")
  cat("  zones: ", max(c(0, x$zones$hot), na.rm = TRUE), " hot, ",
      max(c(0, x$zones$cold), na.rm = TRUE), " cold\n", sep = "")
  invisible(x)
}

# Benchmarks: Getis-Ord Gi* and local Moran (LISA) -------------------------------------

#' Getis-Ord Gi* statistic
#'
#' Standardised \eqn{G_i^*} of Getis and Ord (1992) and Ord and Getis (1995)
#' with binary weights that include the cell itself.
#' @param x numeric values per cell (counts, rates or regression residuals).
#' @param A binary symmetric contiguity matrix with zero diagonal.
#' @return Numeric vector of z-scores: large positive values mark high values
#'   surrounded by high values.
#' @references
#' Getis, A. and Ord, J. K. (1992). The analysis of spatial association by use of
#' distance statistics. *Geographical Analysis*, 24(3), 189--206.
#'
#' Ord, J. K. and Getis, A. (1995). Local spatial autocorrelation statistics:
#' distributional issues and an application. *Geographical Analysis*, 27(4), 286--306.
#' @examples
#' A <- grid_adjacency(10)
#' x <- rnorm(100); x[c(44:46, 54:56, 64:66)] <- x[c(44:46, 54:56, 64:66)] + 3
#' which(local_gi_star(x, A) > qnorm(0.95))
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

#' Local Moran's I with conditional permutation p-values
#'
#' Local indicator of spatial association (Anselin, 1995) with row-standardised
#' weights. The p-values are one-sided (positive association) and use
#' conditional permutations: the neighbours of each cell are replaced by values
#' drawn with replacement from the other cells.
#' @param x numeric values per cell.
#' @param A binary symmetric contiguity matrix with zero diagonal.
#' @param nperm number of permutations (default 999).
#' @return Data frame with columns `I` (local Moran's I), `p` and `quad`, the
#'   quadrant of the Moran scatterplot (`"HH"`, `"LL"`, `"HL"`, `"LH"`).
#' @references Anselin, L. (1995). Local indicators of spatial association--LISA.
#'   *Geographical Analysis*, 27(2), 93--115.
#' @examples
#' set.seed(2)
#' A <- grid_adjacency(10)
#' x <- rnorm(100); x[c(44:46, 54:56, 64:66)] <- x[c(44:46, 54:56, 64:66)] + 3
#' lisa <- local_moran(x, A, nperm = 199)
#' which(lisa$quad == "HH" & lisa$p < 0.05)
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
#' Standardises the covariates or, if `spline_df > 0`, replaces each covariate by
#' a standardised natural spline basis.
#' @param Z covariates (n x k matrix or data frame) or `NULL`.
#' @param spline_df 0: standardised covariates (linear terms); > 0: a natural
#'   spline basis with `spline_df` degrees of freedom for each covariate,
#'   standardised. Columns that are constant or not finite are dropped.
#' @return Numeric matrix, or `NULL` if `Z` is `NULL`.
#' @examples
#' Z <- cbind(a = 1:50, b = rnorm(50))
#' dim(covariate_basis(Z))                  # 50 x 2
#' dim(covariate_basis(Z, spline_df = 3))   # 50 x 6
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
