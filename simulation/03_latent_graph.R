# Latent-graph recovery test ------------------------------------------------------
#
# Question: does the reconstruction recover the latent graph W* that generated
# the data? (The paper's original claim: "pests do not use all connections".)
#
# DGP: queen grid 20 x 20 (A). Y follows a SAR model on W* (a subgraph of A):
#   barriers  river with one bridge + highway on the left half (~85 of 1482 edges cut)
#   random    20% of edges missing at random (hard; negative control)
#   none      W* = A (control: nothing should be cut)
# rho in {0.5, 0.8, 0.95}; Z = (smooth gradient, noise), beta = (1, 0).
#
# Edge-level evaluation (truth: edge in W* = keep, barrier = cut):
#   AUC_*   threshold-free: how well a score separates W* edges from barrier edges
#           (0.5 = no information, 1 = perfect)
#     B_lambda / B_gamma2  b_ij from the FIRST ESF fit (the reconstruction signal)
#     dres                 -|e_i - e_j|, e = OLS residuals of Y on Z (trivial baseline)
#   After the full algorithm: cut_rate (barrier edges removed), keep_rate (W* edges
#   kept), J = cut_rate + keep_rate - 1 (Youden; 0 = chance).
#   J_dres_same: J of the trivial baseline cutting the SAME number of edges.
#
# Run from repo root: Rscript simulation/03_latent_graph.R

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else library(hotspots)
  library(parallel)
})
NREP  <- as.integer(Sys.getenv("NREP", 20))
CORES <- as.integer(Sys.getenv("CORES", detectCores()))

auc <- function(score, keep) {           # P(score of a kept edge > score of a cut edge)
  r <- rank(score); n1 <- sum(keep); n0 <- sum(!keep)
  if (n0 == 0 || n1 == 0) return(NA)
  (sum(r[keep]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

run_one <- function(type, rho, rep) {
  set.seed(round(1e4 * rho) + 1e3 * match(type, c("barriers", "random", "none")) + rep)
  s  <- sim_sar_barrier(k = 20, rho = rho, type = type)
  up <- which(upper.tri(s$A) & s$A == 1)
  keep_true <- s$Wstar[up] == 1
  Ys <- as.numeric(scale(s$Y)); Zs <- scale(s$Z)
  X0 <- cbind(1, Zs)
  e  <- lm.fit(X0, Ys)$residuals
  ij <- arrayInd(up, dim(s$A))
  dres <- -abs(e[ij[, 1]] - e[ij[, 2]])

  # first ESF fit (same lambda as the algorithm, a = 2, all eigenvectors)
  z0  <- moran_resid(e, s$A, X0)$z
  f0  <- esf_fit(Ys, Zs, s$A, s$A, abs(z0)^(-2), mc_min = NULL)
  Bl  <- if (length(f0$S)) reconstruct_B(f0, "lambda")[up] else rep(0, length(up))
  Bg  <- if (length(f0$S)) reconstruct_B(f0, "gamma2")[up] else rep(0, length(up))

  base <- data.frame(type = type, rho = rho, rep = rep,
                     n_barrier = sum(!keep_true), n_eig0 = length(f0$S),
                     AUC_B_lambda = auc(Bl, keep_true), AUC_B_gamma2 = auc(Bg, keep_true),
                     AUC_dres = auc(dres, keep_true))

  variants <- list(
    "all/lambda/q75"  = list(mc = NULL, w = "lambda", r = round_quantile(0.75)),
    "all/lambda/mix"  = list(mc = NULL, w = "lambda", r = round_mixture()),
    "all/gamma2/mix"  = list(mc = NULL, w = "gamma2", r = round_mixture()),
    "g25/lambda/q75"  = list(mc = 0.25, w = "lambda", r = round_quantile(0.75)),
    "g25/lambda/mix"  = list(mc = 0.25, w = "lambda", r = round_mixture()))
  rows <- lapply(names(variants), function(vn) {
    v <- variants[[vn]]
    res <- find_hotspots(s$Y, s$Z, s$A, weight = v$w, rounding = v$r, mc_min = v$mc)
    kept <- res$D[up] == 1
    cut_rate  <- if (any(!keep_true)) mean(!kept[!keep_true]) else NA
    keep_rate <- mean(kept[keep_true])
    n_cut <- sum(!kept)
    # trivial baseline cutting the same number of edges (largest |e_i - e_j|)
    kept_d <- rep(TRUE, length(up)); if (n_cut > 0) kept_d[order(dres)[seq_len(n_cut)]] <- FALSE
    cut_d  <- if (any(!keep_true)) mean(!kept_d[!keep_true]) else NA
    keep_d <- mean(kept_d[keep_true])
    cbind(base, variant = vn, stop = res$stop_reason, edges_kept = mean(kept),
          cut_rate = cut_rate, keep_rate = keep_rate, J = cut_rate + keep_rate - 1,
          J_dres_same = cut_d + keep_d - 1)
  })
  do.call(rbind, rows)
}

jobs <- expand.grid(type = c("barriers", "random", "none"), rho = c(0.5, 0.8, 0.95),
                    rep = seq_len(NREP), stringsAsFactors = FALSE)
t0 <- Sys.time()
out <- mclapply(seq_len(nrow(jobs)), function(i)
  tryCatch(run_one(jobs$type[i], jobs$rho[i], jobs$rep[i]),
           error = function(e) { message("job ", i, ": ", conditionMessage(e)); NULL }),
  mc.cores = CORES, mc.preschedule = FALSE)
res <- do.call(rbind, Filter(is.data.frame, out))
dir.create("simulation/results", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "simulation/results/03_latent_graph.csv", row.names = FALSE)
message("done ", sum(vapply(out, is.data.frame, TRUE)), "/", nrow(jobs), " in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
