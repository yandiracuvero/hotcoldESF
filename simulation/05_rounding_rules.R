# Rounding rules compared ------------------------------------------------------------
#
# Question: does ANY rounding rule make the iterative reconstruction useful?
#
# Everything is fixed as in the 2025 manuscript (Algorithm 1) except the rounding rule:
#   - response: standardised counts; candidate eigenvectors: all; lambda = |z0|^(-2)
#   - B = E_S Delta_S E_S' (eigenvalue weights); nested pruning; Moran stop with A
# Rounding rules (how B is turned into a 0/1 graph at each iteration):
#   none   keep every edge (no reconstruction; control)
#   q50, q75, q80, q90   keep edges with b_ij >= the q-quantile of the positive b_ij
#                        on the current edges (q75 = the manuscript's choice)
#   mixture              2-component Gaussian mixture on those values; cut at the
#                        boundary between components (formal version of the "valley")
#   permutation          Freedman-Lane permutations of the residuals of Y on Z
#                        (30 per iteration); keep edges above the 95th percentile of
#                        b_ij under "no spatial structure"
#
# Scenarios (20 x 20 grid, 400 cells, cell side 50, cluster radius 100):
#   A   3 clusters (500 points each) + 2000 random points; no covariate
#   B   A + 3000 points in an area driven by a covariate (the model receives it)
#   D   SAR process on a grid with a river and a highway (85 of 1482 edges cut),
#       rho = 0.8 and 0.95
#   D0  same SAR process on the full grid (no barriers): nothing should be cut
# 20 replicates per scenario.
#
# Metrics
#   A, B: precision of the top 20 cells (share inside a true cluster) for the four
#         ranking rules; cold cells among the top 20; edges kept; iterations
#   D:    barrier edges cut, real edges kept, Youden J = cut + kept - 1 (0 = chance)
#   D0:   edges kept (should be 100%)
#
# Run from repo root: Rscript simulation/05_rounding_rules.R   (about 30 min on 2 cores)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else
    if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(hotcoldESF)
  library(parallel)
})
NREP  <- as.integer(Sys.getenv("NREP", 20))
CORES <- as.integer(Sys.getenv("CORES", detectCores()))
OUT   <- "simulation/results/05"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rules <- list(none = round_none(), q50 = round_quantile(0.50), q75 = round_quantile(0.75),
              q80 = round_quantile(0.80), q90 = round_quantile(0.90),
              mixture = round_mixture(), permutation = round_permutation(30, 0.95))
rank_names <- c(R0 = "Centrality, whole graph", R1 = "Centrality, hot cells",
                R2 = "Largest excess", R3 = "Hot zones first")

run_job <- function(scn, rep) {
  set.seed(5000 + 100 * match(scn, c("A", "B", "D80", "D95", "D0")) + rep)
  out <- list()
  if (scn %in% c("A", "B")) {
    s <- sim_thomas(k = 20, n_parents = 3, radius = 100, sons = 500, n_background = 2000,
                    n_cov = if (scn == "B") 3000 else 0)
    Z <- if (scn == "B") s$Z else NULL
    for (rn in names(rules)) {
      t0 <- Sys.time()
      r <- find_hotspots(s$Y, Z, s$A, mc_min = NULL, rounding = rules[[rn]])
      for (rk in names(rank_names)) {
        top <- rank_cells(r, rk)[1:20]
        out[[length(out) + 1]] <- data.frame(
          scenario = scn, rep = rep, rule = rn, ranking = rank_names[[rk]],
          precision = mean(top %in% s$footprint), cold = mean(r$labels[top] == "cold"),
          edges_kept = sum(r$D) / sum(r$A), iters = sum(r$trace$accepted) - 1,
          stop = r$stop_reason, J = NA, cut = NA, kept = NA,
          secs = as.numeric(Sys.time() - t0, units = "secs"))
      }
    }
  } else {
    rho  <- switch(scn, D80 = 0.8, D95 = 0.95, D0 = 0.95)
    s    <- sim_sar_barrier(k = 20, rho = rho, type = if (scn == "D0") "none" else "barriers")
    up   <- which(upper.tri(s$A) & s$A == 1)
    real <- s$Wstar[up] == 1
    for (rn in names(rules)) {
      t0 <- Sys.time()
      r <- find_hotspots(s$Y, s$Z, s$A, mc_min = NULL, rounding = rules[[rn]])
      kept <- r$D[up] == 1
      cut_rate  <- if (any(!real)) mean(!kept[!real]) else NA
      keep_rate <- mean(kept[real])
      out[[length(out) + 1]] <- data.frame(
        scenario = scn, rep = rep, rule = rn, ranking = NA, precision = NA, cold = NA,
        edges_kept = mean(kept), iters = sum(r$trace$accepted) - 1, stop = r$stop_reason,
        J = cut_rate + keep_rate - 1, cut = cut_rate, kept = keep_rate,
        secs = as.numeric(Sys.time() - t0, units = "secs"))
    }
  }
  do.call(rbind, out)
}

jobs <- expand.grid(scn = c("A", "B", "D80", "D95", "D0"), rep = seq_len(NREP),
                    stringsAsFactors = FALSE)
jobs$file <- sprintf("%s/%s_r%02d.rds", OUT, jobs$scn, jobs$rep)
todo <- which(!file.exists(jobs$file))
message(length(todo), " of ", nrow(jobs), " datasets to run on ", CORES, " cores")
t0 <- Sys.time()
invisible(mclapply(todo, function(i) {
  r <- tryCatch(run_job(jobs$scn[i], jobs$rep[i]), error = function(e) {
    message("FAILED ", basename(jobs$file[i]), ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, jobs$file[i])
  message(format(Sys.time(), "%H:%M:%S"), "  done ", basename(jobs$file[i]))
}, mc.cores = CORES, mc.preschedule = FALSE))
done <- jobs$file[file.exists(jobs$file)]
res  <- do.call(rbind, lapply(done, readRDS))
write.csv(res, "simulation/results/05_rounding_rules.csv", row.names = FALSE)
message("finished ", length(done), "/", nrow(jobs), " in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
