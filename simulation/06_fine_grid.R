# Fine grid, weak clusters: can rounding help where zone ranking already wins? -----------
#
# Experiment 2 found one setting where ranking by hot zones beat every benchmark:
# a 40 x 40 grid (1600 cells, cell side 25) with weak clusters (50 or 200 points per
# cluster of radius 100). This experiment asks whether the iterative reconstruction
# (rounding + pruning) improves the rankings further in that setting.
#
# Scenario B (clusters + covariate area, covariate given to the model), 40 x 40 grid,
# 50 and 200 points per cluster, 10 replicates each.
# Fixed (baseline configuration): standardised counts, all eigenvectors,
# lambda = |z0|^(-2), eigenvalue weights, nested pruning, Moran stop with A.
# Rounding rules: none, 75% quantile (baseline), 90% quantile, permutation (30).
# Rankings: centrality whole graph (baseline), centrality hot cells, largest excess,
# hot zones first. Benchmarks: top counts, top residuals (counts and log counts).
# Budget: top 80 cells (5%).
#
# Run from repo root: Rscript simulation/06_fine_grid.R   (resumable; ~70 min on 2 cores,
# much faster on a Mac with more cores)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else library(hotspots)
  library(parallel)
})
NREP  <- as.integer(Sys.getenv("NREP", 10))
CORES <- as.integer(Sys.getenv("CORES", detectCores()))
OUT   <- "simulation/results/06"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rules <- list(none = round_none(), q75 = round_quantile(0.75), q90 = round_quantile(0.90),
              permutation = round_permutation(30, 0.95))
rank_names <- c(R0 = "Centrality, whole graph", R1 = "Centrality, hot cells",
                R2 = "Largest excess", R3 = "Hot zones first")
m <- 80

run_job <- function(sons, rep) {
  set.seed(6000 + sons + rep)
  s <- sim_thomas(k = 40, n_parents = 3, radius = 100, sons = sons,
                  n_background = 2000, n_cov = 3000)
  row <- function(rule, ranking, top, r = NULL, secs = NA)
    data.frame(sons = sons, rep = rep, rule = rule, ranking = ranking,
               precision = mean(top %in% s$footprint),
               recall = mean(sapply(s$parent_cells, function(pc) any(top == pc) || any(s$A[pc, top] == 1))),
               cov_fp = mean(top %in% s$cov_cells),
               edges_kept = if (is.null(r)) NA else sum(r$D) / sum(r$A),
               iters = if (is.null(r)) NA else sum(r$trace$accepted) - 1,
               stop = if (is.null(r)) NA else r$stop_reason, secs = secs)
  X0 <- cbind(1, s$Z)
  out <- list(
    row("benchmark", "Top counts", order(-s$Y)[1:m]),
    row("benchmark", "Top residuals", order(-lm.fit(X0, as.numeric(scale(s$Y)))$residuals)[1:m]),
    row("benchmark", "Top residuals (log)", order(-lm.fit(X0, as.numeric(scale(log1p(s$Y))))$residuals)[1:m]))
  for (rn in names(rules)) {
    t0 <- Sys.time()
    r <- find_hotspots(s$Y, s$Z, s$A, mc_min = NULL, rounding = rules[[rn]])
    secs <- as.numeric(Sys.time() - t0, units = "secs")
    for (rk in names(rank_names))
      out[[length(out) + 1]] <- row(rn, rank_names[[rk]], rank_cells(r, rk)[1:m], r, secs)
  }
  do.call(rbind, out)
}

jobs <- expand.grid(sons = c(50, 200), rep = seq_len(NREP))
jobs$file <- sprintf("%s/s%03d_r%02d.rds", OUT, jobs$sons, jobs$rep)
todo <- which(!file.exists(jobs$file))
message(length(todo), " of ", nrow(jobs), " datasets to run on ", CORES, " cores")
t0 <- Sys.time()
invisible(mclapply(todo, function(i) {
  r <- tryCatch(run_job(jobs$sons[i], jobs$rep[i]), error = function(e) {
    message("FAILED ", basename(jobs$file[i]), ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, jobs$file[i])
  message(format(Sys.time(), "%H:%M:%S"), "  done ", basename(jobs$file[i]))
}, mc.cores = CORES, mc.preschedule = FALSE))
done <- jobs$file[file.exists(jobs$file)]
res  <- do.call(rbind, lapply(done, readRDS))
write.csv(res, "simulation/results/06_fine_grid.csv", row.names = FALSE)
message("finished ", length(done), "/", nrow(jobs), " in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
