# Hot/cold map with a permutation band vs Getis-Ord Gi* and LISA ------------------------
#
# The proposed one-pass procedure (hot_cold_map): log(Y + 1), Moran test, MI-Lasso
# ESF fit on the eigenvectors of A (a = 2, all eigenvectors), fhat = spatial
# component; permutation band (Freedman-Lane, same lambda): hot if fhat is above
# the 95th percentile of fhat under no spatial structure, cold if below the 5th.
#   pointwise band: one quantile per cell (about 5% false hot cells expected)
#   max band:       quantile of the map maximum (controls any false hot cell)
#
# Benchmarks (all on the same queen contiguity A):
#   Gi* on counts and on residuals of log(Y + 1) on Z   (hot: z > 1.645)
#   LISA on counts and on residuals of log(Y + 1) on Z  (hot: High-High, p < 0.05)
#   top counts, top residuals of log(Y + 1) on Z        (ranking only)
#
# Scenarios
#   20 x 20 grid (cell side 50), 20 replicates, 200 permutations:
#     A    3 clusters x 500 points + 2000 random points, no covariate
#     B50, B200, B500   3 clusters x 50/200/500 points + 2000 random + 3000 points in a
#                       covariate area (the model receives the covariate)
#     C    covariate area only (no clusters)            -> nothing should be hot
#     N    3000 random points (complete randomness)     -> nothing should be hot
#   40 x 40 grid (cell side 25), 10 replicates, 100 permutations: B50, B200, C
#
# Outputs per dataset and method
#   labels: n_hot, label_precision (share of hot cells inside a cluster disc),
#           sensitivity (share of cluster cells labelled hot), cov_hot (share of
#           hot cells inside the covariate area)
#   ranking (top 5% of cells): precision, recall (clusters touched), cov_fp
#
# Run from the repo root: Rscript simulation/07_hotcold_benchmarks.R  (resumable)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else
    if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(hotcoldESF)
  library(parallel)
})
CORES <- as.integer(Sys.getenv("CORES", detectCores()))
OUT   <- "simulation/results/07"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cfg <- rbind(
  data.frame(k = 20, scn = c("A", "B50", "B200", "B500", "C", "N"), nrep = 20, nperm = 200),
  data.frame(k = 40, scn = c("B50", "B200", "C"), nrep = 10, nperm = 100))
jobs <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i)
  data.frame(cfg[i, c("k", "scn", "nperm")], rep = seq_len(cfg$nrep[i]))))
jobs$file <- sprintf("%s/k%02d_%s_r%02d.rds", OUT, jobs$k, jobs$scn, jobs$rep)
jobs <- jobs[order(-jobs$k), ]                            # long jobs first

simulate <- function(k, scn) {
  switch(substr(scn, 1, 1),
    A = sim_thomas(k, 3, 100, 500, n_background = 2000),
    B = sim_thomas(k, 3, 100, as.numeric(substring(scn, 2)), n_background = 2000, n_cov = 3000),
    C = sim_thomas(k, 0, 100, 0, n_background = 2000, n_cov = 3000),
    N = sim_thomas(k, 0, 100, 0, n_background = 3000))
}

run_job <- function(j) {
  set.seed(7000 + 1000 * j$k + 10 * match(j$scn, c("A", "B50", "B200", "B500", "C", "N")) + j$rep)
  s <- simulate(j$k, j$scn)
  n <- j$k^2; m <- round(0.05 * n)
  Z <- s$Z
  X0 <- if (is.null(Z)) matrix(1, n, 1) else cbind(1, scale(Z))
  ylog <- as.numeric(scale(log1p(s$Y)))
  res  <- lm.fit(X0, ylog)$residuals

  lab_row <- function(method, hot) data.frame(
    k = j$k, scenario = j$scn, rep = j$rep, kind = "labels", method = method,
    n_hot = length(hot),
    label_precision = if (length(hot) && length(s$footprint)) mean(hot %in% s$footprint) else NA,
    sensitivity = if (length(s$footprint)) mean(s$footprint %in% hot) else NA,
    cov_hot = if (length(hot) && length(s$cov_cells)) mean(hot %in% s$cov_cells) else NA,
    precision = NA, recall = NA, cov_fp = NA)
  rank_row <- function(method, ord) { top <- ord[1:m]; data.frame(
    k = j$k, scenario = j$scn, rep = j$rep, kind = "ranking", method = method,
    n_hot = NA, label_precision = NA, sensitivity = NA, cov_hot = NA,
    precision = if (length(s$footprint)) mean(top %in% s$footprint) else NA,
    recall = if (length(s$parent_cells)) mean(sapply(s$parent_cells, function(pc)
      any(top == pc) || any(s$A[pc, top] == 1))) else NA,
    cov_fp = if (length(s$cov_cells)) mean(top %in% s$cov_cells) else NA) }

  t0 <- Sys.time()
  h  <- hot_cold_map(s$Y, Z, s$A, nperm = j$nperm, penalty = "z", mc_min = NULL, fit_type = "lasso")
  esf_secs <- as.numeric(Sys.time() - t0, units = "secs")
  lab_pw <- if (is.null(h$labels_pointwise)) h$labels else h$labels_pointwise
  lab_mx <- if (is.null(h$labels_max)) h$labels else h$labels_max

  gi_c <- local_gi_star(s$Y, s$A); gi_r <- local_gi_star(res, s$A)
  lm_c <- local_moran(s$Y, s$A);   lm_r <- local_moran(res, s$A)
  hh   <- function(l) which(l$quad == "HH" & l$p < 0.05)
  hhord <- function(l) order(-(l$quad == "HH"), -l$I)

  out <- rbind(
    lab_row("ESF band (pointwise)", which(lab_pw == "hot")),
    lab_row("ESF band (max)", which(lab_mx == "hot")),
    lab_row("Gi* counts", which(gi_c > qnorm(0.95))),
    lab_row("Gi* residuals", which(gi_r > qnorm(0.95))),
    lab_row("LISA counts", hh(lm_c)),
    lab_row("LISA residuals", hh(lm_r)),
    rank_row("ESF: largest excess", rank_cells(h, "R2")),
    rank_row("ESF: hot zones first", rank_cells(h, "R3")),
    rank_row("ESF: centrality, hot cells", rank_cells(h, "R1")),
    rank_row("Gi* counts", order(-gi_c)),
    rank_row("Gi* residuals", order(-gi_r)),
    rank_row("LISA counts", hhord(lm_c)),
    rank_row("LISA residuals", hhord(lm_r)),
    rank_row("Top counts", order(-s$Y)),
    rank_row("Top residuals", order(-res)))
  out$esf_stop <- h$stop_reason
  out$esf_secs <- esf_secs
  out$n_cold_pointwise <- sum(lab_pw == "cold")
  out
}

todo <- which(!file.exists(jobs$file))
message(length(todo), " of ", nrow(jobs), " datasets to run on ", CORES, " cores")
t0 <- Sys.time()
invisible(mclapply(todo, function(i) {
  r <- tryCatch(run_job(jobs[i, ]), error = function(e) {
    message("FAILED ", basename(jobs$file[i]), ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, jobs$file[i])
  message(format(Sys.time(), "%H:%M:%S"), "  done ", basename(jobs$file[i]))
}, mc.cores = CORES, mc.preschedule = FALSE))
done <- jobs$file[file.exists(jobs$file)]
res  <- do.call(rbind, lapply(done, readRDS))
write.csv(res, "simulation/results/07_hotcold_benchmarks.csv", row.names = FALSE)
message("finished ", length(done), "/", nrow(jobs), " in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
