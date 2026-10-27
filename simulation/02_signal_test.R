# Step 3b: signal strength x cell size x penalty exponent x response scale ----------
#
# Question: when does ESF beat the simple baselines? Expected: weak clusters spread
# over several cells, where single-cell residuals are noisy.
#
# Factors
#   k      cells per side: 10, 20, 40 (cell side 100, 50, 25; cluster radius 100)
#   sons   points per cluster: 50, 200, 500
#   a      penalty exponent, lambda = |z0|^(-a): 0.5, 1, 2
#   resp   response: standardised counts ("std") or log(Y + 1) ("log")
#   prune  none vs quantile 0.75 (nested pruning)
# Scenarios
#   S2  3 clusters + 2000 uniform + 3000 covariate-driven points (Z given)
#   S3  no clusters: 2000 uniform + 3000 covariate-driven points (nothing is hot)
# Rankings: R0 (published), R2 (hot cells by fhat), R3 (zones first);
# baselines: B0 top counts, B1 top OLS residuals on Z (std and log).
# Budget m = 5% of the cells.
#
# HOW TO RUN (on the Mac, from the repo root, in Terminal):
#     caffeinate -i Rscript simulation/02_signal_test.R
# - uses all cores but one; set CORES=4 to change, NREP=10 replicates by default
# - it is RESUMABLE: every dataset is saved in simulation/results/02/ as it
#   finishes; if you stop it (Ctrl+C) or the Mac sleeps, run it again and it
#   continues where it stopped
# - at the end it writes simulation/results/02_signal_test.csv
# Needs the R packages: glmnet, igraph, devtools (install.packages(...) once).

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else
    if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(hotcoldESF)
  library(parallel)
})

NREP  <- as.integer(Sys.getenv("NREP", 10))
CORES <- as.integer(Sys.getenv("CORES", max(1, detectCores() - 1)))
OUT   <- "simulation/results/02"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

jobs <- rbind(
  expand.grid(scenario = "S2", k = c(10, 20, 40), sons = c(50, 200, 500),
              rep = seq_len(NREP), stringsAsFactors = FALSE),
  expand.grid(scenario = "S3", k = c(10, 20, 40), sons = NA,
              rep = seq_len(NREP), stringsAsFactors = FALSE))
jobs$file <- with(jobs, sprintf("%s/%s_k%02d_s%s_r%02d.rds", OUT, scenario, k,
                                ifelse(is.na(sons), "NA", sons), rep))
# biggest grids first so the long jobs do not end up alone at the end
jobs <- jobs[order(-jobs$k), ]

metrics <- function(top, s, labels = NULL) {
  near <- function(pc) any(top == pc) || any(s$A[pc, top] == 1)
  c(recall    = if (length(s$parent_cells)) mean(sapply(s$parent_cells, near)) else NA,
    precision = if (length(s$footprint)) mean(top %in% s$footprint) else NA,
    cov_fp    = mean(top %in% s$cov_cells),
    cold      = if (is.null(labels)) NA else mean(labels[top] == "cold"))
}

run_job <- function(j) {
  set.seed(7919 * j$k + 101 * ifelse(is.na(j$sons), 0, j$sons) + j$rep +
             ifelse(j$scenario == "S3", 1e5, 0))
  s <- sim_thomas(k = j$k, n_parents = if (j$scenario == "S2") 3 else 0,
                  radius = 100, sons = ifelse(is.na(j$sons), 0, j$sons),
                  n_background = 2000, n_cov = 3000)
  n <- j$k^2; m <- max(3, round(0.05 * n))
  rows <- list()
  add <- function(method, ranking, top, labels = NULL, info = list()) {
    rows[[length(rows) + 1]] <<- data.frame(
      scenario = j$scenario, k = j$k, sons = j$sons, rep = j$rep, m = m,
      method = method, ranking = ranking, t(metrics(top, s, labels)),
      lambda = if (is.null(info$lambda)) NA else info$lambda,
      n_eig0 = if (is.null(info$n_eig0)) NA else info$n_eig0,
      edges_kept = if (is.null(info$edges_kept)) NA else info$edges_kept,
      stop = if (is.null(info$stop)) NA else info$stop,
      secs = if (is.null(info$secs)) NA else info$secs)
  }
  X0 <- cbind(1, s$Z)
  add("B0", "-", order(-s$Y)[1:m])
  add("B1_std", "-", order(-lm.fit(X0, as.numeric(scale(s$Y)))$residuals)[1:m])
  add("B1_log", "-", order(-lm.fit(X0, as.numeric(scale(log1p(s$Y))))$residuals)[1:m])

  for (resp in c("std", "log")) for (a in c(0.5, 1, 2)) for (pr in c("none", "q75")) {
    t0 <- Sys.time()
    res <- find_hotspots(s$Y, s$Z, s$A, a = a, response = resp, mc_min = NULL,
                         rounding = if (pr == "none") round_none() else round_quantile(0.75))
    info <- list(lambda = res$lambda,
                 n_eig0 = if (is.null(res$trace)) 0 else res$trace$n_eig[1],
                 edges_kept = sum(res$D) / sum(res$A), stop = res$stop_reason,
                 secs = as.numeric(Sys.time() - t0, units = "secs"))
    meth <- sprintf("%s_a%s_%s", resp, a, pr)
    for (rk in c("R0", "R2", "R3"))
      add(meth, rk, rank_cells(res, rk)[1:m], res$labels, info)
  }
  do.call(rbind, rows)
}

todo <- which(!file.exists(jobs$file))
message(length(todo), " of ", nrow(jobs), " datasets to run on ", CORES, " cores")
t0 <- Sys.time()
invisible(mclapply(todo, function(i) {
  r <- tryCatch(run_job(jobs[i, ]), error = function(e) {
    message("FAILED ", basename(jobs$file[i]), ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, jobs$file[i])
  message(format(Sys.time(), "%H:%M:%S"), "  done ", basename(jobs$file[i]))
  NULL
}, mc.cores = CORES, mc.preschedule = FALSE))

done <- jobs$file[file.exists(jobs$file)]
res  <- do.call(rbind, lapply(done, readRDS))
write.csv(res, "simulation/results/02_signal_test.csv", row.names = FALSE)
message("finished ", length(done), "/", nrow(jobs), " datasets in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
