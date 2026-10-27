# Step 3: small test ------------------------------------------------------------
#
# Goal: before the full study, decide which variants are worth keeping.
#   - weight:   "lambda" (published) vs "gamma2"
#   - candidates: all eigenvectors (published) vs Griffith's rule (>= 0.25 max)
#   - rounding: none (control: no pruning), quantile 0.50, quantile 0.75
#               (published), mixture. Permutation is left for later (slow).
#   - ranking:  R0 (published), R1, R2, R3, and the baselines
#               B0 = top counts, B1 = top OLS residuals of Y on Z
#
# Scenarios (400 cells, 20 x 20 grid, cell side 50, window 1000 x 1000):
#   S1 clusters      3 Thomas clusters (r = 100, 500 points each) + 2000 uniform
#   S2 clusters+cov  S1 + 3000 points driven by an environmental covariate (bump
#                    far from the parents); the method receives Z = (env, noise)
#   S3 cov only      no clusters: 2000 uniform + 3000 covariate-driven points.
#                    Nothing should be called hot once Z is used.
#
# Metrics for the top m = 20 cells (5% of the area):
#   recall    share of parents whose cell is in the top m or touches it
#   precision share of the top m inside a cluster disc
#   cov_fp    share of the top m in the covariate bump (bump > 0.5)
#   cold      share of the top m labelled cold
#
# Run from the repo root:  Rscript simulation/01_small_test.R
# Output: simulation/results/01_small_test.csv

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else library(hotcoldESF)
  library(parallel)
})

NREP  <- as.integer(Sys.getenv("NREP", 20))
CORES <- as.integer(Sys.getenv("CORES", detectCores()))
m     <- 20

scenarios <- list(
  S1 = list(n_parents = 3, sons = 500, n_background = 2000, n_cov = 0,    useZ = FALSE),
  S2 = list(n_parents = 3, sons = 500, n_background = 2000, n_cov = 3000, useZ = TRUE),
  S3 = list(n_parents = 0, sons = 500, n_background = 2000, n_cov = 3000, useZ = TRUE)
)
rules <- list(none = round_none(), q50 = round_quantile(0.5),
              q75 = round_quantile(0.75), mix = round_mixture())
# candidate eigenvectors: all n (published) or Griffith's rule (eigenvalue >= 0.25 max)
mcs <- list(all = NULL, g25 = 0.25)

metrics <- function(top, s, labels = NULL) {
  near <- function(pc) any(top == pc) || any(s$A[pc, top] == 1)
  c(recall    = if (length(s$parent_cells)) mean(sapply(s$parent_cells, near)) else NA,
    precision = if (length(s$footprint)) mean(top %in% s$footprint) else NA,
    cov_fp    = if (length(s$cov_cells)) mean(top %in% s$cov_cells) else NA,
    cold      = if (is.null(labels)) NA else mean(labels[top] == "cold"))
}

run_one <- function(sc_name, rep) {
  sc <- scenarios[[sc_name]]
  set.seed(1000 * match(sc_name, names(scenarios)) + rep)
  s <- sim_thomas(k = 20, n_parents = sc$n_parents, radius = 100, sons = sc$sons,
                  n_background = sc$n_background, n_cov = sc$n_cov)
  Z <- if (sc$useZ) s$Z else NULL
  rows <- list()
  add <- function(variant, ranking, top, labels, extra) {
    rows[[length(rows) + 1]] <<- data.frame(
      scenario = sc_name, rep = rep, variant = variant, ranking = ranking,
      t(metrics(top, s, labels)), extra, check.names = FALSE)
  }
  none <- data.frame(stop = NA, step0_stop = NA, esf_inadequate = NA,
                     edges_kept = NA, iters = NA, n_hot = NA, hot_zones = NA,
                     n_eig0 = NA)

  # baselines
  add("baseline", "B0", order(-s$Y)[1:m], NULL, none)
  X0 <- if (is.null(Z)) matrix(1, length(s$Y), 1) else cbind(1, Z)
  add("baseline", "B1", order(-lm.fit(X0, as.numeric(scale(s$Y)))$residuals)[1:m],
      NULL, none)

  for (mc in names(mcs)) for (w in c("lambda", "gamma2")) for (rn in names(rules)) {
    if (rn == "none" && w == "gamma2") next          # identical to lambda/none
    res <- find_hotspots(s$Y, Z, s$A, weight = w, rounding = rules[[rn]],
                         mc_min = mcs[[mc]])
    extra <- data.frame(
      stop = res$stop_reason,
      step0_stop = grepl("step 0", res$stop_reason),
      esf_inadequate = isTRUE(res$esf_inadequate),
      edges_kept = sum(res$D) / sum(res$A),
      iters = if (is.null(res$trace)) 0 else sum(res$trace$accepted) - 1,
      n_hot = sum(res$labels == "hot"),
      hot_zones = max(c(0, res$zones$hot), na.rm = TRUE))
    for (rk in c("R0", "R1", "R2", "R3"))
      add(paste(mc, w, rn, sep = "/"), rk, rank_cells(res, rk)[1:m], res$labels,
          cbind(extra, n_eig0 = res$trace$n_eig[1]))
  }
  do.call(rbind, rows)
}

jobs <- expand.grid(sc = names(scenarios), rep = seq_len(NREP), stringsAsFactors = FALSE)
t0 <- Sys.time()
out <- mclapply(seq_len(nrow(jobs)), function(i) {
  tryCatch(run_one(jobs$sc[i], jobs$rep[i]),
           error = function(e) { message("job ", i, ": ", conditionMessage(e)); NULL })
}, mc.cores = CORES, mc.preschedule = FALSE)
bad <- which(!vapply(out, is.data.frame, logical(1)))
if (length(bad)) message("failed jobs: ", paste(bad, collapse = ", "))
res <- do.call(rbind, if (length(bad)) out[-bad] else out)
dir.create("simulation/results", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "simulation/results/01_small_test.csv", row.names = FALSE)
message("done: ", nrow(jobs), " datasets in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
