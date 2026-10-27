# Two follow-up tests on the same 150 datasets as experiments 07 and 10 ----------------------
#
# Test 1. Does the Moran gate alone give Gi* the same false-positive control?
#   Gi* on the residuals of log(Y + 1) on Z, hot if z > 1.645, in four versions:
#     Gi*                    no gate, no correction (as in experiments 07 and 10)
#     Gi* + gate             no hot cells unless the Moran test of the residuals is
#                            significant (5%), the same gate as the ESF procedure
#     Gi* + FDR              Benjamini-Hochberg over the cells (5%), no gate
#     Gi* + gate + FDR       both
#   The ranking (top 5%) is the same for the four versions: cells ordered by Gi* z.
#
# Test 2. Can the ESF fit be improved with weak clusters? Smooth candidates only
#   (eigenvalue >= 0.25 x the largest), penalty chosen by cross-validation (lowest
#   error), in four versions:
#     Lasso, random folds    as in experiment 10 (reference)
#     Lasso, spatial folds   folds = square blocks of 200 x 200 units (4 x 4 cells on
#                            20 x 20, 8 x 8 cells on 40 x 40), assigned at random to 10 folds
#     Ridge, random folds    ridge instead of Lasso: all smooth eigenvectors kept, shrunk
#     Ridge, spatial folds
#
# Run from the repo root: Rscript simulation/13_gate_and_cv.R   (resumable)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else
    if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(hotcoldESF)
  library(parallel)
})
CORES <- as.integer(Sys.getenv("CORES", detectCores()))
OUT   <- "simulation/results/13"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cfg <- rbind(
  data.frame(k = 20, scn = c("A", "B50", "B200", "B500", "C", "N"), nrep = 20, nperm = 200),
  data.frame(k = 40, scn = c("B50", "B200", "C"), nrep = 10, nperm = 100))
jobs <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i)
  data.frame(cfg[i, c("k", "scn", "nperm")], rep = seq_len(cfg$nrep[i]), row.names = NULL)))
jobs$file <- sprintf("%s/k%02d_%s_r%02d.rds", OUT, jobs$k, jobs$scn, jobs$rep)
jobs <- jobs[order(-jobs$k), ]

esf_variants <- data.frame(fit = c("lasso", "lasso", "ridge", "ridge"),
                           folds = c("random", "spatial", "random", "spatial"))

simulate <- function(k, scn) switch(substr(scn, 1, 1),
  A = sim_thomas(k, 3, 100, 500, n_background = 2000),
  B = sim_thomas(k, 3, 100, as.numeric(substring(scn, 2)), n_background = 2000, n_cov = 3000),
  C = sim_thomas(k, 0, 100, 0, n_background = 2000, n_cov = 3000),
  N = sim_thomas(k, 0, 100, 0, n_background = 3000))

run_job <- function(j) {
  seed <- 7000 + 1000 * j$k + 10 * match(j$scn, c("A", "B50", "B200", "B500", "C", "N")) + j$rep
  set.seed(seed)
  s <- simulate(j$k, j$scn)
  n <- j$k^2; m <- round(0.05 * n); foot <- s$footprint
  Zs <- if (is.null(s$Z)) NULL else scale(s$Z)
  X0 <- design0(n, Zs)
  y  <- as.numeric(scale(log1p(s$Y)))
  res <- lm.fit(X0, y)$residuals
  gate_ok <- isTRUE(moran_resid(res, s$A, X0)$p <= 0.05)
  prec <- function(ord) if (length(foot)) mean(ord[1:m] %in% foot) else NA
  row <- function(method, hot, ord, extra = list()) data.frame(
    k = j$k, scenario = j$scn, rep = j$rep, method = method,
    n_hot = length(hot),
    label_precision = if (length(hot) && length(foot)) mean(hot %in% foot) else NA,
    sensitivity = if (length(foot)) mean(foot %in% hot) else NA,
    prec_top = prec(ord), gate_ok = gate_ok,
    lambda = if (is.null(extra$lambda)) NA else extra$lambda,
    n_eig = if (is.null(extra$n_eig)) NA else extra$n_eig,
    secs = if (is.null(extra$secs)) 0 else extra$secs)

  # ---- test 1: Gi* with and without gate / FDR ----
  gi <- local_gi_star(res, s$A)
  p  <- stats::pnorm(gi, lower.tail = FALSE)
  hot_raw <- which(gi > stats::qnorm(0.95))
  hot_fdr <- which(stats::p.adjust(p, "BH") <= 0.05)
  ord_gi  <- order(-gi)
  out <- list(
    row("Gi*", hot_raw, ord_gi),
    row("Gi* + gate", if (gate_ok) hot_raw else integer(0), ord_gi),
    row("Gi* + FDR", hot_fdr, ord_gi),
    row("Gi* + gate + FDR", if (gate_ok) hot_fdr else integer(0), ord_gi))

  # ---- test 2: ESF, smooth candidates, CV penalty ----
  block <- round(200 / s$cell_side)
  for (v in seq_len(nrow(esf_variants))) {
    va <- esf_variants[v, ]
    set.seed(seed + 100 * v)
    fid <- if (va$folds == "spatial") block_folds(j$k, block = block, nfolds = 10) else NULL
    t0 <- Sys.time()
    h <- hot_cold_map(s$Y, s$Z, s$A, nperm = j$nperm, penalty = "cv_min", mc_min = 0.25,
                      fit_type = va$fit, foldid = fid)
    out[[length(out) + 1]] <- row(
      sprintf("ESF %s, %s folds", if (va$fit == "lasso") "Lasso" else "ridge", va$folds),
      which(h$labels == "hot"), rank_cells(h, "R2"),
      list(lambda = h$lambda,
           n_eig = if (is.null(h$fit)) 0 else if (va$fit == "ridge") length(h$fit$S) else length(h$fit$S),
           secs = as.numeric(Sys.time() - t0, units = "secs")))
  }
  do.call(rbind, out)
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
write.csv(res, "simulation/results/13_gate_and_cv.csv", row.names = FALSE)
message("finished ", length(done), "/", nrow(jobs), " in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
