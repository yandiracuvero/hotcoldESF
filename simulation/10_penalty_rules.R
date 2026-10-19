# How the Lasso penalty is chosen: effect on the ESF hot/cold map ------------------------
#
# Same datasets as simulation/07_hotcold_benchmarks.R (same seeds). For each dataset, the
# one-pass procedure (log counts, Moran gate, ESF fit, permutation band) is run with:
#
#   penalty rule (see R/penalty.R)
#     z^-1, z^-2 (default), z^-3, z^-4   lambda = |z0|^(-a)
#     cv_min, cv_1se                      10-fold cross-validation
#     bic                                 BIC along the Lasso path
#     moran                               largest lambda whose residuals pass the Moran test
#   candidate eigenvectors
#     all       all n eigenvectors
#     positive  eigenvalue >= 0.25 x the largest (smooth patterns only, Griffith's rule)
#
# 20 x 20 grid: A, B50, B200, B500, C, N; 20 replicates; 200 permutations.
# 40 x 40 grid: B50, B200, C; 10 replicates; 100 permutations ("moran" with all
#   eigenvectors is skipped on 40 x 40: too slow).
# Gi* on the residuals of log counts on Z is recomputed as the reference.
#
# Outputs per dataset and variant: lambda, number of eigenvectors selected, cells
# labelled hot, label precision, sensitivity, precision of the top 5% of cells
# ranked by the ESF rule "largest excess" and by f-hat alone.
#
# Run from the repo root: Rscript simulation/10_penalty_rules.R   (resumable)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else library(hotspots)
  library(parallel)
})
CORES <- as.integer(Sys.getenv("CORES", detectCores()))-2
OUT   <- "simulation/results/10"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cfg <- rbind(
  data.frame(k = 20, scn = c("A", "B50", "B200", "B500", "C", "N"), nrep = 20, nperm = 200),
  data.frame(k = 40, scn = c("B50", "B200", "C"), nrep = 10, nperm = 100))
jobs <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i)
  data.frame(cfg[i, c("k", "scn", "nperm")], rep = seq_len(cfg$nrep[i]), row.names = NULL)))
jobs$file <- sprintf("%s/k%02d_%s_r%02d.rds", OUT, jobs$k, jobs$scn, jobs$rep)
jobs <- jobs[order(-jobs$k), ]

variants <- rbind(
  data.frame(penalty = "z", a = 1:4, label = paste0("z^-", 1:4)),
  data.frame(penalty = c("cv_min", "cv_1se", "bic", "moran"), a = 2,
             label = c("CV (min)", "CV (1 SE)", "BIC", "Moran stop")))
variants <- merge(variants, data.frame(cand = c("all", "positive")))

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
  X0 <- if (is.null(s$Z)) matrix(1, n, 1) else cbind(1, scale(s$Z))
  res <- lm.fit(X0, as.numeric(scale(log1p(s$Y))))$residuals
  prec <- function(ord) if (length(foot)) mean(ord[1:m] %in% foot) else NA
  gi <- local_gi_star(res, s$A)
  ref <- data.frame(label = "Gi* residuals", cand = "-", lambda = NA, n_eig = NA,
                    n_hot = sum(gi > qnorm(0.95)),
                    label_precision = if (length(foot) && any(gi > qnorm(0.95))) mean(which(gi > qnorm(0.95)) %in% foot) else NA,
                    sensitivity = if (length(foot)) mean(foot %in% which(gi > qnorm(0.95))) else NA,
                    prec_excess = prec(order(-gi)), prec_fhat = prec(order(-gi)),
                    stop = "-", secs = 0)
  rows <- lapply(seq_len(nrow(variants)), function(v) {
    va <- variants[v, ]
    if (j$k == 40 && va$penalty == "moran" && va$cand == "all") return(NULL)
    set.seed(seed + 100 * v)
    t0 <- Sys.time()
    h <- hot_cold_map(s$Y, s$Z, s$A, nperm = j$nperm, a = va$a, penalty = va$penalty,
                      mc_min = if (va$cand == "all") NULL else 0.25)
    hot <- which(h$labels == "hot")
    data.frame(label = va$label, cand = va$cand,
               lambda = if (is.null(h$lambda)) NA else h$lambda,
               n_eig = if (is.null(h$fit)) 0 else length(h$fit$S),
               n_hot = length(hot),
               label_precision = if (length(hot) && length(foot)) mean(hot %in% foot) else NA,
               sensitivity = if (length(foot)) mean(foot %in% hot) else NA,
               prec_excess = prec(rank_cells(h, "R2")),
               prec_fhat = prec(order(-h$fhat, -s$Y)),
               stop = h$stop_reason,
               secs = as.numeric(Sys.time() - t0, units = "secs"))
  })
  cbind(k = j$k, scenario = j$scn, rep = j$rep, rbind(ref, do.call(rbind, rows)))
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
write.csv(res, "simulation/results/10_penalty_rules.csv", row.names = FALSE)
message("finished ", length(done), "/", nrow(jobs), " in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
