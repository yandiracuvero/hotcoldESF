# Fixing the Moran gate ------------------------------------------------------------------
#
# In the full study (17) the gate passed too often with no clusters: 16% for random
# points on 40 x 40 (sparse counts: log(Y + 1) far from normal) and 19-32% for the
# covariate area alone (a linear term does not capture the covariate effect).
#
# Same datasets as experiment 17 (same seeds). Gate variants:
#   normal / permutation      exact-moment normal approximation vs Freedman-Lane
#                             permutation test (999 permutations)
#   linear / spline           covariates as linear terms vs a natural spline with
#                             4 degrees of freedom for each covariate
# For every variant: does the gate pass (5%)?
# For the proposed new setting (permutation + spline): the full ESF procedure and the
# Gi* / LISA benchmarks with the same gate and covariate terms, to compare with 17.
#
# Settings: C, N (should not pass); B weak, S weak, M weak, K strong (should pass).
# 20 x 20: 100 replicates; 40 x 40: 50 replicates.
#
# Run from the repo root: Rscript simulation/19_gate_fixes.R   (resumable)

source("simulation/17_full_study.R")        # simulate(), seeds, scen; does not run the study
OUT19 <- "simulation/results/19"
dir.create(OUT19, showWarnings = FALSE, recursive = TRUE)
keep <- (jobs$scn %in% c("C", "N")) |
        (jobs$scn %in% c("B", "S", "M") & jobs$strength == "weak") |
        (jobs$scn == "K" & jobs$strength == "strong")
jobs19 <- jobs[keep, ]
jobs19$file <- sub("results/17/", "results/19/", jobs19$file)

run_job19 <- function(j) {
  seed <- 20000 + 100000 * (j$k == 40) + 1000 * match(j$scn, names(scen)) +
          500 * (j$strength == "strong") + j$rep
  s <- simulate(j$k, j$scn, j$strength, seed)
  n <- j$k^2; m <- round(0.05 * n); y <- as.numeric(scale(log1p(s$Y)))
  set.seed(seed + 11)
  gates <- do.call(rbind, lapply(c(0, 4), function(df) {
    X0 <- design0(n, covariate_basis(s$Z, df))
    e  <- lm.fit(X0, y)$residuals
    data.frame(covariates = if (df) "spline" else "linear",
               p_normal = moran_resid(e, s$A, X0)$p,
               p_permutation = moran_perm(e, s$A, X0, 999)$p)
  }))
  gates <- cbind(k = j$k, scenario = j$scn, strength = j$strength, rep = j$rep, gates)

  # full procedure with the permutation gate and spline covariates
  set.seed(seed + 7)
  h <- hot_cold_map(s$Y, s$Z, s$A, nperm = NPERM[as.character(j$k)], penalty = "cv_min",
                    mc_min = 0.25, fit_type = "ridge", gate = "permutation", spline_df = 4)
  X0 <- design0(n, covariate_basis(s$Z, 4))
  res <- lm.fit(X0, y)$residuals
  gate_ok <- isTRUE(h$p0 <= 0.05)
  cl <- lapply(seq_along(s$clusters), function(i) union(s$clusters[[i]], s$parent_cells[i]))
  foot <- sort(unique(unlist(cl))); cold_area <- s$cold_cells
  share <- function(a, b) if (length(a)) mean(a %in% b) else NA
  row <- function(method, hot, cold, ord) {
    top <- if (is.null(ord)) NULL else ord[1:m]
    data.frame(k = j$k, scenario = j$scn, strength = j$strength, rep = j$rep, method = method,
      n_hot = length(hot), hot_precision = share(hot, foot),
      hot_sensitivity = if (length(foot)) mean(foot %in% hot) else NA,
      n_cold = length(cold), cold_precision = share(cold, cold_area),
      cold_sensitivity = if (length(cold_area)) mean(cold_area %in% cold) else NA,
      top_precision = if (length(foot)) mean(top %in% foot) else NA,
      clusters_top = if (length(cl)) mean(sapply(cl, function(c) any(c %in% top))) else NA,
      gate_ok = gate_ok)
  }
  none <- integer(0)
  gi <- local_gi_star(res, s$A)
  set.seed(seed + 9); li <- local_moran(res, s$A)
  methods <- rbind(
    row("ESF", which(h$labels == "hot"), which(h$labels == "cold"), rank_cells(h, "R2")),
    row("Gi* + gate", if (gate_ok) which(gi > qnorm(0.95)) else none,
        if (gate_ok) which(gi < qnorm(0.05)) else none, order(-gi)),
    row("LISA + gate", if (gate_ok) which(li$quad == "HH" & li$p < 0.05) else none,
        if (gate_ok) which(li$quad == "LL" & li$p < 0.05) else none, order(-(li$quad == "HH"), -li$I)),
    row("Top residuals", none, none, order(-res)))
  list(gates = gates, methods = methods)
}

todo <- which(!file.exists(jobs19$file))
message(length(todo), " of ", nrow(jobs19), " datasets to run on ", CORES, " cores")
t0 <- Sys.time()
invisible(mclapply(todo, function(i) {
  r <- tryCatch(run_job19(jobs19[i, ]), error = function(e) {
    message("FAILED ", basename(jobs19$file[i]), ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, jobs19$file[i])
  if (i %% 25 == 0) message(format(Sys.time(), "%H:%M:%S"), "  done ", basename(jobs19$file[i]))
}, mc.cores = CORES, mc.preschedule = FALSE))
done <- jobs19$file[file.exists(jobs19$file)]
rr <- lapply(done, readRDS)
write.csv(do.call(rbind, lapply(rr, `[[`, "gates")), "simulation/results/19_gate_pvalues.csv", row.names = FALSE)
write.csv(do.call(rbind, lapply(rr, `[[`, "methods")), "simulation/results/19_gate_methods.csv", row.names = FALSE)
message("finished ", length(done), "/", nrow(jobs19), " in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
