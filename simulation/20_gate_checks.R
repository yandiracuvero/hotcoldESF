# Two quick checks of the Moran gate (gate only, no ESF fit) ---------------------------
#  1. Random points (N) on 40 x 40 with 500 replicates: is the 14-16% pass rate in
#     experiment 19 (50 replicates) real or sampling noise? Normal and permutation tests.
#  2. Spline degrees of freedom 2 and 3 (between linear and the df = 4 of experiment 19):
#     pass rate in C (should be low) and in B weak, S weak (should be high). Same
#     datasets as experiments 17 and 19.
# Run from the repo root: Rscript simulation/20_gate_checks.R
source("simulation/17_full_study.R")
library(parallel)

# ---- check 1 ----
n_rep <- 500
chk1 <- do.call(rbind, mclapply(seq_len(n_rep), function(r) {
  set.seed(500000 + r)
  s <- sim_thomas(40, 0, 100, 0, n_background = 3000)
  y <- as.numeric(scale(log1p(s$Y))); X0 <- matrix(1, 1600, 1); e <- y - mean(y)
  data.frame(rep = r, p_normal = moran_resid(e, s$A, X0)$p, p_permutation = moran_perm(e, s$A, X0, 499)$p)
}, mc.cores = CORES))
write.csv(chk1, "simulation/results/20_gate_N40.csv", row.names = FALSE)

# ---- check 2 ----
sel <- jobs[(jobs$scn == "C") | (jobs$scn %in% c("B", "S") & jobs$strength == "weak"), ]
chk2 <- do.call(rbind, mclapply(seq_len(nrow(sel)), function(i) {
  j <- sel[i, ]
  seed <- 20000 + 100000 * (j$k == 40) + 1000 * match(j$scn, names(scen)) +
          500 * (j$strength == "strong") + j$rep
  s <- simulate(j$k, j$scn, j$strength, seed)
  y <- as.numeric(scale(log1p(s$Y)))
  do.call(rbind, lapply(0:4, function(df) {
    X0 <- design0(j$k^2, covariate_basis(s$Z, df))
    data.frame(k = j$k, scenario = j$scn, strength = j$strength, rep = j$rep, df = df,
               p_normal = moran_resid(lm.fit(X0, y)$residuals, s$A, X0)$p)
  }))
}, mc.cores = CORES, mc.preschedule = FALSE))
write.csv(chk2, "simulation/results/20_gate_df.csv", row.names = FALSE)
message("done")
