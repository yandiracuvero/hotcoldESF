# Smoke test: one Thomas simulation, every weight x rounding option, every ranking.
# Not the simulation study. It only checks that everything runs and gives a first look.
# Run from the repo root after devtools::load_all(), or with the package installed.

library(hotcoldESF)
set.seed(10)
s <- sim_thomas(k = 20, n_parents = 3, radius = 100, sons = 1000, n_background = 1500)
cat("n =", length(s$Y), " parent cells:", s$parent_cells, "\n\n")
rules <- list(q75 = round_quantile(.75), q50 = round_quantile(.5), mix = round_mixture(),
              perm = round_permutation(30))
for (w in c("lambda", "gamma2")) for (nm in names(rules)) {
  t0 <- Sys.time()
  r <- find_hotspots(s$Y, A = s$A, weight = w, rounding = rules[[nm]])
  e0 <- sum(r$A)/2; e1 <- sum(r$D)/2
  m <- 20
  hits <- sapply(c("B0","R0","R1","R2","R3"), function(meth) {
    top <- rank_cells(r, meth)[1:m]
    sum(sapply(s$parent_cells, function(pc) any(top == pc | r$A[pc, top] == 1)))
  })
  cold_in_R0 <- sum(r$labels[rank_cells(r, "R0")[1:m]] == "cold")
  cat(sprintf("%-6s %-4s edges %4d->%4d iters %2d stop=%-12s hot=%3d cold=%3d | parents found in top-20: %s | cold cells in R0 top-20: %d (%.1fs)\n",
      w, nm, e0, e1, sum(r$trace$accepted)-1, r$stop_reason, sum(r$labels=="hot"), sum(r$labels=="cold"),
      paste(names(hits), hits, sep="=", collapse=" "), cold_in_R0, as.numeric(Sys.time()-t0, units="secs")))
}
