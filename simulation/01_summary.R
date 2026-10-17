# Summary tables for simulation/01_small_test.R
# Run from the repo root: Rscript simulation/01_summary.R

x <- read.csv("simulation/results/01_small_test.csv", stringsAsFactors = FALSE)
r2 <- function(v) round(mean(v, na.rm = TRUE), 2)

cat("\n=== A. Ranking quality (top 20 cells), mean over replicates ===\n")
A <- aggregate(cbind(recall, precision, cov_fp, cold) ~ scenario + variant + ranking,
               data = x, FUN = r2, na.action = na.pass)
for (sc in unique(A$scenario)) {
  cat("\n--", sc, "--\n")
  a <- A[A$scenario == sc, -1]
  a <- a[order(-a$precision, a$cov_fp), ]
  print(a, row.names = FALSE)
}

cat("\n=== B. Algorithm behaviour ===\n")
y <- x[x$ranking == "R0" & x$variant != "baseline", ]
B <- aggregate(cbind(edges_kept, iters, n_eig0, esf_inadequate, step0_stop) ~
                 scenario + variant, data = y, FUN = r2, na.action = na.pass)
print(B, row.names = FALSE)

cat("\n=== C. Stop reasons ===\n")
print(with(y, table(paste(scenario, variant), stop)))
