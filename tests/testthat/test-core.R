test_that("queen grid has the right degrees", {
  A <- grid_adjacency(4)
  expect_true(isSymmetric(A))
  expect_equal(sort(unique(rowSums(A))), c(3, 5, 8))   # corner, border, interior
  expect_equal(sum(grid_adjacency(4, queen = FALSE)), 2 * 24)
})

test_that("Moran p-values are roughly uniform under the null", {
  set.seed(1)
  A <- grid_adjacency(10)
  X <- cbind(1, rnorm(100))
  p <- replicate(300, moran_resid(lm.fit(X, rnorm(100))$residuals, A, X)$p)
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.04)
})

test_that("B is invariant to flipping the sign of an eigenvector", {
  set.seed(2)
  E <- qr.Q(qr(matrix(rnorm(25 * 3), 25)))
  fit  <- list(E_S = E, values_S = c(3, 2, 1), gamma_S = c(.5, -.2, .1))
  fit2 <- fit; fit2$E_S[, 2] <- -fit2$E_S[, 2]; fit2$gamma_S[2] <- -fit2$gamma_S[2]
  for (w in c("lambda", "gamma2"))
    expect_equal(reconstruct_B(fit, w), reconstruct_B(fit2, w))
  # ... while fhat = E_S gamma_S is also invariant (and carries the sign)
  expect_equal(as.numeric(fit$E_S %*% fit$gamma_S),
               as.numeric(fit2$E_S %*% fit2$gamma_S))
})

test_that("pruning is nested and the loop terminates", {
  set.seed(3)
  s <- sim_thomas(k = 12, n_parents = 2, radius = 120, sons = 800)
  for (rule in list(round_quantile(0.5), round_mixture())) {
    res <- find_hotspots(s$Y, A = s$A, rounding = rule)
    expect_true(all(res$D <= res$A))                 # subgraph of A
    acc <- res$trace$edges[res$trace$accepted]
    expect_true(all(diff(acc) < 0))                  # strictly shrinking
    expect_false(is.na(res$stop_reason))
  }
})

test_that("CSR gives no structure most of the time", {
  set.seed(4)
  stops <- replicate(20, {
    s <- sim_thomas(k = 10, n_parents = 0, n_background = 2000)
    find_hotspots(s$Y, A = s$A)$stop_reason
  })
  expect_gte(mean(stops == "no residual spatial structure (step 0)"), 0.8)
})

test_that("rank_cells returns a permutation for every method", {
  set.seed(5)
  s <- sim_thomas(k = 12, n_parents = 2, radius = 120, sons = 800)
  res <- find_hotspots(s$Y, A = s$A)
  for (m in c("R0", "R1", "R2", "R3", "B0"))
    expect_equal(sort(rank_cells(res, m)), seq_len(144))
})
