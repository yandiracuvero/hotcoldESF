test_that("Gi* matches the textbook formula on a small grid", {
  set.seed(1); A <- grid_adjacency(5); x <- rpois(25, 3)
  W <- A + diag(25); n <- 25; xb <- mean(x); s <- sqrt(mean(x^2) - xb^2)
  i <- 13
  z <- (sum(W[i, ] * x) - xb * sum(W[i, ])) /
       (s * sqrt((n * sum(W[i, ]^2) - sum(W[i, ])^2) / (n - 1)))
  expect_equal(local_gi_star(x, A)[i], z)
})

test_that("LISA flags a planted high-high block", {
  set.seed(2); A <- grid_adjacency(10); x <- rnorm(100)
  blk <- c(44, 45, 46, 54, 55, 56, 64, 65, 66); x[blk] <- x[blk] + 4
  lm <- local_moran(x, A, nperm = 499)
  expect_true(all(lm$quad[55] == "HH"), lm$p[55] < 0.05)
})

test_that("hot_cold_map finds planted clusters and few hot cells under CSR", {
  set.seed(3)
  s <- sim_thomas(k = 20, n_parents = 3, radius = 100, sons = 500, n_background = 2000)
  h <- hot_cold_map(s$Y, NULL, s$A, nperm = 50)
  hot <- which(h$labels == "hot")
  expect_gt(mean(hot %in% s$footprint), 0.7)
  null_hot <- replicate(5, {
    s0 <- sim_thomas(k = 20, n_parents = 0, n_background = 3000)
    sum(hot_cold_map(s0$Y, NULL, s0$A, nperm = 50)$labels == "hot")
  })
  expect_lt(mean(null_hot), 40)            # pointwise band: about 5% of 400 at most
  expect_equal(sort(rank_cells(h, "R2")), seq_len(400))
})

test_that("permutation gate rejects rarely under exchangeable noise and detects a block", {
  set.seed(11)
  A <- grid_adjacency(15); X <- matrix(1, 225, 1)
  p <- replicate(40, moran_perm(rnorm(225), A, X, nperm = 199)$p)
  expect_lt(mean(p <= 0.05), 0.2)
  e <- rnorm(225); e[1:45] <- e[1:45] + 2
  expect_lt(moran_perm(e - mean(e), A, X, nperm = 199)$p, 0.05)
})

test_that("covariate_basis gives linear or spline columns", {
  Z <- cbind(a = 1:50, b = rnorm(50))
  expect_equal(ncol(covariate_basis(Z)), 2)
  expect_equal(ncol(covariate_basis(Z, 4)), 8)
  expect_null(covariate_basis(NULL))
})
