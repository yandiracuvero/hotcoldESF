# Tests 3 and 4: cluster shapes and sizes that should favour ESF or Gi* -------------------
#
# Gi* looks at a fixed window (the cell and its 8 neighbours). ESF with smooth
# eigenvectors fits patterns at many scales, but cannot follow very small details.
#
# Test 3 (expected to favour ESF)
#   L  large diffuse clusters: 2 clusters, radius 200
#   E  elongated clusters: 3 ellipses with the area of a disc of radius 100,
#      semi-axes 300 x 33, random angle
#   M  mixed sizes: radii 50, 100 and 200, same intensity (points proportional to area)
# Test 4 (expected to favour Gi*)
#   S  many small clusters: 8 clusters, radius 40 (about one 50-unit cell, or 1.6 cells
#      on the 40 x 40 grid)
#
# Each scenario at two strengths (weak, strong); every dataset also has 2000 random
# points and a covariate area of 3000 points, as in scenario B. Window 1000 x 1000.
# 20 x 20 grid: 20 replicates, 200 permutations; 40 x 40: 10 replicates, 100 permutations.
#
# Methods: ESF ridge (smooth eigenvectors, penalty by CV; the proposed version),
# ESF Lasso (smooth, CV), Gi* + Moran gate, Gi* + Moran gate + FDR, top residuals.
# Measures: precision of the top 5% of cells; clusters found (share of clusters with at
# least one cell in the top 5%); label precision, sensitivity, clusters with a hot cell.
#
# Run from the repo root: Rscript simulation/15_cluster_shapes.R   (resumable)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else
    if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(hotcoldESF)
  library(parallel)
})
CORES <- as.integer(Sys.getenv("CORES", detectCores()))
OUT   <- "simulation/results/15"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

scen <- list(
  L = list(n = 2, radius = 200, aspect = 1, weak = 200, strong = 600, min_dist = 400),
  E = list(n = 3, radius = 100, aspect = 9, weak = 100, strong = 300, min_dist = 350),
  M = list(n = 3, radius = c(50, 100, 200), aspect = 1,
           weak = c(24, 94, 377), strong = c(50, 200, 800), min_dist = 350),
  S = list(n = 8, radius = 40, aspect = 1, weak = 15, strong = 40, min_dist = 150))
cfg <- rbind(
  expand.grid(k = 20, scn = names(scen), strength = c("weak", "strong"), nrep = 20, nperm = 200,
              stringsAsFactors = FALSE),
  expand.grid(k = 40, scn = names(scen), strength = c("weak", "strong"), nrep = 10, nperm = 100,
              stringsAsFactors = FALSE))
jobs <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i)
  data.frame(cfg[i, c("k", "scn", "strength", "nperm")], rep = seq_len(cfg$nrep[i]), row.names = NULL)))
jobs$file <- sprintf("%s/k%02d_%s_%s_r%02d.rds", OUT, jobs$k, jobs$scn, jobs$strength, jobs$rep)
jobs <- jobs[order(-jobs$k), ]

simulate <- function(k, scn, strength, seed) {
  p <- scen[[scn]]
  for (t in 0:20) {                                  # rare: no room for the covariate area
    set.seed(seed + 100000 * t)
    s <- tryCatch(sim_thomas(k, p$n, p$radius, p[[strength]], n_background = 2000,
                             n_cov = 3000, min_dist = p$min_dist, aspect = p$aspect),
                  error = function(e) NULL)
    if (!is.null(s)) return(s)
  }
  stop("simulation failed")
}

run_job <- function(j) {
  seed <- 9000 + 1000 * j$k + 100 * match(j$scn, names(scen)) +
          50 * (j$strength == "strong") + j$rep
  s <- simulate(j$k, j$scn, j$strength, seed)
  n <- j$k^2; m <- round(0.05 * n)
  # each cluster = cells whose centre lies inside it, plus the cell of its centre
  # (small clusters can fall between cell centres)
  cl <- lapply(seq_along(s$clusters), function(i) union(s$clusters[[i]], s$parent_cells[i]))
  foot <- sort(unique(unlist(cl)))
  Zs <- scale(s$Z); X0 <- design0(n, Zs)
  res <- lm.fit(X0, as.numeric(scale(log1p(s$Y))))$residuals
  gate_ok <- isTRUE(moran_resid(res, s$A, X0)$p <= 0.05)
  row <- function(method, hot, ord) {
    top <- if (is.null(ord)) NULL else ord[1:m]
    data.frame(k = j$k, scenario = j$scn, strength = j$strength, rep = j$rep, method = method,
      n_foot = length(foot), n_hot = length(hot),
      label_precision = if (length(hot)) mean(hot %in% foot) else NA,
      sensitivity = mean(foot %in% hot),
      clusters_hot = mean(sapply(cl, function(c) any(c %in% hot))),
      prec_top = if (is.null(top)) NA else mean(top %in% foot),
      clusters_top = if (is.null(top)) NA else mean(sapply(cl, function(c) any(c %in% top))),
      found_by_cluster = if (is.null(top)) NA else
        paste(as.integer(sapply(cl, function(c) any(c %in% top))), collapse = ","),
      gate_ok = gate_ok)
  }
  gi <- local_gi_star(res, s$A)
  hot_gi  <- if (gate_ok) which(gi > qnorm(0.95)) else integer(0)
  hot_fdr <- if (gate_ok) which(p.adjust(pnorm(gi, lower.tail = FALSE), "BH") <= 0.05) else integer(0)
  out <- list(row("Gi* + gate", hot_gi, order(-gi)),
              row("Gi* + gate + FDR", hot_fdr, order(-gi)),
              row("Top residuals", integer(0), order(-res)))
  for (ft in c("ridge", "lasso")) {
    set.seed(seed + if (ft == "ridge") 1 else 2)
    h <- hot_cold_map(s$Y, s$Z, s$A, nperm = j$nperm, penalty = "cv_min", mc_min = 0.25,
                      fit_type = ft)
    out[[length(out) + 1]] <- row(paste0("ESF ", ft), which(h$labels == "hot"), rank_cells(h, "R2"))
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
write.csv(res, "simulation/results/15_cluster_shapes.csv", row.names = FALSE)
message("finished ", length(done), "/", nrow(jobs), " in ",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
