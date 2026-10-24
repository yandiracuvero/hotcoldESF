# Full simulation study ------------------------------------------------------------------
#
# Proposed method (ESF): log(Y + 1); covariate adjustment; Moran gate (5%); ESF fit with
# smooth eigenvectors (eigenvalue >= 0.25 x the largest), ridge, penalty by 10-fold CV;
# permutation band (hot above the 95% quantile, cold below the 5%); ranking by largest
# excess (hot cells first, then f-hat).
#
# Benchmarks (all on the residuals of log(Y + 1) on Z, queen contiguity):
#   Gi*                 no gate; hot z > 1.645, cold z < -1.645; ranking by z
#   Gi* + gate          same, labels only if the Moran test is significant
#   Gi* + gate + FDR    same, Benjamini-Hochberg (5%) for hot and for cold
#   LISA + gate         hot High-High, cold Low-Low (p < 0.05, 999 permutations);
#                       ranking High-High first, then local I
#   Top residuals, Top counts   (ranking only)
#
# Scenarios (window 1000 x 1000; all with 2000 random points and a covariate area of
# 3000 points that the model receives, except U):
#   B  3 compact clusters, radius 100              weak 50 / strong 200 points each
#   L  2 large diffuse clusters, radius 200        weak 200 / strong 600
#   E  3 elongated clusters (ellipses 300 x 33)    weak 100 / strong 300
#   M  mixed sizes, radii 50 / 100 / 200           same intensity, weak / strong
#   S  8 small clusters, radius 40                 weak 15 / strong 40
#   K  B (200 points) plus a cold area, radius 150, keeping 50% (weak) / 20% (strong)
#      of its points
#   U  B data (weak 50 / strong 200) but the model does NOT receive the covariate
#   C  covariate area only (no clusters);  N  3000 random points (no clusters)
# Grids: 20 x 20 (cell side 50; 100 replicates) and 40 x 40 (cell side 25; 50 replicates);
#        999 permutations for the band on both grids, as in the application.
#
# Run from the repo root: Rscript simulation/17_full_study.R   (resumable; ~3-4 h on 2 cores)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else library(hotspots)
  library(parallel)
})
CORES <- as.integer(Sys.getenv("CORES", detectCores()))
OUT   <- "simulation/results/17"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NREP  <- c("20" = as.integer(Sys.getenv("NREP20", 100)), "40" = as.integer(Sys.getenv("NREP40", 50)))
NPERM <- c("20" = 999, "40" = 999)

scen <- list(
  B = list(n = 3, radius = 100, aspect = 1, weak = 50, strong = 200, min_dist = 300),
  L = list(n = 2, radius = 200, aspect = 1, weak = 200, strong = 600, min_dist = 400),
  E = list(n = 3, radius = 100, aspect = 9, weak = 100, strong = 300, min_dist = 350),
  M = list(n = 3, radius = c(50, 100, 200), aspect = 1,
           weak = c(24, 94, 377), strong = c(50, 200, 800), min_dist = 350),
  S = list(n = 8, radius = 40, aspect = 1, weak = 15, strong = 40, min_dist = 150),
  K = list(n = 3, radius = 100, aspect = 1, weak = 200, strong = 200, min_dist = 300,
           cold_keep = c(weak = 0.5, strong = 0.2)),
  U = list(n = 3, radius = 100, aspect = 1, weak = 50, strong = 200, min_dist = 300, no_cov = TRUE),
  C = list(n = 0), N = list(n = 0))
settings <- rbind(
  expand.grid(scn = c("B", "L", "E", "M", "S", "K", "U"), strength = c("weak", "strong"),
              stringsAsFactors = FALSE),
  data.frame(scn = c("C", "N"), strength = "none"))
jobs <- do.call(rbind, lapply(c(20, 40), function(k) do.call(rbind, lapply(seq_len(nrow(settings)), function(i)
  data.frame(k = k, scn = settings$scn[i], strength = settings$strength[i], rep = seq_len(NREP[as.character(k)]))))))
jobs$file <- sprintf("%s/k%02d_%s_%s_r%03d.rds", OUT, jobs$k, jobs$scn, jobs$strength, jobs$rep)
jobs <- jobs[order(-jobs$k, jobs$rep), ]                  # long jobs first; replicates interleaved

simulate <- function(k, scn, strength, seed) {
  if (scn == "C") { set.seed(seed); return(sim_thomas(k, 0, 100, 0, n_background = 2000, n_cov = 3000)) }
  if (scn == "N") { set.seed(seed); return(sim_thomas(k, 0, 100, 0, n_background = 3000)) }
  p <- scen[[scn]]
  for (t in 0:20) {
    set.seed(seed + 1e6 * t)
    s <- tryCatch(sim_thomas(k, p$n, p$radius, p[[strength]], n_background = 2000, n_cov = 3000,
                             min_dist = p$min_dist, aspect = p$aspect,
                             n_cold = if (is.null(p$cold_keep)) 0 else 1,
                             cold_keep = if (is.null(p$cold_keep)) 0.2 else p$cold_keep[[strength]]),
                  error = function(e) NULL)
    if (!is.null(s)) { if (isTRUE(p$no_cov)) s$Z_hidden <- s$Z; return(s) }
  }
  stop("simulation failed")
}

run_job <- function(j) {
  seed <- 20000 + 100000 * (j$k == 40) + 1000 * match(j$scn, names(scen)) +
          500 * (j$strength == "strong") + j$rep
  s <- simulate(j$k, j$scn, j$strength, seed)
  Z <- if (j$scn == "U") NULL else s$Z
  n <- j$k^2; m <- round(0.05 * n)
  cl <- lapply(seq_along(s$clusters), function(i) union(s$clusters[[i]], s$parent_cells[i]))
  foot <- sort(unique(unlist(cl))); cold_area <- s$cold_cells; cov_area <- s$cov_cells
  Zs <- if (is.null(Z)) NULL else scale(Z); X0 <- design0(n, Zs)
  res <- lm.fit(X0, as.numeric(scale(log1p(s$Y))))$residuals
  gate <- moran_resid(res, s$A, X0); gate_ok <- isTRUE(gate$p <= 0.05)

  share <- function(a, b) if (length(a)) mean(a %in% b) else NA
  row <- function(method, hot = NULL, cold = NULL, ord = NULL, secs = NA) {
    top <- if (is.null(ord)) NULL else ord[1:m]
    data.frame(k = j$k, scenario = j$scn, strength = j$strength, rep = j$rep, method = method,
      n_foot = length(foot), n_cold_area = length(cold_area),
      n_hot = if (is.null(hot)) NA else length(hot),
      hot_precision = if (is.null(hot)) NA else share(hot, foot),
      hot_sensitivity = if (is.null(hot) || !length(foot)) NA else mean(foot %in% hot),
      clusters_hot = if (is.null(hot) || !length(cl)) NA else mean(sapply(cl, function(c) any(c %in% hot))),
      hot_in_cov = if (is.null(hot)) NA else share(hot, cov_area),
      n_cold = if (is.null(cold)) NA else length(cold),
      cold_precision = if (is.null(cold)) NA else share(cold, cold_area),
      cold_sensitivity = if (is.null(cold) || !length(cold_area)) NA else mean(cold_area %in% cold),
      top_precision = if (is.null(top) || !length(foot)) NA else mean(top %in% foot),
      clusters_top = if (is.null(top) || !length(cl)) NA else mean(sapply(cl, function(c) any(c %in% top))),
      top_in_cov = if (is.null(top)) NA else share(top, cov_area),
      found_by_cluster = if (is.null(top) || !length(cl)) NA else
        paste(as.integer(sapply(cl, function(c) any(c %in% top))), collapse = ","),
      gate_z = gate$z, gate_ok = gate_ok, secs = secs)
  }

  # ESF (proposed)
  set.seed(seed + 7)
  t0 <- Sys.time()
  h <- hot_cold_map(s$Y, Z, s$A, nperm = NPERM[as.character(j$k)], penalty = "cv_min",
                    mc_min = 0.25, fit_type = "ridge")
  esf_secs <- as.numeric(Sys.time() - t0, units = "secs")
  out <- list(row("ESF", which(h$labels == "hot"), which(h$labels == "cold"), rank_cells(h, "R2"), esf_secs))

  # Gi*
  gi <- local_gi_star(res, s$A)
  ph <- pnorm(gi, lower.tail = FALSE); pc <- pnorm(gi)
  g_hot <- which(gi > qnorm(0.95)); g_cold <- which(gi < qnorm(0.05))
  f_hot <- which(p.adjust(ph, "BH") <= 0.05); f_cold <- which(p.adjust(pc, "BH") <= 0.05)
  none <- integer(0)
  out <- c(out, list(
    row("Gi*", g_hot, g_cold, order(-gi)),
    row("Gi* + gate", if (gate_ok) g_hot else none, if (gate_ok) g_cold else none, order(-gi)),
    row("Gi* + gate + FDR", if (gate_ok) f_hot else none, if (gate_ok) f_cold else none, order(-gi))))

  # LISA
  set.seed(seed + 9)
  li <- local_moran(res, s$A)
  l_hot <- which(li$quad == "HH" & li$p < 0.05); l_cold <- which(li$quad == "LL" & li$p < 0.05)
  out <- c(out, list(
    row("LISA + gate", if (gate_ok) l_hot else none, if (gate_ok) l_cold else none,
        order(-(li$quad == "HH"), -li$I)),
    row("Top residuals", ord = order(-res)),
    row("Top counts", ord = order(-s$Y))))
  do.call(rbind, out)
}

if (sys.nframe() == 0) {
  todo <- which(!file.exists(jobs$file))
  message(length(todo), " of ", nrow(jobs), " datasets to run on ", CORES, " cores")
  t0 <- Sys.time()
  invisible(mclapply(todo, function(i) {
    r <- tryCatch(run_job(jobs[i, ]), error = function(e) {
      message("FAILED ", basename(jobs$file[i]), ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) saveRDS(r, jobs$file[i])
    if (i %% 25 == 0) message(format(Sys.time(), "%H:%M:%S"), "  done ", basename(jobs$file[i]))
  }, mc.cores = CORES, mc.preschedule = FALSE))
  done <- jobs$file[file.exists(jobs$file)]
  res  <- do.call(rbind, lapply(done, readRDS))
  write.csv(res, "simulation/results/17_full_study.csv", row.names = FALSE)
  message("finished ", length(done), "/", nrow(jobs), " in ",
          round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
}
