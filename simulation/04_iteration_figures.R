# Figures: what the iterative reconstruction does, and why it is not retained -------
#
# Reads the results of the three experiments (no re-simulation, except figure S2,
# which runs one example in a few seconds) and draws:
#   fig1_cold_zones.png          published ranking picks cold zones (experiment 1)
#   fig2_iteration_no_effect.png iterating does not change precision (experiment 2)
#   fig3_latent_graph.png        reconstruction does not find barriers (experiment 3)
#   figS1_where_esf_helps.png    precision by grid and cluster strength (experiment 2)
#   figS2_example_run.png        one simulated map, step by step
#
# Run from the repo root, in RStudio or the Terminal:
#     Rscript simulation/04_iteration_figures.R
# Figures are written to simulation/figures/. In RStudio, set SHOW <- TRUE to also
# draw them in the Plots pane.

SHOW <- interactive()
suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else library(hotspots)
})
dir.create("simulation/figures", showWarnings = FALSE, recursive = TRUE)

# ---- palette: grey first, one accent -----------------------------------------
GREY   <- "#9AA0A6"; DARK <- "#3C4043"; LIGHT <- "#DADCE0"
ACCENT <- "#1A73E8"; ACC2 <- "#8AB4F8"; WARM <- "#D93025"
mean_na <- function(z) mean(z, na.rm = TRUE)

fig <- function(file, width, height, expr) {
  png(file.path("simulation/figures", file), width = width, height = height,
      res = 150, pointsize = 10)
  on.exit(dev.off())
  par(family = "sans", col.axis = DARK, col.lab = DARK, fg = GREY, las = 1,
      mgp = c(2.4, 0.6, 0), tcl = -0.25)
  eval.parent(substitute(expr))
  if (SHOW) { dev.off(); eval.parent(substitute(expr)); png(tempfile()) }
}
title_left <- function(main, sub = NULL) {
  mtext(main, side = 3, line = if (is.null(sub)) 1 else 2.2, adj = 0, font = 2,
        col = DARK, cex = 1.05)
  if (!is.null(sub)) mtext(sub, side = 3, line = 1, adj = 0, col = GREY, cex = 0.85)
}

x1 <- read.csv("simulation/results/01_small_test.csv", stringsAsFactors = FALSE)
x2 <- read.csv("simulation/results/02_signal_test.csv", stringsAsFactors = FALSE)
x3 <- read.csv("simulation/results/03_latent_graph.csv", stringsAsFactors = FALSE)

# ==== Figure 1: the published ranking picks cold zones (experiment 1) ==========
key <- data.frame(
  variant = c("all/lambda/q75", "baseline", "all/lambda/q75", "all/lambda/q75", "baseline"),
  ranking = c("R0", "B0", "R3", "R2", "B1"),
  label   = c("Published (R0)", "Top counts", "Hot zones first (R3)",
              "Hot cells by f-hat (R2)", "Top residuals on Z"),
  stringsAsFactors = FALSE)
get1 <- function(sc, v, r, col) mean_na(x1[x1$scenario == sc & x1$variant == v & x1$ranking == r, col])
prec <- sapply(c("S1", "S2"), function(sc) mapply(get1, sc, key$variant, key$ranking, "precision"))
cold <- mapply(get1, "S2", key$variant, key$ranking, "cold")

fig("fig1_cold_zones.png", 1500, 700, {
  layout(matrix(1:2, 1), widths = c(1.6, 1))
  par(mar = c(4, 11, 4.5, 1))
  y <- seq_len(nrow(key)); off <- 0.18
  plot(NA, xlim = c(0, 1), ylim = c(0.5, nrow(key) + 0.5), yaxt = "n", xlab = "Precision (share of top-20 cells inside a true cluster)", ylab = "", bty = "n")
  abline(v = seq(0, 1, 0.25), col = LIGHT)
  axis(2, at = y, labels = key$label, tick = FALSE, col.axis = DARK)
  is_pub <- key$ranking == "R0"
  segments(0, y + off, prec[, "S1"], y + off, lwd = 7, lend = 1, col = ifelse(is_pub, ACC2, LIGHT))
  segments(0, y - off, prec[, "S2"], y - off, lwd = 7, lend = 1, col = ifelse(is_pub, ACCENT, GREY))
  text(prec[, "S1"], y + off, sprintf("%.2f", prec[, "S1"]), pos = 4, cex = 0.8, col = DARK, xpd = TRUE)
  text(prec[, "S2"], y - off, sprintf("%.2f", prec[, "S2"]), pos = 4, cex = 0.8, col = DARK, xpd = TRUE)
  legend("bottomright", c("clusters only", "clusters + covariate"), lwd = 7, col = c(LIGHT, GREY),
         bty = "n", cex = 0.8, inset = c(0, -0.02))
  title_left("The published ranking has the lowest precision",
             "Experiment 1: 20 replicates per scenario, 400 cells, top 20 cells")
  par(mar = c(4, 10, 4.5, 3))
  cc <- cold[!is.na(cold)]; lab <- key$label[!is.na(cold)]
  bp <- barplot(100 * cc, horiz = TRUE, xlim = c(0, 100), col = ifelse(grepl("Published", lab), ACCENT, GREY),
                border = NA, names.arg = lab, cex.names = 0.8, xlab = "Cold cells among the top 20 (%)", axes = FALSE)
  axis(1, at = seq(0, 100, 25))
  text(100 * cc, bp, sprintf("%.0f%%", 100 * cc), pos = 4, cex = 0.8, col = DARK, xpd = TRUE)
  title_left("Half of its picks are cold zones", "Scenario with covariate")
})

# ==== Figure 2: iterating does not change precision (experiment 2) =============
s2 <- x2[x2$scenario == "S2" & x2$ranking %in% c("R0", "R2", "R3"), ]
s2$base <- sub("_(none|q75)$", "", s2$method)
s2$prune <- sub(".*_", "", s2$method)
agg <- aggregate(precision ~ base + ranking + k + sons + prune, s2, mean_na)
w <- reshape(agg, idvar = c("base", "ranking", "k", "sons"), timevar = "prune", direction = "wide")
w <- w[complete.cases(w), ]
d <- merge(s2[s2$prune == "q75", c("base", "ranking", "k", "sons", "rep", "precision", "edges_kept")],
           s2[s2$prune == "none", c("base", "ranking", "k", "sons", "rep", "precision")],
           by = c("base", "ranking", "k", "sons", "rep"), suffixes = c("_q75", "_none"))
d$diff <- d$precision_q75 - d$precision_none
maxdiff <- max(abs(w$precision.q75 - w$precision.none))

fig("fig2_iteration_no_effect.png", 1500, 700, {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 4.5, 1))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Precision, single ESF fit (no iteration)",
       ylab = "Precision, iterated reconstruction", bty = "n", asp = 1)
  abline(h = seq(0, 1, 0.25), v = seq(0, 1, 0.25), col = LIGHT)
  abline(0, 1, col = DARK, lty = 2)
  text(0.97, 0.9, "no change", srt = 45, col = DARK, cex = 0.8, adj = 1)
  points(w$precision.none, w$precision.q75, pch = 21, bg = adjustcolor(ACCENT, 0.5), col = "white", cex = 1.4)
  title_left("No systematic gain: settings scatter on both sides",
             sprintf("%d settings (grid x strength x scale x a x ranking), mean of 10 runs each",
                     nrow(w)))
  par(mar = c(4.2, 9, 4.5, 2))
  dp <- d[d$edges_kept < 0.999, ]
  rk <- c(R0 = "Published (R0)", R2 = "Hot cells by f-hat (R2)", R3 = "Hot zones first (R3)")
  st <- t(sapply(names(rk), function(r) { z <- dp$diff[dp$ranking == r]
    c(m = mean(z), lo = mean(z) - 1.96 * sd(z) / sqrt(length(z)), hi = mean(z) + 1.96 * sd(z) / sqrt(length(z)), n = length(z)) }))
  yy <- rev(seq_len(nrow(st)))
  plot(NA, xlim = c(-0.15, 0.15), ylim = c(0.5, nrow(st) + 0.5), yaxt = "n", bty = "n",
       xlab = "Precision change from iterating (mean, 95% CI)", ylab = "")
  abline(v = seq(-0.15, 0.15, 0.05), col = LIGHT); abline(v = 0, col = DARK)
  axis(2, at = yy, labels = rk, tick = FALSE, cex.axis = 0.85)
  segments(st[, "lo"], yy, st[, "hi"], yy, lwd = 2, col = GREY)
  points(st[, "m"], yy, pch = 16, cex = 1.4, col = ifelse(names(rk) == "R0", ACCENT, DARK))
  text(st[, "m"], yy + 0.22, sprintf("%+.2f", st[, "m"]), cex = 0.8, col = DARK)
  text(-0.15, 0.55, "worse", adj = 0, col = GREY, cex = 0.8); text(0.15, 0.55, "better", adj = 1, col = GREY, cex = 0.8)
  title_left("Iterating does not help the rankings that work",
             sprintf("Only the %d of %d runs where iterating removed edges (%.0f%%)",
                     nrow(dp), nrow(d), 100 * nrow(dp) / nrow(d)))
})

# ==== Figure 3: the reconstruction does not find barriers (experiment 3) =======
b <- unique(x3[x3$type == "barriers", c("rho", "rep", "AUC_B_lambda", "AUC_B_gamma2", "AUC_dres")])
au <- aggregate(cbind(AUC_B_lambda, AUC_B_gamma2, AUC_dres) ~ rho, b, mean_na)
ctrl <- aggregate(edges_kept ~ variant, x3[x3$type == "none", ], mean_na)

fig("fig3_latent_graph.png", 1500, 700, {
  layout(matrix(1:2, 1), widths = c(1.3, 1))
  par(mar = c(4.2, 4.2, 4.5, 8))
  plot(NA, xlim = c(0.45, 1), ylim = c(0.4, 1), xlab = expression("Spatial dependence " * rho),
       ylab = "AUC: real edges vs barrier edges", bty = "n", xaxt = "n")
  axis(1, at = au$rho)
  abline(h = seq(0.4, 1, 0.1), col = LIGHT)
  abline(h = 0.5, col = DARK, lty = 2); text(0.46, 0.52, "chance", adj = 0, cex = 0.8, col = DARK)
  lines(au$rho, au$AUC_B_lambda, col = GREY, lwd = 2, type = "o", pch = 16)
  lines(au$rho, au$AUC_B_gamma2, col = GREY, lwd = 2, type = "o", pch = 1, lty = 3)
  lines(au$rho, au$AUC_dres, col = ACCENT, lwd = 2.5, type = "o", pch = 16)
  last <- nrow(au)
  text(1, au$AUC_dres[last], "neighbour difference\nof residuals", pos = 4, col = ACCENT, cex = 0.8, xpd = TRUE)
  text(1, au$AUC_B_gamma2[last] + 0.03, "reconstruction b_ij\n(gamma^2 weights)", pos = 4, col = GREY, cex = 0.8, xpd = TRUE)
  text(1, au$AUC_B_lambda[last] - 0.03, "reconstruction b_ij\n(eigenvalue weights)", pos = 4, col = GREY, cex = 0.8, xpd = TRUE)
  title_left("Reconstructed weights do not detect barriers",
             "SAR on a grid with a river and a highway; 20 replicates per rho")
  par(mar = c(4.2, 9, 4.5, 2))
  ctrl <- ctrl[order(ctrl$edges_kept), ]
  bp <- barplot(100 * ctrl$edges_kept, horiz = TRUE, xlim = c(0, 118), col = GREY, border = NA,
                names.arg = ctrl$variant, cex.names = 0.75, xlab = "Edges kept (%)", axes = FALSE)
  axis(1, at = seq(0, 100, 25))
  abline(v = 100, col = ACCENT, lty = 2)
  text(100, max(bp) + 0.9, "should be 100%", col = ACCENT, cex = 0.8, xpd = TRUE)
  text(103, bp, sprintf("%.0f%%", 100 * ctrl$edges_kept), adj = 0, cex = 0.8, col = DARK, xpd = TRUE)
  title_left("With no barriers, edges are still cut", "Control: data generated on the full grid")
})

# ==== Figure S1: where ESF helps (experiment 2) =================================
pick <- data.frame(v = c("B0|-", "B1_log|-", "std_a2_q75|R0", "log_a2_none|R2", "std_a2_none|R3"),
                   lab = c("Top counts", "Top residuals, log scale", "Published (R0)",
                           "ESF log, hot cells by f-hat", "ESF, hot zones first (R3)"),
                   col = c(LIGHT, DARK, WARM, GREY, ACCENT), lwd = c(2, 2.5, 2, 2, 3),
                   stringsAsFactors = FALSE)
x2$v <- paste(x2$method, x2$ranking, sep = "|")
cfg <- expand.grid(sons = c(50, 200, 500), k = c(10, 20, 40))
fig("figS1_where_esf_helps.png", 1500, 750, {
  par(mar = c(5.5, 4.2, 4.5, 13))
  xx <- seq_len(nrow(cfg)) + rep(c(0, 1, 2), each = 3)
  plot(NA, xlim = range(xx), ylim = c(0, 1), xaxt = "n", xlab = "", ylab = "Precision (top 5% of cells)", bty = "n")
  abline(h = seq(0, 1, 0.25), col = LIGHT)
  axis(1, at = xx, labels = cfg$sons, tick = FALSE, line = -0.5)
  mtext(c("10 x 10 grid", "20 x 20 grid", "40 x 40 grid"), side = 1, line = 2, at = tapply(xx, cfg$k, mean), col = DARK)
  mtext("points per cluster", side = 1, line = 3.5, col = GREY, cex = 0.85)
  ends <- numeric(nrow(pick))
  for (i in seq_len(nrow(pick))) {
    yv <- mapply(function(k, s) mean_na(x2$precision[x2$scenario == "S2" & x2$v == pick$v[i] & x2$k == k & x2$sons == s]),
                 cfg$k, cfg$sons)
    for (g in c(10, 20, 40)) { j <- cfg$k == g; lines(xx[j], yv[j], col = pick$col[i], lwd = pick$lwd[i], type = "o", pch = 16) }
    ends[i] <- yv[length(yv)]
  }
  o <- order(ends); lab_y <- ends[o]                     # keep end labels >= 0.06 apart
  for (j in seq_along(lab_y)[-1]) lab_y[j] <- max(lab_y[j], lab_y[j - 1] + 0.06)
  lab_y <- lab_y - max(0, max(lab_y) - 1.02)
  text(max(xx) + 0.3, lab_y, pick$lab[o], adj = 0, col = ifelse(pick$col[o] == LIGHT, GREY, pick$col[o]), cex = 0.8, xpd = TRUE)
  title_left("ESF zones beat adjusted residuals only on the fine grid with weak clusters",
             "Experiment 2, scenario with covariate; 10 replicates per setting")
})

# ==== Figure S2: one example, step by step ======================================
set.seed(2024)
s <- sim_thomas(k = 20, n_parents = 3, radius = 100, sons = 500, n_background = 2000, n_cov = 3000)
res <- find_hotspots(s$Y, s$Z, s$A, mc_min = NULL, rounding = round_quantile(0.75), keep_graphs = TRUE)
cc <- s$centres
draw_grid <- function(main, value = NULL, D = NULL, top = NULL, pal = NULL) {
  par(mar = c(1, 1, 3, 1))
  plot(NA, xlim = c(0, 1000), ylim = c(0, 1000), asp = 1, axes = FALSE, xlab = "", ylab = "")
  h <- s$cell_side
  if (!is.null(value)) rect(cc[, 1] - h / 2, cc[, 2] - h / 2, cc[, 1] + h / 2, cc[, 2] + h / 2, col = pal(value), border = NA)
  if (!is.null(D)) { e <- which(upper.tri(D) & D == 1, arr.ind = TRUE)
    segments(cc[e[, 1], 1], cc[e[, 1], 2], cc[e[, 2], 1], cc[e[, 2], 2], col = DARK, lwd = 0.6) }
  if (!is.null(top)) points(cc[top, 1], cc[top, 2], pch = 22, cex = 1.6, lwd = 2, col = ACCENT)
  th <- seq(0, 2 * pi, length.out = 60)
  for (i in seq_len(nrow(s$parents))) lines(s$parents[i, 1] + 100 * cos(th), s$parents[i, 2] + 100 * sin(th), col = WARM, lwd = 1.5)
  if (!is.null(s$cov_center)) lines(s$cov_center[1] + 150 * cos(th), s$cov_center[2] + 150 * sin(th), col = DARK, lwd = 1.2, lty = 2)
  rect(0, 0, 1000, 1000, border = GREY)
  mtext(main, side = 3, line = 0.8, adj = 0, cex = 0.8, col = DARK, font = 2)
}
seqpal <- function(v) { r <- (v - min(v)) / diff(range(v)); rgb(colorRamp(c("white", "#1A73E8"))(r), maxColorValue = 255) }
divpal <- function(v) { m <- max(abs(v)); r <- (v / m + 1) / 2; rgb(colorRamp(c("#2166AC", "white", "#D93025"))(r), maxColorValue = 255) }
ng <- length(res$graphs); m <- 20
fig("figS2_example_run.png", 1800, 1250, {
  par(mfrow = c(2, 3), oma = c(2, 0, 3, 0))
  draw_grid("(a) Counts per cell", value = s$Y, pal = seqpal)
  draw_grid(sprintf("(b) Iteration 0: full grid, %d edges", sum(res$graphs[[1]]$D) / 2), D = res$graphs[[1]]$D)
  mid <- min(2, ng)
  draw_grid(sprintf("(c) Iteration %d: %d edges", mid - 1, sum(res$graphs[[mid]]$D) / 2), D = res$graphs[[mid]]$D)
  draw_grid(sprintf("(d) Final graph (iteration %d): %d edges", ng - 1, sum(res$D) / 2), D = res$D)
  draw_grid("(e) Fitted spatial component f-hat (red = hot, blue = cold)", value = res$fhat, pal = divpal)
  draw_grid(sprintf("(f) Top %d cells: published R0 (squares)", m), value = res$fhat, pal = divpal,
            top = rank_cells(res, "R0")[1:m])
  mtext("One simulated map, published settings (all eigenvectors, eigenvalue weights, q75, harmonic centrality)",
        outer = TRUE, side = 3, line = 1, adj = 0.02, font = 2, col = DARK)
  mtext("Red circles: true clusters. Dashed circle: covariate-driven area (should not be hot).",
        outer = TRUE, side = 1, line = 0.5, adj = 0.02, col = GREY, cex = 0.85)
})

message("Figures written to simulation/figures/: ",
        paste(list.files("simulation/figures", pattern = "png$"), collapse = ", "))
