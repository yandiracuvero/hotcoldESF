# Example maps and permutation bands for experiment 07 ------------------------------------
#
# Re-creates three datasets of simulation/07_hotcold_benchmarks.R (same seeds):
#   a) B200: highest ESF ranking precision
#   b) B50:  lowest ESF ranking precision (Gi* on residuals does well)
#   c) B50:  highest ESF ranking precision among the weak-cluster datasets
# and draws, for each:
#   fig8x_maps_<scenario>_r<rep>.png   counts, covariate, residuals, f-hat, labels by
#                                      method, top 5% chosen by ESF and by Gi*
#   fig9x_band_<scenario>_r<rep>.png   the permutation band: f-hat of every cell
#                                      against its 5%-95% null band
#
# Run from the repo root: Rscript simulation/09_example_maps.R   (about 1 min)

suppressMessages(if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else
    if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(hotcoldESF))
dir.create("simulation/figures", showWarnings = FALSE)
GREY <- "#9AA0A6"; DARK <- "#3C4043"; LIGHT <- "#DADCE0"; ACCENT <- "#1A73E8"; WARM <- "#D93025"; COLD <- "#2166AC"

# ---- choose the two datasets from the results --------------------------------------------
x <- read.csv("simulation/results/07_hotcold_benchmarks.csv", stringsAsFactors = FALSE)
r <- x[x$kind == "ranking" & x$k == 20, ]
w <- merge(r[r$method == "ESF: largest excess", c("scenario", "rep", "precision", "esf_stop")],
           r[r$method == "Gi* residuals", c("scenario", "rep", "precision")],
           by = c("scenario", "rep"), suffixes = c("_esf", "_gi"))
w <- w[w$esf_stop == "ok", ]
b <- w[w$scenario == "B200", ]; best  <- b[order(-b$precision_esf, -b$precision_gi), ][1, ]
b <- w[w$scenario == "B50", ];  worst <- b[order(b$precision_esf, -b$precision_gi), ][1, ]
weak  <- b[order(-b$precision_esf), ][1, ]   # weak clusters, highest ESF precision
examples <- rbind(cbind(best, tag = "a"), cbind(worst, tag = "b"), cbind(weak, tag = "c"))
print(examples)

simulate <- function(k, scn) switch(substr(scn, 1, 1),
  A = sim_thomas(k, 3, 100, 500, n_background = 2000),
  B = sim_thomas(k, 3, 100, as.numeric(substring(scn, 2)), n_background = 2000, n_cov = 3000),
  C = sim_thomas(k, 0, 100, 0, n_background = 2000, n_cov = 3000),
  N = sim_thomas(k, 0, 100, 0, n_background = 3000))

for (e in seq_len(nrow(examples))) {
  ex <- examples[e, ]; k <- 20; nperm <- 200
  # same seed and call order as run_job() in 07_hotcold_benchmarks.R
  set.seed(7000 + 1000 * k + 10 * match(ex$scenario, c("A", "B50", "B200", "B500", "C", "N")) + ex$rep)
  s <- simulate(k, ex$scenario)
  n <- k^2; m <- round(0.05 * n)
  X0 <- cbind(1, scale(s$Z))
  res <- lm.fit(X0, as.numeric(scale(log1p(s$Y))))$residuals
  h <- hot_cold_map(s$Y, s$Z, s$A, nperm = nperm, penalty = "z", mc_min = NULL, fit_type = "lasso")
  gi <- local_gi_star(res, s$A); li <- local_moran(res, s$A)
  hot_gi <- which(gi > qnorm(0.95)); hot_li <- which(li$quad == "HH" & li$p < 0.05)
  top_esf <- rank_cells(h, "R2")[1:m]; top_gi <- order(-gi)[1:m]
  foot <- s$footprint; cc <- s$centres; hs <- s$cell_side
  prec <- function(top) mean(top %in% foot)
  sub_title <- sprintf("Scenario %s, replicate %d (20 x 20 grid, 3 clusters of %s points, covariate area of 3000 points)",
                       ex$scenario, ex$rep, substring(ex$scenario, 2))

  # ---- map helpers ----------------------------------------------------------------------
  th <- seq(0, 2 * pi, length.out = 60)
  draw <- function(main, fill = NULL, pts = NULL, pch = 22, pcol = ACCENT) {
    par(mar = c(0.5, 0.5, 2.6, 0.5))
    plot(NA, xlim = c(0, 1000), ylim = c(0, 1000), asp = 1, axes = FALSE, xlab = "", ylab = "")
    if (!is.null(fill)) rect(cc[, 1] - hs / 2, cc[, 2] - hs / 2, cc[, 1] + hs / 2, cc[, 2] + hs / 2, col = fill, border = NA)
    if (!is.null(pts) && length(pts)) points(cc[pts, 1], cc[pts, 2], pch = pch, cex = 1.3, lwd = 1.6, col = pcol)
    for (i in seq_len(nrow(s$parents))) lines(s$parents[i, 1] + 100 * cos(th), s$parents[i, 2] + 100 * sin(th), col = WARM, lwd = 1.5)
    lines(s$cov_center[1] + 150 * cos(th), s$cov_center[2] + 150 * sin(th), col = DARK, lwd = 1.2, lty = 2)
    rect(0, 0, 1000, 1000, border = GREY)
    mtext(main, side = 3, line = 0.6, adj = 0, cex = 0.78, col = DARK, font = 2)
  }
  seqpal <- function(v) { r <- (v - min(v)) / diff(range(v)); rgb(colorRamp(c("white", "#5F6368"))(r), maxColorValue = 255) }
  divpal <- function(v) { mx <- max(abs(v)); if (mx == 0) return(rep("white", length(v)))
    rgb(colorRamp(c(COLD, "white", WARM))((v / mx + 1) / 2), maxColorValue = 255) }
  labfill <- function(hot, cold = integer(0)) { f <- rep("white", n); f[hot] <- "#F4A79F"; f[cold] <- "#9EC3E6"; f }

  # ---- figure 8: maps -------------------------------------------------------------------
  png(sprintf("simulation/figures/fig8%s_maps_%s_r%02d.png", ex$tag, ex$scenario, ex$rep), 2000, 1100, res = 150, pointsize = 10)
  par(mfrow = c(2, 4), oma = c(2.2, 0, 3.2, 0), family = "sans")
  draw("(a) Counts per cell", seqpal(log1p(s$Y)))
  draw("(b) Covariate (environment)", seqpal(s$Z[, "env"]))
  draw("(c) Residuals of log counts on Z", divpal(res))
  draw("(d) ESF spatial component f-hat", divpal(h$fhat))
  draw(sprintf("(e) ESF band: %d hot, %d cold", sum(h$labels == "hot"), sum(h$labels == "cold")),
       labfill(which(h$labels == "hot"), which(h$labels == "cold")))
  draw(sprintf("(f) Gi* on residuals: %d hot (z > 1.645)", length(hot_gi)), labfill(hot_gi))
  draw(sprintf("(g) LISA on residuals: %d High-High", length(hot_li)), labfill(hot_li))
  draw(sprintf("(h) Top %d cells: ESF %.2f, Gi* %.2f", m, prec(top_esf), prec(top_gi)), NULL)
  points(cc[top_gi, 1], cc[top_gi, 2], pch = 1, cex = 1.9, lwd = 1.6, col = DARK)
  points(cc[top_esf, 1], cc[top_esf, 2], pch = 22, cex = 1.2, lwd = 1.6, col = ACCENT)
  legend("bottom", c("ESF: largest excess", "Gi* residuals"), pch = c(22, 1), col = c(ACCENT, DARK),
         pt.lwd = 1.6, bty = "n", cex = 0.95, horiz = TRUE, inset = c(0, 0.01))
  mtext(sub_title, outer = TRUE, side = 3, line = 1.6, adj = 0.01, font = 2, col = DARK)
  mtext("Red circles: true clusters. Dashed circle: covariate area (explained by Z, not a hot spot). (e)-(g): red = hot, blue = cold. (h): precision = share of the top cells inside a red circle.",
        outer = TRUE, side = 3, line = 0.4, adj = 0.01, col = GREY, cex = 0.8)
  mtext(sprintf("Moran test of the residuals: z = %.2f, p = %.3f; penalty lambda = %.3f; %d of %d eigenvectors selected",
                h$z0, h$p0, h$lambda, length(h$fit$S), n),
        outer = TRUE, side = 1, line = 0.8, adj = 0.01, col = GREY, cex = 0.8)
  dev.off()

  # ---- figure 9: permutation band -------------------------------------------------------
  png(sprintf("simulation/figures/fig9%s_band_%s_r%02d.png", ex$tag, ex$scenario, ex$rep), 2000, 900, res = 150, pointsize = 10)
  layout(matrix(1:2, 1), widths = c(1.6, 1))
  par(family = "sans", col.axis = DARK, col.lab = DARK, fg = GREY, las = 1, mgp = c(2.4, 0.6, 0), tcl = -0.25,
      oma = c(2.5, 0, 2.6, 0))
  incl <- seq_len(n) %in% foot
  lab_col <- ifelse(h$labels == "hot", WARM, ifelse(h$labels == "cold", COLD, GREY))
  band_panel <- function(idx, main, sub) {
    o <- idx[order(h$fhat[idx])]; xs <- seq_along(o)
    yl <- range(c(h$fhat[o], h$upper[o], h$lower[o]))
    plot(NA, xlim = c(1, length(o)), ylim = yl, bty = "n", xlab = "Cells, ordered by f-hat", ylab = "f-hat")
    abline(h = 0, col = LIGHT)
    polygon(c(xs, rev(xs)), c(h$upper[o], rev(h$lower[o])), col = "#E8EAED", border = NA)
    lines(xs, h$upper[o], col = GREY, lwd = 0.8); lines(xs, h$lower[o], col = GREY, lwd = 0.8)
    points(xs, h$fhat[o], pch = ifelse(incl[o], 17, 16), cex = ifelse(incl[o], 0.9, 0.55), col = lab_col[o])
    mtext(main, side = 3, line = 1.6, adj = 0, font = 2, col = DARK, cex = 0.95)
    mtext(sub, side = 3, line = 0.5, adj = 0, col = GREY, cex = 0.8)
  }
  par(mar = c(4, 4, 3.6, 1))
  band_panel(seq_len(n), "All cells", "Grey area: 5%-95% band of f-hat under no spatial structure (one band per cell)")
  par(xpd = NA)
  legend(x = par("usr")[1], y = par("usr")[3] - 0.15 * diff(par("usr")[3:4]),
         c("hot (above the band)", "neutral", "cold (below the band)", "triangle: cell inside a true cluster"),
         pch = c(16, 16, 16, 17), col = c(WARM, GREY, COLD, DARK), bty = "n", cex = 0.85, horiz = TRUE)
  par(xpd = FALSE)
  band_panel(foot, sprintf("Cells inside a true cluster (%d)", length(foot)),
             sprintf("%d of %d labelled hot", sum(h$labels[foot] == "hot"), length(foot)))
  mtext(paste0("Permutation band, ", sub_title), outer = TRUE, side = 3, line = 0.8, adj = 0.01, font = 2, col = DARK)
  dev.off()
}
message("written: fig8*/fig9* in simulation/figures/")
