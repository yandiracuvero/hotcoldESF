# Figure for simulation/15_cluster_shapes.R ----------------------------------------------
# Run from the repo root: Rscript simulation/16_cluster_shapes_figure.R
# Writes simulation/figures/fig14_cluster_shapes.png
# Panels: rows = weak / strong clusters, columns = scenarios L, E, M (test 3) and S (test 4).
# Each panel: methods on the x axis, one line per measure (mean over the 20 x 20 and
# 40 x 40 grids, each grid averaged over its replicates).

x <- read.csv("simulation/results/15_cluster_shapes.csv", stringsAsFactors = FALSE)
mean_na <- function(z) if (all(is.na(z))) NA else mean(z, na.rm = TRUE)
DARK <- "#3C4043"; GREY <- "#9AA0A6"; LIGHT <- "#E8EAED"
meas <- c(prec_top = "Precision of the top 5% of cells", clusters_top = "Clusters with a cell in the top 5%",
          label_precision = "Hot cells inside a true cluster", sensitivity = "Cluster cells labelled hot")
COL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100"); PCH <- c(16, 15, 17, 18)
meth <- c("ESF ridge", "ESF lasso", "Gi* + gate", "Gi* + gate + FDR", "Top residuals")
mlab <- c("ESF ridge", "ESF Lasso", "Gi* + gate", "Gi* + gate\n+ FDR", "Top\nresiduals")
scn  <- c(L = "L: large diffuse clusters", E = "E: elongated clusters",
          M = "M: mixed sizes (radius 50, 100, 200)", S = "S: many small clusters")

val <- function(sc, st, mm, col) mean_na(sapply(c(20, 40), function(k)
  mean_na(x[[col]][x$scenario == sc & x$strength == st & x$method == mm & x$k == k])))

png("simulation/figures/fig14_cluster_shapes.png", 2600, 1500, res = 150, pointsize = 10)
layout(rbind(1:4, 5:8, 9), heights = c(1, 1, 0.22))
par(family = "sans", col.axis = DARK, col.lab = DARK, fg = GREY, las = 1, tcl = -0.25,
    mgp = c(2.4, 0.6, 0), oma = c(0, 0, 4, 0))
for (st in c("weak", "strong")) for (sc in names(scn)) {
  par(mar = c(4.5, 4, 3.2, 1))
  plot(NA, xlim = c(0.7, length(meth) + 0.3), ylim = c(0, 1), xaxt = "n", bty = "n",
       xlab = "", ylab = if (sc == "L") "Share" else "")
  abline(h = seq(0, 1, 0.25), col = LIGHT)
  axis(1, at = seq_along(meth), labels = FALSE)
  text(seq_along(meth), -0.12, mlab, xpd = NA, cex = 0.8, col = DARK, adj = c(0.5, 1))
  for (i in seq_along(meas)) {
    y <- sapply(meth, val, sc = sc, st = st, col = names(meas)[i])
    ok <- which(!is.na(y))
    lines(ok, y[ok], col = COL[i], lwd = 2)
    points(ok, y[ok], col = COL[i], pch = PCH[i], cex = 1.3)
  }
  mtext(scn[sc], side = 3, line = 1.4, adj = 0, font = 2, col = DARK, cex = 0.85)
  mtext(paste(if (st == "weak") "Weak" else "Strong", "clusters;",
              if (sc == "S") "test 4 (expected to favour Gi*)" else "test 3 (expected to favour ESF)"),
        side = 3, line = 0.4, adj = 0, col = GREY, cex = 0.75)
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", meas, col = COL, pch = PCH, lwd = 2, horiz = TRUE, bty = "n", cex = 1, text.col = DARK)
mtext("Cluster shapes and sizes: ESF (smooth eigenvectors, penalty by CV) vs Gi* on residuals with the Moran gate",
      outer = TRUE, side = 3, line = 2.2, adj = 0.01, font = 2, col = DARK, cex = 1.05)
mtext("Mean of the 20 x 20 (20 replicates) and 40 x 40 (10 replicates) grids. All datasets also have 2000 random points and a covariate area of 3000 points. Higher is better.",
      outer = TRUE, side = 3, line = 0.9, adj = 0.01, col = GREY, cex = 0.8)
dev.off()
message("written: fig14_cluster_shapes.png")
