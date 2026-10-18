# Figures for simulation/07_hotcold_benchmarks.R ----------------------------------------
# Run from the repo root: Rscript simulation/08_benchmark_figures.R
# Writes simulation/figures/fig6_labels.png and fig7_ranking.png

suppressMessages(if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE))
x <- read.csv("simulation/results/07_hotcold_benchmarks.csv", stringsAsFactors = FALSE)
dir.create("simulation/figures", showWarnings = FALSE)
GREY <- "#9AA0A6"; DARK <- "#3C4043"; LIGHT <- "#DADCE0"; ACCENT <- "#1A73E8"; WARM <- "#D93025"
mean_na <- function(z) mean(z, na.rm = TRUE)
x$cfg <- paste0(x$k, "x", x$k, " ", x$scenario)
cfgs <- c("20x20 A", "20x20 B500", "20x20 B200", "20x20 B50", "40x40 B200", "40x40 B50", "20x20 C", "20x20 N", "40x40 C")
clab <- c("A\n500", "B\n500", "B\n200", "B\n50", "B\n200", "B\n50", "C", "N", "C")
grid_lab <- c(rep("20 x 20", 4), rep("40 x 40", 2), "20 x 20", "20 x 20", "40 x 40")
title_left <- function(main, sub) {
  mtext(main, side = 3, line = 2.2, adj = 0, font = 2, col = DARK, cex = 1.05)
  mtext(sub, side = 3, line = 1, adj = 0, col = GREY, cex = 0.85)
}
dotplot <- function(tab, methods, cols, ylab, ylim, which_cfg, main, sub, ref = NULL) {
  xs <- seq_along(which_cfg)
  plot(NA, xlim = c(0.5, length(xs) + 0.5), ylim = ylim, xaxt = "n", bty = "n", xlab = "", ylab = ylab, las = 1)
  abline(h = pretty(ylim), col = LIGHT)
  if (!is.null(ref)) abline(h = ref, col = DARK, lty = 2)
  axis(1, at = xs, labels = FALSE); mtext(clab[match(which_cfg, cfgs)], side = 1, at = xs, line = 1.6, cex = 0.75, col = DARK)
  mtext(grid_lab[match(which_cfg, cfgs)], side = 1, at = xs, line = 3.6, cex = 0.65, col = GREY)
  off <- seq(-0.25, 0.25, length.out = length(methods))
  for (j in seq_along(methods)) {
    v <- sapply(which_cfg, function(cf) mean_na(tab[[cf]][[methods[j]]]))
    points(xs + off[j], v, pch = 16, cex = 1.3, col = cols[j])
  }
  title_left(main, sub)
}
fig <- function(file, w, h, expr) { png(file.path("simulation/figures", file), w, h, res = 150, pointsize = 10)
  par(family = "sans", col.axis = DARK, col.lab = DARK, fg = GREY, mgp = c(2.6, 0.6, 0), tcl = -0.25)
  on.exit(dev.off()); eval.parent(substitute(expr)) }

# ---- labels (hot / not hot) --------------------------------------------------------------
L <- x[x$kind == "labels", ]
lm <- c("ESF band (pointwise)", "Gi* residuals", "LISA residuals", "Gi* counts", "LISA counts")
lc <- c(ACCENT, DARK, GREY, WARM, "#F29C94")
tabL <- function(col) lapply(setNames(cfgs, cfgs), function(cf) lapply(setNames(lm, lm), function(mm) L[[col]][L$cfg == cf & L$method == mm]))
fig("fig6_labels.png", 1900, 1300, {
  par(mfrow = c(2, 2), mar = c(5, 4.4, 4.5, 1), oma = c(0, 0, 0, 0))
  sc <- cfgs[1:6]
  dotplot(tabL("label_precision"), lm, lc, "Share of hot cells inside a cluster", c(0, 1), sc,
          "Hot cells inside a true cluster", "Label precision; scenarios with clusters")
  dotplot(tabL("sensitivity"), lm, lc, "Share of cluster cells labelled hot", c(0, 1), sc,
          "Cluster cells labelled hot", "Sensitivity; scenarios with clusters")
  dotplot(tabL("n_hot"), lm, lc, "Cells labelled hot", c(0, 260), cfgs[7:9],
          "Cells labelled hot when there are no clusters", "C: covariate area only; N: complete randomness")
  legend("topleft", lm, pch = 16, col = lc, bty = "n", cex = 0.85, inset = c(0.02, 0.02))
  dotplot(tabL("cov_hot"), lm, lc, "Share of hot cells in the covariate area", c(0, 1), cfgs[c(2:6, 7, 9)],
          "Hot cells inside the covariate area", "Scenarios B and C")
})

# ---- ranking (top 5% of cells) ---------------------------------------------------------
R <- x[x$kind == "ranking", ]
rm_ <- c("ESF: largest excess", "ESF: hot zones first", "Gi* residuals", "LISA residuals", "Top residuals", "Top counts")
rc_ <- c(ACCENT, "#8AB4F8", DARK, GREY, "#7E57C2", WARM)
tabR <- lapply(setNames(cfgs, cfgs), function(cf) lapply(setNames(rm_, rm_), function(mm) R$precision[R$cfg == cf & R$method == mm]))
fig("fig7_ranking.png", 1700, 750, {
  par(mar = c(5, 4.4, 4.5, 13))
  dotplot(tabR, rm_, rc_, "Precision (top 5% of cells)", c(0, 1), cfgs[1:6],
          "Precision of the top 5% of cells by method", "Scenarios A and B; 20 replicates (20 x 20) or 10 (40 x 40)")
  legend(par("usr")[2] + 0.1, 1, rm_, pch = 16, col = rc_, bty = "n", cex = 0.8, xpd = TRUE)
})
message("written: fig6_labels.png, fig7_ranking.png")
