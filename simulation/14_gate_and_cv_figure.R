# Figure for simulation/13_gate_and_cv.R (same layout as fig12_penalty_choice) -------------
# Run from the repo root: Rscript simulation/14_gate_and_cv_figure.R
# Writes simulation/figures/fig13_gate_and_cv.png
#
# The current ESF rule (|z|^-2, all eigenvectors) is added from experiment 10 for reference.

x13 <- read.csv("simulation/results/13_gate_and_cv.csv", stringsAsFactors = FALSE)
x10 <- read.csv("simulation/results/10_penalty_rules.csv", stringsAsFactors = FALSE)
x10 <- x10[x10$label == "z^-2" & x10$cand == "all", ]
x10 <- data.frame(k = x10$k, scenario = x10$scenario, rep = x10$rep,
                  method = "ESF current rule (|z|^-2, all)", n_hot = x10$n_hot,
                  label_precision = x10$label_precision, sensitivity = x10$sensitivity,
                  prec_top = x10$prec_excess)
x <- rbind(x13[, names(x10)], x10)
x$cfg <- paste0(x$k, "x", x$k, " ", x$scenario)
mean_na <- function(z) if (all(is.na(z))) NA else mean(z, na.rm = TRUE)
DARK <- "#3C4043"; GREY <- "#9AA0A6"; LIGHT <- "#E8EAED"
COL <- c(top = "#2a78d6", inside = "#eb6834", found = "#1baf7a")
PCH <- c(top = 16, inside = 15, found = 17)

meth <- c("Gi*", "Gi* + gate", "Gi* + FDR", "Gi* + gate + FDR",
          "ESF current rule (|z|^-2, all)",
          "ESF Lasso, random folds", "ESF Lasso, spatial folds",
          "ESF ridge, random folds", "ESF ridge, spatial folds")
mlab <- c("Gi*", "Gi* + Moran gate", "Gi* + FDR", "Gi* + Moran gate + FDR",
          "ESF, current rule", "ESF Lasso, CV random folds", "ESF Lasso, CV spatial folds",
          "ESF ridge, CV random folds", "ESF ridge, CV spatial folds")
blocks <- list(1:4, 5, 6:9)
strong <- c("20x20 A", "20x20 B500", "20x20 B200", "40x40 B200")
weak   <- c("20x20 B50", "40x40 B50")
nulls  <- c("20x20 C", "20x20 N", "40x40 C")

val <- function(col, cfgs) sapply(meth, function(mm) mean_na(sapply(cfgs, function(cf)
  mean_na(x[[col]][x$cfg == cf & x$method == mm]))))
M <- list(
  top    = list(s = val("prec_top", strong), w = val("prec_top", weak), lab = "Precision of the top 5% of cells"),
  inside = list(s = val("label_precision", strong), w = val("label_precision", weak), lab = "Hot cells inside a true cluster"),
  found  = list(s = val("sensitivity", strong), w = val("sensitivity", weak), lab = "Cluster cells labelled hot"))
fp <- sapply(meth, function(mm) mean_na(sapply(nulls, function(cf) {
  d <- x[x$cfg == cf & x$method == mm, ]; 100 * mean(d$n_hot / d$k^2) })))

png("simulation/figures/fig13_gate_and_cv.png", 2200, 1450, res = 150, pointsize = 10)
layout(matrix(1:2, 2), heights = c(3, 1.3))
par(family = "sans", col.axis = DARK, col.lab = DARK, fg = GREY, las = 1, tcl = -0.25,
    mgp = c(2.6, 0.6, 0), oma = c(0, 0, 0, 19))
xs <- seq_along(meth); sep <- c(4.5, 5.5)

par(mar = c(1, 4.5, 5, 1))
plot(NA, xlim = c(0.6, length(meth) + 0.4), ylim = c(0, 1), xaxt = "n", bty = "n",
     xlab = "", ylab = "Share (mean over scenarios)")
abline(h = seq(0, 1, 0.25), col = LIGHT); abline(v = sep, col = GREY, lty = 3)
for (j in names(M)) for (part in c("w", "s")) {
  y <- M[[j]][[part]]; lt <- if (part == "s") 1 else 2; ok <- !is.na(y)
  for (blk in blocks) { b <- blk[ok[blk]]; if (length(b) > 1) lines(b, y[b], col = COL[j], lwd = 2, lty = lt) }
  points(xs[ok], y[ok], pch = PCH[j], col = COL[j], cex = 1.2)
}
text(c(2.5, 5, 7.5), 1.06, c("Gi* on residuals", "Reference", "ESF, smooth eigenvectors, penalty by CV"),
     col = GREY, cex = 0.85, xpd = NA)
ends <- sapply(names(M), function(j) c(M[[j]]$s[length(meth)], M[[j]]$w[length(meth)]))
labs <- c(t(outer(sapply(M, `[[`, "lab"), c(", strong clusters", ", weak clusters"), paste0)))
yv <- c(ends); o <- order(yv); yl <- yv[o]
for (i in seq_along(yl)[-1]) yl[i] <- max(yl[i], yl[i - 1] + 0.065)
yl <- yl - max(0, max(yl) - 1.02)
jj <- rep(names(M), each = 2)[o]; ss <- rep(c(1, 2), 3)[o]
xe <- length(meth) + 0.6
par(xpd = NA)
for (i in seq_along(yl)) {
  segments(xe, yl[i], xe + 0.4, yl[i], col = COL[jj[i]], lwd = 2, lty = ss[i])
  points(xe + 0.2, yl[i], pch = PCH[jj[i]], col = COL[jj[i]], cex = 1)
  text(xe + 0.55, yl[i], labs[o][i], adj = 0, cex = 0.8, col = DARK)
}
par(xpd = FALSE)
mtext("Moran gate and FDR for Gi*; spatial folds and ridge for ESF", side = 3, line = 3, adj = 0, font = 2, col = DARK, cex = 1.1)
mtext("Same 150 datasets. Solid: strong clusters (A, B500, B200 on 20x20; B200 on 40x40). Dashed: weak clusters (B50 on 20x20 and 40x40). Higher is better.",
      side = 3, line = 1.9, adj = 0, col = GREY, cex = 0.85)

par(mar = c(9, 4.5, 3, 1))
plot(NA, xlim = c(0.6, length(meth) + 0.4), ylim = c(0, 6), xaxt = "n", bty = "n", xlab = "", ylab = "% of cells")
abline(h = seq(0, 6, 2), col = LIGHT); abline(v = sep, col = GREY, lty = 3)
rect(xs - 0.3, 0, xs + 0.3, fp, col = DARK, border = NA)
text(xs, fp, sprintf("%.1f", fp), pos = 3, cex = 0.8, col = DARK, xpd = NA)
axis(1, at = xs, labels = FALSE)
text(xs, -0.9, mlab, srt = 30, adj = 1, xpd = NA, cex = 0.85, col = DARK)
mtext("Cells labelled hot when there are no clusters", side = 3, line = 1.4, adj = 0, font = 2, col = DARK)
mtext("Mean of scenarios C (20x20), N (20x20) and C (40x40). Lower is better.", side = 3, line = 0.4, adj = 0, col = GREY, cex = 0.85)
dev.off()
message("written: fig13_gate_and_cv.png")
