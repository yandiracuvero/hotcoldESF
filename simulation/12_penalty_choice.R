# One figure to choose the penalty rule (results of simulation/10_penalty_rules.R) ---------
# Run from the repo root: Rscript simulation/12_penalty_choice.R
# Writes simulation/figures/fig12_penalty_choice.png
#
# x axis: the 17 variants (Gi* benchmark + 8 penalty rules x 2 candidate sets)
# top panel (0-1): precision of the top 5%, hot cells inside a true cluster, cluster
#   cells labelled hot; solid = strong clusters (20x20 A, B500, B200; 40x40 B200),
#   dashed = weak clusters (20x20 and 40x40 B50); mean of the scenario means
# bottom panel: cells labelled hot when there are no clusters (% of the map; mean of
#   20x20 C, 20x20 N, 40x40 C). Different unit, so a separate panel (no second axis).

x <- read.csv("simulation/results/10_penalty_rules.csv", stringsAsFactors = FALSE)
x$cfg <- paste0(x$k, "x", x$k, " ", x$scenario)
mean_na <- function(z) if (all(is.na(z))) NA else mean(z, na.rm = TRUE)
DARK <- "#3C4043"; GREY <- "#9AA0A6"; LIGHT <- "#E8EAED"
COL <- c(top = "#2a78d6", inside = "#eb6834", found = "#1baf7a")
PCH <- c(top = 16, inside = 15, found = 17)

rules <- c("z^-1", "z^-2", "z^-3", "z^-4", "CV (min)", "CV (1 SE)", "BIC", "Moran stop")
rlab  <- c("|z|^-1", "|z|^-2 (default)", "|z|^-3", "|z|^-4", "CV, lowest error", "CV, 1 SE", "BIC", "Moran stop")
v <- rbind(data.frame(label = "Gi* residuals", cand = "-", txt = "Benchmark: Gi* on residuals"),
           data.frame(label = rules, cand = "all", txt = paste0(rlab, ", all")),
           data.frame(label = rules, cand = "positive", txt = paste0(rlab, ", smooth")))
strong <- c("20x20 A", "20x20 B500", "20x20 B200", "40x40 B200")
weak   <- c("20x20 B50", "40x40 B50")
nulls  <- c("20x20 C", "20x20 N", "40x40 C")

val <- function(col, cfgs) sapply(seq_len(nrow(v)), function(i) {
  m <- sapply(cfgs, function(cf) {
    d <- x[x$cfg == cf & x$label == v$label[i] & x$cand == v$cand[i], ]
    if (nrow(d) == 0) NA else mean_na(d[[col]])
  })
  mean_na(m)
})
M <- list(
  top    = list(s = val("prec_excess", strong), w = val("prec_excess", weak),
                lab = "Precision of the top 5% of cells"),
  inside = list(s = val("label_precision", strong), w = val("label_precision", weak),
                lab = "Hot cells inside a true cluster"),
  found  = list(s = val("sensitivity", strong), w = val("sensitivity", weak),
                lab = "Cluster cells labelled hot"))
fp <- sapply(seq_len(nrow(v)), function(i) mean_na(sapply(nulls, function(cf) {
  d <- x[x$cfg == cf & x$label == v$label[i] & x$cand == v$cand[i], ]
  if (nrow(d) == 0) NA else 100 * mean(d$n_hot / d$k^2) })))

png("simulation/figures/fig12_penalty_choice.png", 2400, 1500, res = 150, pointsize = 10)
layout(matrix(1:2, 2), heights = c(3, 1.25))
par(family = "sans", col.axis = DARK, col.lab = DARK, fg = GREY, las = 1, tcl = -0.25,
    mgp = c(2.6, 0.6, 0), oma = c(0, 0, 0, 19))
xs <- seq_len(nrow(v)); sep <- c(1.5, 9.5)

# ---- top panel -------------------------------------------------------------------------
par(mar = c(1, 4.5, 5, 1))
plot(NA, xlim = c(0.6, nrow(v) + 0.4), ylim = c(0, 1), xaxt = "n", bty = "n",
     xlab = "", ylab = "Share (mean over scenarios)")
abline(h = seq(0, 1, 0.25), col = LIGHT)
abline(v = sep, col = GREY, lty = 3)
for (j in names(M)) {
  for (part in c("w", "s")) {
    y <- M[[j]][[part]]; lt <- if (part == "s") 1 else 2
    ok <- !is.na(y)
    for (blk in list(1, 2:9, 10:17)) {                   # no line across group borders
      b <- blk[ok[blk]]
      if (length(b) > 1) lines(b, y[b], col = COL[j], lwd = 2, lty = lt)
    }
    points(xs[ok], y[ok], pch = PCH[j], col = COL[j], cex = 1.2)
  }
}
text(c(1, 5.5, 13.5), 1.06, c("Benchmark", "All eigenvectors", "Smooth eigenvectors only"),
     col = GREY, cex = 0.85, xpd = NA)
# direct labels at the right end (text in ink, a coloured marker beside it)
ends <- sapply(names(M), function(j) c(M[[j]]$s[17], M[[j]]$w[17]))
labs <- c(t(outer(sapply(M, `[[`, "lab"), c(", strong clusters", ", weak clusters"), paste0)))
yv <- c(ends); o <- order(yv); yl <- yv[o]
for (i in seq_along(yl)[-1]) yl[i] <- max(yl[i], yl[i - 1] + 0.065)
yl <- yl - max(0, max(yl) - 1.02)
jj <- rep(names(M), each = 2)[o]; ss <- rep(c(1, 2), 3)[o]
par(xpd = NA)
for (i in seq_along(yl)) {
  segments(17.6, yl[i], 18.2, yl[i], col = COL[jj[i]], lwd = 2, lty = ss[i])
  points(17.9, yl[i], pch = PCH[jj[i]], col = COL[jj[i]], cex = 1)
  text(18.4, yl[i], labs[o][i], adj = 0, cex = 0.8, col = DARK)
}
par(xpd = FALSE)
mtext("Penalty rules compared on the same 150 datasets", side = 3, line = 3, adj = 0, font = 2, col = DARK, cex = 1.1)
mtext("Solid: strong clusters (A, B500, B200 on 20x20; B200 on 40x40). Dashed: weak clusters (B50 on 20x20 and 40x40). Higher is better.",
      side = 3, line = 1.9, adj = 0, col = GREY, cex = 0.85)

# ---- bottom panel ----------------------------------------------------------------------
par(mar = c(9.5, 4.5, 3, 1))
plot(NA, xlim = c(0.6, nrow(v) + 0.4), ylim = c(0, 6), xaxt = "n", bty = "n",
     xlab = "", ylab = "% of cells")
abline(h = seq(0, 6, 2), col = LIGHT)
abline(v = sep, col = GREY, lty = 3)
ok <- !is.na(fp)
rect(xs[ok] - 0.3, 0, xs[ok] + 0.3, fp[ok], col = DARK, border = NA)
text(xs[ok], fp[ok], sprintf("%.1f", fp[ok]), pos = 3, cex = 0.75, col = DARK, xpd = NA)
text(xs[!ok], 0.3, "not run", cex = 0.7, col = GREY)
axis(1, at = xs, labels = FALSE)
text(xs, -0.9, v$txt, srt = 40, adj = 1, xpd = NA, cex = 0.8, col = DARK)
mtext("Cells labelled hot when there are no clusters", side = 3, line = 1.4, adj = 0, font = 2, col = DARK)
mtext("Mean of scenarios C (20x20), N (20x20) and C (40x40). Lower is better.", side = 3, line = 0.4, adj = 0, col = GREY, cex = 0.85)
dev.off()
message("written: fig12_penalty_choice.png")
