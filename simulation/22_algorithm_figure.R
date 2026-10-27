# Figure: the procedure on one simulated dataset (Section 2, Algorithm 1) -------------------
# Scenario K (strong): three clusters of 200 points, a cold area keeping 20% of its points,
# 2000 uniform points and a covariate area of 3000 points; 20 x 20 grid; B = 999.
# Seed 10 was chosen as a typical dataset: its label precisions (hot 0.78, cold 0.80) are close
# to the medians over seeds 1-30 (0.78 and 0.80).
# Run from the repo root: Rscript simulation/22_algorithm_figure.R
# Writes simulation/figures/fig0_algorithm.png

suppressMessages(if (file.exists("DESCRIPTION")) pkgload::load_all(".", quiet = TRUE) else library(hotcoldESF))
SEED <- 10
set.seed(SEED)
s <- sim_thomas(20, 3, 100, 200, n_background = 2000, n_cov = 3000, min_dist = 300, n_cold = 1, cold_keep = 0.2)
set.seed(SEED + 7)
h <- hot_cold_map(s$Y, s$Z, s$A, nperm = 999, penalty = "cv_min", mc_min = 0.25, fit_type = "ridge")

# the null distributions of fhat for three cells (same permutations as the band)
k <- s$k; n <- k^2; side <- s$cell_side
y  <- as.numeric(scale(log1p(s$Y))); Zs <- scale(s$Z)
ols <- lm.fit(design0(n, Zs), y)
Sm <- ridge_smoother(y, Zs, eigen(s$A, symmetric = TRUE), h$lambda, mc_min = 0.25)
set.seed(SEED + 7); invisible(select_lambda(y, Zs, s$A, eigen(s$A, symmetric = TRUE), mc_min = 0.25,
                                            rule = "cv_min", z0 = h$z0, mix = 0))   # same RNG stream as h
null <- Sm %*% vapply(seq_len(999), function(b) sample(ols$residuals), numeric(n))
f <- h$fhat; lab <- h$labels
stopifnot(all.equal(apply(null, 1, quantile, 0.95, names = FALSE), h$upper))   # same band as hot_cold_map
pick <- c(hot = which.max(f), neutral = which(lab == "neutral")[which.max(f[lab == "neutral"])],
          cold = which.min(f))

INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; AXIS <- "#c3c2b7"
HOT <- "#e34948"; COLD <- "#2a78d6"; NEUT <- "#f0efec"
GREENS <- colorRampPalette(c("#f2f7ee", "#a8d08d", "#4d9a3a", "#1f5f1a"))(100)
DIV <- colorRampPalette(c("#184f95", "#6da7ec", "#f0efec", "#ec8a89", "#a32727"))(101)

cellxy <- function(i) c((i - 1) %% k, (i - 1) %/% k) * side
blank <- function() plot(NA, xlim = c(0, 1000), ylim = c(0, 1000), asp = 1, axes = FALSE, xlab = "", ylab = "",
                         xaxs = "i", yaxs = "i")
grid_img <- function(v, cols, zlim) {
  m <- matrix(v, k, k)                                     # cell i: x = (i-1) %% k, y = (i-1) %/% k
  blank()
  image(seq(0, 1000, length.out = k + 1), seq(0, 1000, length.out = k + 1), m, col = cols, zlim = zlim,
        add = TRUE, useRaster = TRUE)
}
truth <- function(col = INK) {
  th <- seq(0, 2 * pi, length.out = 100)
  for (i in seq_len(nrow(s$parents))) lines(s$parents[i, 1] + s$radius[i] * cos(th), s$parents[i, 2] + s$radius[i] * sin(th),
                                            lty = 2, col = col, lwd = 0.8)
  cc <- s$centres[s$cold_cells, , drop = FALSE]; ctr <- colMeans(cc)
  lines(ctr[1] + 150 * cos(th), ctr[2] + 150 * sin(th), lty = 3, col = col, lwd = 1)
}
frame <- function() rect(0, 0, 1000, 1000, border = AXIS, lwd = 0.8)
title_p <- function(letter, main, sub, sub2 = NULL) {
  mtext(if (is.character(main)) paste0(letter, "  ", main) else main, side = 3, line = 2.3, adj = 0, font = 2, col = INK, cex = 1)
  if (!is.character(main)) mtext(letter, side = 3, line = 2.3, adj = 0, font = 2, col = INK, cex = 1, at = par("usr")[1] - 0.06 * diff(par("usr")[1:2]))
  mtext(sub, side = 3, line = 1.3, adj = 0, col = INK2, cex = 0.9)
  if (!is.null(sub2)) mtext(sub2, side = 3, line = 0.3, adj = 0, col = INK2, cex = 0.9)
}
cbar <- function(cols, labs, at) {                         # small horizontal colour bar under a map
  u <- par("usr"); y0 <- u[3] - 0.07 * diff(u[3:4]); y1 <- y0 - 0.035 * diff(u[3:4])
  xs <- seq(u[1] + 0.15 * diff(u[1:2]), u[2] - 0.15 * diff(u[1:2]), length.out = length(cols) + 1)
  rect(xs[-length(xs)], y1, xs[-1], y0, col = cols, border = NA, xpd = NA)
  text(xs[1] + at * (xs[length(xs)] - xs[1]), y1, labs, pos = 1, cex = 0.85, col = INK2, xpd = NA, offset = 0.25)
}

png("simulation/figures/fig0_algorithm.png", width = 6.69, height = 2.45, units = "in", res = 300, pointsize = 8, type = "cairo")
layout(matrix(1:5, 1), widths = c(1, 1, 1, 1.05, 1))
par(family = "sans", oma = c(0.6, 0.2, 0.1, 0.2), cex = 1)

# (a) points and covariate area
par(mar = c(2.6, 0.5, 3.6, 0.5))
blank()
bx <- seq(0, 1000, 10); bz <- outer(bx, bx, function(x, yy) exp(-((x - s$cov_center[1])^2 + (yy - s$cov_center[2])^2) / (2 * 150^2)))
.filled <- function() for (lv in c(0.25, 0.5, 0.75)) polygon(contourLines(bx, bx, bz, levels = lv)[[1]][c("x", "y")],
                                                           col = adjustcolor("#898781", 0.18), border = NA)
.filled()
points(s$points, pch = 16, cex = 0.06, col = adjustcolor(INK, 0.55))
truth(); frame()
title_p("a", "Data", "points and covariate")

# (b) counts per cell
grid_img(log1p(s$Y), GREENS, range(log1p(s$Y))); truth(); frame()
title_p("b", "Counts", "lines 2-3", bquote("Moran " * z[0] == .(round(h$z0, 1))))
cbar(GREENS, c(0, max(s$Y)), c(0, 1))

# (c) spatial component
zl <- max(abs(f)) * c(-1, 1)
grid_img(f, DIV, zl); truth(INK2); frame()
for (nm in names(pick)) { p <- cellxy(pick[nm]) + side / 2
  points(p[1], p[2], pch = 21, bg = "white", col = INK, cex = 1.7, lwd = 0.8)
  text(p[1], p[2], c(hot = "1", neutral = "2", cold = "3")[nm], cex = 0.75, font = 2, col = INK) }
title_p("c", "Spatial component", "lines 4-5")
cbar(DIV, formatC(c(zl[1], 0, zl[2]), format = "f", digits = 2), c(0, 0.5, 1))

# (d) permutation band for the three marked cells
par(mar = c(2.6, 0.6, 3.6, 0.6))
xr <- range(c(null[pick, ], f[pick])); xr <- xr + c(-0.05, 0.05) * diff(xr)
plot(NA, xlim = xr, ylim = c(0.4, 3.6), axes = FALSE, xlab = "", ylab = "")
abline(v = 0, col = AXIS, lwd = 0.8)
for (j in 1:3) {
  i <- pick[j]; d <- density(null[i, ]); yy <- 4 - j
  polygon(c(d$x, rev(d$x)), yy + c(d$y, -rev(d$y)) / max(d$y) * 0.32, col = "#e1e0d9", border = MUTED, lwd = 0.6)
  q <- quantile(null[i, ], c(0.05, 0.95))
  segments(q, yy - 0.36, q, yy + 0.36, col = INK2, lwd = 1)
  colv <- c(hot = HOT, neutral = MUTED, cold = COLD)[as.character(lab[i])]
  points(f[i], yy, pch = 18, cex = 1.6, col = colv)
  text(xr[1], yy + 0.5, paste0(j, ": ", lab[i]), adj = c(0, 0.5), cex = 0.9, col = INK, font = 2)
}
axis(1, col = AXIS, col.axis = INK2, cex.axis = 0.85, tcl = -0.2, mgp = c(2, 0.35, 0))
mtext(expression(hat(f)[i]), side = 1, line = 1.5, col = INK2, cex = 0.9)
title_p("d", "Permutation band", "lines 6-9", "cells 1-3 of (c)")


# (e) labels
par(mar = c(2.6, 0.5, 3.6, 0.5))
lv <- c(hot = 3, neutral = 2, cold = 1)[as.character(lab)]
grid_img(lv, c(COLD, NEUT, HOT), c(1, 3)); truth(); frame()
title_p("e", "Hot and cold cells", "lines 9-10")
u <- par("usr")
legend(mean(u[1:2]), u[3] - 0.03 * diff(u[3:4]), c("hot", "neutral", "cold"), fill = c(HOT, NEUT, COLD), xjust = 0.5,
       border = AXIS, horiz = TRUE, bty = "n", cex = 0.85, text.col = INK2, xpd = NA, x.intersp = 0.3,
       text.width = strwidth(c("hot  ", "neutral  ", "cold"), cex = 0.85))
dev.off()
message("written: fig0_algorithm.png")
