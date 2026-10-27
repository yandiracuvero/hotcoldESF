# Figure for Section 3.2, choice 7 (ridge with a cross-validated penalty) ------------------
# Data: simulation/results/10_penalty_rules.csv (eight Lasso penalty rules x two candidate
# sets) and simulation/results/13_gate_and_cv.csv (ridge, penalty by 10-fold CV, smooth
# candidates, random folds). Both use the same datasets: 20 replicates on 20x20, 10 on 40x40.
# Run from the repo root: Rscript simulation/21_penalty_choice_figure.R
# Writes simulation/figures/fig17_penalty_choice.png

x <- read.csv("simulation/results/10_penalty_rules.csv", stringsAsFactors = FALSE)
y <- read.csv("simulation/results/13_gate_and_cv.csv", stringsAsFactors = FALSE)
x$set <- paste0(x$scenario, "_", x$k); y$set <- paste0(y$scenario, "_", y$k)
x$val <- x$prec_excess; y$val <- y$prec_top

INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
SMOOTH <- "#2a78d6"; ALL <- "#898781"; REF <- "#eb6834"

rules <- c("z^-1", "z^-2", "z^-3", "z^-4", "CV (min)", "CV (1 SE)", "BIC", "Moran stop")
rlab  <- c(expression("Lasso, " * group("|", z[0], "|")^-1),
           expression("Lasso, " * group("|", z[0], "|")^-2 * " (MI-Lasso)"),
           expression("Lasso, " * group("|", z[0], "|")^-3),
           expression("Lasso, " * group("|", z[0], "|")^-4),
           "Lasso, CV lowest error", "Lasso, CV 1 SE", "Lasso, BIC", "Lasso, Moran stop")
strong <- c("A_20", "B500_20", "B200_20", "B200_40")

stat <- function(d, sets, worst = FALSE) {
  if (!worst) { v <- d$val[d$set == sets]; return(c(m = mean(v), se = sd(v) / sqrt(length(v)))) }
  ms <- sapply(sets, function(s) mean(d$val[d$set == s])); if (any(is.na(ms))) return(c(m = NA, se = NA))
  s <- sets[which.min(ms)]; v <- d$val[d$set == s]; c(m = mean(v), se = sd(v) / sqrt(length(v)))
}
panels <- list(list(sets = "B50_20", worst = FALSE, main = expression(bold("Weak clusters, 20" %*% "20"))),
               list(sets = "B50_40", worst = FALSE, main = expression(bold("Weak clusters, 40" %*% "40"))),
               list(sets = strong,   worst = TRUE,  main = expression(bold("Strong clusters (worst)"))))

nr <- length(rules) + 1; ypos <- c(nr + 0.6 - seq_along(rules), 0.6)   # ridge row set apart at the bottom
png("simulation/figures/fig17_penalty_choice.png", width = 6.69, height = 3.4, units = "in", res = 300, pointsize = 8,
    type = "cairo")
layout(matrix(c(1, 2, 3, 4, 5, 5, 5, 5), 2, byrow = TRUE), widths = c(1.05, 1, 1, 1), heights = c(1, 0.12))
par(family = "sans", col.axis = INK2, fg = AXIS, las = 1, tcl = -0.2, mgp = c(2, 0.5, 0), cex = 1, oma = c(0, 0, 0.4, 0))
# row labels
par(mar = c(3.2, 0, 1.8, 0.3)); plot(NA, xlim = c(0, 1), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "")
for (i in seq_along(rules)) text(1, ypos[i], rlab[i], adj = c(1, 0.5), col = if (rules[i] == "CV (min)") INK else INK2,
                                 font = if (rules[i] == "CV (min)") 2 else 1, cex = 0.9)
text(1, ypos[nr], "Ridge, CV lowest error", adj = c(1, 0.5), col = INK, font = 2, cex = 0.9)
for (p in panels) {
  par(mar = c(3.2, 0.8, 1.8, 0.8))
  plot(NA, xlim = c(0, 1), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "")
  abline(v = seq(0, 1, 0.25), col = GRID, lwd = 0.8); abline(h = 1.1, col = AXIS, lwd = 0.8)
  axis(1, at = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"), col = AXIS, cex.axis = 0.85, mgp = c(2, 0.4, 0), tcl = -0.2)
  mtext("Precision of the top 5%", side = 1, line = 1.8, col = INK2, cex = 0.9)
  mtext(p$main, side = 3, line = 0.5, adj = 0, col = INK, cex = 1)
  if (!p$worst) {                                    # reference: Gi* on residuals (same data)
    g <- stat(x[x$label == "Gi* residuals", ], p$sets)["m"]
    abline(v = g, col = REF, lty = 2, lwd = 1.2)
    text(g, nr + 0.05, expression(G[i]^"*"), col = REF, adj = c(-0.25, 0.5), cex = 0.9, xpd = NA)
  }
  for (i in seq_along(rules)) for (cd in c("all", "positive")) {
    s <- stat(x[x$label == rules[i] & x$cand == cd, ], p$sets, p$worst); if (is.na(s["m"])) next
    off <- if (cd == "all") 0.17 else -0.17; col <- if (cd == "all") ALL else SMOOTH
    segments(max(0, s["m"] - 2 * s["se"]), ypos[i] + off, min(1, s["m"] + 2 * s["se"]), ypos[i] + off, col = col, lwd = 1)
    points(s["m"], ypos[i] + off, pch = if (cd == "all") 1 else 16, col = col, cex = 1, lwd = 1.2)
  }
  s <- stat(y[y$method == "ESF ridge, random folds", ], p$sets, p$worst)
  segments(max(0, s["m"] - 2 * s["se"]), ypos[nr], min(1, s["m"] + 2 * s["se"]), ypos[nr], col = SMOOTH, lwd = 1)
  points(s["m"], ypos[nr], pch = 18, col = SMOOTH, cex = 1.5)
  text(s["m"], ypos[nr], formatC(s["m"], format = "f", digits = 2), pos = 2,
       col = INK, cex = 0.85, offset = 0.6)
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", c("Lasso, all eigenvectors", "Lasso, smooth candidates", "Ridge, smooth candidates",
                   expression(G[i]^"*" * " on the residuals")),
       pch = c(1, 16, 18, NA), lty = c(NA, NA, NA, 2), col = c(ALL, SMOOTH, SMOOTH, REF), pt.cex = c(1, 1, 1.5, 1),
       pt.lwd = 1.2, horiz = TRUE, bty = "n", text.col = INK2, cex = 0.9, x.intersp = 0.5, seg.len = 1.4)
dev.off()
message("written: fig17_penalty_choice.png")
