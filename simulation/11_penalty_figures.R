# Figures for simulation/10_penalty_rules.R ---------------------------------------------
# Run from the repo root: Rscript simulation/11_penalty_figures.R
# Writes simulation/figures/fig10_penalty_ranking.png and fig11_penalty_labels.png

x <- read.csv("simulation/results/10_penalty_rules.csv", stringsAsFactors = FALSE)
dir.create("simulation/figures", showWarnings = FALSE)
DARK <- "#3C4043"; GREY <- "#9AA0A6"
mean_na <- function(z) if (all(is.na(z))) NA else mean(z, na.rm = TRUE)
x$cfg <- paste0(x$k, "x", x$k, " ", x$scenario)
rules <- c("z^-1", "z^-2", "z^-3", "z^-4", "CV (min)", "CV (1 SE)", "BIC", "Moran stop")
rlab  <- c("|z|^-1", "|z|^-2 (default)", "|z|^-3", "|z|^-4", "CV, lowest error", "CV, 1 SE", "BIC", "Moran stop")
rows  <- rbind(data.frame(label = "Gi* residuals", cand = "-", txt = "Benchmark: Gi* on residuals"),
               data.frame(label = rules, cand = "all", txt = paste0(rlab, ", all")),
               data.frame(label = rules, cand = "positive", txt = paste0(rlab, ", smooth")))

table_plot <- function(cols, clab, value, main, sub, digits = 2, zlim = c(0, 1), rev_col = FALSE) {
  M <- sapply(cols, function(cf) mapply(function(l, c) value(x[x$cfg == cf & x$label == l & x$cand == c, ]),
                                        rows$label, rows$cand))
  nr <- nrow(M); nc <- ncol(M)
  plot(NA, xlim = c(0, nc), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "")
  pal <- colorRamp(c("white", "#1A73E8"))
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    v <- M[i, j]; yy <- nr - i
    r <- if (is.na(v)) 0 else min(1, max(0, (v - zlim[1]) / diff(zlim)))
    if (rev_col) r <- 1 - r
    rect(j - 1, yy, j, yy + 1, col = if (is.na(v)) "#F1F3F4" else rgb(pal(r * 0.8), maxColorValue = 255), border = "white")
    text(j - 0.5, yy + 0.5, if (is.na(v)) "-" else formatC(v, format = "f", digits = digits),
         cex = 0.75, col = if (r > 0.55) "white" else DARK)
  }
  axis(2, at = nr - seq_len(nr) + 0.5, labels = rows$txt, las = 1, tick = FALSE, cex.axis = 0.75, col.axis = DARK)
  mtext(clab, side = 3, at = seq_len(nc) - 0.5, line = 0.2, cex = 0.7, col = DARK)
  abline(h = c(nr - 1, nr - 1 - length(rules)), col = DARK, lwd = 1)
  mtext(main, side = 3, line = 4, adj = 0, font = 2, col = DARK, cex = 0.95)
  mtext(sub, side = 3, line = 3, adj = 0, col = GREY, cex = 0.75)
}
withc <- c("20x20 A", "20x20 B500", "20x20 B200", "20x20 B50", "40x40 B200", "40x40 B50")
withl <- c("A\n20x20", "B500\n20x20", "B200\n20x20", "B50\n20x20", "B200\n40x40", "B50\n40x40")
nulc  <- c("20x20 C", "20x20 N", "40x40 C"); null_l <- c("C\n20x20", "N\n20x20", "C\n40x40")

png("simulation/figures/fig10_penalty_ranking.png", 2100, 1150, res = 150, pointsize = 10)
par(mfrow = c(1, 2), mar = c(1, 12.5, 6, 1), family = "sans")
table_plot(withc, withl, function(d) mean_na(d$prec_excess), "Precision of the top 5% of cells",
           "Ranking: largest excess. Mean over replicates")
table_plot(withc, withl, function(d) stats::median(d$n_eig), "Eigenvectors selected (median)",
           "All: 400 (20x20) or 1600 (40x40) candidates; smooth: about 21% of them",
           digits = 0, zlim = c(0, 400))
dev.off()

png("simulation/figures/fig11_penalty_labels.png", 2000, 2300, res = 150, pointsize = 10)
layout(matrix(c(1, 2, 3, 0), 2, byrow = TRUE))
par(mar = c(1, 12.5, 6, 1), family = "sans")
table_plot(withc, withl, function(d) mean_na(d$label_precision), "Hot cells inside a true cluster",
           "Label precision; mean over replicates with at least one hot cell")
table_plot(withc, withl, function(d) mean_na(d$sensitivity), "Cluster cells labelled hot", "Sensitivity")
table_plot(nulc, null_l, function(d) mean_na(d$n_hot), "Cells labelled hot with no clusters",
           "C: covariate area only; N: complete randomness. Mean", digits = 0, zlim = c(0, 100), rev_col = FALSE)
dev.off()
message("written: fig10_penalty_ranking.png, fig11_penalty_labels.png")
