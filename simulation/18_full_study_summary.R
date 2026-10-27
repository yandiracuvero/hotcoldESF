# Tables and figures for simulation/17_full_study.R ---------------------------------------
# Run from the repo root: Rscript simulation/18_full_study_summary.R
# Writes
#   simulation/results/17_summary.csv            mean and Monte Carlo SE per setting x method
#   paper/tables/sim_ranking.tex, sim_labels.tex  LaTeX tables for the paper
#   simulation/figures/fig15_full_ranking.png     ranking (top 5% of cells)
#   simulation/figures/fig16_full_labels.png      hot/cold labels

x <- read.csv("simulation/results/17_full_study.csv", stringsAsFactors = FALSE)
dir.create("paper/tables", showWarnings = FALSE, recursive = TRUE)
mean_na <- function(z) if (all(is.na(z))) NA else mean(z, na.rm = TRUE)
se_na <- function(z) { z <- z[!is.na(z)]; if (length(z) < 2) NA else sd(z) / sqrt(length(z)) }
DARK <- "#3C4043"; GREY <- "#9AA0A6"; LIGHT <- "#E8EAED"

scn_lab <- c(B = "B: compact clusters", L = "L: large diffuse clusters", E = "E: elongated clusters",
             M = "M: mixed sizes", S = "S: many small clusters", K = "K: clusters + cold area",
             U = "U: covariate not in the model", C = "C: covariate area only", N = "N: random points")
setting <- function(d) paste(d$scenario, d$strength, d$k)

# ---- summary table (all measures) ----------------------------------------------------------
measures <- c("top_precision", "clusters_top", "top_in_cov", "n_hot", "hot_precision", "hot_sensitivity",
              "clusters_hot", "hot_in_cov", "n_cold", "cold_precision", "cold_sensitivity")
key <- unique(x[, c("k", "scenario", "strength", "method")])
summ <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
  d <- x[x$k == key$k[i] & x$scenario == key$scenario[i] & x$strength == key$strength[i] &
         x$method == key$method[i], ]
  cbind(key[i, ], n_rep = nrow(d), gate_pass = mean(d$gate_ok),
        as.data.frame(as.list(c(sapply(measures, function(mm) mean_na(d[[mm]])),
                                setNames(sapply(measures, function(mm) se_na(d[[mm]])), paste0(measures, "_se"))))))
}))
write.csv(summ, "simulation/results/17_summary.csv", row.names = FALSE)
get <- function(k, sc, st, me, col) {
  v <- summ[summ$k == k & summ$scenario == sc & summ$strength == st & summ$method == me, col]
  if (length(v)) v else NA
}

# ---- LaTeX tables ------------------------------------------------------------------------
fmt <- function(v) ifelse(is.na(v), "--", formatC(v, format = "f", digits = 2))
rank_m <- c("ESF", "Gi*", "LISA + gate", "Top residuals", "Top counts")
lab_m  <- c("ESF", "Gi*", "Gi* + gate", "Gi* + gate + FDR", "LISA + gate")
tex_m  <- function(m) gsub("\\*", "$^*$", m)
sc_cl  <- c("B", "L", "E", "M", "S", "K", "U")
lines_r <- c("\\begin{tabular}{ll", paste(rep("r", length(rank_m)), collapse = ""), "}", "\\toprule",
             paste0("Scenario & Grid / strength & ", paste(tex_m(rank_m), collapse = " & "), " \\\\"), "\\midrule")
lines_r <- c(paste0(lines_r[1:3], collapse = ""), lines_r[-(1:3)])
for (sc in sc_cl) {
  for (k in c(20, 40)) for (st in c("weak", "strong")) {
    v <- sapply(rank_m, function(me) get(k, sc, st, me, "top_precision"))
    b <- which(v == max(v, na.rm = TRUE))
    s <- fmt(v); s[b] <- paste0("\\textbf{", s[b], "}")
    lines_r <- c(lines_r, paste0(if (k == 20 && st == "weak") sc else "", " & ", k, "$\\times$", k, " ", st,
                                 " & ", paste(s, collapse = " & "), " \\\\"))
  }
  if (sc != tail(sc_cl, 1)) lines_r <- c(lines_r, "\\addlinespace")
}
lines_r <- c(lines_r, "\\bottomrule", "\\end{tabular}")
writeLines(lines_r, "paper/tables/sim_ranking.tex")

lines_l <- c(paste0("\\begin{tabular}{ll", paste(rep("r", length(lab_m)), collapse = ""), "}"), "\\toprule",
             paste0("Measure & Setting & ", paste(tex_m(lab_m), collapse = " & "), " \\\\"), "\\midrule")
add <- function(name, col, sets, digits = 2, pct = FALSE) {
  for (i in seq_along(sets)) {
    st <- sets[[i]]
    v <- sapply(lab_m, function(me) mean_na(sapply(st$k, function(k) get(k, st$sc, st$st, me, col)) /
                                             if (pct) (st$k^2 / 100) else 1))
    lines_l <<- c(lines_l, paste0(if (i == 1) name else "", " & ", st$lab, " & ",
                                  paste(ifelse(is.na(v), "--", formatC(v, format = "f", digits = digits)), collapse = " & "), " \\\\"))
  }
  lines_l <<- c(lines_l, "\\addlinespace")
}
avg_sets <- function(st) lapply(c(20, 40), function(k) list(k = k, sc = NA, st = st))
mk <- function(sc, st, k, lab) list(sc = sc, st = st, k = k, lab = lab)
add("Hot cells inside a cluster", "hot_precision",
    list(mk("B", "strong", 20, "B strong 20x20"), mk("B", "weak", 20, "B weak 20x20"),
         mk("M", "strong", 40, "M strong 40x40"), mk("S", "weak", 40, "S weak 40x40")))
add("Cluster cells labelled hot", "hot_sensitivity",
    list(mk("B", "strong", 20, "B strong 20x20"), mk("B", "weak", 20, "B weak 20x20"),
         mk("M", "strong", 40, "M strong 40x40"), mk("S", "weak", 40, "S weak 40x40")))
add("\\% of cells hot, no clusters", "n_hot",
    list(mk("C", "none", 20, "C 20x20"), mk("N", "none", 20, "N 20x20"),
         mk("C", "none", 40, "C 40x40"), mk("N", "none", 40, "N 40x40")), digits = 1, pct = TRUE)
add("Cold cells inside the cold area", "cold_precision",
    list(mk("K", "strong", 20, "K strong 20x20"), mk("K", "strong", 40, "K strong 40x40")))
add("Cold-area cells labelled cold", "cold_sensitivity",
    list(mk("K", "strong", 20, "K strong 20x20"), mk("K", "strong", 40, "K strong 40x40")))
add("Hot cells in the covariate area", "hot_in_cov",
    list(mk("U", "strong", 20, "U strong 20x20"), mk("U", "strong", 40, "U strong 40x40")))
lines_l <- c(head(lines_l, -1), "\\bottomrule", "\\end{tabular}")
lines_l <- gsub("x20", "$\\\\times$20", gsub("x40", "$\\\\times$40", lines_l, fixed = FALSE))
writeLines(lines_l, "paper/tables/sim_labels.tex")

# ---- figures: drawn at their printed size (17 cm wide, 8 pt text) -------------------------
FIG_W <- 6.69; PT <- 8                                  # \textwidth of the paper, in inches
fig_png <- function(file, h) png(file, width = FIG_W, height = h, units = "in", res = 300, pointsize = PT, type = "cairo")
fig_par <- function() par(family = "sans", col.axis = DARK, col.lab = DARK, fg = GREY, las = 1, tcl = -0.2,
                          mgp = c(2, 0.45, 0), cex.axis = 0.85, cex.lab = 0.9, mar = c(3, 3.1, 1.6, 0.4), oma = c(0, 0, 0.3, 0), cex = 1)
short <- c(B = "B: compact clusters", L = "L: large clusters", E = "E: elongated clusters",
           M = "M: mixed sizes", S = "S: small clusters", K = "K: cold area", U = "U: covariate missing",
           C = "C: covariate only", N = "N: random points")

# ---- figure 15: ranking ---------------------------------------------------------------------
RCOL <- c("ESF" = "#2a78d6", "Gi*" = "#eb6834", "LISA + gate" = "#1baf7a",
          "Top residuals" = "#eda100", "Top counts" = "#e87ba4")
RPCH <- c(16, 15, 17, 18, 4)
xset <- data.frame(k = c(20, 20, 40, 40), st = c("weak", "strong", "weak", "strong"),
                   lab = c("weak", "strong", "weak", "strong"))
panel_dots <- function(sc, col, methods, cols, pchs, ylab, ylim = c(0, 1), main = short[sc]) {
  plot(NA, xlim = c(0.6, 4.4), ylim = ylim, xaxt = "n", bty = "n", xlab = "", ylab = ylab)
  abline(h = pretty(ylim), col = LIGHT, lwd = 0.6); abline(v = 2.5, col = GREY, lty = 3, lwd = 0.6)
  axis(1, at = 1:4, labels = FALSE)
  text(1:4, ylim[1] - 0.09 * diff(ylim), xset$lab, xpd = NA, cex = 0.85, col = DARK, adj = c(0.5, 1))
  text(c(1.5, 3.5), ylim[1] - 0.2 * diff(ylim), c("20\u00d720", "40\u00d740"), xpd = NA, cex = 0.85, col = DARK, adj = c(0.5, 1))
  off <- seq(-0.26, 0.26, length.out = length(methods))
  for (i in seq_along(methods)) {
    v  <- mapply(function(k, st) get(k, sc, st, methods[i], col), xset$k, xset$st)
    se <- mapply(function(k, st) get(k, sc, st, methods[i], paste0(col, "_se")), xset$k, xset$st)
    segments(1:4 + off[i], v - 2 * se, 1:4 + off[i], v + 2 * se, col = cols[i], lwd = 0.8)
    points(1:4 + off[i], v, pch = pchs[i], col = cols[i], cex = 0.9)
  }
  mtext(main, side = 3, line = 0.4, adj = 0, font = 2, col = DARK, cex = 1)
}
legend_row <- function(cols, pchs) {
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", names(cols), col = cols, pch = pchs, horiz = TRUE, bty = "n", text.col = DARK, cex = 0.95,
         x.intersp = 0.6, text.width = NA)
}
fig_png("simulation/figures/fig15_full_ranking.png", 4.3)
layout(rbind(1:4, 5:8, 9), heights = c(1, 1, 0.16)); fig_par()
for (sc in sc_cl) panel_dots(sc, "top_precision", names(RCOL), RCOL, RPCH, "Precision (top 5%)")
# 8th panel: scenario S, share of clusters with a cell in the top 5%
panel_dots("S", "clusters_top", names(RCOL), RCOL, RPCH, "Share of clusters", main = "S: clusters found")
legend_row(RCOL, RPCH)
dev.off()

# ---- figure 16: labels -------------------------------------------------------------------
LCOL <- c("ESF" = "#2a78d6", "Gi*" = "#e87ba4", "Gi* + gate" = "#eb6834",
          "Gi* + gate + FDR" = "#eda100", "LISA + gate" = "#1baf7a")
LPCH <- c(16, 4, 15, 18, 17)
fig_png("simulation/figures/fig16_full_labels.png", 6.2)
layout(rbind(1:4, 5:8, 9:12, 13), heights = c(1, 1, 1, 0.14)); fig_par()
for (sc in c("B", "L", "M", "S")) panel_dots(sc, "hot_precision", names(LCOL), LCOL, LPCH, "Hot cells in a cluster")
for (sc in c("B", "L", "M", "S")) panel_dots(sc, "hot_sensitivity", names(LCOL), LCOL, LPCH, "Cluster cells hot")
# false hot cells with no clusters (% of cells)
nul <- data.frame(sc = c("C", "N", "C", "N"), k = c(20, 20, 40, 40),
                  lab = c("C", "N", "C", "N"))
plot(NA, xlim = c(0.6, 4.4), ylim = c(0, 6), xaxt = "n", bty = "n", xlab = "", ylab = "% of map hot")
abline(h = 0:6, col = LIGHT, lwd = 0.6); axis(1, at = 1:4, labels = FALSE)
text(1:4, -0.54, nul$lab, xpd = NA, cex = 0.85, col = DARK, adj = c(0.5, 1))
text(c(1.5, 3.5), -1.2, c("20\u00d720", "40\u00d740"), xpd = NA, cex = 0.85, col = DARK, adj = c(0.5, 1))
off <- seq(-0.26, 0.26, length.out = length(LCOL))
for (i in seq_along(LCOL)) points(1:4 + off[i], mapply(function(sc, k) 100 * get(k, sc, "none", names(LCOL)[i], "n_hot") / k^2,
                                                      nul$sc, nul$k), pch = LPCH[i], col = LCOL[i], cex = 0.9)
mtext("C, N: no clusters", side = 3, line = 0.4, adj = 0, font = 2, col = DARK, cex = 1)
panel_dots("K", "cold_precision", names(LCOL), LCOL, LPCH, "Cold cells in cold area")
panel_dots("K", "cold_sensitivity", names(LCOL), LCOL, LPCH, "Cold area labelled cold")
panel_dots("U", "hot_in_cov", names(LCOL), LCOL, LPCH, "Hot cells in covariate area")
legend_row(LCOL, LPCH)
dev.off()
message("written: 17_summary.csv, sim_ranking.tex, sim_labels.tex, fig15, fig16")
