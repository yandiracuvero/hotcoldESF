# Madrid application: hot/cold maps, persistence, out-of-sample trap placement ----------------
#
# Proposed method (as in the simulation study): log(Y + 1); the 19 covariates as linear
# terms; Moran gate (5%); ESF with smooth eigenvectors, ridge, penalty by 10-fold CV;
# permutation band (999 permutations; hot above the 95% quantile, cold below the 5%).
# Benchmarks: Gi* and LISA on the residuals of log(Y + 1) on Z, with the same gate.
# Sensitivity: the same procedure with a natural spline (df 4) for each covariate.
#
# 1. Hot/cold map per summer (2010-2013)
# 2. Persistence: number of summers each cell is hot
# 3. Out-of-sample: choose m cells with the data of summer t; evaluate with summer t + 1
#    (2010->11, 11->12, 12->13). m = 25 cells (5%) and 50 cells (10%).
#      sightings captured   share of the t + 1 sightings that fall in the m cells
#      excess captured      share of the t + 1 excess (observed minus expected from a
#                           Poisson regression of the t + 1 counts on Z, positive part
#                           summed over all cells) that falls in the m cells
#    Interventions between summers are unknown; they affect all methods alike.
#
# Run from the repo root: Rscript application/02_madrid_analysis.R  (about 2 min)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else
    if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(hotcoldESF)
  library(sf)
})
d  <- readRDS("application/data/madrid_grid.rds")
A  <- d$A; Z <- d$Z; n <- nrow(A); years <- d$years
OUT <- "application/results"; FIG <- "application/figures"
dir.create(OUT, showWarnings = FALSE); dir.create(FIG, showWarnings = FALSE)
NPERM <- 999

# ---- fits per summer ---------------------------------------------------------------------
fit_year <- function(Y, spline_df = 0, seed = 1) {
  set.seed(seed)
  h  <- hot_cold_map(Y, Z, A, nperm = NPERM, penalty = "cv_min", mc_min = 0.25,
                     fit_type = "ridge", spline_df = spline_df)
  X0 <- design0(n, covariate_basis(Z, spline_df))
  res <- lm.fit(X0, as.numeric(scale(log1p(Y))))$residuals
  gate_ok <- isTRUE(h$p0 <= 0.05)
  gi <- local_gi_star(res, A)
  li <- local_moran(res, A)
  lab <- function(hot, cold) factor(ifelse(seq_len(n) %in% hot, "hot",
                                    ifelse(seq_len(n) %in% cold, "cold", "neutral")),
                                    levels = c("hot", "neutral", "cold"))
  none <- integer(0)
  list(h = h, res = res, gi = gi, li = li, gate_ok = gate_ok,
       labels = list(
         ESF = h$labels,
         `Gi*` = if (gate_ok) lab(which(gi > qnorm(0.95)), which(gi < qnorm(0.05))) else lab(none, none),
         LISA = if (gate_ok) lab(which(li$quad == "HH" & li$p < 0.05), which(li$quad == "LL" & li$p < 0.05)) else lab(none, none)),
       rank = list(
         ESF = rank_cells(h, "R2"), `Gi*` = order(-gi), LISA = order(-(li$quad == "HH"), -li$I),
         `Top residuals` = order(-res), `Top counts` = order(-Y)))
}
fits  <- lapply(seq_along(years), function(t) fit_year(d$Y[, t], 0, seed = years[t]))
fitsS <- lapply(seq_along(years), function(t) fit_year(d$Y[, t], 4, seed = years[t]))
names(fits) <- names(fitsS) <- years

# ---- 1. summary per summer -----------------------------------------------------------------
per_year <- do.call(rbind, lapply(c("linear", "spline"), function(cv) {
  ff <- if (cv == "linear") fits else fitsS
  do.call(rbind, lapply(seq_along(years), function(t) {
    f <- ff[[t]]
    do.call(rbind, lapply(names(f$labels), function(me) data.frame(
      covariates = cv, year = years[t], method = me, sightings = sum(d$Y[, t]),
      moran_z = f$h$z0, gate_p = f$h$p0,
      n_hot = sum(f$labels[[me]] == "hot"), n_cold = sum(f$labels[[me]] == "cold"),
      sightings_in_hot = sum(d$Y[f$labels[[me]] == "hot", t]) / sum(d$Y[, t]))))
  }))
}))
write.csv(per_year, file.path(OUT, "madrid_per_summer.csv"), row.names = FALSE)

# ---- 2. persistence ---------------------------------------------------------------------
persist <- sapply(c("ESF", "Gi*", "LISA"), function(me)
  rowSums(sapply(fits, function(f) f$labels[[me]] == "hot")))
persist_cold <- sapply(c("ESF", "Gi*", "LISA"), function(me)
  rowSums(sapply(fits, function(f) f$labels[[me]] == "cold")))
write.csv(data.frame(cellid = d$grid$cellid, hot_ = persist, cold_ = persist_cold),
          file.path(OUT, "madrid_persistence.csv"), row.names = FALSE)

# agreement between consecutive summers (share of hot cells in t that are hot in t + 1)
agree <- do.call(rbind, lapply(1:3, function(t) data.frame(
  from = years[t], to = years[t + 1], method = c("ESF", "Gi*", "LISA"),
  hot_again = sapply(c("ESF", "Gi*", "LISA"), function(me) {
    a <- fits[[t]]$labels[[me]] == "hot"; b <- fits[[t + 1]]$labels[[me]] == "hot"
    if (sum(a)) sum(a & b) / sum(a) else NA }))))
write.csv(agree, file.path(OUT, "madrid_hot_again.csv"), row.names = FALSE)

# ---- 3. out-of-sample ---------------------------------------------------------------------
excess_next <- function(Ynext) {
  g <- glm(Ynext ~ ., data = data.frame(Ynext = Ynext, Z), family = poisson())
  Ynext - fitted(g)
}
oos <- do.call(rbind, lapply(c("linear", "spline"), function(cv) {
  ff <- if (cv == "linear") fits else fitsS
  do.call(rbind, lapply(1:3, function(t) {
    Yn <- d$Y[, t + 1]; ex <- excess_next(Yn); pos <- sum(pmax(ex, 0))
    do.call(rbind, lapply(c(25, 50), function(m) {
      r <- ff[[t]]$rank
      rows <- do.call(rbind, lapply(names(r), function(me) {
        top <- r[[me]][1:m]
        data.frame(covariates = cv, from = years[t], to = years[t + 1], m = m, method = me,
                   sightings_captured = sum(Yn[top]) / sum(Yn),
                   excess_captured = sum(ex[top]) / pos)
      }))
      rbind(rows, data.frame(covariates = cv, from = years[t], to = years[t + 1], m = m,
                             method = "Random cells", sightings_captured = m / n,
                             excess_captured = 0))
    }))
  }))
}))
write.csv(oos, file.path(OUT, "madrid_out_of_sample.csv"), row.names = FALSE)

# LaTeX table: linear covariates, mean over the three pairs, plus each pair
meths <- c("ESF", "Gi*", "LISA", "Top residuals", "Top counts", "Random cells")
o <- oos[oos$covariates == "linear", ]
fmt <- function(v) formatC(100 * v + 0, format = "f", digits = 1)
tab <- c("\\begin{tabular}{lrrrrrr}", "\\toprule",
         " & \\multicolumn{3}{c}{Sightings captured (\\%)} & \\multicolumn{3}{c}{Excess captured (\\%)} \\\\",
         "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
         "Method & 2010$\\to$11 & 2011$\\to$12 & 2012$\\to$13 & 2010$\\to$11 & 2011$\\to$12 & 2012$\\to$13 \\\\")
for (m in c(25, 50)) {
  tab <- c(tab, "\\midrule", paste0("\\multicolumn{7}{l}{\\emph{", m, " cells (", round(100 * m / n), "\\% of the city)}} \\\\"))
  for (me in meths) {
    v <- o[o$m == m & o$method == me, ]
    v <- v[order(v$from), ]
    tab <- c(tab, paste0(gsub("\\*", "$^*$", me), " & ", paste(fmt(v$sightings_captured), collapse = " & "),
                         " & ", paste(fmt(v$excess_captured), collapse = " & "), " \\\\"))
  }
}
tab <- c(tab, "\\bottomrule", "\\end{tabular}")
dir.create("paper/tables", showWarnings = FALSE, recursive = TRUE)
writeLines(tab, "paper/tables/madrid_out_of_sample.tex")

# ---- figures -----------------------------------------------------------------------------
DARK <- "#3C4043"; GREY <- "#9AA0A6"
HOT <- "#D93025"; COLD <- "#2166AC"
geom <- st_geometry(d$grid)
seqcol <- function(v, top = "#3C4043") {
  r <- log1p(v) / log1p(max(v)); rgb(colorRamp(c("white", top))(r), maxColorValue = 255) }
labcol <- function(l) ifelse(l == "hot", "#F4A79F", ifelse(l == "cold", "#9EC3E6", "white"))
map <- function(col, main, sub = NULL) {
  par(mar = c(0.3, 0.3, 2.8, 0.3))
  plot(geom, col = col, border = "#DADCE0", lwd = 0.3)
  plot(st_geometry(d$window), add = TRUE, border = GREY, lwd = 0.8)
  mtext(main, side = 3, line = 1.3, adj = 0, font = 2, col = DARK, cex = 0.85)
  if (!is.null(sub)) mtext(sub, side = 3, line = 0.3, adj = 0, col = GREY, cex = 0.7)
}

png(file.path(FIG, "madrid_maps.png"), 2400, 2200, res = 150, pointsize = 10)
par(mfrow = c(4, 4), oma = c(0, 0, 3.5, 0), family = "sans")
for (t in seq_along(years)) map(seqcol(d$Y[, t]), paste("Sightings", years[t]),
                                paste0(sum(d$Y[, t]), " sightings; max ", max(d$Y[, t]), " per cell"))
for (t in seq_along(years)) {
  f <- fits[[t]]$h$fhat; mx <- max(abs(f))
  map(if (mx > 0) rgb(colorRamp(c(COLD, "white", HOT))((f / mx + 1) / 2), maxColorValue = 255) else "white",
      paste("ESF spatial component", years[t]),
      sprintf("Moran z = %.1f; red: more sightings than the covariates explain", fits[[t]]$h$z0))
}
for (t in seq_along(years)) map(labcol(fits[[t]]$labels$ESF), paste("ESF hot/cold", years[t]),
  sprintf("%d hot, %d cold cells", sum(fits[[t]]$labels$ESF == "hot"), sum(fits[[t]]$labels$ESF == "cold")))
for (t in seq_along(years)) map(labcol(fits[[t]]$labels$`Gi*`), paste("Gi* hot/cold", years[t]),
  sprintf("%d hot, %d cold cells", sum(fits[[t]]$labels$`Gi*` == "hot"), sum(fits[[t]]$labels$`Gi*` == "cold")))
mtext("Cockroach sightings in Madrid, April-September, 1 km cells", outer = TRUE, side = 3, line = 2,
      adj = 0.01, font = 2, col = DARK, cex = 1.1)
mtext("Rows: sightings; ESF spatial component (sightings not explained by the 19 covariates); ESF labels; Gi* labels (same residuals and gate)",
      outer = TRUE, side = 3, line = 0.6, adj = 0.01, col = GREY, cex = 0.85)
dev.off()

pcol <- c("white", "#FCE8E6", "#F4A79F", "#E8665A", "#A50E0E")
png(file.path(FIG, "madrid_persistence.png"), 2200, 900, res = 150, pointsize = 10)
par(mfrow = c(1, 3), oma = c(2.5, 0, 3.5, 0), family = "sans")
for (me in c("ESF", "Gi*", "LISA")) map(pcol[persist[, me] + 1], me,
  sprintf("%d cells hot in all 4 summers; %d in at least 2", sum(persist[, me] == 4), sum(persist[, me] >= 2)))
par(fig = c(0, 1, 0, 0.1), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = paste(0:4, c("summers", "summer", "summers", "summers", "summers")),
       fill = pcol, border = "#DADCE0", horiz = TRUE, bty = "n", cex = 1, text.col = DARK)
mtext("Persistence of hot spots: number of summers (2010-2013) each cell is labelled hot",
      outer = TRUE, side = 3, line = 2, adj = 0.01, font = 2, col = DARK, cex = 1.05)
dev.off()

OC <- c("ESF" = "#2a78d6", "Gi*" = "#eb6834", "LISA" = "#1baf7a", "Top residuals" = "#eda100",
        "Top counts" = "#e87ba4", "Random cells" = GREY)
png(file.path(FIG, "madrid_out_of_sample.png"), 2200, 900, res = 150, pointsize = 10)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 4, 1), oma = c(0, 0, 0, 11), family = "sans",
    col.axis = DARK, col.lab = DARK, fg = GREY, las = 1)
for (col in c("sightings_captured", "excess_captured")) {
  plot(NA, xlim = c(0.6, 3.4), ylim = c(0, if (col == "sightings_captured") 0.5 else 1), xaxt = "n", bty = "n",
       xlab = "Cells chosen with summer t, evaluated in summer t + 1", ylab = "Share captured (25 cells)")
  abline(h = pretty(c(0, 1)), col = "#E8EAED")
  axis(1, at = 1:3, labels = c("2010 -> 11", "2011 -> 12", "2012 -> 13"))
  off <- seq(-0.25, 0.25, length.out = length(OC))
  for (i in seq_along(OC)) {
    v <- o[o$m == 25 & o$method == names(OC)[i], ]; v <- v[order(v$from), col]
    points(1:3 + off[i], v, pch = 16, col = OC[i], cex = 1.4)
  }
  mtext(if (col == "sightings_captured") "Sightings captured" else "Excess captured (sightings above what the covariates predict)",
        side = 3, line = 1.5, adj = 0, font = 2, col = DARK)
}
par(xpd = NA)
legend(par("usr")[2] + 0.1, par("usr")[4], names(OC), pch = 16, col = OC, bty = "n", text.col = DARK)
dev.off()

message("done. Gate p-values: ", paste(sprintf("%.3g", sapply(fits, function(f) f$h$p0)), collapse = ", "))
