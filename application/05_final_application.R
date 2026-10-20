# Madrid: final application with the final covariate set --------------------------------------
#
# Final models (cockroaches April-September; rats all months):
#   M1 exposure      log(1 + residential points), Densidad, log cell area
#   M2 environment   M1 + dist_verd, AreaVer, YCons, dist_agua, dist_fuent, D_GatoOK, Desempleo
#   M3 rat reports   M2 + log(1 + rat reports in the same period)
# Method: log counts, Moran gate, ESF smooth eigenvectors + ridge + CV, 999 permutations.
#
# 1. 1 km grid, each summer 2010-2013: hot/cold maps (M1, M2, M3), persistence
# 2. 1 km grid, out-of-sample: cells chosen with summer t, evaluated with summer t + 1
#    (and 2010-2012 pooled -> 2013); ESF M1-M3, Gi* and LISA on M2 residuals (with the
#    gate), top residuals (M2), top counts, random cells
# 3. Grid size: 500 m, 1 km, 2 km; 2010-2012 pooled; hot cells compared on a common
#    250 m reference lattice (a 250 m cell is hot for a grid if its centre lies in a hot cell)
#
# Run from the repo root: Rscript application/05_final_application.R   (about 30 min)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else library(hotspots)
  library(sf); library(spdep)
})
OUT <- "application/results"; FIG <- "application/figures"; TAB <- "paper/tables"
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
NPERM <- 999
base <- readRDS("application/data/madrid_grid.rds")
window <- base$window

env_vars <- c("dist_verd", "AreaVer", "YCons", "dist_agua", "dist_fuent", "D_GatoOK", "Desempleo")
cv  <- readRDS("application/data/covariates_clean.rds")
cpt <- st_as_sf(cv[, c("x", "y")], coords = c("x", "y"), crs = st_crs(window))
p   <- read.csv("application/data/RC_Madrid_Year_Month_Week.csv", sep = ";")
ppt <- st_as_sf(p, coords = c("X", "Y"), crs = st_crs(window))
summer <- p$Month %in% 4:9

# ---- build a grid with its adjacency, covariates and counts --------------------------------
make_grid <- function(size) {
  bb <- st_as_sfc(st_bbox(base$grid))
  g <- st_sf(geometry = st_intersection(st_make_grid(bb, cellsize = c(size, size)), window))
  g <- g[as.numeric(st_area(g)) >= 0.0075 * size^2, ]
  n <- nrow(g)
  A <- nb2mat(poly2nb(g, queen = TRUE), style = "B", zero.policy = TRUE); dimnames(A) <- NULL
  first <- function(ii) vapply(ii, function(i) if (length(i)) i[1] else NA_integer_, 1L)
  cid <- first(st_intersects(cpt, g))
  cellmean <- function(v) {
    m <- as.numeric(tapply(cv[[v]], factor(cid, levels = seq_len(n)), mean, na.rm = TRUE))
    repeat { miss <- is.na(m); if (!any(miss)) break
      k <- as.numeric(A %*% !miss); fill <- miss & k > 0; if (!any(fill)) break
      m[fill] <- (as.numeric(A %*% ifelse(miss, 0, m)) / pmax(k, 1))[fill] }
    m[is.na(m)] <- mean(m, na.rm = TRUE)
    m
  }
  EXP <- cbind(log_residential = log1p(tabulate(cid[!is.na(cid)], nbins = n)),
               Densidad = cellmean("Densidad"), log_area = log(as.numeric(st_area(g))))
  ENV <- sapply(env_vars, cellmean)
  pc <- first(st_intersects(ppt, g))
  cnt <- function(sel) tabulate(pc[sel & !is.na(pc)], nbins = n)
  list(size = size, g = g, n = n, A = A, EXP = EXP, ENV = ENV,
       cock = function(yrs) cnt(p$Vector == "Cucaracha" & summer & p$Year %in% yrs),
       rat  = function(yrs) cnt(p$Vector == "Rata" & p$Year %in% yrs))
}
design <- function(G, model, yrs) {
  Z <- cbind(G$EXP, if (model != "M1") G$ENV)
  if (model == "M3") Z <- cbind(Z, log_rat_reports = log1p(G$rat(yrs)))
  Z
}
fit <- function(G, model, yrs, seed) {
  set.seed(seed)
  hot_cold_map(G$cock(yrs), design(G, model, yrs), G$A, nperm = NPERM,
               penalty = "cv_min", mc_min = 0.25, fit_type = "ridge")
}
models <- c("M1", "M2", "M3")
years <- 2010:2013

# ---- 1. 1 km grid: each summer, pooled, persistence --------------------------------------------
G1 <- make_grid(1000)
per <- list()
for (m in models) for (yr in years) per[[paste(m, yr)]] <- fit(G1, m, yr, seed = yr)
pooled <- lapply(setNames(models, models), function(m) fit(G1, m, 2010:2012, seed = 2026))
lab <- function(h, what = "hot") h$labels == what

summ <- do.call(rbind, lapply(names(per), function(k) { h <- per[[k]]; mm <- sub(" .*", "", k); yr <- as.integer(sub(".* ", "", k))
  Y <- G1$cock(yr)
  data.frame(model = mm, period = yr, reports = sum(Y), moran_z = h$z0, gate_p = h$p0,
             n_hot = sum(lab(h)), n_cold = sum(lab(h, "cold")), share_reports_in_hot = sum(Y[lab(h)]) / sum(Y)) }))
summ <- rbind(summ, do.call(rbind, lapply(models, function(m) { h <- pooled[[m]]; Y <- G1$cock(2010:2012)
  data.frame(model = m, period = "2010-2012", reports = sum(Y), moran_z = h$z0, gate_p = h$p0,
             n_hot = sum(lab(h)), n_cold = sum(lab(h, "cold")), share_reports_in_hot = sum(Y[lab(h)]) / sum(Y)) })))
write.csv(summ, file.path(OUT, "final_maps_summary.csv"), row.names = FALSE)

persist <- sapply(models, function(m) rowSums(sapply(years, function(yr) lab(per[[paste(m, yr)]]))))
hot_again <- do.call(rbind, lapply(models, function(m) do.call(rbind, lapply(1:3, function(t) {
  a <- lab(per[[paste(m, years[t])]]); b <- lab(per[[paste(m, years[t + 1])]])
  data.frame(model = m, from = years[t], to = years[t + 1], hot_again = if (sum(a)) sum(a & b) / sum(a) else NA) }))))
write.csv(hot_again, file.path(OUT, "final_hot_again.csv"), row.names = FALSE)
write.csv(data.frame(cellid = seq_len(G1$n), persist), file.path(OUT, "final_persistence.csv"), row.names = FALSE)

# ---- 2. out-of-sample ------------------------------------------------------------------------
excess <- function(Y) Y - fitted(glm(Y ~ ., data = data.frame(Y, G1$EXP, G1$ENV), family = poisson()))
rankings <- function(train_yrs, fits3) {
  Y <- G1$cock(train_yrs); Z2 <- design(G1, "M2", train_yrs)
  res <- lm.fit(design0(G1$n, scale(Z2)), as.numeric(scale(log1p(Y))))$residuals
  gi <- local_gi_star(res, G1$A); set.seed(1); li <- local_moran(res, G1$A)
  c(lapply(setNames(fits3, paste("ESF", models)), function(h) rank_cells(h, "R2")),
    list("Gi* (M2 residuals)" = order(-gi), "LISA (M2 residuals)" = order(-(li$quad == "HH"), -li$I),
         "Top residuals (M2)" = order(-res), "Top counts" = order(-Y)))
}
pairs <- list(list(train = 2010, test = 2011), list(train = 2011, test = 2012),
              list(train = 2012, test = 2013), list(train = 2010:2012, test = 2013))
oos <- do.call(rbind, lapply(pairs, function(pr) {
  fits3 <- if (length(pr$train) == 1) lapply(models, function(m) per[[paste(m, pr$train)]]) else pooled
  r <- rankings(pr$train, fits3)
  Yn <- G1$cock(pr$test); ex <- excess(Yn); pos <- sum(pmax(ex, 0))
  do.call(rbind, lapply(c(25, 50), function(m) rbind(
    do.call(rbind, lapply(names(r), function(nm) { top <- r[[nm]][1:m]
      data.frame(train = paste(range(pr$train), collapse = "-"), test = pr$test, m = m, ranking = nm,
                 sightings = sum(Yn[top]) / sum(Yn), excess = sum(ex[top]) / pos) })),
    data.frame(train = paste(range(pr$train), collapse = "-"), test = pr$test, m = m, ranking = "Random cells",
               sightings = m / G1$n, excess = 0))))
}))
oos$train <- sub("^(\\d+)-\\1$", "\\1", oos$train)
write.csv(oos, file.path(OUT, "final_out_of_sample.csv"), row.names = FALSE)

rk <- c("ESF M1", "ESF M2", "ESF M3", "Gi* (M2 residuals)", "LISA (M2 residuals)",
        "Top residuals (M2)", "Top counts", "Random cells")
rk_tex <- c("ESF, M1 exposure", "ESF, M2 environment", "ESF, M3 rat reports", "Gi$^*$ (M2 residuals)",
            "LISA (M2 residuals)", "Top residuals (M2)", "Top counts", "Random cells")
cols <- c("2010/2011", "2011/2012", "2012/2013", "2010-2012/2013")
f1 <- function(v) formatC(100 * v + 0, format = "f", digits = 1)
tab <- c("\\begin{tabular}{lrrrrrrrr}", "\\toprule",
         " & \\multicolumn{4}{c}{Sightings captured (\\%)} & \\multicolumn{4}{c}{Excess captured (\\%)} \\\\",
         "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
         paste0("Ranking & ", paste(rep(gsub("/", "$\\\\to$", cols), 2), collapse = " & "), " \\\\"))
for (m in c(25, 50)) {
  tab <- c(tab, "\\midrule", paste0("\\multicolumn{9}{l}{\\emph{", m, " cells (", round(100 * m / G1$n), "\\% of the city)}} \\\\"))
  for (i in seq_along(rk)) {
    v <- oos[oos$m == m & oos$ranking == rk[i], ]
    v <- v[match(c("2010", "2011", "2012", "2010-2012"), v$train), ]
    tab <- c(tab, paste0(rk_tex[i], " & ", paste(f1(v$sightings), collapse = " & "), " & ",
                         paste(f1(v$excess), collapse = " & "), " \\\\"))
  }
}
writeLines(c(tab, "\\bottomrule", "\\end{tabular}"), file.path(TAB, "madrid_final_out_of_sample.tex"))

# ---- 3. grid size ----------------------------------------------------------------------------
grids <- list("500 m" = make_grid(500), "1 km" = G1, "2 km" = make_grid(2000))
gfit <- lapply(grids, function(G) lapply(setNames(models, models), function(m)
  if (G$size == 1000) pooled[[m]] else fit(G, m, 2010:2012, seed = 2026)))
ref <- st_sf(geometry = st_intersection(st_make_grid(st_as_sfc(st_bbox(base$grid)), cellsize = c(250, 250)), window))
ref <- ref[as.numeric(st_area(ref)) >= 0.5 * 250^2, ]
refc <- suppressWarnings(st_centroid(ref))
on_ref <- function(G, h) { i <- vapply(st_intersects(refc, G$g), function(x) if (length(x)) x[1] else NA_integer_, 1L)
  out <- h$labels[i] == "hot"; out[is.na(out)] <- FALSE; out }
hotref <- lapply(models, function(m) sapply(names(grids), function(gn) on_ref(grids[[gn]], gfit[[gn]][[m]])))
names(hotref) <- models
jac <- function(a, b) sum(a & b) / max(1, sum(a | b))
grid_cmp <- do.call(rbind, lapply(models, function(m) {
  H <- hotref[[m]]
  data.frame(model = m,
             share_city_hot_500m = mean(H[, 1]), share_city_hot_1km = mean(H[, 2]), share_city_hot_2km = mean(H[, 3]),
             overlap_500m_1km = jac(H[, 1], H[, 2]), overlap_1km_2km = jac(H[, 2], H[, 3]),
             overlap_500m_2km = jac(H[, 1], H[, 3]), hot_in_all_three = mean(rowSums(H) == 3) / max(mean(rowSums(H) > 0), 1e-9))
}))
grid_stats <- do.call(rbind, lapply(names(grids), function(gn) { G <- grids[[gn]]
  do.call(rbind, lapply(models, function(m) { h <- gfit[[gn]][[m]]; Y <- G$cock(2010:2012)
    data.frame(grid = gn, cells = G$n, zero_cells = mean(Y == 0), model = m, moran_z = h$z0, gate_p = h$p0,
               n_hot = sum(lab(h)), share_reports_in_hot = sum(Y[lab(h)]) / sum(Y)) })) }))
write.csv(grid_cmp, file.path(OUT, "final_grid_overlap.csv"), row.names = FALSE)
write.csv(grid_stats, file.path(OUT, "final_grid_stats.csv"), row.names = FALSE)

# ---- figures ---------------------------------------------------------------------------------
DARK <- "#3C4043"; GREY <- "#9AA0A6"
labcol <- function(l) ifelse(l == "hot", "#F4A79F", ifelse(l == "cold", "#9EC3E6", "white"))
map <- function(geom, col, main, sub = NULL) {
  par(mar = c(0.3, 0.3, 3, 0.3))
  plot(geom, col = col, border = "#DADCE0", lwd = 0.25)
  plot(st_geometry(window), add = TRUE, border = GREY, lwd = 0.8)
  mtext(main, side = 3, line = 1.4, adj = 0, font = 2, col = DARK, cex = 0.8)
  if (!is.null(sub)) mtext(sub, side = 3, line = 0.3, adj = 0, col = GREY, cex = 0.68)
}
g1 <- st_geometry(G1$g)
mlab <- c(M1 = "M1 exposure", M2 = "M2 environment", M3 = "M3 rat reports")

# per-summer maps and persistence (M2, M3)
pcol <- c("white", "#FCE8E6", "#F4A79F", "#E8665A", "#A50E0E")
png(file.path(FIG, "final_summers.png"), 2600, 1500, res = 150, pointsize = 10)
par(mfrow = c(2, 5), oma = c(2.5, 0, 3.2, 0), family = "sans")
for (m in c("M2", "M3")) {
  for (yr in years) { h <- per[[paste(m, yr)]]
    map(g1, labcol(h$labels), paste(mlab[m], yr), sprintf("%d hot, %d cold", sum(lab(h)), sum(lab(h, "cold")))) }
  map(g1, pcol[persist[, m] + 1], paste(mlab[m], ": persistence"),
      sprintf("%d cells hot in all 4 summers; %d in 2 or more", sum(persist[, m] == 4), sum(persist[, m] >= 2)))
}
par(fig = c(0, 1, 0, 0.06), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
legend("center", c("hot", "cold", paste(1:4, "summers hot")), fill = c("#F4A79F", "#9EC3E6", pcol[2:5]),
       border = "#DADCE0", horiz = TRUE, bty = "n", text.col = DARK)
mtext("Cockroach hot and cold cells per summer (1 km cells), and persistence over 2010-2013",
      outer = TRUE, side = 3, line = 1.5, adj = 0.01, font = 2, col = DARK)
dev.off()

# grid sizes (pooled 2010-2012)
png(file.path(FIG, "final_grid_sizes.png"), 2200, 2000, res = 150, pointsize = 10)
par(mfrow = c(3, 3), oma = c(0, 0, 3.2, 0), family = "sans")
for (m in models) for (gn in names(grids)) { h <- gfit[[gn]][[m]]; G <- grids[[gn]]
  st <- grid_stats[grid_stats$grid == gn & grid_stats$model == m, ]
  map(st_geometry(G$g), labcol(h$labels), paste0(mlab[m], ", ", gn),
      sprintf("%d cells (%.0f%% without reports); Moran z = %.1f; %d hot", G$n, 100 * st$zero_cells, h$z0, sum(lab(h)))) }
mtext("Hot and cold cells by grid size; cockroaches April-September 2010-2012",
      outer = TRUE, side = 3, line = 1.5, adj = 0.01, font = 2, col = DARK)
dev.off()

print(summ, digits = 3); print(hot_again, digits = 2); print(grid_stats, digits = 3); print(grid_cmp, digits = 2)
print(aggregate(cbind(sightings, excess) ~ ranking + m, oos, mean), digits = 2)
