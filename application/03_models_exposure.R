# Madrid: nested models (exposure, + environment, + rat reports) ----------------------------
#
# Cockroach sightings, April-September. Main map: 2010-2012 pooled; 2013 kept for an
# out-of-sample check. Same method as the simulation study (log counts, Moran gate,
# ESF smooth eigenvectors + ridge + CV, 999 permutations).
#
#   M1 exposure        where people live and can report:
#                        log(1 + residential points in the cell)  [TODO_MADRID119m rows]
#                        population density (Densidad), log cell area (edge cells are cut)
#   M2 + environment   M1 + risk factors chosen a priori (Tamayo-Uria et al. 2013;
#                        dictionary section 10), PROVISIONAL:
#                        D_AlcSuc, D_AlcDef, D_RsuDefi, D_Mercados, D_GatoOK, dist_agua,
#                        dist_verd, AreaVer, YCons, Euros, Desempleo
#   M3 + rat reports   M2 + log(1 + rat reports in the cell, all months, same years):
#                        hot spots specific to cockroaches
#   Rats (comparison)  rat reports, all months, 2010-2012, model M2
#
# Covariates are cell means of the point data. Cells with no points (parks, empty land)
# get the mean of their neighbouring cells (repeated until all cells have a value).
#
# Run from the repo root: Rscript application/03_models_exposure.R   (about 5 min)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else
    if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(hotcoldESF)
  library(sf)
})
d <- readRDS("application/data/madrid_grid.rds")
grid <- d$grid; A <- d$A; n <- nrow(A)
OUT <- "application/results"; FIG <- "application/figures"
NPERM <- 999

# ---- covariates ------------------------------------------------------------------------
env_vars <- c("D_AlcSuc", "D_AlcDef", "D_RsuDefi", "D_Mercados", "D_GatoOK", "dist_agua",
              "dist_verd", "AreaVer", "YCons", "Euros", "Desempleo")
todo <- read.csv("application/data/TODO_MADRID119m.txt", sep = ";", dec = ",")
pts  <- st_as_sf(todo[, c("x", "y", "Densidad", env_vars)], coords = c("x", "y"), crs = st_crs(grid))
cid  <- vapply(st_intersects(pts, grid), function(i) if (length(i)) i[1] else NA_integer_, 1L)
n_res <- tabulate(cid[!is.na(cid)], nbins = n)
cellmean <- function(v) {
  m <- tapply(todo[[v]], factor(cid, levels = seq_len(n)), mean, na.rm = TRUE)
  m <- as.numeric(m)
  while (any(is.na(m))) {                               # neighbour mean for empty cells
    nb_mean <- as.numeric(A %*% ifelse(is.na(m), 0, m)) / pmax(as.numeric(A %*% !is.na(m)), 1)
    fill <- is.na(m) & as.numeric(A %*% !is.na(m)) > 0
    if (!any(fill)) break
    m[fill] <- nb_mean[fill]
  }
  m
}
EXP <- cbind(log_residential = log1p(n_res), Densidad = cellmean("Densidad"),
             log_area = log(as.numeric(st_area(grid))))
ENV <- sapply(env_vars, cellmean)

# ---- responses ----------------------------------------------------------------------------
p <- read.csv("application/data/RC_Madrid_Year_Month_Week.csv", sep = ";")
pc <- vapply(st_intersects(st_as_sf(p, coords = c("X", "Y"), crs = st_crs(grid)), grid),
             function(i) if (length(i)) i[1] else NA_integer_, 1L)
cnt <- function(sel) tabulate(pc[sel & !is.na(pc)], nbins = n)
summer <- p$Month %in% 4:9
cock_train <- cnt(p$Vector == "Cucaracha" & summer & p$Year %in% 2010:2012)
cock_test  <- cnt(p$Vector == "Cucaracha" & summer & p$Year == 2013)
rat_train  <- cnt(p$Vector == "Rata" & p$Year %in% 2010:2012)
RAT <- cbind(log_rat_reports = log1p(rat_train))

models <- list(
  "M1 exposure"           = list(Y = cock_train, Z = EXP),
  "M2 + environment"      = list(Y = cock_train, Z = cbind(EXP, ENV)),
  "M3 + rat reports"      = list(Y = cock_train, Z = cbind(EXP, ENV, RAT)),
  "Rats, M2 (comparison)" = list(Y = rat_train,  Z = cbind(EXP, ENV)))

fits <- lapply(names(models), function(nm) {
  m <- models[[nm]]
  set.seed(2026)
  h <- hot_cold_map(m$Y, m$Z, A, nperm = NPERM, penalty = "cv_min", mc_min = 0.25, fit_type = "ridge")
  X0 <- design0(n, scale(m$Z)); y <- as.numeric(scale(log1p(m$Y)))
  ols <- lm.fit(X0, y)
  r2 <- 1 - sum(ols$residuals^2) / sum((y - mean(y))^2)
  list(h = h, r2 = r2, coef = setNames(ols$coefficients[-1], colnames(m$Z)))
})
names(fits) <- names(models)

# ---- summaries ---------------------------------------------------------------------------
summ <- do.call(rbind, lapply(names(fits), function(nm) { f <- fits[[nm]]; Y <- models[[nm]]$Y
  data.frame(model = nm, reports = sum(Y), r2_covariates = f$r2, moran_z = f$h$z0, gate_p = f$h$p0,
             n_hot = sum(f$h$labels == "hot"), n_cold = sum(f$h$labels == "cold"),
             share_reports_in_hot = sum(Y[f$h$labels == "hot"]) / sum(Y)) }))
write.csv(summ, file.path(OUT, "madrid_models_summary.csv"), row.names = FALSE)

coefs <- do.call(rbind, lapply(names(fits), function(nm)
  data.frame(model = nm, covariate = names(fits[[nm]]$coef), std_coef = fits[[nm]]$coef)))
write.csv(coefs, file.path(OUT, "madrid_models_coefficients.csv"), row.names = FALSE)

hot <- sapply(fits, function(f) f$h$labels == "hot")
jacc <- outer(seq_len(ncol(hot)), seq_len(ncol(hot)), Vectorize(function(i, j)
  sum(hot[, i] & hot[, j]) / max(1, sum(hot[, i] | hot[, j]))))
dimnames(jacc) <- list(colnames(hot), colnames(hot))
write.csv(round(jacc, 2), file.path(OUT, "madrid_models_overlap.csv"))

# out-of-sample: cells ranked with 2010-2012, cockroach summer 2013
ex13 <- cock_test - fitted(glm(cock_test ~ ., data = data.frame(cock_test, EXP, ENV), family = poisson()))
oos <- do.call(rbind, lapply(c(25, 50), function(m) {
  r <- c(lapply(fits[1:3], function(f) rank_cells(f$h, "R2")),
         list("Top counts 2010-12" = order(-cock_train)))
  do.call(rbind, lapply(names(r), function(nm) { top <- r[[nm]][1:m]
    data.frame(m = m, ranking = nm, sightings_2013 = sum(cock_test[top]) / sum(cock_test),
               excess_2013 = sum(ex13[top]) / sum(pmax(ex13, 0))) }))
}))
write.csv(oos, file.path(OUT, "madrid_models_out_of_sample.csv"), row.names = FALSE)

# ---- figure ------------------------------------------------------------------------------
DARK <- "#3C4043"; GREY <- "#9AA0A6"; HOT <- "#D93025"; COLD <- "#2166AC"
geom <- st_geometry(grid)
map <- function(col, main, sub = NULL) {
  par(mar = c(0.3, 0.3, 3, 0.3))
  plot(geom, col = col, border = "#DADCE0", lwd = 0.3)
  plot(st_geometry(d$window), add = TRUE, border = GREY, lwd = 0.8)
  mtext(main, side = 3, line = 1.4, adj = 0, font = 2, col = DARK, cex = 0.8)
  if (!is.null(sub)) mtext(sub, side = 3, line = 0.3, adj = 0, col = GREY, cex = 0.68)
}
divcol <- function(f) { mx <- max(abs(f)); if (mx == 0) return("white")
  rgb(colorRamp(c(COLD, "white", HOT))((f / mx + 1) / 2), maxColorValue = 255) }
labcol <- function(l) ifelse(l == "hot", "#F4A79F", ifelse(l == "cold", "#9EC3E6", "white"))
seqcol <- function(v) rgb(colorRamp(c("white", DARK))(log1p(v) / log1p(max(v))), maxColorValue = 255)

png(file.path(FIG, "madrid_models.png"), 2600, 1500, res = 150, pointsize = 10)
par(mfrow = c(2, 5), oma = c(0, 0, 3.5, 0), family = "sans")
map(seqcol(cock_train), "Cockroach reports", sprintf("Apr-Sep 2010-2012: %d", sum(cock_train)))
for (nm in names(fits)[1:3]) { f <- fits[[nm]]
  map(divcol(f$h$fhat), paste0(nm, ": spatial component"),
      sprintf("covariates explain %.0f%%; Moran z = %.1f", 100 * f$r2, f$h$z0)) }
map(divcol(fits[[4]]$h$fhat), "Rats, M2: spatial component",
    sprintf("covariates explain %.0f%%; Moran z = %.1f", 100 * fits[[4]]$r2, fits[[4]]$h$z0))
map(seqcol(rat_train), "Rat reports", sprintf("all months 2010-2012: %d", sum(rat_train)))
for (nm in names(fits)) { f <- fits[[nm]]
  map(labcol(f$h$labels), paste0(nm, ": hot/cold"),
      sprintf("%d hot, %d cold; hot cells hold %.0f%% of reports", sum(f$h$labels == "hot"),
              sum(f$h$labels == "cold"), 100 * sum(models[[nm]]$Y[f$h$labels == "hot"]) / sum(models[[nm]]$Y))) }
mtext("Cockroach hot spots under three nested models, and the rat map for comparison",
      outer = TRUE, side = 3, line = 2, adj = 0.01, font = 2, col = DARK, cex = 1.05)
mtext("M1: exposure (residential points, density, cell area). M2: M1 + environmental risk factors (provisional set). M3: M2 + rat reports. Red: more reports than the model explains.",
      outer = TRUE, side = 3, line = 0.7, adj = 0.01, col = GREY, cex = 0.8)
dev.off()

print(summ, digits = 3); print(round(jacc, 2)); print(oos, digits = 2)
print(reshape(coefs, idvar = "covariate", timevar = "model", direction = "wide"), digits = 2)
