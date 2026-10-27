# Madrid: nested models with the documented covariates, and their relevance ----------------
#
# Cockroach sightings, April-September; main map 2010-2012 pooled; 2013 out-of-sample.
# Same method as the simulation study (log counts, Moran gate, ESF smooth eigenvectors +
# ridge + CV, 999 permutations).
#
#   M1 exposure        residential points in the cell (log), Densidad, log cell area
#   M2 + environment   M1 + dist_verd, AreaVer, YCons, dist_agua, dist_fuent,
#                      D_GatoOK (cat feeding points; D_Gatos as a variant), Desempleo
#                      (D_Mercados and Euros removed after the relevance check)
#   M3 + rat reports   M2 + log(1 + rat reports, all months, 2010-2012)
#   M2 + field         sensitivity: M2 + variables recorded in the vector-control field
#                      campaigns (D_AlcSuc, D_AlcDef, D_AlcRes, D_DefCons, D_Solar, D_Edif_rui)
#   Rats, M2           comparison
#
# Relevance of each environmental covariate (M2, cockroaches):
#   ESF-adjusted coefficient: OLS of log counts on the covariates plus the smooth
#     eigenvectors selected by Lasso-CV (the eigenvectors absorb the spatial
#     autocorrelation, so the standard errors are usable; plain OLS shown for contrast)
#   stability: the same coefficient in each summer 2010-2013 (sign and size)
#   drop-one: hot cells of M2 without the covariate vs full M2 (overlap of hot cells)
#
# Run from the repo root: Rscript application/04_models_covariates.R   (about 5 min)

suppressMessages({
  if (file.exists("DESCRIPTION")) devtools::load_all(".", quiet = TRUE) else
    if (file.exists("../DESCRIPTION")) devtools::load_all("..", quiet = TRUE) else library(hotcoldESF)
  library(sf); library(glmnet)
})
d <- readRDS("application/data/madrid_grid.rds")
grid <- d$grid; A <- d$A; n <- nrow(A)
OUT <- "application/results"; FIG <- "application/figures"
NPERM <- 999

# ---- covariates per cell ------------------------------------------------------------------
# D_Mercados and Euros were removed after the relevance check (step 12): not significant,
# sign not stable across summers, and the hot map barely changes without them.
env_vars   <- c("dist_verd", "AreaVer", "YCons", "dist_agua", "dist_fuent", "D_GatoOK", "Desempleo")
field_vars <- c("D_AlcSuc", "D_AlcDef", "D_AlcRes", "D_DefCons", "D_Solar", "D_Edif_rui")
cv <- readRDS("application/data/covariates_clean.rds")
pts <- st_as_sf(cv[, c("x", "y")], coords = c("x", "y"), crs = st_crs(grid))
cid <- vapply(st_intersects(pts, grid), function(i) if (length(i)) i[1] else NA_integer_, 1L)
n_res <- tabulate(cid[!is.na(cid)], nbins = n)
cellmean <- function(v) {
  m <- as.numeric(tapply(cv[[v]], factor(cid, levels = seq_len(n)), mean, na.rm = TRUE))
  repeat {                                              # empty cells: mean of neighbours
    miss <- is.na(m); if (!any(miss)) break
    k <- as.numeric(A %*% !miss); fill <- miss & k > 0
    if (!any(fill)) break
    m[fill] <- (as.numeric(A %*% ifelse(miss, 0, m)) / pmax(k, 1))[fill]
  }
  m
}
EXP <- cbind(log_residential = log1p(n_res), Densidad = cellmean("Densidad"),
             log_area = log(as.numeric(st_area(grid))))
ENV <- sapply(env_vars, cellmean)
ENV_gatos <- ENV; ENV_gatos[, "D_GatoOK"] <- cellmean("D_Gatos")
colnames(ENV_gatos)[colnames(ENV_gatos) == "D_GatoOK"] <- "D_Gatos"
FIELD <- sapply(field_vars, cellmean)

# ---- responses ----------------------------------------------------------------------------
p <- read.csv("application/data/RC_Madrid_Year_Month_Week.csv", sep = ";")
pc <- vapply(st_intersects(st_as_sf(p, coords = c("X", "Y"), crs = st_crs(grid)), grid),
             function(i) if (length(i)) i[1] else NA_integer_, 1L)
cnt <- function(sel) tabulate(pc[sel & !is.na(pc)], nbins = n)
summer <- p$Month %in% 4:9
cock <- function(yrs) cnt(p$Vector == "Cucaracha" & summer & p$Year %in% yrs)
cock_train <- cock(2010:2012); cock_test <- cock(2013)
rat_train  <- cnt(p$Vector == "Rata" & p$Year %in% 2010:2012)
RAT <- cbind(log_rat_reports = log1p(rat_train))

models <- list(
  "M1 exposure"            = list(Y = cock_train, Z = EXP),
  "M2 + environment"       = list(Y = cock_train, Z = cbind(EXP, ENV)),
  "M2 (D_Gatos variant)"   = list(Y = cock_train, Z = cbind(EXP, ENV_gatos)),
  "M3 + rat reports"       = list(Y = cock_train, Z = cbind(EXP, ENV, RAT)),
  "M2 + field variables"   = list(Y = cock_train, Z = cbind(EXP, ENV, FIELD)),
  "Rats, M2 (comparison)"  = list(Y = rat_train,  Z = cbind(EXP, ENV)))

run <- function(Y, Z, seed = 2026) {
  set.seed(seed)
  hot_cold_map(Y, Z, A, nperm = NPERM, penalty = "cv_min", mc_min = 0.25, fit_type = "ridge")
}
r2_of <- function(Y, Z) { y <- as.numeric(scale(log1p(Y))); X0 <- design0(n, scale(Z))
  e <- lm.fit(X0, y)$residuals; 1 - sum(e^2) / sum((y - mean(y))^2) }
fits <- lapply(models, function(m) run(m$Y, m$Z))

summ <- do.call(rbind, lapply(names(fits), function(nm) { h <- fits[[nm]]; Y <- models[[nm]]$Y
  data.frame(model = nm, n_covariates = ncol(models[[nm]]$Z), r2_covariates = r2_of(Y, models[[nm]]$Z),
             moran_z = h$z0, gate_p = h$p0, n_hot = sum(h$labels == "hot"), n_cold = sum(h$labels == "cold"),
             share_reports_in_hot = sum(Y[h$labels == "hot"]) / sum(Y)) }))
write.csv(summ, file.path(OUT, "madrid_cov_models_summary.csv"), row.names = FALSE)

hot <- sapply(fits, function(h) h$labels == "hot")
jacc <- function(a, b) sum(a & b) / max(1, sum(a | b))
overlap <- outer(seq_len(ncol(hot)), seq_len(ncol(hot)), Vectorize(function(i, j) jacc(hot[, i], hot[, j])))
dimnames(overlap) <- list(colnames(hot), colnames(hot))
write.csv(round(overlap, 2), file.path(OUT, "madrid_cov_models_overlap.csv"))

# ---- relevance of the environmental covariates (M2, cockroaches) ---------------------------
eig <- eigen(A, symmetric = TRUE)
smooth <- which(eig$values >= 0.25 * eig$values[1])
esf_coef <- function(Y, Z, seed = 1) {
  y <- as.numeric(scale(log1p(Y))); Zs <- scale(Z); E <- eig$vectors[, smooth]
  set.seed(seed)
  cvf <- cv.glmnet(cbind(Zs, E), y, penalty.factor = c(rep(0, ncol(Zs)), rep(1, ncol(E))))
  b <- as.numeric(coef(cvf, s = "lambda.min"))[-1]
  S <- which(b[-(1:ncol(Zs))] != 0)
  fit_esf <- lm(y ~ Zs + E[, S, drop = FALSE])
  fit_ols <- lm(y ~ Zs)
  ce <- summary(fit_esf)$coefficients[1 + seq_len(ncol(Zs)), , drop = FALSE]
  co <- summary(fit_ols)$coefficients[1 + seq_len(ncol(Zs)), , drop = FALSE]
  data.frame(covariate = colnames(Z), coef_esf = ce[, 1], se_esf = ce[, 2], p_esf = ce[, 4],
             coef_ols = co[, 1], p_ols = co[, 4], n_eigenvectors = length(S))
}
Z2 <- cbind(EXP, ENV)
rel <- esf_coef(cock_train, Z2)
by_year <- sapply(2010:2013, function(yr) esf_coef(cock(yr), Z2, seed = yr)$coef_esf)
colnames(by_year) <- paste0("coef_", 2010:2013)
rel <- cbind(rel, by_year,
             same_sign_all_years = apply(by_year, 1, function(b) abs(sum(sign(b))) == 4))
rel_rat <- esf_coef(rat_train, Z2)
rel$coef_esf_rats <- rel_rat$coef_esf; rel$p_esf_rats <- rel_rat$p_esf

# drop-one: hot cells without each environmental covariate
full <- fits[["M2 + environment"]]$labels == "hot"
drop1 <- do.call(rbind, lapply(env_vars, function(v) {
  h <- run(cock_train, Z2[, colnames(Z2) != v])
  data.frame(covariate = v, n_hot_without = sum(h$labels == "hot"),
             overlap_with_full = jacc(h$labels == "hot", full))
}))
rel <- merge(rel, drop1, by = "covariate", all.x = TRUE, sort = FALSE)
write.csv(rel, file.path(OUT, "madrid_cov_relevance.csv"), row.names = FALSE)

# out-of-sample (2013)
ex13 <- cock_test - fitted(glm(cock_test ~ ., data = data.frame(cock_test, EXP, ENV), family = poisson()))
oos <- do.call(rbind, lapply(c(25, 50), function(m) {
  r <- c(lapply(fits[c("M1 exposure", "M2 + environment", "M3 + rat reports", "M2 + field variables")],
                function(h) rank_cells(h, "R2")), list("Top counts 2010-12" = order(-cock_train)))
  do.call(rbind, lapply(names(r), function(nm) { top <- r[[nm]][1:m]
    data.frame(m = m, ranking = nm, sightings_2013 = sum(cock_test[top]) / sum(cock_test),
               excess_2013 = sum(ex13[top]) / sum(pmax(ex13, 0))) }))
}))
write.csv(oos, file.path(OUT, "madrid_cov_out_of_sample.csv"), row.names = FALSE)

# ---- figures ------------------------------------------------------------------------------
DARK <- "#3C4043"; GREY <- "#9AA0A6"; HOT <- "#D93025"; COLD <- "#2166AC"; LIGHT <- "#E8EAED"
geom <- st_geometry(grid)
map <- function(col, main, sub = NULL) {
  nl <- if (is.null(sub)) 0 else length(strsplit(sub, "\n")[[1]])
  par(mar = c(0.2, 0.2, 1.5 + nl * 1.05, 0.2))
  plot(geom, col = col, border = "#DADCE0", lwd = 0.2)
  plot(st_geometry(d$window), add = TRUE, border = GREY, lwd = 0.6)
  mtext(main, side = 3, line = 0.3 + nl * 1.05, adj = 0, font = 2, col = DARK, cex = 1)
  if (!is.null(sub)) mtext(sub, side = 3, line = 0.3, adj = 0, col = DARK, cex = 0.9, padj = 0)
}
labcol <- function(l) ifelse(l == "hot", "#F4A79F", ifelse(l == "cold", "#9EC3E6", "white"))
show <- c("M1 exposure", "M2 + environment", "M3 + rat reports", "M2 + field variables", "Rats, M2 (comparison)")
# figures are drawn at their printed size (17 cm wide, 8 pt text)
png(file.path(FIG, "madrid_cov_models.png"), width = 6.69, height = 2.2, units = "in", res = 300, pointsize = 8, type = "cairo")
par(mfrow = c(1, 5), oma = c(0, 0, 0, 0), family = "sans"); par(cex = 1)
for (nm in show) { h <- fits[[nm]]; s <- summ[summ$model == nm, ]
  map(labcol(h$labels), sub(" (comparison)", "", nm, fixed = TRUE), sprintf("explained: %.0f%%\n%d hot, %d cold\nhot hold %.0f%% of reports",
      100 * s$r2_covariates, s$n_hot, s$n_cold, 100 * s$share_reports_in_hot)) }
dev.off()

# coefficient plot: ESF-adjusted coefficient +/- 2 SE, cockroaches; by-year points
png(file.path(FIG, "madrid_cov_relevance.png"), width = 5.69, height = 3.3, units = "in", res = 300, pointsize = 8, type = "cairo")
par(mar = c(3.4, 7, 3.4, 1), family = "sans", col.axis = DARK, col.lab = DARK, fg = GREY, las = 1,
    mgp = c(2.1, 0.5, 0), tcl = -0.2, cex.axis = 0.9)
r <- rel[match(colnames(Z2), rel$covariate), ]; k <- nrow(r); yy <- rev(seq_len(k))
xl <- range(c(r$coef_esf - 2 * r$se_esf, r$coef_esf + 2 * r$se_esf, as.matrix(r[, paste0("coef_", 2010:2013)])))
plot(NA, xlim = xl, ylim = c(0.5, k + 0.5), yaxt = "n", bty = "n",
     xlab = "Standardised coefficient (log reports per SD of the covariate)", ylab = "")
abline(v = 0, col = DARK); abline(h = yy, col = LIGHT)
axis(2, at = yy, labels = r$covariate, tick = FALSE)
for (yr in 2010:2013) points(r[[paste0("coef_", yr)]], yy, pch = 1, col = GREY, cex = 0.9)
segments(r$coef_esf - 2 * r$se_esf, yy, r$coef_esf + 2 * r$se_esf, yy, col = "#2a78d6", lwd = 1.5)
points(r$coef_esf, yy, pch = 16, col = "#2a78d6", cex = 1.1)
mtext("Covariate effects on cockroach reports (model M2)", side = 3, line = 1.9, adj = 0, font = 2, col = DARK)
mtext("Blue: 2010-2012 pooled, ESF-adjusted, +/- 2 SE. Grey circles: each summer separately.",
      side = 3, line = 0.7, adj = 0, col = DARK, cex = 0.9)
dev.off()

print(summ, digits = 3); print(round(overlap, 2))
print(rel[, c("covariate", "coef_esf", "se_esf", "p_esf", "p_ols", paste0("coef_", 2010:2013),
              "same_sign_all_years", "coef_esf_rats", "p_esf_rats", "n_hot_without", "overlap_with_full")], digits = 2)
print(oos, digits = 2)
