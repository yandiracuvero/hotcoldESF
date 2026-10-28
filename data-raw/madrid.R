# Build the datasets madrid, madrid_A and madrid_boundary ----------------------------------
# Run from the root of the paper repository (it reads application/data/), with the
# package source as working copy:  Rscript data-raw/madrid.R
# Needs sf and spdep (not package dependencies). Same 1 km grid, counts and cell
# covariates as application/05_final_application.R.
suppressMessages({ library(sf); library(spdep) })
APP <- Sys.getenv("APP_DATA", "application/data")
base   <- readRDS(file.path(APP, "madrid_grid.rds")); window <- base$window
cv     <- readRDS(file.path(APP, "covariates_clean.rds"))
cpt    <- st_as_sf(cv[, c("x", "y")], coords = c("x", "y"), crs = st_crs(window))
p      <- read.csv(file.path(APP, "RC_Madrid_Year_Month_Week.csv"), sep = ";")
ppt    <- st_as_sf(p, coords = c("X", "Y"), crs = st_crs(window))
size <- 1000
bb <- st_as_sfc(st_bbox(base$grid))
sq <- st_make_grid(bb, cellsize = c(size, size))
g  <- st_intersection(st_sf(sq_id = seq_along(sq), geometry = sq), window)
g  <- g[as.numeric(st_area(g)) >= 0.0075 * size^2, ]
n  <- nrow(g)
A  <- nb2mat(poly2nb(g, queen = TRUE), style = "B", zero.policy = TRUE); dimnames(A) <- NULL
first <- function(ii) vapply(ii, function(i) if (length(i)) i[1] else NA_integer_, 1L)
cid <- first(st_intersects(cpt, g))
cellmean <- function(v) {                       # mean over address points; empty cells from neighbours
  m <- as.numeric(tapply(cv[[v]], factor(cid, levels = seq_len(n)), mean, na.rm = TRUE))
  repeat { miss <- is.na(m); if (!any(miss)) break
    k <- as.numeric(A %*% !miss); fill <- miss & k > 0; if (!any(fill)) break
    m[fill] <- (as.numeric(A %*% ifelse(miss, 0, m)) / pmax(k, 1))[fill] }
  m[is.na(m)] <- mean(m, na.rm = TRUE)
  m
}
pc  <- first(st_intersects(ppt, g))
cnt <- function(sel) tabulate(pc[sel & !is.na(pc)], nbins = n)
summer <- p$Month %in% 4:9
centre <- st_coordinates(st_centroid(sq[g$sq_id]))      # centre of the full 1 km square
madrid <- data.frame(
  cell = seq_len(n), x = centre[, 1], y = centre[, 2],
  area = as.numeric(st_area(g)) / 1e6,
  cockroaches_2010 = cnt(p$Vector == "Cucaracha" & summer & p$Year == 2010),
  cockroaches_2011 = cnt(p$Vector == "Cucaracha" & summer & p$Year == 2011),
  cockroaches_2012 = cnt(p$Vector == "Cucaracha" & summer & p$Year == 2012),
  cockroaches_2013 = cnt(p$Vector == "Cucaracha" & summer & p$Year == 2013),
  rats_2010 = cnt(p$Vector == "Rata" & p$Year == 2010),
  rats_2011 = cnt(p$Vector == "Rata" & p$Year == 2011),
  rats_2012 = cnt(p$Vector == "Rata" & p$Year == 2012),
  rats_2013 = cnt(p$Vector == "Rata" & p$Year == 2013),
  residential = tabulate(cid[!is.na(cid)], nbins = n),
  density = cellmean("Densidad"),
  dist_green = cellmean("dist_verd"), green_area = cellmean("AreaVer"),
  year_built = cellmean("YCons"), dist_water = cellmean("dist_agua"),
  dist_fountain = cellmean("dist_fuent"), dist_cat = cellmean("D_GatoOK"),
  unemployment = cellmean("Desempleo"))
madrid_A <- A
bd <- st_coordinates(st_simplify(st_geometry(window), dTolerance = 50))
madrid_boundary <- data.frame(x = bd[, "X"], y = bd[, "Y"],
                              ring = as.integer(interaction(bd[, "L1"], bd[, ncol(bd)], drop = TRUE)))
save(madrid, file = "data/madrid.rda", compress = "xz")
save(madrid_A, file = "data/madrid_A.rda", compress = "xz")
save(madrid_boundary, file = "data/madrid_boundary.rda", compress = "xz")
message("n = ", n, "; reports (cockroaches, summers 2010-2012) = ",
        sum(madrid$cockroaches_2010 + madrid$cockroaches_2011 + madrid$cockroaches_2012))
