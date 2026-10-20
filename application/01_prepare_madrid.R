# Madrid cockroach data: grid, adjacency, covariates and counts per summer ------------------
#
# Reproduces the data preparation of the original application (Codigo/50Location.Rmd):
#   - cockroach sightings ("Cucaracha"), April to September, 2010-2013
#   - window: union of the Madrid districts (small holes closed with a +/- 500 m buffer),
#     cut to the bounding box of the sightings plus 1 km
#   - square grid of 1000 m cells, cut to the window; cells with less than 7500 m2 dropped
#   - queen contiguity between cells (spdep::poly2nb)
#   - 19 covariates from TODO_MADRID119m.txt (point data), averaged per cell
#
# Input (application/data/): RC_Madrid_Year_Month_Week.csv, Distritos_20210712.*,
#   TODO_MADRID119m.txt (107 MB; not kept in git: copy it from 02SepTest/Codigo).
# Output: application/data/madrid_grid.rds
#
# Run from the repo root: Rscript application/01_prepare_madrid.R

suppressMessages({ library(sf); library(spdep) })
dd <- "application/data"

pts <- read.csv(file.path(dd, "RC_Madrid_Year_Month_Week.csv"), sep = ";")
pts <- pts[pts$Vector == "Cucaracha" & pts$Month %in% 4:9, ]

mot <- st_union(st_geometry(st_read(file.path(dd, "Distritos_20210712.shp"), quiet = TRUE)))
outer <- st_sfc(st_polygon(list(st_geometry(mot)[[1]][[1]])), crs = st_crs(mot))
mot <- st_buffer(st_buffer(outer, 500), -500)

r10 <- function(x) round(x / 10) * 10
bb  <- st_as_sfc(st_bbox(c(xmin = r10(min(pts$X)) - 1000, ymin = r10(min(pts$Y)) - 1000,
                           xmax = r10(max(pts$X)) + 1000, ymax = r10(max(pts$Y)) + 1000),
                         crs = st_crs(mot)))
window <- st_intersection(bb, mot)

grid <- st_sf(geometry = st_intersection(st_make_grid(bb, cellsize = c(1000, 1000)), window))
grid <- grid[as.numeric(st_area(grid)) >= 7500, ]
grid$cellid <- seq_len(nrow(grid))

nb <- poly2nb(grid, queen = TRUE)
A  <- nb2mat(nb, style = "B", zero.policy = TRUE); dimnames(A) <- NULL

# covariates: same 19 variables and the same filling as the original code
todo <- read.csv(file.path(dd, "TODO_MADRID119m.txt"), sep = ";", header = TRUE, dec = ",")
vars <- c("AreaVer", "D_AlcDef", "D_AlcRes", "D_CC", "D_Conten", "D_Desc", "D_Edif_rui",
          "D_Gatos", "D_JarBas", "D_JarMadr", "D_Mercados", "D_Perros", "D_RsuDefi", "D_Solar",
          "RPD", "YCons", "dist_agua", "dist_fuent", "dist_verd")
cov <- todo[, c("x", "y", vars)]
for (v in vars) {                                   # tidyr::fill(.direction = "down")
  x <- cov[[v]]; ok <- !is.na(x)
  cov[[v]] <- if (any(ok)) x[c(NA, which(ok))[cumsum(ok) + 1]] else x
}
cov_sf <- st_as_sf(cov, coords = c("x", "y"), crs = st_crs(grid))
cell_of_cov <- st_intersects(cov_sf, grid)
cid <- vapply(cell_of_cov, function(i) if (length(i)) i[1] else NA_integer_, 1L)
Z <- sapply(vars, function(v) tapply(cov[[v]], factor(cid, levels = grid$cellid), mean, na.rm = TRUE))
for (v in vars) {                                   # fill(.direction = "downup") for empty cells
  x <- Z[, v]; ok <- !is.na(x)
  if (any(!ok)) {
    dn <- x[c(NA, which(ok))[cumsum(ok) + 1]]
    up <- rev({ r <- rev(x); okr <- !is.na(r); r[c(NA, which(okr))[cumsum(okr) + 1]] })
    Z[, v] <- ifelse(is.na(dn), up, dn)
  }
}
Z <- Z[, apply(Z, 2, function(z) sd(z, na.rm = TRUE) > 0), drop = FALSE]

# counts per cell and year
pts_sf <- st_as_sf(pts, coords = c("X", "Y"), crs = st_crs(grid))
cell_of_pt <- vapply(st_intersects(pts_sf, grid), function(i) if (length(i)) i[1] else NA_integer_, 1L)
years <- 2010:2013
Y <- sapply(years, function(yr) tabulate(cell_of_pt[pts$Year == yr], nbins = nrow(grid)))
colnames(Y) <- years

saveRDS(list(grid = grid, window = window, A = A, Z = Z, Y = Y, years = years,
             n_points = table(pts$Year), n_outside = tapply(is.na(cell_of_pt), pts$Year, sum),
             vars = colnames(Z)),
        file.path(dd, "madrid_grid.rds"))
message("cells: ", nrow(grid), "; edges: ", sum(A) / 2, "; covariates: ", ncol(Z),
        "; sightings per summer: ", paste(colSums(Y), collapse = ", "),
        "; outside the grid: ", paste(tapply(is.na(cell_of_pt), pts$Year, sum), collapse = ", "))
