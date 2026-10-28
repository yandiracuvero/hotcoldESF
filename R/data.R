#' Cockroach and rat reports in Madrid, on a 1 km grid
#'
#' Confirmed sightings of cockroaches and rats reported by citizens to the
#' Technical Unit for Vector Control of the municipality of Madrid (Spain),
#' 2010--2013, counted on a grid of 1 km cells, with the covariates of the
#' application in the package's paper. Only counts per cell are included, not
#' the locations of the reports.
#'
#' The study region is the union of the Madrid districts, cut to the bounding
#' box of the sightings plus 1 km, divided into 1 km squares; cells with less
#' than 0.75% of their area inside the region are dropped. Cockroach counts
#' cover April to September of each year (87% of the cockroach reports); rat
#' counts cover the whole year.
#'
#' The covariates come from a database that records, for each residential
#' address point of the city (about 120,000), environmental, building and census
#' variables at that location (Tamayo-Uría et al., 2013, 2014). The value of a
#' covariate in a cell is the mean over its address points. The 162 cells
#' without residential address points take the mean of their neighbours
#' (repeatedly, until every cell has a value) and have `residential = 0`.
#'
#' @format A data frame with 493 rows (cells) and 21 columns:
#' \describe{
#'   \item{cell}{cell number; rows and columns of [madrid_A] follow this order.}
#'   \item{x, y}{centre of the 1 km square, in metres (ETRS89 / UTM zone 30N, EPSG:25830).}
#'   \item{area}{area of the cell inside the study region, in km^2 (border cells are smaller).}
#'   \item{cockroaches_2010, cockroaches_2011, cockroaches_2012, cockroaches_2013}{cockroach
#'     reports in the cell, April--September of each year.}
#'   \item{rats_2010, rats_2011, rats_2012, rats_2013}{rat reports in the cell, January--December.}
#'   \item{residential}{number of residential address points in the cell.}
#'   \item{density}{population density of the census areas of the address points.}
#'   \item{dist_green}{distance to the nearest green area, in metres.}
#'   \item{green_area}{surface of the nearest green area, in m^2.}
#'   \item{year_built}{year of construction of the building.}
#'   \item{dist_water}{distance to the nearest body of water, in metres.}
#'   \item{dist_fountain}{distance to the nearest drinking fountain, in metres.}
#'   \item{dist_cat}{distance to the nearest cat-feeding point, in metres.}
#'   \item{unemployment}{unemployment rate of the census area, in percent.}
#' }
#' @source Reports: Technical Unit for Vector Control, municipality of Madrid.
#'   Covariates: the database of Tamayo-Uría et al. (2013). Prepared by the
#'   authors of the package (script `data-raw/madrid.R` in the source repository).
#' @references
#' Tamayo-Uría, I., Cámara, R. and Escobar, F. (2013). Identificación de factores
#' de riesgo de infestación de plagas urbanas. *Revista de Salud Ambiental*,
#' 13(2), 148--157.
#'
#' Tamayo-Uría, I., Mateu, J. and Mughini-Gras, L. (2014). Risk factors and
#' spatial distribution of urban rat infestations. *Journal of Pest Science*,
#' 87, 107--115.
#' @seealso [madrid_A] for the contiguity matrix and [madrid_boundary] for the
#'   outline of the region. `vignette("madrid", package = "hotcoldESF")` reproduces
#'   the main map of the paper.
#' @examples
#' str(madrid)
#' # cockroach reports, summers 2010-2012, and the covariates of model M2
#' Y <- with(madrid, cockroaches_2010 + cockroaches_2011 + cockroaches_2012)
#' Z <- with(madrid, cbind(log_residential = log1p(residential), density,
#'                         log_area = log(area), dist_green, green_area, year_built,
#'                         dist_water, dist_fountain, dist_cat, unemployment))
#' \donttest{
#' set.seed(2026)
#' h <- hot_cold_map(Y, Z, madrid_A)
#' h
#' }
"madrid"

#' Contiguity matrix of the Madrid grid
#'
#' Binary queen contiguity matrix of the 493 cells of [madrid]: entry (i, j) is
#' 1 when cells i and j share a border or a corner. Every cell has at least one
#' neighbour.
#' @format A 493 x 493 numeric matrix with zero diagonal (1798 edges).
#' @seealso [madrid]
#' @examples
#' table(rowSums(madrid_A))
"madrid_A"

#' Outline of the Madrid study region
#'
#' Boundary of the study region of [madrid] (the union of the Madrid districts),
#' simplified to a tolerance of 50 m, for plotting.
#' @format A data frame with columns `x`, `y` (ETRS89 / UTM zone 30N, metres) and
#'   `ring` (polygon ring; draw each ring with [graphics::polygon()]).
#' @seealso [madrid]
#' @examples
#' plot(madrid$x, madrid$y, pch = 15, cex = 0.6, col = "grey80", asp = 1,
#'      xlab = "", ylab = "")
#' for (r in unique(madrid_boundary$ring))
#'   with(subset(madrid_boundary, ring == r), polygon(x, y, border = "grey30"))
"madrid_boundary"
