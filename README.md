# hotcoldESF

<!-- badges: start -->
<!-- badges: end -->

**hotcoldESF** maps covariate-adjusted hot and cold spots of area-level counts.
A hot (cold) spot is an area with more (fewer) events than its covariates predict,
and with neighbours that share this excess. The procedure

1. tests the residuals of a regression on the covariates for spatial autocorrelation
   (Moran's I with its exact moments) and stops if there is none;
2. estimates the spatial component by eigenvector spatial filtering, with the smooth
   eigenvectors of the contiguity matrix and a ridge penalty chosen by cross-validation;
3. labels each area hot, cold or neutral with a permutation band;
4. groups hot and cold areas into zones and orders the areas into a priority list.

This repository also contains the simulation study, the application to cockroach
reports in Madrid, and the LaTeX source of the paper that introduces the method.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yandiracuvero/hotcoldESF", build_vignettes = TRUE)
```

## Example

```r
library(hotcoldESF)

# cockroach reports in Madrid, summers 2010-2012, 1 km cells
Y <- with(madrid, cockroaches_2010 + cockroaches_2011 + cockroaches_2012)
Z <- with(madrid, cbind(log_residential = log1p(residential), density,
                        log_area = log(area), dist_green, green_area, year_built,
                        dist_water, dist_fountain, dist_cat, unemployment))
set.seed(2026)
h <- hot_cold_map(Y, Z, madrid_A)
h                      # 68 hot and 61 cold cells
rank_cells(h)[1:25]    # priority list for a budget of 25 cells
```

Two vignettes walk through the method: `vignette("hotcoldESF")` on a simulated
dataset and `vignette("madrid")` on the Madrid data.

## Repository contents

| Folder | Contents |
|---|---|
| `R/`, `man/`, `tests/`, `vignettes/`, `data/` | the R package |
| `data-raw/` | how the package datasets are built from `application/data/` |
| `simulation/` | simulation study: scripts, results and figures |
| `application/` | Madrid application: scripts, results and figures |
| `paper/` | LaTeX source of the paper |
| `docs/` | notes, including `cran-submission.md` |

The large raw covariate file (`application/data/TODO_MADRID119m.txt`, 103 MB) is not
in the repository; `application/00_clean_covariates.R` describes how it is processed.

## Data

The Madrid reports were collected by the Technical Unit for Vector Control of the
municipality of Madrid; the covariates come from the database of Tamayo-Uría et al.
(2013). See `?madrid` for the sources and the processing.
<!-- TODO: data licence (to agree with the co-authors) -->

## Licence

Code: MIT (see `LICENSE`).
<!-- TODO: data licence -->

## Authors

Yandira Cuvero, Rowan Cherodian, Jorge Mateu, Carlos Almeida.
