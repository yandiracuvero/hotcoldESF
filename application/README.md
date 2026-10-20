# application

Madrid cockroach sightings (April-September, 2010-2013), 1 km cells.

- `01_prepare_madrid.R`: grid, queen adjacency, 19 covariates, counts per summer -> `data/madrid_grid.rds`.
  Needs `data/TODO_MADRID119m.txt` (107 MB, not in git; copy it from 02SepTest/Codigo).
- `02_madrid_analysis.R`: hot/cold maps per summer (ESF, Gi*, LISA), persistence, out-of-sample table
  -> `results/`, `figures/`, `paper/tables/madrid_out_of_sample.tex`.

Run both from the repo root.
