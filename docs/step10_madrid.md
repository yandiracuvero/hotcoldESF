# Step 10: Madrid application (20 Oct 2026)

Scripts: `application/01_prepare_madrid.R` (data), `application/02_madrid_analysis.R` (analysis, about 4 min).
Results: `application/results/` (per-summer summary, persistence, hot-again, out-of-sample), figures `application/figures/madrid_maps.png`, `madrid_persistence.png`, `madrid_out_of_sample.png`, table `paper/tables/madrid_out_of_sample.tex`.

## Data (reproduces the original preparation in 50Location.Rmd)

- Cockroach sightings, April–September: 1983, 1852, 1960, 1318 (2010–2013); 3 sightings in 2013 fall outside the grid.
- 493 cells of 1 km (cells under 7500 m² dropped), queen contiguity (1798 edges); 53% of cells have no sightings in any summer.
- 19 covariates averaged per cell (same variables and filling as the original code).
- Method: as in the simulation study (log counts, linear covariates, Moran gate, ESF smooth eigenvectors + ridge + CV, 999 permutations). Sensitivity: splines (df 4) for each covariate.

## Results

- The gate passes clearly every summer (Moran z 13.7–15.3; with splines 6.2–8.4).
- ESF labels 92–106 cells hot and 83–104 cold per summer; the hot cells hold 62–67% of the sightings. Gi*: 79–87 hot; LISA: 60–68 hot.
- Hot areas are stable: 79–92% of ESF hot cells are hot again the next summer (Gi* 78–87%, LISA 72–89%). 66 cells are hot in all four summers with ESF (Gi* 59, LISA 39).
- The hot areas form a large central-southern band. Cold areas are the periphery.
- With splines, fewer cells are hot (69–89) and they hold 45–52% of sightings; the main zones remain.

## Out-of-sample (cells chosen with summer t, evaluated in summer t + 1), 25 cells (5%)

| | Sightings captured | Excess captured |
|---|---|---|
| ESF | 23–25% | 25–33% |
| Gi* | 21–23% | 22–31% |
| LISA | 25–27% | 34–41% |
| Top residuals | 19–27% | 26–45% |
| Top counts | 28–30% | 32–42% |
| Random | 5% | 0% |

With 50 cells: ESF 38–39%, Gi* 37%, LISA 37–40%, top counts 48–51% of sightings.

## Points to discuss

1. **Top counts predicts next summer's sightings best.** Sightings are very persistent, so last year's counts are a strong predictor of raw sightings. The hot-spot methods target a different quantity: excess beyond what the covariates explain. The paper should say which goal trap placement serves: catching the most cockroaches (counts) or finding places where infestation exceeds what the environment predicts (hot spots, e.g. to investigate causes or to target interventions).
2. **Sightings are citizen reports.** The covariates include no population or housing density, so the central hot band may partly reflect where more people live and report. Is there a population variable (e.g. RPD?) or should one be added?
3. ESF beats Gi* on sightings captured in all three pairs, but LISA and top counts are higher.
