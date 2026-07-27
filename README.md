# Simulating Fine-Scale Population Flow for Infectious Disease Transmission Modelling

> **Read the full thesis online:** [https://iprincegh.github.io/HighRes-InfectiousDisease-Modelling/](https://iprincegh.github.io/HighRes-InfectiousDisease-Modelling/)


## Abstract

Infectious disease spreads where people take it, yet most surveillance systems map risk at administrative boundaries drawn without transmission in mind. Predicting where and when transmission concentrates requires models that fuse fine-scale mobility data with epidemiological surveillance while quantifying uncertainty. This thesis develops a Bayesian spatiotemporal model that maps COVID-19 transmission risk at 100 m × 100 m grid resolution and daily temporal resolution by combining anonymised mobile-phone-derived mobility data with case surveillance, census demographics, and land use records.

The model integrates an origin–destination infection pressure kernel, derived from anonymised visitor flows, with spatially and temporally structured random effects under a Negative Binomial likelihood fitted via Integrated Nested Laplace Approximations (INLA). Six model variants of increasing complexity are compared, from a spatial-only baseline to a full Knorr-Held Type IV spatiotemporal interaction model in which each area's epidemic trajectory is allowed to differ.

The approach is tested on two contrasting study areas, New York City (NYC) and Utah, across three pandemic waves spanning distinct epidemiological regimes: Wave 1 (March–May 2020, initial lockdown), Wave 2 (December 2020–February 2021, Alpha-dominant winter surge), and Wave 3 (December 2021–February 2022, Omicron). The Type IV interaction model returns the best or statistically equivalent fit across all six study configurations. Out-of-sample evaluation through both temporal forward and spatial blocked cross-validation shows strong generalisation, with NYC Wave 1 reaching area-level R² = 0.80 under spatial holdout. Thirteen covariates, including building density, infection pressure, point-of-interest activity density, and contact intensity, are consistently significant in New York City. Building density (relative risk 2.20–2.47) and origin–destination infection pressure (relative risk 1.47–1.54) are the dominant spatial risk factors. Mobility covariates are most informative during intermediate pandemic phases, when movement patterns vary markedly but case volumes have not yet saturated all spatial units. A counterfactual decomposition of the linear predictor quantifies the Proportion of Cases Attributable to Mobility (PCAtM): in NYC, the six mobility covariates collectively account for 14–40% of predicted cases depending on the wave, providing a single policy-translatable metric that bridges Bayesian model output and public health decision-making.

**Keywords:** Bayesian inference, spatiotemporal modelling, COVID-19, human mobility, INLA, disease mapping, change of support

## Repository Structure

```
helpers/                        # Shared R utility functions
  _utils.R                      # Core helpers (I/O, model wrappers, scoring, VIF)
  _adjacency.R                  # 10-stage augmented adjacency graph builder
  _cv_helpers.R                 # Cross-validation scaffolding and transport globals

pipeline/                       # Numbered analysis pipeline (run in order)
  00_setup.Rmd                  # Package installation and CONFIG object
  01_spatial_support.Rmd        # Build spatial support (100 m grid for NYC, counties for Utah)
  02a_population_offset.Rmd     # Dasymetric population disaggregation → E offsets
  02b_case_data.Rmd             # Process NYC MODZCTA + Utah NYT/CDC case counts
  03_advan_subset.Rmd           # Process Advan monthly patterns → POI, OD, chain
  03b_diurnal_dwell_extract.Rmd # Extracts and processes daily dwell-time
  04_mobility_diagnostics.Rmd   # QC of mobility indicators
  05_temporal_intensity.Rmd     # Baseline temporal intensity μ₀(g,t)
  06_od_kernel.Rmd              # OD flow kernel ρ(g,t) and ABM adjacency pairs
  07a_spatial_covariates.Rmd    # Land use, transport, building density covariates
  07b_demographics.Rmd          # Age structure and income/employment covariates
  08a_build_model_data.Rmd      # Assemble grid×time panel with all covariates
  08b_disaggregation_viz.Rmd    # Disaggregation conservation checks
  09_fit_model.Rmd              # Fit NegBin BYM2 + RW1/RW2 INLA models, select by WAIC
  09b_unified_component.Rmd     # Unified component extraction
  10a_temporal_cv.Rmd           # CPO/PIT + rolling-origin temporal CV + MODZCTA CV
  10b_spatial_cv.Rmd            # Area-level spatial CV (K-means NYC, LOCO Utah)
  11_compute_pcatm.Rmd          # PCAtM analysis (post-CV, pre-visualisation)
  12_visualization.Rmd          # Spatial/temporal random effects, fitted vs observed
  12b_visualization.Rmd         # Additional visualisations
  run_pipeline.Rmd              # Master runner (knits all stages)

generate_figures/                # Publication figure generators
  generate_publication_figures.R    # Main thesis figures → figures_clean/
  generate_cv_figures.R             # CV diagnostic figures (CPO, PIT, coverage)
  generate_cross_wave_figures.R     # Cross-wave facet panels
  generate_nonlinearity_diagnostics.R  # Nonlinearity scatter diagnostics
```

## Requirements

- **R ≥ 4.3** with packages (installed by `00_setup.Rmd`):
  `data.table`, `sf`, `spdep`, `INLA`, `ggplot2`, `patchwork`, `arrow`, `terra`, `viridis`, `scales`, `jsonlite`, `duckdb`, `lubridate`
- ~16 GB RAM (INLA fits on large grids)
- Data files in `data/` (not tracked — see structure below)

## How to Run

> **Important:** All scripts assume the **project root** as the R working directory.
> Helper scripts are sourced as `source("helpers/_utils.R")` with explicit folder prefixes.

### Option B — Run individual stages

```r
setwd("<project-root>")
rmarkdown::render("pipeline/02b_case_data.Rmd")
```

### Generate figures

```r
setwd("<project-root>")
source("generate_figures/generate_publication_figures.R")
source("generate_figures/generate_cv_figures.R")
source("generate_figures/generate_cross_wave_figures.R")
source("generate_figures/generate_nonlinearity_diagnostics.R")
source("generate_figures/generate_pcatm_oof_figures.R")
source("generate_figures/generate_cpo_pit_plots.R")
```

## Data Structure

- `data/advan_monthly_patterns/` — [Advan Research Monthly Patterns](https://docs.deweydata.io/docs/advan-research-monthly-patterns) (monthly mobility files)
- `data/building_footprints/` — [Microsoft US Building Footprints](https://github.com/microsoft/USBuildingFootprints)
- `data/cbg_shapefiles/` — [US Census TIGER/Line Block Groups](https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html)
- `data/census_tract/` — [US Census Tract Data](https://www.census.gov/geographies/reference-files/time-series/geo/tallies.html)
- `data/demographics/` — [ACS S0101](https://data.census.gov/table/ACSST5Y2020.S0101) and [ACS DP03](https://data.census.gov/table/ACSDP5Y2020.DP03)
- `data/landuse/` — [NYC MapPLUTO 25v3](https://www.nyc.gov/site/planning/data-maps/open-data/dwn-pluto-mappluto.page) and [Utah AGRC GIS Open Data](https://gis.utah.gov/)
- `data/MODZCTA/` — [NYC DOHMH Coronavirus Trends / MODZCTA resources](https://github.com/nychealth/coronavirus-data/tree/master/trends)
- `data/population/` — [WorldPop 100m Population Estimates](https://doi.org/10.5258/SOTON/WP00839)
- `data/transport/` — [NYC DOT Bridge Ratings](https://data.cityofnewyork.us/Transportation/Bridge-Ratings/4yue-vjfc), [NYU Subway Stations](https://geo.nyu.edu/catalog/nyu-2451-34758), and [TIGER/Line NY Roads](https://catalog.data.gov/dataset/tiger-line-shapefile-2023-state-new-york-primary-and-secondary-roads)
- `data/us_CovidCases/` — [NYC DOHMH Coronavirus Trends](https://github.com/nychealth/coronavirus-data/tree/master/trends) and [NYT COVID-19 Data](https://github.com/nytimes/covid-19-data)

## Outputs

All outputs are written to `outputs/spatiotemporal_risk/{state}_wave{1,2,3}/` with subdirectories `models/`, `validation/`, `figures/`, and `tables/`.
