# Clustering trajectories of modelled estimates

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23202766.svg)](https://doi.org/10.5281/zenodo.23202766)

Code and data for *Clustering trajectories of modelled estimates: an uncertainty aware framework, with an application to global cardiovascular disease burden, 1990 to 2023* (Afrifa-Yamoah, draft v1.0, October 2026).

The framework clusters country trajectories of modelled health estimates while carrying the published uncertainty through to cluster membership and to the choice of the number of clusters. Every table and figure in the manuscript and its supplement is produced by two R scripts from the data files in `data/`; the LaTeX tabulars written to `output/tables/` are the ones the manuscript inputs.

## Quick start

```bash
git clone https://github.com/EAfrifa-Yamoah/Clustering_trajectories_of_modelled_estimates.git
cd Clustering_trajectories_of_modelled_estimates
./run_all.sh          # about 2 minutes; builds output/
```

or step by step:

```bash
Rscript analysis/01_run_analysis.R        # every table and figure, output/final.RData
Rscript analysis/02_supplement_tables.R   # Tables S1 to S7
Rscript tests/test_reproduction.R         # checks the run against the archived partition
```

Set `GBD_B=200` in the environment for a quick run with fewer perturbation draws (the regression test will then fail, as intended).

## Requirements

R 4.3 or later with the packages `cluster`, `nnet` and `maps` (all on CRAN; on Ubuntu `apt install r-cran-cluster r-cran-nnet r-cran-maps`). No other packages. `sessionInfo.txt` records the versions used for the archived run.

## Layout

```
R/framework.R                 the method: equations (1) to (12) as documented functions
analysis/01_run_analysis.R    primary analysis, deaths, 2000 window, bootstrap, narrow interval subset,
                              attribution, tables 1 to 5, figures 1 to 8
analysis/02_supplement_tables.R   supplementary tables from output/final.RData
tests/test_reproduction.R     regression check against data/reference_membership.csv
tools/unwrap.R, unwrap_csv.R  convert GBD connector responses to panel CSVs (extraction only)
data/                         input panels and lookups (see below)
output/                       figures/, tables/ (LaTeX tabulars), final.RData, final_membership.csv (generated)
```

## Data

All inputs are public GBD 2023 outputs of the Institute for Health Metrics and Evaluation, retrieved through the IHME GBD 2023 Model Context Protocol connector in September and October 2026. They are redistributed here under IHME's free of charge non commercial user agreement; see https://www.healthdata.org/data-tools-practices/data-practices/ihme-free-charge-non-commercial-user-agreement.

| File | Content |
|---|---|
| `data/panel_cvd.csv` | cardiovascular diseases (GBD level 2), both sexes, age standardised and all ages DALY rate per 100,000 with 95% interval on the all ages rate; 203 countries, 1990 to 2023 |
| `data/panel_cvd_deaths.csv` | the same for death rates |
| `data/sdi_2023.csv` | Socio-demographic Index, 2023, 204 locations |
| `data/superregion.csv` | GBD super region per location |
| `data/income_group.csv` | World Bank income group used to partition the extraction |
| `data/analysis_set.csv` | the 201 countries in the analysis set with region and income group |
| `data/reference_membership.csv` | the archived partition and membership probabilities (seed 20261003, B = 1000) |

The extraction schedule and the connector's response caps are documented in Section S1 of the supplement. The `tools/` scripts are kept for anyone repeating the extraction; they are not needed to reproduce the analysis.

## Method in brief

For country *i* and year *t*, with *y* the log age standardised rate and *s* the log scale standard deviation implied by the published 95% interval:

1. **Level and shape.** `shape_of()` subtracts each country's mean log rate, so clusters describe trajectory form, not burden magnitude.
2. **Clustering.** Ward (ward.D2) on Euclidean distance between shape vectors; one tree, cut at each *k*.
3. **Membership probability.** `stability(mode = "value")` perturbs every value within its published SD, reclusters, aligns labels by majority overlap, and records the share of draws in which each country returns to its reference cluster. `mode = "country"` runs Hennig's country bootstrap instead for comparison.
4. **Number of clusters.** The largest *k* at which every cluster's median membership probability is at least 0.80. Silhouette and the gap statistic are reported alongside.
5. **Assignment.** Countries with membership probability below 0.70 are reported as unassigned, not forced into a cluster.
6. **Attribution.** Moran's I on nearest neighbour weights; logistic and multinomial models on SDI, super region and interval width.

`framework()` runs steps 1 to 5 and returns everything the scripts need.

## Citation

If you use the code or data, please cite the archived release:

Afrifa-Yamoah E. (2026). *Clustering trajectories of modelled estimates: uncertainty aware trajectory clustering code and data* (Version 1.0.0) [Computer software]. Zenodo. [https://doi.org/10.5281/zenodo.23202766](https://doi.org/10.5281/zenodo.23202766)

and the accompanying manuscript:

Afrifa-Yamoah E. Clustering trajectories of modelled estimates: an uncertainty aware framework, with an application to global cardiovascular disease burden, 1990 to 2023. Draft, 2026.

Citation metadata are in `CITATION.cff`; GitHub's "Cite this repository" button reads from it.

## Licence

Code: MIT (see `LICENSE`). Data: IHME free of charge non commercial user agreement (see above). The manuscript and supplement are not part of this repository.
