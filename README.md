# Clustering trajectories of modelled estimates

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23202765.svg)](https://doi.org/10.5281/zenodo.23202765)

Code and data for *Clustering trajectories of modelled estimates: an uncertainty aware framework, with an application to global cardiovascular disease burden, 1990 to 2023* (Afrifa-Yamoah, draft v1.2, 2026).

The framework clusters country trajectories of modelled health estimates while carrying the published uncertainty through to cluster membership and to the choice of the number of clusters. Every table and figure of the primary analysis is produced by two R scripts from the data files in `data/`. Two further scripts examine how far the result depends on choices the primary analysis fixes (independence of the perturbations across years, the clustering algorithm, the scale of the perturbation) and test the selection rule on simulated data with a known number of trajectory types. The LaTeX tabulars written to `output/tables/` are the ones the manuscript inputs.

The repository is at version 1.1.1. Version 1.1.0 added the sensitivity and simulation scripts to the primary analysis of v1.0.0; version 1.1.1 changes documentation and citation metadata only. See [Citation](#citation) for the archived releases.

## Quick start

```bash
git clone https://github.com/EAfrifa-Yamoah/Clustering_trajectories_of_modelled_estimates.git
cd Clustering_trajectories_of_modelled_estimates
./run_all.sh          # primary analysis, sensitivity analyses and regression check; builds output/
./run_all.sh --sim    # the same, plus the simulation with known truth
```

or step by step:

```bash
Rscript analysis/01_run_analysis.R                       # primary tables and figures, output/final.RData
Rscript analysis/02_supplement_tables.R                  # Tables S1 to S7
Rscript analysis/03_sensitivity_dependence_algorithm.R   # Table 6, Tables S8 to S11, Figure S3
Rscript analysis/04_simulation_known_truth.R             # Table S12 (optional, slow)
Rscript tests/test_reproduction.R                        # checks the run against the archived results
```

The scripts must run in this order: `03` reads `output/final.RData` written by `01`, and `04` reads both that file and `output/sensitivity.RData` written by `03`.

Runtimes measured on a four core Linux machine (Ubuntu 24.04, R 4.3.3) at the default settings:

| Step | Runtime |
|---|---|
| `01_run_analysis.R` | about 1 minute |
| `02_supplement_tables.R` | under 1 second |
| `03_sensitivity_dependence_algorithm.R` | about 3 minutes |
| `04_simulation_known_truth.R` | about 11 minutes |
| `tests/test_reproduction.R` | under 1 second |

`./run_all.sh` therefore takes about 4 minutes and `./run_all.sh --sim` about 15 minutes on such a machine. `03` and `04` run their jobs in parallel, so both take longer on fewer cores; `01` does not use more than one.

### Run settings

| Environment variable | Default | Used by | Meaning |
|---|---|---|---|
| `GBD_B` | 1000 | `01`, `03` | perturbation draws per stability run. `GBD_B=200` gives a quick run; the regression test then fails, as intended |
| `SIM_CORES` | all detected cores | `03`, `04` | worker processes for `parallel::mclapply` |
| `SIM_R` | 50 | `04` | simulated data sets per design cell |
| `SIM_B` | 100 | `04` | perturbation draws per simulated data set |
| `SIM_REUSE` | unset | `04` | `SIM_REUSE=1` rebuilds the summary and Table S12 from the archived `output/tables/sim_replicates.csv` without rerunning the simulation |

Every run in `03` and every data set and analysis in `04` sets its own seed (all derived from 20261003), so the results do not depend on the number of cores and any single row can be reproduced alone.

## Requirements

R 4.3 or later with the packages `cluster`, `nnet` and `maps` (all on CRAN; on Ubuntu `apt install r-cran-cluster r-cran-nnet r-cran-maps`), plus `parallel`, which ships with R. No other packages. `sessionInfo.txt` records the versions used for the archived run.

`03` and `04` use forked workers through `parallel::mclapply`, which R supports on Linux and macOS only. On Windows set `SIM_CORES=1`; the scripts then run serially and give the same results, more slowly.

## Layout

```
R/framework.R                 the method: equations (1) to (13) as documented functions (numbering of
                              manuscript draft v1.2)
analysis/01_run_analysis.R    primary analysis, deaths, 2000 window, bootstrap, narrow interval subset,
                              attribution, tables 1 to 5, figures 1 to 8
analysis/02_supplement_tables.R   supplementary tables from output/final.RData
analysis/03_sensitivity_dependence_algorithm.R
                              sensitivity of the primary result to dependence between years, clustering
                              algorithm and perturbation scale (Table 6, Tables S8 to S11, Figure S3)
analysis/04_simulation_known_truth.R
                              simulation with one, two or three known components (Table S12)
tests/test_reproduction.R     regression check against data/reference_membership.csv and, when
                              output/sensitivity.RData exists, against the archived sensitivity results
tools/unwrap.R, unwrap_csv.R  convert GBD connector responses to panel CSVs (extraction only)
data/                         input panels and lookups (see below)
output/                       figures/, tables/ (LaTeX tabulars and CSVs); final.RData, final_membership.csv,
                              sensitivity.RData and simulation.RData are generated and not tracked
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

## Sensitivity analyses

The primary analysis makes two choices that the published estimates cannot settle. It perturbs each country's values independently across years, although the errors of a modelled series are likely to be correlated over time, and it clusters with Ward's criterion. `framework()` and `stability()` therefore take two further arguments.

- `rho` is the lag one correlation of the perturbations across years. `ar1_noise()` draws a stationary AR(1) series with unit marginal variance for each country, so every value keeps the standard deviation implied by its published interval whatever `rho` is. `rho = 0` (the default) gives the independent draws of the primary analysis; `rho = 1` gives one common shock per country.
- `method` (in `framework()`) or `cl` (in `stability()`) selects the clustering algorithm through `clusterer()`: `"ward"` (default), `"pam"` (partitioning around medoids), `"kmeans"` (Hartigan and Wong, 25 random starts) or `"average"` (average linkage).

With the defaults the functions consume the random number stream as before and reproduce the archived primary run exactly; the regression test checks this. `framework()` now returns `k = NA` when no partition meets the stability rule, and then reports cluster labels and membership probabilities at *k* = 2.

`analysis/03_sensitivity_dependence_algorithm.R` applies these to the primary panel (DALYs, 1990 to 2023) in four parts.

| Part | Question | Output in `output/tables/` |
|---|---|---|
| (a) | Does AR(1) dependence of the perturbations (`rho` = 0.5 to 1, Ward) change the selected *k* or the assignments? A second `rho = 0` run with a different seed gauges Monte Carlo error. The shape perturbation variance is also decomposed by `rho` | `sens_dependence.csv`, `sens_dependence_variance.csv`, `tabS8.tex` |
| (b) | Do PAM, k-means and average linkage select the same *k* as Ward, alone and combined with `rho` = 0.7, 0.9 and 0.95? | `sens_algorithm.csv`, `tab6.tex`, `tabS9.tex` |
| (c) | Do Ward, PAM and k-means give a country the same label at *k* = 2, and which countries differ? How do their three cluster partitions compare? | `sens_agreement_k2.csv`, `sens_three_clusters.csv`, `tabS10.tex` |
| (d) | How does the weakest cluster's median membership at *k* = 2, 3 and 4 respond when every published SD is multiplied by 0.1, 0.25, 0.5 or 2? | `sens_scale.csv`, `tabS11.tex` |

Figure S3 (`output/figures/figS3.png`) plots the weakest cluster median membership against *k* by `rho` and by algorithm, and the *k* = 2 membership probabilities at `rho = 0.9` against those at `rho = 0`.

### What the archived run shows

These statements describe the tracked CSV files (seed 20261003, B = 1000) and are specific to the cardiovascular case study.

- **Dependence, Ward.** The stability rule selects two clusters at every `rho` examined. For `rho` up to 0.99 between 167 and 172 of the 201 countries are assigned, and 1 to 4 countries change assignment status relative to the primary run. An independent replicate at `rho = 0` with a different seed changes 3, so differences of this size are within Monte Carlo error. At `rho = 1` the picture differs: 187 countries are assigned and 17 change status, since a common shock per country leaves little perturbation of shape (about 3% of the variance at `rho = 0`).
- **Algorithm.** Silhouette selects 2 and the gap statistic 3 under Ward, PAM and k-means alike. The stability rule does not agree across algorithms. With independent perturbations it selects 2 under Ward and 4 under both PAM and k-means; with AR(1) perturbations it selects 2 under Ward, 2 or 3 under PAM and 3 under k-means. Under average linkage no *k* meets the rule, and its two cluster partition isolates 2 countries.
- **Labels at *k* = 2.** Ward, PAM and k-means give the same label to 170 of 201 countries. Of the 31 that differ, 29 are unassigned in the primary analysis. Among countries assigned both in the primary analysis and under the alternative algorithm, labels coincide for 170 of 170 (PAM) and 166 of 168 (k-means) with independent perturbations. The reference partitions themselves differ more than this suggests (adjusted Rand index with Ward 0.72 for PAM and 0.48 for k-means); the disagreement falls almost entirely on countries the primary analysis leaves unassigned.
- **Perturbation scale.** The selected *k* depends on the size of the perturbation. With the published SDs multiplied by 0.25, Ward's weakest cluster median membership is 0.85 at *k* = 3 and 0.81 at *k* = 4, both above the 0.80 threshold, against 0.70 and 0.53 at the published scale.

Taken together, the two cluster assignment of the countries the primary analysis assigns is robust to these choices, whereas the number of clusters is not.

## Simulation with known truth

`analysis/04_simulation_known_truth.R` asks what the stability rule selects, and how accurate the assigned labels are, when the number of trajectory types is known.

- **Truth.** *K* = 1, 2 or 3 Gaussian components in shape space, with centroids, sizes and pooled within component covariance taken from the observed shapes (*K* = 2 from the Ward two cluster partition, *K* = 3 from the k-means three cluster partition).
- **Estimates.** Truth plus error with the observed log scale SDs, resampled within the matching observed cluster, and AR(1) correlation 0 or 0.9 across years.
- **Analysis.** The framework under Ward, PAM and k-means, with perturbations that either match the error correlation or wrongly assume independence. Three arms (error correlation, assumed correlation): (0, 0), (0.9, 0), (0.9, 0.9).
- **Outcomes.** The *k* selected by the stability rule, and at the true *K* the share of units assigned and the share of assigned and of unassigned units whose label matches the true component.

It writes `output/tables/sim_replicates.csv` (one row per data set, arm and algorithm), `sim_summary.csv` and `tabS12.tex`.

### What the archived run shows

In this simulation setting (50 data sets per cell, 100 perturbation draws each):

- **One true component.** The rule rarely reports that no partition is stable. It does so in 18% to 34% of data sets under Ward and in none under PAM or k-means; in all other data sets it selects two or more clusters.
- **Two true components.** The rule selects 2 in 62% to 80% of data sets under Ward, 2% to 30% under PAM and 0% to 8% under k-means.
- **Three true components.** The rule selects 3 in 58% to 86% of data sets under Ward, 56% to 96% under PAM and 26% to 98% under k-means. PAM and k-means recover 3 most often when the assumed correlation matches the error correlation of 0.9 (96% and 98%).
- **Labels.** At the true *K*, 87% to 95% of assigned units carry the correct label, against 44% to 67% of unassigned units.

Several features limit what can be read from these figures. With 50 data sets per cell a proportion has a Monte Carlo standard error of up to about 7 percentage points. The components are Gaussian and calibrated to the observed partitions, with the two component truth taken from Ward and the three component truth from k-means, which may favour the algorithm that defined each truth. The simulation uses 100 draws per data set where the case study uses 1000. The results describe the behaviour of the rule in this design and are not evidence about its behaviour on other panels.

## Regression test

`tests/test_reproduction.R` stops with an error unless the run reproduces the archived partition and membership probabilities in `data/reference_membership.csv` to within 1e-9, together with the headline results of the primary analysis. When `output/sensitivity.RData` is present and B = 1000 it also checks that Ward selects two clusters at every `rho`, that PAM and k-means select four with independent perturbations, and that the label agreement counts (170 the same, 29 of the differing countries unassigned) match the archived run. The simulation is not covered by the test. A rerun at the default settings on Ubuntu 24.04 with R 4.3.3 and `cluster` 2.1.6 reproduced the tracked `sim_replicates.csv`, `sim_summary.csv` and `tabS12.tex` byte for byte; on another platform compare a rerun against the tracked `sim_summary.csv`.

## Citation

If you use the code or data, please cite the archived release:

Afrifa-Yamoah E. (2026). *Clustering trajectories of modelled estimates: uncertainty aware trajectory clustering code and data* (Version 1.1.1) [Computer software]. Zenodo. [https://doi.org/10.5281/zenodo.23202765](https://doi.org/10.5281/zenodo.23202765)

and the accompanying manuscript:

Afrifa-Yamoah E. Clustering trajectories of modelled estimates: an uncertainty aware framework, with an application to global cardiovascular disease burden, 1990 to 2023. Draft v1.2, 2026.

The DOI in this citation, [10.5281/zenodo.23202765](https://doi.org/10.5281/zenodo.23202765), is the Zenodo concept DOI: it covers all versions and resolves to the latest archived release. A DOI specific to a version exists only once Zenodo has archived that release, so it cannot be written into the files of the release itself.

| Version | DOI | Content |
|---|---|---|
| 1.1.1 | assigned by Zenodo when the release is archived; listed on the concept DOI page | as 1.1.0, with measured runtimes in the README and corrected version and DOI metadata in `CITATION.cff` |
| 1.1.0 | [10.5281/zenodo.24372735](https://doi.org/10.5281/zenodo.24372735) | primary analysis, sensitivity analyses (`analysis/03_*.R`), simulation (`analysis/04_*.R`), the `rho` and `method` arguments in `R/framework.R` |
| 1.0.0 | [10.5281/zenodo.23202766](https://doi.org/10.5281/zenodo.23202766) | primary analysis only |

Where a DOI specific to the version used is available, cite that one. Code, data and output files are identical in versions 1.1.0 and 1.1.1, so results do not depend on which of the two is used.

The v1.1.0 archive is a snapshot of commit `984db47`. Its code, data and output files are identical to those in this repository; its copies of `README.md` and `CITATION.cff` predate the runtime table and the version 1.1.0 metadata, so `CITATION.cff` inside the archive still reads version 1.0.0.

Citation metadata are in `CITATION.cff`, which records version 1.1.1 and the concept DOI; GitHub's "Cite this repository" button reads from it.

## Licence

Code: MIT (see `LICENSE`). Data: IHME free of charge non commercial user agreement (see above). The manuscript and supplement are not part of this repository.
