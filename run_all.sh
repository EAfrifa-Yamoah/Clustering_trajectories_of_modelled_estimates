#!/usr/bin/env bash
# Reproduce every table and figure, render the supplement tables, run the sensitivity analyses and the
# regression check. The simulation with known truth takes about twenty minutes on two cores and runs
# only when asked for:  ./run_all.sh --sim
set -euo pipefail
cd "$(dirname "$0")"
Rscript analysis/01_run_analysis.R
Rscript analysis/02_supplement_tables.R
Rscript analysis/03_sensitivity_dependence_algorithm.R
if [[ "${1:-}" == "--sim" ]]; then Rscript analysis/04_simulation_known_truth.R; fi
Rscript tests/test_reproduction.R
