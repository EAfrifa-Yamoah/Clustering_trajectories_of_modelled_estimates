#!/usr/bin/env bash
# Reproduce every table and figure, render the supplement tables, run the regression check, build the PDFs.
set -euo pipefail
cd "$(dirname "$0")"
Rscript analysis/01_run_analysis.R
Rscript analysis/02_supplement_tables.R
Rscript tests/test_reproduction.R
