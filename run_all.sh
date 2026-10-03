#!/usr/bin/env bash
# Reproduce every table and figure, render the supplement tables, run the regression check, build the PDFs.
set -euo pipefail
cd "$(dirname "$0")"
Rscript analysis/01_run_analysis.R
Rscript analysis/02_supplement_tables.R
Rscript tests/test_reproduction.R
if command -v pdflatex >/dev/null; then
  mkdir -p manuscript/figures manuscript/tables
  cp output/figures/fig*.png manuscript/figures/; cp output/tables/tab*.tex manuscript/tables/
  (cd manuscript && pdflatex -interaction=nonstopmode manuscript.tex >/dev/null && bibtex manuscript >/dev/null && \
   pdflatex -interaction=nonstopmode manuscript.tex >/dev/null && pdflatex -interaction=nonstopmode manuscript.tex >/dev/null && \
   pdflatex -interaction=nonstopmode supplement.tex >/dev/null && pdflatex -interaction=nonstopmode supplement.tex >/dev/null)
  echo "Built manuscript/manuscript.pdf and manuscript/supplement.pdf"
fi
