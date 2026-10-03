#!/usr/bin/env Rscript
# Hyphenate conventional compound modifiers in manuscript.tex (and the supplementary file if present).
# Only compounds that conventionally take a hyphen when used attributively are changed; nominal uses stay open.
files <- commandArgs(trailingOnly = TRUE); if (!length(files)) files <- "manuscript.tex"
rules <- c(
  "age standardised"              = "age-standardised",
  "high income"                   = "high-income",
  "non fatal"                     = "non-fatal",
  "non zero"                      = "non-zero",
  "mid 1990s"                     = "mid-1990s",
  "log normal"                    = "log-normal",
  "log scale standard deviation"  = "log-scale standard deviation",
  "log scale SD"                  = "log-scale SD",
  "within cluster"                = "within-cluster",
  "row standardised"              = "row-standardised",
  "well measured"                 = "well-measured",
  "poorly measured"               = "poorly-measured",
  "best measured"                 = "best-measured",
  "well regarded"                 = "well-regarded",
  "late stalling"                 = "late-stalling",
  "two cluster solution"          = "two-cluster solution",
  "three cluster solution"        = "three-cluster solution",
  "three shape structure"         = "three-shape structure",
  "draw level"                    = "draw-level",
  "all ages rate"                 = "all-ages rate",
  "all ages interval"             = "all-ages interval",
  "zero to five star scale"       = "zero-to-five-star scale",
  "one standard error bars"       = "one-standard-error bars",
  "low order polynomial"          = "low-order polynomial",
  "level based classification"    = "level-based classification",
  "per observation standard deviation" = "per-observation standard deviation",
  "great circle distance"         = "great-circle distance",
  "covariate driven"              = "covariate-driven",
  "wide interval series"          = "wide-interval series",
  "six nearest neighbour weights" = "six-nearest-neighbour weights",
  "nearest neighbour weight matrix" = "nearest-neighbour weight matrix",
  "cluster wise"                  = "cluster-wise",
  "mean centred"                  = "mean-centred",
  "cause of death data"           = "cause-of-death data",
  "point estimates alone"         = "point estimates alone",
  "\\$k\\$ cluster partition"     = "$k$-cluster partition",
  "the \\$k\\$ cluster"           = "the $k$-cluster",
  "a lower order"                 = "a lower-order")
for (f in files) {
  s <- readLines(f, warn = FALSE); n0 <- 0
  for (k in names(rules)) { hits <- sum(grepl(k, s, fixed = !grepl("\\\\", k))); n0 <- n0 + hits
    s <- if (grepl("\\\\", k)) gsub(k, rules[[k]], s) else gsub(k, rules[[k]], s, fixed = TRUE) }
  writeLines(s, f); cat(sprintf("%s: %d replacement lines touched\n", f, n0)) }
