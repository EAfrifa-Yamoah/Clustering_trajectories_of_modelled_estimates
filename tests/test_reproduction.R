# test_reproduction.R
# Regression check: the primary partition reproduced from output/final.RData must match the
# membership file archived with the manuscript (data/reference_membership.csv) exactly.
root <- if (file.exists("R/framework.R")) "." else ".."; setwd(root); source("R/framework.R")
stopifnot(file.exists("output/final.RData")); load("output/final.RData")
ref <- read.csv("data/reference_membership.csv")
stopifnot(identical(as.integer(ref$assigned), as.integer(meta$assigned)))
stopifnot(all(abs(ref$p - meta$p) < 1e-9))
stopifnot(rD$k == 2, rD$k_gap == 3, rD$k_sil == 2, sum(meta$assigned == 0) == 31)
stopifnot(abs(ari(rD$assigned[rD$assigned > 0 & rM$assigned > 0], rM$assigned[rD$assigned > 0 & rM$assigned > 0]) - 1) < 1e-9)
cat("Reproduction check passed: partition, membership probabilities and headline results match the archived run.\n")
