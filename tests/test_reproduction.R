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
if (file.exists("output/sensitivity.RData") && B == 1000) { load("output/sensitivity.RData")
  stopifnot(all(dep$k_stab == 2))                                            # Ward selects two clusters at every rho
  stopifnot(identical(alg$k_stab[alg$method == "pam" & alg$rho == 0], 4L), identical(alg$k_stab[alg$method == "kmeans" & alg$rho == 0], 4L))
  stopifnot(sum(same3) == 170, sum(!same3 & rD$assigned == 0) == 29)         # label agreement across algorithms at k = 2
  cat("Sensitivity check passed: selected k by dependence and algorithm, and label agreement, match the archived run.\n") }
cat("Reproduction check passed: partition, membership probabilities and headline results match the archived run.\n")
