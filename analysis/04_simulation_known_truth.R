# 04_simulation_known_truth.R
# Small simulation with known truth: what does the stability rule select, and how accurate are the
# assigned labels, when the number of trajectory types is known?
#
# Design (stylised, calibrated to the cardiovascular panel):
#   truth     K = 1, 2 or 3 Gaussian components in shape space. Centroids, sizes and the pooled
#             within component covariance are taken from the observed shapes (K = 1: one component;
#             K = 2: the Ward two cluster partition; K = 3: the k-means three cluster partition).
#   estimate  y = truth + s * u, where s is a row of the observed log scale SD matrix drawn from the
#             countries of the matching observed cluster, and u is AR(1) with lag one correlation
#             rho_true (0 or 0.9) and unit marginal variance.
#   analysis  the framework applied to y and s with perturbations of correlation rho_an, under Ward,
#             PAM and k-means. Arms (rho_true, rho_an): (0, 0), (0.9, 0), (0.9, 0.9).
#   outcome   the k selected by the stability rule (pi = 0.80), and at k = K the share of units
#             assigned (tau = 0.70) and the share of assigned and of unassigned units whose reference
#             label matches the true component.
# Each data set and each analysis has its own seed, so results do not depend on the number of cores.
# Run from the repository root after 01 and 03:  Rscript analysis/04_simulation_known_truth.R
# Runtime about 30 minutes on two cores at the defaults (SIM_R = 50 replicates, SIM_B = 100 draws).
root <- if (file.exists("R/framework.R")) "." else if (file.exists("../R/framework.R")) ".." else stop("run from the repository root")
setwd(root); source("R/framework.R"); suppressPackageStartupMessages(library(parallel))
stopifnot(file.exists("output/final.RData"), file.exists("output/sensitivity.RData"))
load("output/final.RData"); load("output/sensitivity.RData")
R <- as.integer(Sys.getenv("SIM_R", "50")); B <- as.integer(Sys.getenv("SIM_B", "100"))
nc <- as.integer(Sys.getenv("SIM_CORES", max(1, detectCores()))); seed0 <- 20261003
n <- nrow(Y); nT <- ncol(Y); sh <- rD$shape; tau <- 0.70; pi_ <- 0.80; ks <- 2:5

# ---- truth: components estimated from the observed shapes -----------------------------------------
truth_def <- function(g) { g <- as.integer(factor(g)); K <- max(g)
  cen <- t(sapply(1:K, function(c) colMeans(sh[g == c, , drop = FALSE])))
  res <- sh - cen[g, , drop = FALSE]; e <- eigen(crossprod(res) / (n - K), symmetric = TRUE)
  list(K = K, g = g, cen = cen, root = e$vectors %*% diag(sqrt(pmax(e$values, 0)))) }
truths <- list(`1` = truth_def(rep(1, n)), `2` = truth_def(rD$stab[["2"]]$ref), `3` = truth_def(ra$kmeans$stab[["3"]]$ref))

gen <- function(tr, rho_true) {
  X <- tr$cen[tr$g, , drop = FALSE] + matrix(rnorm(n * nT), n) %*% t(tr$root)       # true shapes
  rows <- integer(n); for (c in 1:tr$K) { i <- which(tr$g == c); rows[i] <- i[sample.int(length(i), length(i), replace = TRUE)] }
  Ss <- S[rows, ]; list(Y = X + Ss * ar1_noise(n, nT, rho_true), S = Ss, g = tr$g) }

analyse <- function(d, method, rho_an) {
  cl <- clusterer(method); shp <- shape_of(d$Y)
  st <- lapply(ks, function(k) stability(d$Y, d$S, NULL, k, B, "value", rho = rho_an, cl = cl, ref = cl(shp, k))); names(st) <- ks
  weak <- sapply(st, function(s) min(tapply(s$p, s$ref, median))); ok <- ks[weak >= pi_]
  out <- c(setNames(weak, paste0("weak", ks)), k = if (length(ok)) max(ok) else NA, assigned = NA, acc_assigned = NA, acc_unassigned = NA, acc_all = NA)
  K <- max(d$g)
  if (K >= 2) { s <- st[[as.character(K)]]; map <- apply(table(s$ref, d$g), 1, which.max); hit <- map[s$ref] == d$g; a <- s$p >= tau
    out[c("assigned", "acc_assigned", "acc_unassigned", "acc_all")] <- c(mean(a), mean(hit[a]), if (any(!a)) mean(hit[!a]) else NA, mean(hit)) }
  out }

arms <- data.frame(rho_true = c(0, 0.9, 0.9), rho_an = c(0, 0, 0.9)); meths <- c("ward", "pam", "kmeans")
jobs <- expand.grid(rep = seq_len(R), rho_true = c(0, 0.9), K = 1:3)
if (Sys.getenv("SIM_REUSE") == "1" && file.exists("output/tables/sim_replicates.csv")) sim <- read.csv("output/tables/sim_replicates.csv") else {
res <- mclapply(seq_len(nrow(jobs)), function(j) {
  set.seed(seed0 + 1000 + j); d <- gen(truths[[as.character(jobs$K[j])]], jobs$rho_true[j])
  a <- arms[arms$rho_true == jobs$rho_true[j], ]
  do.call(rbind, lapply(seq_len(nrow(a)), function(q) do.call(rbind, lapply(seq_along(meths), function(m) {
    set.seed(seed0 + 500000 + 100 * j + 10 * q + m)
    data.frame(K = jobs$K[j], rho_true = a$rho_true[q], rho_an = a$rho_an[q], method = meths[m], rep = jobs$rep[j], t(analyse(d, meths[m], a$rho_an[q]))) })))) },
  mc.cores = nc, mc.preschedule = FALSE)
sim <- do.call(rbind, res); write.csv(sim, "output/tables/sim_replicates.csv", row.names = FALSE) }

# ---- summary --------------------------------------------------------------------------------------
sim$method <- factor(sim$method, levels = meths)
cells <- split(sim, list(sim$method, sim$rho_an, sim$rho_true, sim$K), drop = TRUE)
sumry <- do.call(rbind, lapply(cells, function(d) data.frame(
  K = d$K[1], rho_true = d$rho_true[1], rho_an = d$rho_an[1], method = as.character(d$method[1]), reps = nrow(d),
  none = mean(is.na(d$k)), k2 = mean(!is.na(d$k) & d$k == 2), k3 = mean(!is.na(d$k) & d$k == 3), k4 = mean(!is.na(d$k) & d$k == 4), k5 = mean(!is.na(d$k) & d$k == 5),
  weak_at_K = if (d$K[1] >= 2) mean(d[[paste0("weak", d$K[1])]]) else NA,
  assigned = mean(d$assigned), acc_assigned = mean(d$acc_assigned), acc_unassigned = mean(d$acc_unassigned, na.rm = TRUE), acc_all = mean(d$acc_all))))
rownames(sumry) <- NULL; num <- sapply(sumry, is.numeric); sumry[num] <- round(sumry[num], 3)
cat(sprintf("Simulation with known truth: %d replicates per cell, B = %d draws\n", max(sumry$reps), B)); print(sumry, row.names = FALSE)
write.csv(sumry, "output/tables/sim_summary.csv", row.names = FALSE)

# Table S12 for the supplement
pct <- function(x) ifelse(is.na(x) | is.nan(x), "", sprintf("%.0f", 100 * x)); mname <- c(ward = "Ward", pam = "PAM", kmeans = "$k$-means")
first <- !duplicated(sumry[, c("K", "rho_true", "rho_an")])
s12 <- data.frame(`True $K$` = ifelse(!duplicated(sumry$K), sumry$K, ""), `$\\rho$ of errors` = ifelse(first, sprintf("%.1f", sumry$rho_true), ""),
  `$\\rho$ assumed` = ifelse(first, sprintf("%.1f", sumry$rho_an), ""), Algorithm = mname[sumry$method],
  None = pct(sumry$none), `$k = 2$` = pct(sumry$k2), `$k = 3$` = pct(sumry$k3), `$k = 4$` = pct(sumry$k4), `$k = 5$` = pct(sumry$k5),
  `Assigned` = pct(sumry$assigned), `Correct, assigned` = pct(sumry$acc_assigned), `Correct, unassigned` = pct(sumry$acc_unassigned), check.names = FALSE)
tex_table(s12, "output/tables/tabS12.tex", "llllrrrrrrrr")
save(sim, sumry, truths, R, B, file = "output/simulation.RData"); cat("\nDONE\n")
