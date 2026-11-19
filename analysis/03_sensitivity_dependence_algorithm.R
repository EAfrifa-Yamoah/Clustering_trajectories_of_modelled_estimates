# 03_sensitivity_dependence_algorithm.R
# Sensitivity of the primary result (DALYs, 1990 to 2023) to two choices, run after 01_run_analysis.R:
#   (a) temporal dependence of the perturbations: AR(1) on the log scale with lag one correlation rho,
#       every value keeping the marginal standard deviation implied by its published interval;
#   (b) clustering algorithm: the same perturbation and stability rule with partitioning around medoids
#       (PAM), k-means and average linkage in place of Ward's criterion, alone and combined with (a);
#   (c) agreement of the k = 2 labels across algorithms, set against the primary assignment;
#   (d) scale of the perturbation: every published SD multiplied by a constant.
# Every run sets its own seed, so any row can be reproduced alone and the output does not depend on
# the number of cores. Writes output/tables/sens_*.csv, tab6.tex, tabS8 to tabS11, output/figures/figS3.png
# and output/sensitivity.RData.
# Run from the repository root:  Rscript analysis/03_sensitivity_dependence_algorithm.R
# Runtime about 8 minutes on two cores at B = 1000. Requires R >= 4.3 with package cluster.
root <- if (file.exists("R/framework.R")) "." else if (file.exists("../R/framework.R")) ".." else stop("run from the repository root")
setwd(root); source("R/framework.R"); suppressPackageStartupMessages(library(parallel))
stopifnot(file.exists("output/final.RData")); load("output/final.RData")   # rD (primary run), meta, Y, S
B <- as.integer(Sys.getenv("GBD_B", "1000")); seed0 <- 20261003; nc <- as.integer(Sys.getenv("SIM_CORES", max(1, detectCores())))
tab <- function(f) file.path("output/tables", f); fig <- function(f) file.path("output/figures", f)
n <- nrow(Y); nT <- ncol(Y); tau <- 0.70; pi_ <- 0.80
pc <- 100 * (exp(Y[, "2023"] - Y[, "1990"]) - 1)
run <- function(seed, ...) { set.seed(seed); framework(Y, S, B = B, ...) }
# labels ordered by median change so that 1 is the least decline, whatever the algorithm calls it
ordlab <- function(g) match(g, order(tapply(pc, g, median), decreasing = TRUE))

# ---- (a) temporal dependence, Ward ------------------------------------------------------------
# rho = 0 is the primary run. A second rho = 0 run with a different seed shows the Monte Carlo error.
rhos <- c(0.5, 0.7, 0.9, 0.95, 0.99, 1)
rr <- c(list(rD), mclapply(0:length(rhos), function(j) if (j == 0) run(seed0 + 100) else run(seed0 + j, rho = rhos[j]), mc.cores = nc, mc.preschedule = FALSE))
names(rr) <- c("0", "0 (replicate)", format(rhos))
stopifnot(all(sapply(rr, function(r) identical(r$stab[["2"]]$ref, rD$stab[["2"]]$ref))))   # one Ward tree, shared references
ref2 <- rD$cluster                                         # 1 stalled, 2 sustained decline
p2 <- sapply(rr, function(r) r$stab[["2"]]$p)              # n x runs membership probability at k = 2
as2 <- apply(p2, 2, function(p) ifelse(p >= tau, ref2, 0))
dep <- data.frame(rho = names(rr), k_stab = sapply(rr, `[[`, "k"), t(sapply(rr, function(r) round(r$weakest, 2))),
  assigned = colSums(as2 > 0), stalled = colSums(as2 == 1), decline = colSums(as2 == 2), unassigned = colSums(as2 == 0),
  median_p = round(apply(p2, 2, median), 3), mean_p = round(colMeans(p2), 3),
  spearman_with_rho0 = round(apply(p2, 2, function(p) cor(p, p2[, 1], method = "spearman")), 3),
  status_changed = colSums((as2 > 0) != (as2[, 1] > 0)), check.names = FALSE, row.names = NULL)
names(dep)[3:7] <- paste0("weakest_k", 2:6)
cat("\n(a) Temporal dependence of the perturbations, Ward, B =", B, "\n"); print(dep, row.names = FALSE)
write.csv(dep, tab("sens_dependence.csv"), row.names = FALSE)

# Where the perturbation variance goes. For country i the perturbation of the log series has covariance
# D_i R(rho) D_i with D_i = diag(s_i) and R the AR(1) correlation matrix. Centring removes the level
# component. Reported per rho, averaged over countries:
#   total = trace of the centred covariance (all of the shape perturbation), and
#   along = variance of the projection on the unit vector joining the two cluster centroids (the
#           direction in which a perturbation moves a country towards the other cluster's mean).
Cn <- diag(nT) - 1 / nT; sh <- rD$shape
dvec <- colMeans(sh[ref2 == 2, ]) - colMeans(sh[ref2 == 1, ]); dvec <- dvec / sqrt(sum(dvec^2))
vdec <- as.data.frame(t(sapply(c(0, rhos), function(rho) { Rm <- rho^abs(outer(1:nT, 1:nT, "-")); if (rho == 1) Rm[] <- 1
  v <- t(sapply(1:n, function(i) { Sg <- (S[i, ] %o% S[i, ]) * Rm; c(sum(diag(Cn %*% Sg %*% Cn)), drop(dvec %*% Sg %*% dvec)) }))
  c(rho = rho, total = mean(v[, 1]), along = mean(v[, 2])) })))
vdec$total_rel <- vdec$total / vdec$total[1]; vdec$along_rel <- vdec$along / vdec$along[1]
cat("\nShape perturbation variance by rho (mean over countries; _rel is relative to rho = 0)\n"); print(round(vdec, 4), row.names = FALSE)
write.csv(vdec, tab("sens_dependence_variance.csv"), row.names = FALSE)

# ---- (b) clustering algorithm, alone and with AR(1) perturbations ----------------------------------
meth <- c("pam", "kmeans", "average"); rho_b <- c(0.7, 0.9, 0.95)
jb <- rbind(data.frame(method = meth, rho = 0, seed = seed0 + 200 + seq_along(meth)),
            expand.grid(method = c("pam", "kmeans"), rho = rho_b, seed = NA, stringsAsFactors = FALSE))
jb$seed[is.na(jb$seed)] <- seed0 + 300 + seq_len(sum(is.na(jb$seed)))
rb <- mclapply(seq_len(nrow(jb)), function(j) run(jb$seed[j], method = jb$method[j], rho = jb$rho[j]), mc.cores = nc, mc.preschedule = FALSE)
ra <- c(list(ward = rD), setNames(rb[jb$rho == 0], jb$method[jb$rho == 0]))              # independent perturbations
allruns <- c(list(list(method = "ward", rho = 0, r = rD)),
             lapply(which(names(rr) %in% format(rho_b)), function(i) list(method = "ward", rho = as.numeric(names(rr)[i]), r = rr[[i]])),
             lapply(seq_len(nrow(jb)), function(j) list(method = jb$method[j], rho = jb$rho[j], r = rb[[j]])))
alg <- do.call(rbind, lapply(allruns, function(a) { r <- a$r; g2 <- ordlab(r$stab[["2"]]$ref); a2 <- ifelse(r$stab[["2"]]$p >= tau, g2, 0)
  both <- a2 > 0 & rD$assigned > 0
  data.frame(method = a$method, rho = a$rho, k_sil = r$k_sil, k_gap = r$k_gap, k_stab = r$k, t(round(r$weakest, 2)),
    smallest_k2 = min(table(g2)), smallest_k3 = min(table(r$stab[["3"]]$ref)),
    assigned_k2 = sum(a2 > 0), unassigned_k2 = sum(a2 == 0), ari_reference_k2 = round(ari(g2, rD$cluster), 3),
    n_both = sum(both), same_label = sum(a2[both] == rD$assigned[both]), ari_both = round(ari(a2[both], rD$assigned[both]), 3), check.names = FALSE) }))
names(alg)[6:10] <- paste0("weakest_k", 2:6); alg <- alg[order(match(alg$method, c("ward", meth)), alg$rho), ]; rownames(alg) <- NULL
cat("\n(b) Clustering algorithm by perturbation dependence, B =", B, "\n"); print(alg, row.names = FALSE)
write.csv(alg, tab("sens_algorithm.csv"), row.names = FALSE)

# ---- (c) do the three algorithms give the same k = 2 label, and whom does that concern? ----------
gW <- rD$cluster; gP <- ordlab(ra$pam$stab[["2"]]$ref); gK <- ordlab(ra$kmeans$stab[["2"]]$ref)
same3 <- gW == gP & gW == gK
agree <- as.data.frame.matrix(table(factor(rD$assigned, levels = c(1, 2, 0), labels = c("Assigned stalled", "Assigned decline", "Unassigned")),
                                    factor(same3, levels = c(TRUE, FALSE), labels = c("same", "differs"))))
agree$pam_stalled <- tapply(gP == 1, factor(rD$assigned, levels = c(1, 2, 0)), sum); agree$kmeans_stalled <- tapply(gK == 1, factor(rD$assigned, levels = c(1, 2, 0)), sum)
agree$median_p <- round(tapply(rD$p, factor(rD$assigned, levels = c(1, 2, 0)), median), 2)
agree$median_change <- round(tapply(pc, factor(rD$assigned, levels = c(1, 2, 0)), median), 1)
cat("\n(c) Agreement of the k = 2 label across Ward, PAM and k-means, by primary assignment\n"); print(agree)
cat(sprintf("same label under all three: %d of %d; of the %d that differ, %d are unassigned in the primary analysis; assigned but differing: %s\n",
    sum(same3), n, sum(!same3), sum(!same3 & rD$assigned == 0), paste(sub(" \\(.*", "", meta$name[!same3 & rD$assigned > 0]), collapse = ", ")))
cat(sprintf("median primary membership probability: same label %.3f, differing %.3f\n", median(rD$p[same3]), median(rD$p[!same3])))
w3 <- ordlab(rD$stab[["3"]]$ref); p3 <- ordlab(ra$pam$stab[["3"]]$ref); k3 <- ordlab(ra$kmeans$stab[["3"]]$ref)
k3tab <- data.frame(algorithm = c("Ward", "PAM", "k-means"), rbind(table(w3), table(p3), table(k3)),
  rbind(round(tapply(pc, w3, median)), round(tapply(pc, p3, median)), round(tapply(pc, k3, median))),
  rbind(round(tapply(rD$stab[["3"]]$p, w3, median), 2), round(tapply(ra$pam$stab[["3"]]$p, p3, median), 2), round(tapply(ra$kmeans$stab[["3"]]$p, k3, median), 2)))
names(k3tab) <- c("algorithm", paste0("n_", 1:3), paste0("change_", 1:3), paste0("median_p_", 1:3))
k3tab$high_income_in_3 <- c(sum(meta$region[w3 == 3] == "HI"), sum(meta$region[p3 == 3] == "HI"), sum(meta$region[k3 == 3] == "HI"))
cat("\nThree cluster partitions (clusters ordered from least to most decline)\n"); print(k3tab, row.names = FALSE)
cat(sprintf("fastest decline cluster: in all three %d, in at least one %d; Ward's third cluster inside PAM's %d/%d, inside k-means' %d/%d\n",
    sum(w3 == 3 & p3 == 3 & k3 == 3), sum(w3 == 3 | p3 == 3 | k3 == 3), sum(w3 == 3 & p3 == 3), sum(w3 == 3), sum(w3 == 3 & k3 == 3), sum(w3 == 3)))
cat(sprintf("ARI at k = 3: Ward-PAM %.3f, Ward-kmeans %.3f, PAM-kmeans %.3f; at k = 2: %.3f, %.3f, %.3f\n",
    ari(w3, p3), ari(w3, k3), ari(p3, k3), ari(gW, gP), ari(gW, gK), ari(gP, gK)))
cat(sprintf("mean interval width (log scale SD) in Ward's three clusters: %s\n", paste(round(tapply(rowMeans(S), w3, median), 3), collapse = ", ")))
write.csv(agree, tab("sens_agreement_k2.csv")); write.csv(k3tab, tab("sens_three_clusters.csv"), row.names = FALSE)

# ---- (d) scale of the perturbation ----------------------------------------------------------------
# Weakest cluster median membership at k = 2, 3, 4 when every published SD is multiplied by cfac.
jd <- expand.grid(cfac = c(0.1, 0.25, 0.5, 2), method = c("ward", "pam", "kmeans"), stringsAsFactors = FALSE); jd$seed <- seed0 + 400 + seq_len(nrow(jd))
rd <- mclapply(seq_len(nrow(jd)), function(j) { set.seed(jd$seed[j]); cl <- clusterer(jd$method[j])
  sapply(2:4, function(k) { st <- stability(Y, S * jd$cfac[j], rD$hc, k, B, "value", cl = cl, ref = cl(sh, k)); min(tapply(st$p, st$ref, median)) }) }, mc.cores = nc, mc.preschedule = FALSE)
scl <- cbind(jd[, c("method", "cfac")], round(do.call(rbind, rd), 2)); names(scl)[3:5] <- paste0("weakest_k", 2:4)
base <- rbind(rD$weakest[1:3], ra$pam$weakest[1:3], ra$kmeans$weakest[1:3])
scl <- rbind(scl, data.frame(method = c("ward", "pam", "kmeans"), cfac = 1, weakest_k2 = round(base[, 1], 2), weakest_k3 = round(base[, 2], 2), weakest_k4 = round(base[, 3], 2)))
scl <- scl[order(match(scl$method, c("ward", "pam", "kmeans")), scl$cfac), ]; rownames(scl) <- NULL
cat("\n(d) Weakest cluster median membership when published SDs are multiplied by cfac\n"); print(scl, row.names = FALSE)
write.csv(scl, tab("sens_scale.csv"), row.names = FALSE)

# ---- tables for the manuscript and supplement --------------------------------------------------
f2 <- function(x) sprintf("%.2f", x); mname <- c(ward = "Ward", pam = "PAM", kmeans = "$k$-means", average = "Average linkage")
kk <- function(k) ifelse(is.na(k), "none", as.character(k))
# Table 6 (main text): number of clusters and assignment by dependence and algorithm
t6 <- alg[alg$method != "average", ]
t6 <- data.frame(Algorithm = ifelse(duplicated(t6$method), "", mname[t6$method]), `$\\rho$` = sprintf("%.2f", t6$rho), `$k$ (gap)` = t6$k_gap, `$k$ (stability)` = kk(t6$k_stab),
  `$k = 2$` = f2(t6$weakest_k2), `$k = 3$` = f2(t6$weakest_k3), `$k = 4$` = f2(t6$weakest_k4), `Assigned at $k = 2$` = t6$assigned_k2,
  `Same label as primary` = sprintf("%d of %d", t6$same_label, t6$n_both), check.names = FALSE)
tex_table(t6, tab("tab6.tex"), "lrrrrrrrr")
# Table S8: dependence, Ward, all rho
s8 <- data.frame(`$\\rho$` = dep$rho, `Selected $k$` = dep$k_stab, `$k = 2$` = f2(dep$weakest_k2), `$k = 3$` = f2(dep$weakest_k3), `$k = 4$` = f2(dep$weakest_k4),
  `$k = 5$` = f2(dep$weakest_k5), `$k = 6$` = f2(dep$weakest_k6), Stalled = dep$stalled, Decline = dep$decline, Unassigned = dep$unassigned,
  `Status changed` = dep$status_changed, `Spearman` = sprintf("%.2f", dep$spearman_with_rho0), check.names = FALSE)
v8 <- vdec[match(suppressWarnings(as.numeric(sub(" .*", "", dep$rho))), vdec$rho), ]
s8$`Total variance` <- sprintf("%.2f", v8$total_rel); s8$`Variance along $\\mathbf{d}$` <- sprintf("%.2f", v8$along_rel)
tex_table(s8, tab("tabS8.tex"), "lrrrrrrrrrrrrr")
# Table S9: algorithm by dependence, all k
s9 <- data.frame(Algorithm = ifelse(duplicated(alg$method), "", mname[alg$method]), `$\\rho$` = sprintf("%.2f", alg$rho), `$k$ (sil.)` = alg$k_sil, `$k$ (gap)` = alg$k_gap,
  `$k$ (stability)` = kk(alg$k_stab), `$k = 2$` = f2(alg$weakest_k2), `$k = 3$` = f2(alg$weakest_k3), `$k = 4$` = f2(alg$weakest_k4), `$k = 5$` = f2(alg$weakest_k5),
  `$k = 6$` = f2(alg$weakest_k6), `Smallest cluster, $k = 2$` = alg$smallest_k2, `ARI with Ward, $k = 2$` = sprintf("%.2f", alg$ari_reference_k2), check.names = FALSE)
tex_table(s9, tab("tabS9.tex"), "lrrrrrrrrrrr")
# Table S10: agreement of the k = 2 label across algorithms by primary assignment
s10 <- data.frame(`Primary assignment` = rownames(agree), Countries = agree$same + agree$differs, `Same under all three` = agree$same, Differs = agree$differs,
  `PAM: stalled` = agree$pam_stalled, `$k$-means: stalled` = agree$kmeans_stalled, `Median $p_i$` = f2(agree$median_p),
  `Median change (%)` = sprintf("%.1f", agree$median_change), check.names = FALSE)
tex_table(s10, tab("tabS10.tex"), "lrrrrrrr")
# Table S11: scale of the perturbation
s11 <- data.frame(Algorithm = ifelse(duplicated(scl$method), "", mname[scl$method]), `Multiplier $c$` = sprintf("%.2f", scl$cfac),
  `$k = 2$` = f2(scl$weakest_k2), `$k = 3$` = f2(scl$weakest_k3), `$k = 4$` = f2(scl$weakest_k4), check.names = FALSE)
tex_table(s11, tab("tabS11.tex"), "lrrrr")

# ---- Figure S3 -----------------------------------------------------------------------------------
# (a) weakest cluster median by k for each rho, Ward; (b) the same for each algorithm at rho = 0 and 0.9;
# (c) membership probability at k = 2 under rho = 0.9 against rho = 0, coloured by the primary cluster.
png(fig("figS3.png"), width = 2400, height = 860, res = 220); layout(matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE), heights = c(1, 0.2)); par(mar = c(4, 4.4, 1, 1))
sel <- c("0", "0.70", "0.90", "0.95", "1.00"); gr <- c("black", "grey25", "grey45", "grey60", "grey75"); pch_a <- c(16, 17, 15, 18, 1)
plot(NA, xlim = c(2, 6), ylim = c(0, 1), xlab = "Number of clusters k", ylab = "Weakest cluster median membership", las = 1); abline(h = pi_, lty = 2)
for (i in seq_along(sel)) lines(2:6, rr[[sel[i]]]$weakest, type = "b", pch = pch_a[i], col = gr[i], lwd = 1.3)
plot(NA, xlim = c(2, 6), ylim = c(0, 1), xlab = "Number of clusters k", ylab = "Weakest cluster median membership", las = 1); abline(h = pi_, lty = 2)
mb <- c("ward", "pam", "kmeans"); pch_b <- c(16, 17, 15)
for (i in seq_along(mb)) for (q in c(0, 0.9)) { a <- allruns[[which(sapply(allruns, function(a) a$method == mb[i] && a$rho == q))]]
  lines(2:6, a$r$weakest, type = "b", pch = pch_b[i], lty = if (q == 0) 1 else 3, col = if (q == 0) "black" else "grey45", lwd = 1.3) }
colr <- ifelse(rD$assigned == 0, gU, pal[rD$cluster])
plot(p2[, "0"], p2[, "0.90"], col = adjustcolor(colr, 0.8), pch = 16, las = 1, xlab = expression("Membership probability at k = 2, " * rho == 0), ylab = expression("Membership probability at k = 2, " * rho == 0.9))
abline(0, 1, col = "grey70", lty = 3); abline(h = tau, v = tau, lty = 2, col = "grey60")
par(mar = c(0, 0, 0, 0)); plot.new(); legend("center", ncol = 3, bty = "n", pch = pch_a, col = gr, lty = 1, cex = 0.85, legend = parse(text = paste0("rho == ", c("0", "0.7", "0.9", "0.95", "1"))))
plot.new(); legend("center", ncol = 2, bty = "n", cex = 0.85, pch = c(pch_b, NA, NA), lty = c(NA, NA, NA, 1, 3), col = c(rep("black", 4), "grey45"),
                   legend = c("Ward", "PAM", "k-means", expression(rho == 0), expression(rho == 0.9)))
plot.new(); legend("center", ncol = 1, bty = "n", cex = 0.85, pch = 16, col = c(pal[1:2], gU), legend = c(lab_short, "Unassigned")); dev.off()

save(rr, ra, rb, jb, allruns, p2, as2, dep, vdec, alg, agree, k3tab, scl, rhos, meth, same3, file = "output/sensitivity.RData")
cat("\nDONE\n")
