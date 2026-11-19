# framework.R
# Uncertainty aware trajectory clustering of modelled estimates.
# Functions implementing equations (1) to (13) of the manuscript (v1.2 numbering). Base R plus the cluster package.
# Source this file; nothing is executed on load.

suppressPackageStartupMessages(library(cluster))

#' Load a GBD panel into log rate and log scale SD matrices.
#' @param f path to a panel CSV with columns location_id, year, rate_age_std_combined,
#'   rate_combined, rate_lower_combined, rate_upper_combined
#' @param yrs integer vector of years to keep
#' @param ids character vector of location_id values fixing the row order
#' @return list(Y = n x T matrix of log age standardised rates,
#'              S = n x T matrix of log scale SD from the published interval, equation (1))
load_panel <- function(f, yrs, ids) {
  d <- read.csv(f); d <- d[d$year %in% yrs, ]
  W <- function(v) {
    m <- reshape(d[, c("location_id", "year", v)], idvar = "location_id", timevar = "year", direction = "wide")
    rn <- m$location_id; m <- as.matrix(m[, -1]); rownames(m) <- rn
    colnames(m) <- sub(paste0(v, "."), "", colnames(m)); m[ids, order(as.integer(colnames(m)))] }
  Y <- log(W("rate_age_std_combined"))
  S <- (log(W("rate_upper_combined")) - log(W("rate_lower_combined"))) / (2 * 1.96)
  list(Y = Y, S = S)
}

#' Adjusted Rand index (Hubert and Arabie 1985).
ari <- function(a, b) {
  t <- table(a, b); nn <- sum(t); c2 <- function(x) x * (x - 1) / 2
  sij <- sum(c2(t)); si <- sum(c2(rowSums(t))); sj <- sum(c2(colSums(t)))
  e <- si * sj / c2(nn); (sij - e) / ((si + sj) / 2 - e)
}

#' Ward (ward.D2) hierarchical clustering cut at k clusters.
hcl <- function(X, k) cutree(hclust(dist(X), "ward.D2"), k)

#' Level and shape decomposition, equation (3). Returns the shape matrix (rows sum to zero).
shape_of <- function(Y) sweep(Y, 1, rowMeans(Y))

#' Standard normal perturbations with AR(1) dependence across years within a country.
#' Each row is a stationary AR(1) series with unit marginal variance and lag one correlation rho,
#' so every value keeps the marginal standard deviation implied by its published interval whatever
#' rho is. rho = 0 gives the independent draws of equation (6) and consumes the random number
#' stream exactly as matrix(rnorm(n * T), n) does; rho = 1 gives one common shock per country.
#' @param n,T number of countries and years; rho lag one correlation in [0, 1]
ar1_noise <- function(n, T, rho = 0) {
  Z <- matrix(rnorm(n * T), n); if (rho == 0) return(Z)
  U <- Z; for (t in 2:T) U[, t] <- rho * U[, t - 1] + sqrt(1 - rho^2) * Z[, t]
  U
}

#' Clustering algorithms available to the framework. Each returns function(X, k) -> integer labels.
#' ward: Ward's criterion (ward.D2). pam: partitioning around medoids on Euclidean distance.
#' kmeans: Hartigan and Wong k-means, 25 random starts. average: average linkage (UPGMA).
clusterer <- function(method = c("ward", "pam", "kmeans", "average")) {
  switch(match.arg(method),
    ward    = hcl,
    pam     = function(X, k) cluster::pam(X, k, cluster.only = TRUE),
    kmeans  = function(X, k) kmeans(X, k, nstart = 25, iter.max = 100)$cluster,
    average = function(X, k) cutree(hclust(dist(X), "average"), k))
}

#' Membership probabilities under value perturbation (equations 6 to 8) or country bootstrap.
#' @param Y,S log rate and log SD matrices
#' @param hc reference hclust object on shape_of(Y); ignored when ref is supplied
#' @param k number of clusters
#' @param B number of draws
#' @param mode "value" perturbs every value within its published SD; "country" resamples countries
#' @param rho lag one correlation of the perturbations across years (value mode only); 0 = independent
#' @param cl clustering function(X, k) used to recluster each draw; default Ward
#' @param ref reference labels; default cutree(hc, k)
#' @return list(ref = reference labels, p = membership probability per country)
stability <- function(Y, S, hc, k, B, mode = c("value", "country"), rho = 0, cl = hcl, ref = NULL) {
  mode <- match.arg(mode); if (is.null(ref)) ref <- cutree(hc, k); shape <- shape_of(Y)
  agree <- numeric(nrow(Y)); drawn <- numeric(nrow(Y))
  for (b in seq_len(B)) {
    if (mode == "value") {
      Yb <- Y + ar1_noise(nrow(Y), ncol(Y), rho) * S; idx <- seq_len(nrow(Y)); clb <- cl(shape_of(Yb), k)
    } else {
      idx <- sort(unique(sample(nrow(Y), replace = TRUE))); clb <- cl(shape[idx, ], k) }
    lab <- apply(table(clb, ref[idx]), 1, which.max)            # equation (7), majority overlap alignment
    agree[idx] <- agree[idx] + (lab[clb] == ref[idx]); drawn[idx] <- drawn[idx] + 1 }
  list(ref = ref, p = agree / drawn)
}

#' The full framework: shape clustering, three k selection criteria, membership and assignment.
#' @param tau assignment threshold, equation (9); pi_ stability rule, equation (10)
#' @param kmax largest k examined; gapB reference sets for the gap statistic, equation (11)
#' @param rho lag one correlation of the value perturbations across years; 0 reproduces equation (6)
#' @param method clustering algorithm, see clusterer(); "ward" is the manuscript's primary choice
#' @return list with shape, hc, sil, gap, gap_se, k_sil, k_gap, k, weakest, unas, cluster, p, assigned, stab.
#'   Clusters are relabelled so that cluster 1 has the least decline over the window. k is NA when no
#'   partition meets the stability rule; cluster, p and assigned are then reported at k = 2.
framework <- function(Y, S, B = 1000, tau = 0.70, pi_ = 0.80, kmax = 6, mode = "value", gapB = 100,
                      rho = 0, method = "ward") {
  cl <- clusterer(method); shape <- shape_of(Y); hc <- hclust(dist(shape), "ward.D2"); dsh <- dist(shape)
  parts <- lapply(2:kmax, function(k) if (method == "ward") cutree(hc, k) else cl(shape, k)); names(parts) <- 2:kmax
  sil <- sapply(parts, function(g) mean(silhouette(g, dsh)[, 3]))
  gap <- clusGap(shape, function(x, k) list(cluster = cl(x, k)), K.max = kmax, B = gapB, verbose = FALSE)
  k_gap <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], "Tibs2001SEmax")
  stab <- lapply(2:kmax, function(k) stability(Y, S, hc, k, B, mode, rho = rho, cl = cl, ref = parts[[as.character(k)]]))
  names(stab) <- 2:kmax
  med <- sapply(stab, function(s) min(tapply(s$p, s$ref, median)))
  unas <- sapply(stab, function(s) mean(s$p < tau))
  ok <- as.integer(names(stab))[med >= pi_]; k <- if (length(ok)) max(ok) else NA_integer_; kr <- if (is.na(k)) 2L else k
  ref <- stab[[as.character(kr)]]$ref; p <- stab[[as.character(kr)]]$p
  pc <- Y[, ncol(Y)] - Y[, 1]; o <- order(tapply(pc, ref, median), decreasing = TRUE)
  map <- setNames(seq_along(o), o); ref <- unname(map[as.character(ref)])
  list(shape = shape, hc = hc, sil = unname(sil), gap = gap$Tab[, "gap"], gap_se = gap$Tab[, "SE.sim"],
       k_sil = unname(which.max(sil)) + 1, k_gap = k_gap, k = k, weakest = med, unas = unas,
       cluster = ref, p = p, assigned = ifelse(p >= tau, ref, 0), stab = stab, rho = rho, method = method)
}

#' Moran's I with a supplied row standardised weight matrix, equation (13).
moranI <- function(z, W) { z <- z - mean(z); (length(z) / sum(W)) * as.numeric(t(z) %*% W %*% z) / sum(z^2) }

#' k nearest neighbour row standardised weights from a distance matrix.
knn_weights <- function(D, kn) t(apply(D, 1, function(x) { w <- numeric(length(x)); w[order(x)[2:(kn + 1)]] <- 1 / kn; w }))

#' Great circle distance matrix (km) from a two column matrix of lon, lat.
gc_dist <- function(C) { r <- pi / 180
  outer(seq_len(nrow(C)), seq_len(nrow(C)), Vectorize(function(i, j) {
    a <- sin((C[j, 2] - C[i, 2]) * r / 2)^2 + cos(C[i, 2] * r) * cos(C[j, 2] * r) * sin((C[j, 1] - C[i, 1]) * r / 2)^2
    2 * 6371 * asin(sqrt(a)) })) }

#' Country centroids from the maps package world database (largest polygon per country).
country_centroids <- function(names, xw) {
  mp <- maps::map("world", plot = FALSE, fill = TRUE); nm <- sub(":.*", "", mp$names)
  grp <- cumsum(c(TRUE, is.na(mp$x[-length(mp$x)])))
  t(sapply(names, function(g) { q <- if (g %in% names(xw)) xw[[g]] else g; if (is.na(q)) return(c(NA, NA))
    i <- which(nm == q); if (!length(i)) return(c(NA, NA)); ok <- grp %in% i & !is.na(mp$x)
    big <- as.integer(names(which.max(table(grp[ok])))); s <- grp == big & !is.na(mp$x)
    c(mean(mp$x[s]), mean(mp$y[s])) })) }

#' Median (IQR) formatter and a minimal LaTeX tabular writer.
q0 <- function(x, d = 0) sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), median(x), quantile(x, .25), quantile(x, .75))
tex_table <- function(df, file, align = NULL) {
  esc <- function(x) gsub("%", "\\\\%", gsub("&", "\\\\&", as.character(x)))
  al <- if (is.null(align)) paste0("l", strrep("r", ncol(df) - 1)) else align
  out <- c(paste0("\\begin{tabular}{", al, "}"), "\\toprule", paste0(paste(esc(names(df)), collapse = " & "), " \\\\"), "\\midrule",
           apply(df, 1, function(x) paste0(paste(esc(x), collapse = " & "), " \\\\")), "\\bottomrule", "\\end{tabular}")
  writeLines(out, file) }

# GBD name to maps::world name crosswalk, and labels
gbd_to_maps <- list("Antigua and Barbuda"="Antigua","Brunei Darussalam"="Brunei","Czechia"="Czech Republic","Republic of Korea"="South Korea",
  "Democratic People's Republic of Korea"="North Korea","Russian Federation"="Russia","Saint Kitts and Nevis"="Saint Kitts",
  "Trinidad and Tobago"="Trinidad","United Kingdom"="UK","United States Virgin Islands"="Virgin Islands, US","United States of America"="USA",
  "Bolivia (Plurinational State of)"="Bolivia","Cabo Verde"="Cape Verde","Congo"="Republic of Congo","Côte d'Ivoire"="Ivory Coast",
  "Eswatini"="Swaziland","Lao People's Democratic Republic"="Laos","Micronesia (Federated States of)"="Micronesia",
  "United Republic of Tanzania"="Tanzania","Viet Nam"="Vietnam","Syrian Arab Republic"="Syria","Iran (Islamic Republic of)"="Iran",
  "Republic of Moldova"="Moldova","Saint Vincent and the Grenadines"="Saint Vincent","Tuvalu"=NA,"Türkiye"="Turkey",
  "Venezuela (Bolivarian Republic of)"="Venezuela")
regname <- c(CEECA="Central Europe, Eastern Europe and Central Asia", HI="High income", LAC="Latin America and Caribbean",
             NAME="North Africa and Middle East", SA="South Asia", SEAO="Southeast Asia, East Asia and Oceania", SSA="Sub-Saharan Africa")
incname <- c(low="Low", lmi="Lower middle", umi="Upper middle", hi="High")
pal <- c("#E69F00", "#56B4E9", "#009E73"); gU <- "grey40"; gN <- "grey88"
lab_short <- c("Cluster 1: stalled", "Cluster 2: sustained decline")
