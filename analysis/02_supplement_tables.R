# 02_supplement_tables.R
# Renders Tables S1 to S7 of the supplement from output/final.RData. Run after 01_run_analysis.R.
# Build the supplementary tables as LaTeX from the final analysis outputs.
root <- if (file.exists("R/framework.R")) "." else ".."; setwd(root); source("R/framework.R"); load("output/final.RData")
esc <- function(x) gsub("%", "\\\\%", gsub("&", "\\\\&", as.character(x)))
tex <- function(df, file, align, long = FALSE, caption = NULL, label = NULL) {
  hdr <- paste0(paste(esc(names(df)), collapse = " & "), " \\\\")
  rows <- apply(df, 1, function(x) paste0(paste(esc(x), collapse = " & "), " \\\\"))
  out <- if (long) c(paste0("\\begin{longtable}{", align, "}"), paste0("\\caption{", caption, "}\\label{", label, "}\\\\"), "\\toprule", hdr, "\\midrule", "\\endfirsthead",
      paste0("\\multicolumn{", ncol(df), "}{l}{\\small\\textit{Table \\thetable{} continued}}\\\\"), "\\toprule", hdr, "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot", rows, "\\end{longtable}")
    else c(paste0("\\begin{tabular}{", align, "}"), "\\toprule", hdr, "\\midrule", rows, "\\bottomrule", "\\end{tabular}")
  writeLines(out, file) }

# S1: all countries
s1 <- read.csv("output/tables/tableS1_country_membership.csv", check.names = FALSE)
s1$Region <- c("Central Europe, Eastern Europe and Central Asia"="CEECA", "High income"="HI", "Latin America and Caribbean"="LAC", "North Africa and Middle East"="NAME", "South Asia"="SA", "Southeast Asia, East Asia and Oceania"="SEAO", "Sub-Saharan Africa"="SSA")[s1$Region]
s1$Country <- gsub(" \\(.*\\)", "", s1$Country); s1$Country <- sub("Democratic People's Republic of Korea", "Korea, DPR", s1$Country)
s1$Income <- c(Low="L", `Lower middle`="LM", `Upper middle`="UM", High="H")[s1$Income]; s1 <- s1[, c("Country","Region","Income","SDI","Change DALYs (%)","Change deaths (%)","Peak year","Interval width","Cluster (DALYs)","p (DALYs)","Cluster (deaths)","p (deaths)","Cluster (2000 window)","p (2000 window)")]
names(s1) <- c("Country","Region","Inc.","SDI","$\\Delta$DALY","$\\Delta$deaths","Peak","$w_i$","$a_i$","$p_i$","$a_i^{\\mathrm{d}}$","$p_i^{\\mathrm{d}}$","$a_i^{\\mathrm{w}}$","$p_i^{\\mathrm{w}}$")
tex(s1, "output/tables/tabS1.tex", "llrrrrrrrrrrrr", long = TRUE,
    caption = "Country level results. Region codes: CEECA, Central Europe, Eastern Europe and Central Asia; HI, high income; LAC, Latin America and Caribbean; NAME, North Africa and Middle East; SA, South Asia; SEAO, Southeast Asia, East Asia and Oceania; SSA, sub-Saharan Africa. Inc.: World Bank income group (L low, LM lower middle, UM upper middle, H high). $\\Delta$DALY and $\\Delta$deaths: percentage change in the age-standardised rate, 1990 to 2023. Peak: year of maximum centred DALY trajectory. $w_i$: mean log-scale standard deviation of the published interval. $a_i$: assigned cluster (1 stalled, 2 sustained decline, 0 unassigned) and $p_i$: membership probability, for the primary DALY analysis; superscript d, the deaths analysis; superscript w, the 2000 to 2023 window.", label = "S1")

# S2: Moran's I by neighbours
s2 <- read.csv("output/tables/tableS2_moran_neighbours.csv"); names(s2) <- c("Neighbours $\\kappa$", "Moran's $I$", "Permutation $p$"); s2[,2] <- round(s2[,2],3)
tex(s2, "output/tables/tabS2.tex", "rrr")

# S3: threshold sensitivity
s3a <- data.frame(`$\\pi$` = c(0.70,0.75,0.80,0.85,0.90), check.names = FALSE); s3a$`Selected $k$` <- sapply(s3a[[1]], function(p) max((2:6)[rD$weakest >= p]))
s3b <- data.frame(`$\\tau$` = c(0.60,0.70,0.80,0.90), check.names = FALSE); s3b$`Unassigned` <- sapply(s3b[[1]], function(t) sum(rD$p < t)); s3b$`Assigned stalled` <- sapply(s3b[[1]], function(t) sum(rD$p >= t & rD$cluster==1)); s3b$`Assigned decline` <- sapply(s3b[[1]], function(t) sum(rD$p >= t & rD$cluster==2))
tex(s3a, "output/tables/tabS3a.tex", "rr"); tex(s3b, "output/tables/tabS3b.tex", "rrrr")

# S4: k selection diagnostics for every robustness run
runs <- list(`DALYs 1990 to 2023 (primary)`=rD, `Deaths 1990 to 2023`=rM, `DALYs 2000 to 2023`=r2, `DALYs, narrowest half (n = 101)`=rs, `Country bootstrap`=rB)
s4 <- do.call(rbind, lapply(names(runs), function(nm) { r <- runs[[nm]]; kk <- 2:(length(r$sil)+1)
  data.frame(Analysis = c(nm, rep("", length(kk)-1)), k = kk, Silhouette = round(r$sil,3), Gap = sprintf("%.3f (%.3f)", r$gap[kk], r$gap_se[kk]),
             `Weakest median $p$` = round(r$weakest,2), `Share unassigned` = round(r$unas,2), check.names = FALSE) }))
tex(s4, "output/tables/tabS4.tex", "lrrrrr")

# S5: narrow half three cluster solution by region
s5 <- as.data.frame.matrix(table(regname[meta$region[keep]], factor(rs$assigned, levels=c(1,2,3,0), labels=c("Cluster 1 (stalled)","Cluster 2 (rise and fall)","Cluster 3 (steady decline)","Unassigned")))); s5$Total <- rowSums(s5); s5 <- rbind(s5, Total=colSums(s5)); s5 <- cbind(`Super region`=rownames(s5), s5)
tex(s5, "output/tables/tabS5.tex", "lrrrrr")
s5b <- data.frame(Cluster = c("Cluster 1 (stalled)","Cluster 2 (rise and fall)","Cluster 3 (steady decline)"), Members = sapply(1:3, function(c) paste(gsub(" \\(.*\\)","",meta$name[keep][rs$assigned==c]), collapse="; ")))
writeLines(c("\\begin{description}", paste0("\\item[", s5b$Cluster, "] ", s5b$Members, "."), "\\end{description}"), "output/tables/tabS5b.tex")

# S6: unassigned countries in the primary analysis
u <- meta[meta$assigned==0,]; u <- u[order(u$p),]
s6 <- data.frame(Country = gsub(" \\(.*\\)","",u$name), Region = regname[u$region], SDI = u$sdi, `Change (%)` = round(pc[meta$assigned==0][order(meta$p[meta$assigned==0])],1), `Reference cluster` = u$cluster, `$p_i$` = round(u$p,2), check.names = FALSE)
tex(s6, "output/tables/tabS6.tex", "llrrrr")
cat("wrote S1 to S6\n")
