# 01_run_analysis.R
# Primary analysis and all robustness runs. One seed, B = 1000. Produces every table and figure in
# output/. Run from the repository root:  Rscript analysis/01_run_analysis.R
# Runtime about 75 s. Requires R >= 4.3 with packages cluster, nnet, maps.
root <- if (file.exists("R/framework.R")) "." else if (file.exists("../R/framework.R")) ".." else stop("run from the repository root")
setwd(root); source("R/framework.R"); suppressPackageStartupMessages({library(nnet); library(maps)})
set.seed(20261003); B <- as.integer(Sys.getenv("GBD_B", "1000"))   # GBD_B=200 for a quick check
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE); dir.create("output/tables", showWarnings = FALSE)
fig <- function(f) file.path("output/figures", f); tab <- function(f) file.path("output/tables", f)

meta0 <- read.csv("data/analysis_set.csv"); sdi <- read.csv("data/sdi_2023.csv")
meta0$sdi <- sdi$sdi[match(meta0$location_id, sdi$location_id)]
ids <- as.character(meta0$location_id); n <- length(ids); xw <- gbd_to_maps
# ---- 1. primary: DALYs 1990 to 2023 ----------------------------------------------------------
yrs <- 1990:2023; D <- load_panel("data/panel_cvd.csv", yrs, ids); Y <- D$Y; S <- D$S
rD <- framework(Y, S, B=B); meta <- meta0; meta$cluster <- rD$cluster; meta$p <- rD$p; meta$assigned <- rD$assigned
level <- rowMeans(Y); shape <- rD$shape; pc <- 100*(exp(Y[,"2023"]-Y[,"1990"])-1); peak <- yrs[apply(shape,1,which.max)]
cat(sprintf("PRIMARY: k sil %d, gap %d, ua %d | assigned %d, unassigned %d | R2 level %.3f, ARI level tertiles %.3f\n",
    rD$k_sil, rD$k_gap, rD$k, sum(meta$assigned>0), sum(meta$assigned==0),
    summary(lm(level~factor(meta$cluster)))$r.squared, ari(meta$cluster, cut(level, quantile(level,0:3/3), include.lowest=TRUE))))
cat("Weakest medians k=2..6:", round(rD$weakest,2), " unassigned share:", round(rD$unas,2), "\n")
cat("Silhouette:", round(rD$sil,3), " Gap:", round(rD$gap,3), "\n")
cat("Unassigned:", paste(meta$name[meta$assigned==0], collapse="; "), "\n")
pam2 <- pam(shape, 2)$clustering; cat(sprintf("ARI(Ward, PAM) at k=2 = %.3f\n", ari(rD$cluster, pam2)))

# ---- 2. deaths, 2000 window, bootstrap -----------------------------------------------------
Mt <- load_panel("data/panel_cvd_deaths.csv", yrs, ids); rM <- framework(Mt$Y, Mt$S, B=B); pcM <- 100*(exp(Mt$Y[,"2023"]-Mt$Y[,"1990"])-1)
D2 <- load_panel("data/panel_cvd.csv", 2000:2023, ids); r2 <- framework(D2$Y, D2$S, B=B); pc2 <- 100*(exp(D2$Y[,"2023"]-D2$Y[,"2000"])-1)
rB <- framework(Y, S, B=B, mode="country")
cat(sprintf("DEATHS: k sil %d gap %d ua %d; weakest %s; assigned %d; ARI vs DALY (assigned both) %.3f; switchers %d\n", rM$k_sil, rM$k_gap, rM$k,
    paste(round(rM$weakest,2),collapse=" "), sum(rM$assigned>0), {b<-rD$assigned>0&rM$assigned>0; ari(rD$assigned[b],rM$assigned[b])}, sum((rD$assigned==1&rM$assigned==2)|(rD$assigned==2&rM$assigned==1))))
print(table(DALY=rD$assigned, deaths=rM$assigned))
cat(sprintf("WINDOW 2000: k sil %d gap %d ua %d; weakest %s; assigned %d; HI in decline %d/36\n", r2$k_sil, r2$k_gap, r2$k, paste(round(r2$weakest,2),collapse=" "), sum(r2$assigned>0), sum(r2$assigned[meta$region=="HI"]==2)))
print(table(full=rD$assigned, window2000=r2$assigned)); print(round(tapply(pc2, r2$assigned, median),1))
sl <- function(a,b) (Y[,as.character(b)]-Y[,as.character(a)])/(b-a)*100; hi <- meta$region=="HI"
cat(sprintf("HI annual %% change: 2000-2010 %.2f, 2013-2023 %.2f (median)\n", median(sl(2000,2010)[hi]), median(sl(2013,2023)[hi])))
cat(sprintf("BOOTSTRAP: weakest %s -> k %d; VALUE: weakest %s -> k %d; rho(p) at k=3 %.3f\n", paste(round(rB$weakest,2),collapse=" "), rB$k,
    paste(round(rD$weakest,2),collapse=" "), rD$k, cor(rB$stab[["3"]]$p, rD$stab[["3"]]$p, method="spearman")))

# ---- 3. interval width -----------------------------------------------------------------------
w <- rowMeans(S); cat("Width by cluster:", round(tapply(w, meta$assigned, median),3), " by region:", paste(names(sort(tapply(w,meta$region,median))), round(sort(tapply(w,meta$region,median)),3)), "\n")
keep <- w <= median(w); rs <- framework(Y[keep,], S[keep,], B=B, kmax=5);
cat(sprintf("NARROW HALF (n=%d; SSA %d): k ua %d, weakest %s; ARI vs reference %.3f\n", sum(keep), sum(meta$region[keep]=="SSA"), rs$k, paste(round(rs$weakest,2),collapse=" "), ari(rs$assigned, rD$assigned[keep])))
print(table(narrow=rs$assigned, reference=rD$assigned[keep])); print(table(narrow=rs$assigned, region=meta$region[keep]))
a <- data.frame(y=factor(rD$assigned), sdi=meta$sdi, w=w, region=meta$region)[rD$assigned>0,]
g0 <- glm(y~1,binomial,a); g1 <- glm(y~sdi,binomial,a); g2 <- glm(y~sdi+region,binomial,a); g3 <- glm(y~sdi+w,binomial,a); g4 <- glm(y~sdi+region+w,binomial,a)
cat(sprintf("OR decline per 0.1 SDI %.2f; width LR p (vs SDI) %.3g, OR per 0.1 width %.2f; width LR p (vs SDI+region) %.3g\n",
    exp(coef(g1)["sdi"]*0.1), pchisq(deviance(g1)-deviance(g3),1,lower.tail=FALSE), exp(coef(g3)["w"]*0.1), pchisq(deviance(g2)-deviance(g4),1,lower.tail=FALSE)))

# ---- 4. spatial and attribution --------------------------------------------------------------
cent <- country_centroids(meta$name, xw)
ok <- complete.cases(cent) & meta$assigned>0; C <- cent[ok,]; r <- pi/180
Dm <- gc_dist(C)
mor <- t(sapply(c(4,5,6,8,10), function(kn) { W <- knn_weights(Dm, kn)
  z <- as.numeric(meta$assigned[ok]==1); I <- moranI(z,W); nul <- replicate(999, moranI(sample(z),W)); c(knn=kn, I=I, p=mean(nul>=I)) }))
print(round(mor,3))
m <- data.frame(y=factor(meta$assigned, levels=c(1,2,0), labels=c("stalled","decline","unassigned")), sdi=meta$sdi, region=meta$region)
m0 <- multinom(y~1,m,trace=FALSE); m1 <- multinom(y~sdi,m,trace=FALSE); m2 <- multinom(y~sdi+region,m,trace=FALSE)
lrm <- as.numeric(2*(logLik(m2)-logLik(m1))); dfm <- attr(logLik(m2),"df")-attr(logLik(m1),"df")
prm <- predict(m1, data.frame(sdi=seq(0.3,0.9,0.1)), type="probs")
cat(sprintf("MULTINOM: R2 SDI %.3f, +region %.3f, LR %.1f on %d df p %.3g; P(unassigned) max %.2f at SDI %.1f\n",
    1-as.numeric(logLik(m1))/as.numeric(logLik(m0)), 1-as.numeric(logLik(m2))/as.numeric(logLik(m0)), lrm, dfm, pchisq(lrm,dfm,lower.tail=FALSE), max(prm[,"unassigned"]), seq(0.3,0.9,0.1)[which.max(prm[,"unassigned"])]))

# ---- tables (CSV and LaTeX) -------------------------------------------------------------------
grp <- factor(meta$assigned, levels=c(1,2,0), labels=c(lab_short,"Unassigned"))
t1 <- data.frame(Group=levels(grp), n=as.vector(table(grp)), `Rate 1990`=tapply(exp(Y[,"1990"]),grp,q0), `Rate 2023`=tapply(exp(Y[,"2023"]),grp,q0),
  `Change (%)`=tapply(pc,grp,q0,1), `Peak year`=tapply(peak,grp,q0), `SDI 2023`=tapply(meta$sdi,grp,q0,2), `Interval width`=tapply(w,grp,q0,3),
  `Membership probability`=tapply(meta$p,grp,q0,2), check.names=FALSE)
t2 <- as.data.frame.matrix(table(regname[meta$region], grp)); t2$Total <- rowSums(t2); t2 <- rbind(t2, Total=colSums(t2)); t2 <- cbind(`Super region`=rownames(t2), t2)
t3 <- data.frame(k=2:6, Silhouette=round(rD$sil,3), Gap=sprintf("%.3f (%.3f)", rD$gap[2:6], rD$gap_se[2:6]), `Weakest cluster median p (value)`=round(rD$weakest,2),
  `Share unassigned (value)`=round(rD$unas,2), `Weakest cluster median p (country bootstrap)`=round(rB$weakest,2), check.names=FALSE)
t4 <- data.frame(Model=c("Intercept only","SDI","SDI + super region","SDI + interval width","SDI + super region + interval width"),
  Deviance=round(c(deviance(g0),deviance(g1),deviance(g2),deviance(g3),deviance(g4)),1), df=c(g0$df.residual,g1$df.residual,g2$df.residual,g3$df.residual,g4$df.residual),
  AIC=round(c(AIC(g0),AIC(g1),AIC(g2),AIC(g3),AIC(g4)),1), `McFadden R2`=round(1-c(deviance(g0),deviance(g1),deviance(g2),deviance(g3),deviance(g4))/deviance(g0),3), check.names=FALSE)
t5 <- data.frame(Analysis=c("DALYs 1990 to 2023 (primary)","Deaths 1990 to 2023","DALYs 2000 to 2023","DALYs, narrowest half of intervals (n = 101)","Country bootstrap instead of value perturbation"),
  `k (silhouette)`=c(rD$k_sil,rM$k_sil,r2$k_sil,rs$k_sil,rB$k_sil), `k (gap)`=c(rD$k_gap,rM$k_gap,r2$k_gap,rs$k_gap,rB$k_gap), `k (stability)`=c(rD$k,rM$k,r2$k,rs$k,rB$k),
  `Weakest median p at k = 3`=round(c(rD$weakest[2],rM$weakest[2],r2$weakest[2],rs$weakest[2],rB$weakest[2]),2),
  Assigned=c(sum(rD$assigned>0),sum(rM$assigned>0),sum(r2$assigned>0),sum(rs$assigned>0),sum(rB$assigned>0)),
  `ARI with primary (jointly assigned)`=round(c(1, {b<-rD$assigned>0&rM$assigned>0; ari(rD$assigned[b],rM$assigned[b])}, {b<-rD$assigned>0&r2$assigned>0; ari(rD$assigned[b],r2$assigned[b])},
     {b<-rs$assigned>0&rD$assigned[keep]>0; ari(rs$assigned[b],rD$assigned[keep][b])}, {b<-rD$assigned>0&rB$assigned>0; ari(rD$assigned[b],rB$assigned[b])}),3), check.names=FALSE)
for (nmt in c("t1","t2","t3","t4","t5")) write.csv(get(nmt), tab(paste0("table_", nmt, ".csv")), row.names=FALSE)
tex_table(t1, tab("tab1.tex"), "lrrrrrrrr"); tex_table(t2, tab("tab2.tex")); tex_table(t3, tab("tab3.tex")); tex_table(t4, tab("tab4.tex")); tex_table(t5, tab("tab5.tex"), "lrrrrrr")
s1 <- data.frame(Country=meta$name, Region=regname[meta$region], Income=incname[meta$income], SDI=meta$sdi, `Change DALYs (%)`=round(pc,1), `Change deaths (%)`=round(pcM,1),
  `Peak year`=peak, `Interval width`=round(w,3), `Cluster (DALYs)`=meta$assigned, `p (DALYs)`=round(meta$p,2), `Cluster (deaths)`=rM$assigned, `p (deaths)`=round(rM$p,2),
  `Cluster (2000 window)`=r2$assigned, `p (2000 window)`=round(r2$p,2), check.names=FALSE); s1 <- s1[order(s1$Region, s1$Country),]
write.csv(s1, tab("tableS1_country_membership.csv"), row.names=FALSE)
write.csv(data.frame(mor), tab("tableS2_moran_neighbours.csv"), row.names=FALSE)

# ---- figures -----------------------------------------------------------------------------------
colr <- ifelse(meta$assigned==0, gU, pal[pmax(meta$assigned,1)]); nA <- table(factor(meta$assigned,levels=c(1,2,0)))
# Fig 1: trajectories + membership probability
png(fig("fig1.png"), width=2200, height=950, res=220); layout(matrix(c(1,2,3,3),2,2,byrow=TRUE), heights=c(1,0.16), widths=c(1.4,1)); par(mar=c(4,4.2,1,1))
plot(NA, xlim=range(yrs), ylim=range(shape), xlab="Year", ylab="Cardiovascular DALY rate, log ratio to country mean", las=1); abline(h=0,col="grey70",lty=3)
for (i in which(meta$assigned==0)) lines(yrs, shape[i,], col=adjustcolor(gU,0.25))
for (c in 1:2) { s <- shape[meta$assigned==c,,drop=FALSE]; for (i in seq_len(nrow(s))) lines(yrs, s[i,], col=adjustcolor(pal[c],0.18)); lines(yrs, colMeans(s), col=pal[c], lwd=3) }
boxplot(p~cluster, meta, col=adjustcolor(pal[1:2],0.6), xlab="Cluster", ylab="Membership probability", las=1, pch=16, cex=0.6); abline(h=0.7,lty=2)
par(mar=c(0,0,0,0)); plot.new(); legend("center", horiz=TRUE, bty="n", lwd=3, col=c(pal[1:2],gU), legend=c(paste0(lab_short, " (n = ", nA[1:2], ")"), paste0("Unassigned (n = ", nA[3], ")")))
dev.off()
# Fig 2: map
wn <- map("world", plot=FALSE, ylim=c(-58,84))$names; base <- sub(":.*","",wn)
mapname <- setNames(as.list(meta$name), meta$name); for (g in names(xw)) mapname[[g]] <- xw[[g]]
key <- do.call(rbind, lapply(names(mapname), function(g) if (!is.na(mapname[[g]][1])) data.frame(gbd=g, mapn=mapname[[g]])))
colv <- rep(gN, length(wn)); hit <- match(base, key$mapn); okh <- !is.na(hit); cbg <- setNames(colr, meta$name); colv[okh] <- cbg[key$gbd[hit[okh]]]
small <- meta$name[meta$name %in% c("Singapore","Bahrain","Malta","Maldives","Seychelles","Comoros","Mauritius","Barbados","Grenada","Dominica","Saint Lucia","Antigua and Barbuda","Saint Kitts and Nevis","Saint Vincent and the Grenadines","Nauru","Kiribati","Marshall Islands","Palau","Tonga","Samoa","American Samoa","Micronesia (Federated States of)","Northern Mariana Islands","Guam","Sao Tome and Principe","Cabo Verde","Andorra","Monaco","San Marino","Bermuda")]
png(fig("fig2.png"), width=2400, height=1350, res=230); layout(matrix(1:2,2,1), heights=c(1,0.2))
map("world", fill=TRUE, col=colv, border="white", lwd=0.2, mar=c(0,0,0,0), ylim=c(-58,84))
for (s in small) { rg <- tryCatch(map("world", mapname[[s]][1], plot=FALSE)$range, error=function(e) NULL); if (!is.null(rg)) points(mean(rg[1:2]), mean(rg[3:4]), pch=21, bg=cbg[s], col="black", cex=0.7, lwd=0.3) }
par(mar=c(0,0,0,0)); plot.new()
legend("center", ncol=2, bty="n", pch=22, pt.cex=2.2, pt.bg=c(pal[1:2],gU,gN), cex=0.95, legend=c(
  sprintf("Cluster 1: stalled (median change %.0f%%), n = %d", median(pc[meta$assigned==1]), nA[1]), sprintf("Cluster 2: sustained decline (median change %.0f%%), n = %d", median(pc[meta$assigned==2]), nA[2]),
  sprintf("Unassigned (membership probability < 0.70), n = %d", nA[3]), "Not in analysis set"))
dev.off()
# Fig 3: k selection, now four panels including country bootstrap
png(fig("fig3.png"), width=2400, height=700, res=220); par(mfrow=c(1,4), mar=c(4,4.4,1,1))
plot(2:6, rD$sil, type="b", pch=16, xlab="Number of clusters k", ylab="Mean silhouette width", las=1); abline(v=rD$k_sil,col="grey60",lty=3)
plot(1:6, rD$gap, type="b", pch=16, xlab="Number of clusters k", ylab="Gap statistic", las=1, ylim=range(c(rD$gap-rD$gap_se, rD$gap+rD$gap_se)))
arrows(1:6, rD$gap-rD$gap_se, 1:6, rD$gap+rD$gap_se, angle=90, code=3, length=0.03); abline(v=rD$k_gap,col="grey60",lty=3)
plot(2:6, rD$weakest, type="b", pch=16, xlab="Number of clusters k", ylab="Weakest cluster median membership probability", las=1, ylim=c(0,1)); abline(h=0.8,lty=2); abline(v=rD$k,col="grey60",lty=3)
points(2:6, 1-rD$unas, type="b", pch=1, lty=3); legend("bottomleft", bty="n", pch=c(16,1), lty=c(1,3), legend=c("Weakest cluster median","Share assigned"), cex=0.85)
plot(2:6, rB$weakest, type="b", pch=16, xlab="Number of clusters k", ylab="Weakest cluster median (country bootstrap)", las=1, ylim=c(0,1)); abline(h=0.8,lty=2); abline(v=rB$k,col="grey60",lty=3)
dev.off()
# Fig 4: level vs shape, SDI
png(fig("fig4.png"), width=2400, height=900, res=230); layout(matrix(c(1,2,3,3),2,2,byrow=TRUE), heights=c(1,0.14)); par(mar=c(4,4.4,1,1))
plot(exp(Y[,"1990"]), pc, col=adjustcolor(colr,0.8), pch=16, log="x", las=1, xlab="Age standardised CVD DALY rate, 1990 (per 100,000, log scale)", ylab="Change 1990 to 2023 (%)"); abline(h=0,col="grey70",lty=3)
plot(meta$sdi, pc, col=adjustcolor(colr,0.8), pch=16, las=1, xlab="Socio-demographic Index, 2023", ylab="Change 1990 to 2023 (%)"); abline(h=0,col="grey70",lty=3)
par(mar=c(0,0,0,0)); plot.new(); legend("center", horiz=TRUE, bty="n", pch=16, col=c(pal[1:2],gU), legend=c(lab_short,"Unassigned")); dev.off()
# Fig 5: composition
regshort <- c(CEECA="C. and E. Europe,\nCentral Asia", HI="High income", LAC="Latin America,\nCaribbean", NAME="N. Africa,\nMiddle East", SA="South Asia", SEAO="SE Asia, E. Asia,\nOceania", SSA="Sub-Saharan\nAfrica")
png(fig("fig5.png"), width=2400, height=1000, res=230); layout(matrix(c(1,2,3,3),2,2,byrow=TRUE), heights=c(1,0.12), widths=c(1.6,1)); par(mar=c(5,4.4,1,1))
tr <- table(factor(meta$assigned,levels=c(1,2,0)), meta$region); tr <- tr[, order(-tr[1,]/colSums(tr))]
bp <- barplot(prop.table(tr,2)*100, col=c(pal[1:2],gU), border=NA, las=1, names.arg=rep("",ncol(tr)), ylab="Countries (%)"); text(bp,-6,regshort[colnames(tr)],xpd=TRUE,cex=0.62,adj=c(0.5,1))
ti <- table(factor(meta$assigned,levels=c(1,2,0)), factor(meta$income,levels=c("low","lmi","umi","hi")))
bp2 <- barplot(prop.table(ti,2)*100, col=c(pal[1:2],gU), border=NA, las=1, names.arg=rep("",ncol(ti)), ylab="Countries (%)"); text(bp2,-6,paste0(incname[colnames(ti)],"\nincome"),xpd=TRUE,cex=0.7,adj=c(0.5,1))
par(mar=c(0,0,0,0)); plot.new(); legend("center", horiz=TRUE, bty="n", pch=22, pt.cex=2, pt.bg=c(pal[1:2],gU), legend=c(lab_short,"Unassigned")); dev.off()
# Fig 6: robustness. (a) deaths vs DALYs change coloured by DALY cluster; (b) 2000 window trajectories; (c) membership p value vs bootstrap
png(fig("fig6.png"), width=2400, height=800, res=220); layout(matrix(c(1,2,3,4,4,4),2,3,byrow=TRUE), heights=c(1,0.14)); par(mar=c(4,4.4,1,1))
plot(pc, pcM, col=adjustcolor(colr,0.8), pch=16, las=1, xlab="Change in DALY rate 1990 to 2023 (%)", ylab="Change in death rate 1990 to 2023 (%)"); abline(0,1,col="grey70",lty=3)
sh2 <- r2$shape; col2 <- ifelse(r2$assigned==0, gU, pal[pmax(r2$assigned,1)])
plot(NA, xlim=c(2000,2023), ylim=range(sh2), xlab="Year", ylab="CVD DALY rate, log ratio to 2000 to 2023 mean", las=1); abline(h=0,col="grey70",lty=3)
for (i in which(r2$assigned==0)) lines(2000:2023, sh2[i,], col=adjustcolor(gU,0.25))
for (c in 1:2) { s <- sh2[r2$assigned==c,,drop=FALSE]; for (i in seq_len(nrow(s))) lines(2000:2023, s[i,], col=adjustcolor(pal[c],0.18)); lines(2000:2023, colMeans(s), col=pal[c], lwd=3) }
plot(rD$stab[["3"]]$p, rB$stab[["3"]]$p, col=adjustcolor(colr,0.8), pch=16, las=1, xlab="Membership probability at k = 3, value perturbation", ylab="Membership probability at k = 3, country bootstrap"); abline(0,1,col="grey70",lty=3); abline(h=0.7,v=0.7,lty=2,col="grey60")
par(mar=c(0,0,0,0)); plot.new(); legend("center", horiz=TRUE, bty="n", pch=16, col=c(pal[1:2],gU), legend=c(lab_short,"Unassigned")); dev.off()
# Fig 7: interval width. (a) width by region; (b) width vs membership probability; (c) narrow half trajectories
png(fig("fig7.png"), width=2400, height=800, res=220); layout(matrix(c(1,2,3,4,4,4),2,3,byrow=TRUE), heights=c(1,0.14)); par(mar=c(4.5,4.4,1,1))
ow <- names(sort(tapply(w, meta$region, median))); boxplot(w~factor(meta$region,levels=ow), col="grey85", las=2, names=rep("",7), xlab="", ylab="Mean log scale SD of published interval", pch=16, cex=0.6)
text(1:7, par("usr")[3]-0.012, regshort[ow], xpd=TRUE, cex=0.55, adj=c(0.5,1))
plot(w, meta$p, col=adjustcolor(colr,0.8), pch=16, las=1, xlab="Mean log scale SD of published interval", ylab="Membership probability (primary)"); abline(h=0.7,lty=2,col="grey60")
shs <- rs$shape; cols <- ifelse(rs$assigned==0, gU, pal[pmax(rs$assigned,1)])
plot(NA, xlim=range(yrs), ylim=range(shs), xlab="Year", ylab="Narrowest half: log ratio to country mean", las=1); abline(h=0,col="grey70",lty=3)
for (i in which(rs$assigned==0)) lines(yrs, shs[i,], col=adjustcolor(gU,0.25))
for (c in sort(unique(rs$assigned[rs$assigned>0]))) { s <- shs[rs$assigned==c,,drop=FALSE]; for (i in seq_len(nrow(s))) lines(yrs, s[i,], col=adjustcolor(pal[c],0.18)); lines(yrs, colMeans(s), col=pal[c], lwd=3) }
par(mar=c(0,0,0,0)); plot.new(); legend("center", horiz=TRUE, bty="n", pch=16, col=c(pal[1:3],gU), legend=c("Cluster 1 (least decline)","Cluster 2","Cluster 3 (narrow half only)","Unassigned")); dev.off()
# Fig 8: multinomial predicted probabilities across SDI
png(fig("fig8.png"), width=1600, height=900, res=220); par(mar=c(4,4.4,1,1)); sg <- seq(0.2,0.95,0.01); pr <- predict(m1, data.frame(sdi=sg), type="probs")
plot(NA, xlim=range(sg), ylim=c(0,1), xlab="Socio-demographic Index, 2023", ylab="Predicted probability", las=1)
lines(sg, pr[,"stalled"], col=pal[1], lwd=3); lines(sg, pr[,"decline"], col=pal[2], lwd=3); lines(sg, pr[,"unassigned"], col=gU, lwd=3)
rug(meta$sdi[meta$assigned==0], col=gU); legend("top", horiz=TRUE, bty="n", lwd=3, col=c(pal[1:2],gU), legend=c("Stalled","Sustained decline","Unassigned")); dev.off()

write.csv(meta, "output/final_membership.csv", row.names=FALSE); save.image("output/final.RData"); cat("DONE\n")
