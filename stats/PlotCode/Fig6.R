###CODE for FIG6, number of DD genes identified by different methods

NM = c("GSE45719","GSE63818","GSE48968","GSE60749",
"GSE74596", "GSE71585", "EMTAB2805","GSE79102",
"GSE64016","GSE75748", "GSE84465")
total_g = c(45686,65218,27723,22443,23337,24057,38390,19974,19084,19097,23368)
nDD_mst = c(6061, 668,  1502, 5948, 282,  2259, 2202, 311,  486,  5497, 12725)
nDD_des = c(6590, 4667, 1541, 7878, 2561, 2002, 2062, 541,  1091, 5951, 11217)
nDD_sc =  c(7013, 2199, 1496, 10493,821,  2698, 2468, 1402, 1191, 7106, 10598)
nDD_scb = c(9898, 5054, 372,  11954,5418, 7062, 767,  3107, 3441, 6415, 13221)
#total_g = rep(1,13)
TOT = (nDD_mst + nDD_des + nDD_sc + nDD_scb) / total_g
ORD = order(TOT)
pdf("DD95.pdf", height = 6, width = 10)
par(mar=c(7,5,4,1)+.1)
plot((nDD_scb[ORD] / total_g[ORD])[1:11], type = "b", lwd = 2, col = "green",
ylab = "", xaxt = 'n', xlab = "")
mtext("Proportion of DD genes", side=2, line=2.2, cex=1.2)
lines((nDD_des[ORD]/ total_g[ORD])[1:11] , type = "b", lwd = 2, col = "red")
lines( (nDD_sc[ORD]/total_g[ORD])[1:11] , type = "b", lwd = 2, col = "blue")
lines( (nDD_mst[ORD] / total_g[ORD])[1:11] , type = "b", lwd = 2)
axis(1, at=1:11, labels=NM[ORD][1:11], cex.axis= 1.2, las = 2)
legend("topleft", legend=c("MAST", "DESeq2", "scDD", "scDDboost"),
col=c("black", "red", "blue","green"),lty = 1, cex = 1.2,lwd = 2)

dev.off()
