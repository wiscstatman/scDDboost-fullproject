###CODE for FIG6, number of DD genes identified by different methods
NM = c("GSE45719","GSE63818","GSE48968","GSE94383",
"GSE74596", "GSE71585", "EMTAB2805","GSE79102",
"GSE64016","GSE75748", "GSE84465","GSE52529")
total_g = c(45686,65218,27723,21383,23337,24057,45686,19974,19084,19097,23368,47192)
nDD_mst = c(6061, 668,  1502, 1104, 546,  2259, 1112, 311,  486,  5560, 13054,986)
nDD_des = c(6590, 4667, 1541, 1652, 2561, 2002, 1365, 541,  1091, 6026, 8116,1365)
nDD_sc =  c(7013, 2199, 1496, 873,  2176,  2698, 1484, 1402, 1191, 5323, 11781,2438)
nDD_scb = c(9898, 5054, 372,  927,  4882, 6776, 5447, 3107, 3395, 6094, 16644,2725)
#total_g = rep(1,13)
TOT = (nDD_mst + nDD_des + nDD_sc + nDD_scb) / total_g
ORD = order(TOT)
ORD = order(nDD_scb / total_g)
pdf("DD95.pdf", height = 6, width = 10)
par(mar=c(7,5,4,1)+.1)
plot((nDD_scb[ORD] / total_g[ORD])[1:12], type = "b", lwd = 2, col = "green",
ylab = "", xaxt = 'n', xlab = "")
mtext("Proportion of DD genes", side=2, line=2.2, cex=1.2)
lines((nDD_des[ORD]/ total_g[ORD])[1:12] , type = "b", lwd = 2, col = "red")
lines( (nDD_sc[ORD]/total_g[ORD])[1:12] , type = "b", lwd = 2, col = "blue")
lines( (nDD_mst[ORD] / total_g[ORD])[1:12] , type = "b", lwd = 2)
axis(1, at=1:12, labels=NM[ORD][1:12], cex.axis= 1.2, las = 2)
legend("topleft", legend=c("MAST", "DESeq2", "scDD", "scDDboost"),
col=c("black", "red", "blue","green"),lty = 1, cex = 1.2,lwd = 2)

dev.off()
