NM = c("GSE45719","GSE63818","GSE48968","GSE60749",
       "GSE74596", "GSE71585", "EMTAB2805","GSE79102",
       "GSE64016","GSE75748", "GSE84465","GSE94383")
total_g = c(45686,65218,27723,22443,23337,24057,38390,19974,19084,19097,23368,21383)
nDD_mst = c(9437,4565,3449,8044,1783,2179,5609,1506,2714,7170,12725,2932)
nDD_des = c(8620,4667,1541,7878,2561,4173,2062,541,1191,5951,11217,1652)
nDD_sc = c(7259,2199,1724,9184,1154,1570,2468,1010,1091,7106,10598,873)
nDD_scb = c(9392,5023,434,11954,7740,3051,3324,3953,3631,5057,13221,1824)
#total_g = rep(1,13)
TOT = (nDD_mst + nDD_des + nDD_sc + nDD_scb) / total_g
ORD = order(TOT)
pdf("DD95.pdf", height = 6, width = 10)
par(mar=c(7,5,4,1)+.1)
plot((nDD_scb[ORD] / total_g[ORD]), type = "b", lwd = 2, col = "green", 
     ylab = "", xaxt = 'n', xlab = "")
mtext("Proportion of DD genes", side=2, line=2.2, cex=1.2)
lines((nDD_des[ORD]/ total_g[ORD]) , type = "b", lwd = 2, col = "red")
lines( (nDD_sc[ORD]/total_g[ORD]) , type = "b", lwd = 2, col = "blue")
lines( (nDD_mst[ORD] / total_g[ORD]) , type = "b", lwd = 2)
axis(1, at=1:12, labels=NM[ORD], cex.axis= 1.2, las = 2)
legend("topleft", legend=c("MAST", "DESeq2", "scDD", "scDDboost"),
       col=c("black", "red", "blue","green"),lty = 1, cex = 1.2,lwd = 2)

dev.off()
