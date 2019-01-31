###CODE for FIG6, number of DD genes identified by different methods

NM = c("GSE45719","GSE63818","GSE48968","GSE52529","GSE60749",
       "GSE74596", "GSE71585", "EMTAB2805","GSE79102","GSE57872",
       "GSE64016","GSE75748", "GSE84465")
total_g = c(45686,65218,27723,47192,22443,23337,24057,38390,19974,65218,19084,19097,23368)
nDD_mst = c(9437,4565,3449,4135,8044,1783,2179,5609,1506,31925,2714,7170,12725)
nDD_des = c(8620,4667,1541,1828,7878,2561,4173,2062,541,8976,1191,5951,11217)
nDD_sc = c(7259,2199,1724,1410,9184,1154,1570,2468,1010,11866,1091,7106,10598)
nDD_scb = c(15481,9306,4663,4860,10080,4014,4994,5433,4437,33281,663,5953,13221)
total_g = rep(1,13)
TOT = (nDD_mst + nDD_des + nDD_sc + nDD_scb) / total_g
ORD = order(TOT)
pdf("DD_0.95.pdf", height = 5, width = 10)

plot((nDD_scb[ORD] / total_g[ORD])[1:11], type = "b", lwd = 2, col = "green", 
     ylab = "proportion of DD genes" ,xaxt = "n")

lines((nDD_des[ORD]/ total_g[ORD])[1:11] , type = "b", lwd = 2, col = "red")
lines( (nDD_sc[ORD]/total_g[ORD])[1:11] , type = "b", lwd = 2, col = "blue")
lines( (nDD_mst[ORD] / total_g[ORD])[1:11] , type = "b", lwd = 2)
axis(1, at=1:11, labels=NM[ORD][1:11],cex.axis=0.5) 
legend(1, 15000, legend=c("MAST", "Deseq2", "scDD", "scDDb"),
       col=c("black", "red", "blue","green"),lty = 1, cex=0.8)

dev.off()