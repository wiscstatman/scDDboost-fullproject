###plot that scDDboost may lose FDR contorl when over estimating K

load("caseForlossControl.RData")

intra = rep(0, 6)
inter = rep(0, 6)
for (i in 2:7) {
    clusRes = pam(D_c, i, diss = T)
    intra[i - 1] = as.numeric(clusRes$objective[1])
    x = clusRes$id.med
    inter[i - 1] = mean(D_c[x, x])
}
s = intra/inter

s.fdr = rep(0,6)
ngenes = nrow(data_null)
for(i in 1:6){
    s.fdr[i] = listsize(pdd[[i + 1]], 0.05) / ngenes
}
pdf("breakFDR.pdf")
plot(s, type = "l", ylim = c(0,1.6), lwd = 2, ylab = "", xlab = "K",main = "EMTAB2805",
xaxt = "n",yaxt = "n")
lines(s.fdr,lwd = 2, lty = 3)
abline(h = 0.05, col = "red", lwd = 2)
axis(1, at=1:6, labels=2:7, cex.axis= 1.2, las = 1)
axis(2, at=c(0,0.05,0.5,1,1.5), labels = c(0,0.05,0.5,1,1.5), cex.axis = 1.2, las = 1)
legend("topright", legend=c("Validity Score", "False Discovery Rate"),
lty = c(1,3), cex = 1.2,lwd = 2)
dev.off()

