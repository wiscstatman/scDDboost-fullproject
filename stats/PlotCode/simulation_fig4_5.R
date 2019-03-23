library(ggplot2)

#True positive rate, False positive rate, FDR order from K = 3, 7, 15 and (0.1,0.4), (-0.1, 0.3), (0.3, 0.5) , (-0.1,1)

# MAST
TP_mst = c(0.0928,0.0559,0.165,0.31,0.0008,0,0.006,0.041,0,0,0.00029,0.006)
FP_mst = c(0.00215,0.002,0.005,0.007,0,0,0,0.001,0,0,0,0.00014)
FDR_mst = c(0.0226,0.037,0.03,0.023,0,0,0,0.023,0,0,0,0.023)

# scDD
TP_sc = c(0.077,0.0434,0.148,0.283,0,0,0.003,0.028,0,0,0,0.00029)
FP_sc = c(0.0019,0.0017,0.0036,0.0045,0,0,0,0.0004,0,0,0,0)
FDR_sc = c(0.024,0.038,0.024,0.015,0,0,0,0.015,0,0,0,0)

# DESeq2
TP_des = c(0.123,0.083,0.191,0.33,0.0096,0.0027,0.0232,0.085,0.00094,0.057,0.00326,0.057)
FP_des = c(0.0049,0.0025,0.008,0.013,0.0003,0,0.0005,0.0025,0,0.0008,0.00007,0.0008)
FDR_des = c(0.038,0.03,0.04,0.039,0.03,0,0.023,0.028,0,0.013,0.021,0.013)

# scDDboost
TP_scb = c(0.0966,0.0288,0.34,0.52,0.0316,0.0225,0.039,0.23,0,0.006,0.00065,0.027)
FP_scb = c(0,0,0.0053,0.022,0.0003,0,0,0.002,0,0,0,0)
FDR_scb = c(0,0,0.0155,0.021,0.01,0,0,0.009,0,0,0,0)


pdf("simuTPR.pdf", height = 6, width = 10)
par(mar=c(7,5,4,1)+.1)
plot(TP_scb, type = "b", lwd = 2, col = "green", 
     ylab = "", xaxt = 'n', xlab = "K", ylim = c(-0.025,0.5))
mtext("TPR", side=2, line=2.2, cex=1.2, yaxt = 'n')
lines(TP_des , type = "b", lwd = 2, col = "red")
lines( TP_sc , type = "b", lwd = 2, col = "blue")
lines( TP_mst , type = "b", lwd = 2)
#axis(1, at=1:12, cex.axis= 1.2, las = 2) 
legend("topleft", legend=c("MAST", "DESeq2", "scDD", "scDDboost"),
       col=c("black", "red", "blue","green"),lty = 1, cex = 1.2,lwd = 2)
polygon( c(1,1, 4.5 ,4.5) , c(-0.025,-0.005,-0.005,-0.025) , col="lightgrey",
         border=FALSE )
text( 3,  -0.0125, "4" )

polygon( c(4.5,4.5, 8.5 ,8.5) , c(-0.025,-0.005,-0.005,-0.025) , col="magenta",
         border=FALSE )
text( 6.5,  -0.0125, "7" )

polygon( c(8.5,8.5, 12 ,12) , c(-0.025,-0.005,-0.005,-0.025) , col="yellow",
         border=FALSE )
text( 10.5,  -0.0125, "15" )
dev.off()


pdf("simuFPR.pdf", height = 6, width = 10)
par(mar=c(7,5,4,1)+.1)
plot(FP_scb, type = "b", lwd = 2, col = "green", 
     ylab = "", xaxt = 'n', xlab = "K", ylim = c(-0.025,0.5))
mtext("FPR", side=2, line=2.2, cex=1.2, yaxt = 'n')
lines(FP_des , type = "b", lwd = 2, col = "red")
lines( FP_sc , type = "b", lwd = 2, col = "blue")
lines( FP_mst , type = "b", lwd = 2)
#axis(1, at=1:12, cex.axis= 1.2, las = 2) 
legend("topleft", legend=c("MAST", "DESeq2", "scDD", "scDDboost"),
       col=c("black", "red", "blue","green"),lty = 1, cex = 1.2,lwd = 2)
polygon( c(1,1, 4.5 ,4.5) , c(-0.025,-0.005,-0.005,-0.025) , col="lightgrey",
         border=FALSE )
text( 3,  -0.0125, "4" )

polygon( c(4.5,4.5, 8.5 ,8.5) , c(-0.025,-0.005,-0.005,-0.025) , col="magenta",
         border=FALSE )
text( 6.5,  -0.0125, "7" )

polygon( c(8.5,8.5, 12 ,12) , c(-0.025,-0.005,-0.005,-0.025) , col="yellow",
         border=FALSE )
text( 10.5,  -0.0125, "15" )
dev.off()

pdf("simuFDR.pdf", height = 6, width = 10)
par(mar=c(7,5,4,1)+.1)
plot(FDR_scb, type = "b", lwd = 2, col = "green",
ylab = "", xaxt = 'n', xlab = "K", ylim = c(-0.025,0.8))
mtext("TPR", side=2, line=2.2, cex=1.2, yaxt = 'n')
lines(FDR_des , type = "b", lwd = 2, col = "red")
lines( FDR_sc , type = "b", lwd = 2, col = "blue")
lines( FDR_mst , type = "b", lwd = 2)
#axis(1, at=1:12, cex.axis= 1.2, las = 2)
legend("topleft", legend=c("MAST", "DESeq2", "scDD", "scDDboost"),
col=c("black", "red", "blue","green"),lty = 1, cex = 1.2,lwd = 2)
polygon( c(1,1, 4.5 ,4.5) , c(-0.025,-0.005,-0.005,-0.025) , col="lightgrey",
border=FALSE )
text( 3,  -0.0125, "4" )

polygon( c(4.5,4.5, 8.5 ,8.5) , c(-0.025,-0.005,-0.005,-0.025) , col="magenta",
border=FALSE )
text( 6.5,  -0.0125, "7" )

polygon( c(8.5,8.5, 12 ,12) , c(-0.025,-0.005,-0.005,-0.025) , col="yellow",
border=FALSE )
text( 10.5,  -0.0125, "15" )
dev.off()
