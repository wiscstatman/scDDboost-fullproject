library(ggplot2)

#True positive rate order from K = 4, 7, 15 and (0.1,0.4), (-0.1, 0.3), (0.3, 0.5) , (-0.1,1)
TP_mst = c(0.106,0.0856,0.145,0.268,0.071,0.06,0.0946,0.157,0.051,0.048,0.06,0.086)
TP_sc = c(0.042,0.028,0.066,0.0286,0.041,0.031,0.05,0.097,0.018,0.018,0.024,0.038)
TP_des = c(0.02,0.01,0.049,0.172,0.0096,0.0027,0.023,0.085,0.0009,0.0007,0.003,0.035)
TP_scb = c(0.102,0.138,0.276,0.361,0.031,0.32,0.047,0.21,0.158,0.078,0.39,0.205)

#False positive rate
FP_mst = c(0.085,0.0806,0.083,0.0826,0.039,0.039,0.037,0.039,0.01,0.01,0.01,0.009)
FP_sc = c(0.028,0.026,0.024,0,0.025,0.023,0.022,0.02,0.0049,0.004,0.004,0.004)
FP_des = c(0.0008,0.0005,0.0016,0.004,0.0003,0,0.0005,0.002,0,0,0,0.0005)
FP_scb = c(0.0006,0.04,0.004,0.022,0.0003,0.31,0.0002,0.001,0.047,0.022,0.1,0.0002)

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
