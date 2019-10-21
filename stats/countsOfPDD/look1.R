
#ss <- read.csv("DD_0.99.csv", row.names=1)
ss <- read.csv("DD_0.95.csv", row.names=1)

tt <- ss[-6,1:3]  ## dropping one extreme DD case

oo <- order( -tt[,1] )
oo <- order( -rowSums(tt) )
#oo <- order( -apply(tt,1,var) )
###oo <- c(6,4,5, 3,7, 8,2,1 )
qq <- tt[oo,]

pdf( file="conquer95.pdf", height=5, width=8 )
cls <- c("blue", "green", "darkorange" )
matplot( qq, type="b", axes=FALSE, xlab="",ylab="", col=cls, lwd=4,
	 lty=c(5,2,1), cex=2 )
axis( side=1, at=1:8, labels=rownames(qq), cex.axis=.85 )
axis( side=2, las=1, cex.axis=.9 )
box(bty="l")
dev.off()
