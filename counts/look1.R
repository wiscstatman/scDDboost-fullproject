
ss <- read.csv("DD_0.99.csv", row.names=1)

tt <- ss[-6,1:3]  ## dropping one extreme DD case
###oo <- order( rowSums(tt) )
oo <- c(6,4,5, 3,7, 8,2,1 )
qq <- tt[oo,]

#matplot( qq, type="b", axes=FALSE, xlab="",ylab="" )
#axis( side=1, at=1:8, labels=rownames(qq), las=-1 )


#x <- c( rep(1,8), rep(2,8), rep(3,8) )
#x <- c( qq[,1], qq[,2], qq[,3] )
#y <- rep( (8:1), 3 )
#plot( x, y)

