###polygon plot
library(RColorBrewer)
library(MCMCpack)
### {1,2} {3,4,5} {6,7}
set.seed(33)
K = 7
##get proportion p and q such that p1 + p2 = q1 + q2, p3 + p4 + p5 = q3 + q4 + q5 and p6 + p7 = q6 + q7
whole_prp = rdirichlet(1,c(2,3,2))

p = c(rdirichlet(1,c(1,1))*whole_prp[1],
rdirichlet(1,c(1,1,1))*whole_prp[2],
rdirichlet(1,c(1,1))*whole_prp[3])

q = c(rdirichlet(1,c(1,1))*whole_prp[1],
rdirichlet(1,c(1,1,1))*whole_prp[2],
rdirichlet(1,c(1,1))*whole_prp[3])

p12 = q12 = sum(p[1:2])
p345 = q345 = sum(p[3:5])
p67 = q67 = sum(p[6:7])

pdf("prop.pdf")
scols <- brewer.pal(7,"Set2")
x<- 1:7
types <- sample(x)
plot( x, rep(1,7), axes=FALSE, xlab="",ylab="", type="n", xlim=c(0,8),
ylim=c(0,5)  )
##position for polygon
u <- 1.8; v <- .4
u2 <- 1.1; v2 <- .4

###parameter to scale up the length
scale_ = 8

###first two
polygon( c(0, p12 *scale_,  p12 * scale_,0) , c( rep(u,2), rep(v,2) ), col="lightgrey",
border=FALSE )
text( p12*scale_ / 2,  (3*u+v)/4, expression(Psi[12]) )
polygon( c(0, p[1] * scale_, p[1] * scale_,0) , c( rep(u2,2), rep(v2,2) ), col=scols[1],
border=FALSE )
text( p[1] * scale_ / 2, (u+3*v)/4, expression(psi[1]))
polygon( c( p[1] * scale_, p12 *scale_, p12 *scale_, p[1] * scale_),
c( rep(u2,2), rep(v2,2) ), col=scols[2],
border=FALSE )
text( (p[1]  + p[2] / 2) * scale_ , (u+3*v)/4, expression(psi[2]) )

###second three
polygon( c( p12 *scale_,  (p12 + p345) * scale_ , (p12 + p345) * scale_,  p12 *scale_) ,
c( rep(u,2), rep(v,2) ), col="green",
border=FALSE )
text( (p12 + p345/2) * scale_, (3*u+v)/4, expression(Psi[345]) )
polygon( c( p12 *scale_, (p12 + p[3]) *scale_ ,(p12 + p[3]) *scale_, p12 *scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[3],
border=FALSE )
text( (p12 + p[3]/2) * scale_, (u+3*v)/4, expression(psi[3]) )
polygon( c( (p12 + p[3]) *scale_,(p12 + p[3] + p[4]) *scale_,(p12 + p[3] + p[4]) *scale_,(p12 + p[3]) *scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[4],
border=FALSE )
text( (p12 + p[3] + p[4] / 2) * scale_, (u+3*v)/4, expression(psi[4]) )
polygon( c((p12 + p[3] + p[4]) *scale_,(p12 + p345) *scale_,(p12 + p345) *scale_,(p12 + p[3] + p[4]) *scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[5],
border=FALSE )
text( (p12 + p[3] + p[4] + p[5] / 2) * scale_, (u+3*v)/4, expression(psi[5]) )


###last two
polygon( c((p12 + p345) * scale_,  scale_,scale_, (p12 + p345) * scale_) , c( rep(u,2), rep(v,2) ), col="magenta",
border=FALSE )
text( (p12 + p345 + p67/2) * scale_, (3*u+v)/4, expression(Psi[67]) )
polygon( c((p12 + p345) * scale_,(p12 + p345 + p[6]) * scale_,(p12 + p345 + p[6]) * scale_,(p12 + p345 ) * scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[6],
border=FALSE )
text( (p12 + p345 + p[6]/2) * scale_, (u+3*v)/4, expression(psi[6]) )
polygon( c((p12 + p345 + p[6]) * scale_,scale_,scale_,(p12 + p345 + p[6]) * scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[7],
border=FALSE )
text( (p12 + p345 + p[6] + p[7]/2) * scale_, (u+3*v)/4, expression(psi[7]) )






u <- u + 2
v <- v + 2
u2 <- u2 + 2; v2 <- v2+2

###
###first two
polygon( c(0, q12 *scale_,  q12 * scale_,0) , c( rep(u,2), rep(v,2) ), col="lightgrey",
border=FALSE )
text( q12*scale_ / 2,  (3*u+v)/4, expression(Psi[12]) )
polygon( c(0, q[1] * scale_, q[1] * scale_,0) , c( rep(u2,2), rep(v2,2) ), col=scols[1],
border=FALSE )
text( q[1] * scale_ / 2, (u+3*v)/4, expression(psi[1]))
polygon( c( q[1] * scale_, q12 *scale_, q12 *scale_, q[1] * scale_),
c( rep(u2,2), rep(v2,2) ), col=scols[2],
border=FALSE )
text( (q[1]  + q[2] / 2) * scale_ , (u+3*v)/4, expression(psi[2]) )

###second three
polygon( c( q12 *scale_,  (q12 + q345) * scale_ , (q12 + q345) * scale_,  q12 *scale_) ,
c( rep(u,2), rep(v,2) ), col="green",
border=FALSE )
text( (q12 + q345/2) * scale_, (3*u+v)/4, expression(Psi[345]) )
polygon( c( q12 *scale_, (q12 + q[3]) *scale_ ,(q12 + q[3]) *scale_, q12 *scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[3],
border=FALSE )
text( (q12 + q[3]/2) * scale_, (u+3*v)/4, expression(psi[3]) )
polygon( c( (q12 + q[3]) *scale_,(q12 + q[3] + q[4]) *scale_,(q12 + q[3] + q[4]) *scale_,(q12 + q[3]) *scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[4],
border=FALSE )
text( (q12 + q[3] + q[4] / 2) * scale_, (u+3*v)/4, expression(psi[4]) )
polygon( c((q12 + q[3] + q[4]) *scale_,(q12 + q345) *scale_,(q12 + q345) *scale_,(q12 + q[3] + q[4]) *scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[5],
border=FALSE )
text( (q12 + q[3] + q[4] + q[5] / 2) * scale_, (u+3*v)/4, expression(psi[5]) )


###last two
polygon( c((q12 + q345) * scale_,  scale_,scale_, (q12 + q345) * scale_) , c( rep(u,2), rep(v,2) ), col="magenta",
border=FALSE )
text( (q12 + q345 + q67/2) * scale_, (3*u+v)/4, expression(Psi[67]) )
polygon( c((q12 + q345) * scale_,(q12 + q345 + q[6]) * scale_,(q12 + q345 + q[6]) * scale_,(q12 + q345 ) * scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[6],
border=FALSE )
text( (q12 + q345 + q[6]/2) * scale_, (u+3*v)/4, expression(psi[6]) )
polygon( c((q12 + q345 + q[6]) * scale_,scale_,scale_,(q12 + q345 + q[6]) * scale_) , c( rep(u2,2), rep(v2,2) ), col=scols[7],
border=FALSE )
text( (q12 + q345 + q[6] + q[7]/2) * scale_, (u+3*v)/4, expression(psi[7]) )

dev.off()



