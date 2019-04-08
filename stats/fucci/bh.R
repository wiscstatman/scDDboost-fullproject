
# MAN, April 8, 2019

# check t-stat BH method


load("results.RData")

tdat <- log2( data_count + 1/2 )

X <- tdat[,cd==1]
Y <- tdat[,cd==2]

pv.t <- rep(1,nrow(X))
pv.w <- rep(1,nrow(X))

for( i in 1:nrow(X) )
 {
  if( var( c(X[i,], Y[i,] ) ) > 0 )
   {
     pv.t[i]  <- ( t.test( X[i,], Y[i,] ) )$p.value
     pv.w[i]  <- ( wilcox.test( X[i,], Y[i,] ) )$p.value
   }
 }

qv.t <- p.adjust( pv.t, method="BH" )  ## looks like a lot 
qv.w <- p.adjust( pv.w, method="BH" )  ## l

sum( qv.t <= 0.05 )  ## I find 3361...not sure how this compares to other methods
sum( qv.w <= 0.05 )  ## I find 3496 ...not sure how this compares to other methods

