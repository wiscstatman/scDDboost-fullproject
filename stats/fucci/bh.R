
# MAN, April 8, 2019

# check t-stat BH method


load("results.RData")

tdat <- log2( data_count + 1/2 )

X <- tdat[,cd==1]
Y <- tdat[,cd==2]

pv <- rep(1,nrow(X))

for( i in 1:nrow(X) )
 {
  if( var( c(X[i,], Y[i,] ) ) > 0 )
   {
     pv[i]  <- ( t.test( X[i,], Y[i,] ) )$p.value
   }
 }

qv <- p.adjust( pv, method="BH" )  ## looks like a lot 

sum( qv <= 0.05 )  ## I find 3361...not sure how this compares to other methods
