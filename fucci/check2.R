
load("results.RData")


listsize <- function(pDD, FDR=0.01)
 {
  ee <- 1-pDD
  oe <- sort(ee)
  ff <- cumsum(oe)/(1:length(oe))
  return( sum( ff <= FDR ) )
 }


