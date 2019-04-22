
#' index of DD genes under FDR control
#'
#' @param pDD probability of genes being DD
#' @param FDR fdr to be controlled
#' @return index of positive genes


lsz = function(pDD, FDR=0.01)

{
    
    ee <- 1-pDD
    
    oe <- sort(ee)
    
    or = order(ee)
    
    ff <- cumsum(oe)/(1:length(oe))
    
    return(or[which(ff < FDR)])
    
}
