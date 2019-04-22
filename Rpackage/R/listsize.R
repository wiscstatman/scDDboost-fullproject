
#' number of DD genes under FDR control
#'
#' @param pDD estimated probability of being DD
#' @param FDR fdr to be controlled
#' @return number of positive genes

listsize <- function(pDD, FDR=0.01)
{
    
    ee <- 1-pDD
    
    oe <- sort(ee)
    
    ff <- cumsum(oe)/(1:length(oe))
    
    return( sum( ff <= FDR ) )
    
}
