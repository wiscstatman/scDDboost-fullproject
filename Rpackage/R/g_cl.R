#' gene_level cluster
#'
#' @param data transcripts
#' @param ncores number of cores for parallel computing
#' @return return a matrix whose row represent gene specific cluster
#' @export

g_cl = function(data, ncores){
    nr = nrow(data)
    Phi_mdf <- rep(1, nr)
    bt <- rep(1, nr)
    tryCatch(
    {
        MV<-CalcMV(data,Plot = F)
        Phi_mdf<-MV$Phi_mdf
        Q_mdf<-MV$Q_mdf
        bt<-1/Q_mdf-1
    },
    error=function(w) {
     message("estimation of hyper parameter failed, try naively assigned parameters")
        }
    ,finally = {
    bp <- BiocParallel::MulticoreParam(ncores)
    clus = bplapply(1:nr,function(i) MCP(data[i,],1,c(Phi_mdf[i],bt[i]))+1,BPPARAM = bp)
    cl = matrix(0, nrow = nrow(data), ncol = ncol(data))
    for(i in 1:length(clus)){
        cl[i,] = clus[[i]]
    }
    return(cl)
    })
}


