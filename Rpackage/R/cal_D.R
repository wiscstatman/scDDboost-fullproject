
#' calculate distance matrix
#'
#' @param data transcripts
#' @param ncores number of cores for parallel computing
#' @return distance matrix
#' @example
#' data(sim_dat)
#' data_counts = assays(sim_dat)$count
#' D_c = cal_D(data_counts,4)

#' @export



cal_D = function(data,ncores){
    nc = ncol(data)
    nr = nrow(data)
    
    geneMax = apply(data,1,max)
    tmp = which(geneMax > 0)
    cl = g_cl(data[tmp,], ncores)
    
    geneK = apply(cl,1,max)
    tmp = which(geneK > 1)
    f_cl = cl[tmp,]
    geneK = geneK[tmp]
    
    
    ng = nrow(f_cl)
    
    
    D_E = as.matrix(dist(t(f_cl),method = "manhattan")) / nrow(f_cl)
    
    m_ = max(D_E)
    
    D_E = D_E / m_
    
    D_cor = (1 - cor(data)) / 2
   
    w1 = 1 / sd(D_E)
    
    w2 = 1 / sd(D_cor) 
    
    total = w1 + w2
    
    w1 = w1 / total
    
    w2 = w2 / total
    
    D_c = w2*D_cor + w1*D_E
    
    return(D_c)
}
