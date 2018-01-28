
#' calculate distance matrix
#'
#' @param data transcripts
#' @param ncores number of cores for parallel computing
#' @return distance matrix




cal_D = function(data,ncores){
    nc = ncol(data)
    nr = nrow(data)
    cl = g_cl(data + 1, ncores)
    D_cor = matrix(0, nrow = nc, ncol = nc)
    
    D_cor = sapply(1:nc, function(x) sapply(1:nc, function(y) {if(x < y){
        #tmp = which(cl[,x] != cl[,y])
        #return(1 - cor(data[tmp, x], data[tmp, y]))
        return(1 - cor(data[,x], data[,y]))
    }
    else{
        return(0)
    }
    }))
    
    D_E = sapply(1:nc, function(i) sapply(1:nc, function(j)
    if(i<j){length(which(cl[,i] != cl[,j]))/nr}else{0}))
    # for(i in 1:nc){
    #  for(j in 1:nc){
    #  if(i<j){
    #     tmp = which(cl[,i] != cl[,j])
    #      D_cor[i,j] = 1 - cor(data[tmp, i], data[tmp, j])
    #      }
    #
    #          }
    #        }
    D_cor = D_cor + t(D_cor)
    
    D_E = D_E + t(D_E)
    
    w1 = sd(D_E) / mean(D_E)
    
    w2 = sd(D_cor) / mean(D_cor)
    
    total = w1 + w2
    
    w1 = w1 / total
    
    w2 = w2 / total
    
    D_c = w2*D_cor + w1*D_E
    
    return(D_c)
}
