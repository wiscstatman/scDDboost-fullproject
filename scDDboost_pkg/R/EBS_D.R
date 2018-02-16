
#' adjusted distance matrix by EBSeq
#'
#' @param ccl cluster of cells
#' @param part partition patterns
#' @param Y map pattern
#' @return adjusted distance matrix by EBSeq


EBS_D = function(ccl, part, Y){
    nY = length(Y)
    nc = length(ccl)
    Ebs = t(sapply(1:nY, function(i) sapply(ccl,function(j) part[Y[i],j])))
    L = nrow(Ebs)
    ebs_D = sapply(1:nc, function(i) sapply(1:nc, function(j) if(i<j){
            return(1 - length(which(Ebs[,i] == Ebs[, j]))/L)}else{return(0)}))
    # ebs_D = matrix(0, nrow = nc, ncol = nc)
    # for(i in 1:nc){
    #    for(j in 1:nc){
    #        if(i < j){
    #            tmp = which(Ebs[,i] != Ebs[,j])
    #            ebs_D[i, j] = 1 - cor(Ebs[tmp,i], Ebs[tmp,j])
    #        }
    #    }
    #}
    ebs_D = ebs_D + t(ebs_D)
    
    #ebs_D = new_D(ccl, part, Y)
    
    
    
    return(ebs_D)
}
