detK = function(D,epi = 0.15)
{
     intra = rep(0,8)
     inter = rep(0,8)
    
     for(i in 2:9){
         clusRes = pam(D,i,diss = T)
         intra[i - 1] = as.numeric(clusRes$objective[1])
         x = clusRes$id.med
         inter[i - 1] = sum(D[x,x]) / (i * (i - 1))
     }
    
    s = intra / inter
    mins = min(s)
    ss = s/mins - 1
    
    if(min(ss) < epi){
    K = which(ss < epi)[1] + 1
        }else{
        K = 9
        }
    return(K)
}
