detK = function(D,epi = 0.15)
{
     intra = rep(0,8)
     inter = rep(0,8)
    
     for(i in 2:9){
         clusRes = pam(D,i)
         intra[i - 1] = as.numeric(clusRes$objective[1])
         x = clusRes$id.med
         inter[i - 1] = sum(D[x,x]) / (i * (i - 1))
     }
    
    s = intra / inter
    
    ss = abs(1 - s[2:8] / s[1:7])
    
    if(min(ss) < epi){
    K = which.max(ss < epi)[1] + 2
        }else{
        K = 9
        }
    return(K)
}