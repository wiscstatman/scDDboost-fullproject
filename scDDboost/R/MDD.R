


lpzgt <- function(z,pp)
{
    # log prob
    # z a vector of counts
    # pp a partition
    # p(z|t) assuming alpha_k=1;
    tt <- tapply(z,pp,sum) ## sufficient stats over blocks
    nn <- table(pp)  ## counts block sizes
    res <- sum(lgamma(tt+1)) + sum(lgamma(nn)) - sum(lgamma(tt+nn))
    res
}

lpt1t2 <- function(z1,z2,pp)
{
    K <- length(z1)
    n1 <- sum(z1); n2 <- sum(z2)
    t1 <- tapply(z1,pp,sum) ## sufficient stats over blocks
    t2 <- tapply(z2,pp,sum) ## sufficient stats over blocks
    nn <- table(pp)  ## counts block sizes
    tmp1 <- lgamma(n1+1) + lgamma(n2+1) - sum(lgamma(t1+1)) - sum(lgamma(t2+1))
    tmp2 <- lgamma(K) - sum( lgamma(nn) )
    tmp3 <- sum(lgamma(t1+t2+nn))-lgamma(n1+n2+K)
    res <- tmp1+tmp2+tmp3
    res
}

#' posterior of proportion change given mixture double dirichlet prior
#'
#' @param z1 counts of each group in condition 1
#' @param z2 counts of each group in condition 2
#' @param pat partition patterns
#' @return posterior of proportion change

MDD = function(z1, z2, pat){
    np = nrow(pat)
    lpz <- numeric( np )
    for( j in 1:nrow(pat) )
    {
        pp <- pat[j,]
        lpz[j] <- lpzgt(z1,pp)+lpzgt(z2,pp)+lpt1t2(z1,z2,pp)
    }
    post <- exp(lpz-max(lpz))
    post <- post/sum(post)
    return(post)
}
