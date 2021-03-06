
# like try3, but dividing distances by random exponential weights 
# [started 11/26/18]

library(DPpackage)  ## for MCMC

set.seed(75751)

# a toy data set
sz = 1
mu <- 2*c( rep(-3,15 * sz), rep(-1,20 * sz), rep(0,20 * sz), rep(1,30 * sz), rep(5,15 * sz) )
sig <- 1
true.clust <- match(mu, unique(mu))

n <- length(mu)
y <- rnorm(n,mean=mu, sd=sig)
f <- as.factor(1:n)


# try DPlmm

prior <- list( alpha=1, nu0=1, tau1=1, tau2=1, tinv=diag(1), mub=0, Sb=diag(1) ) ## try7-b.pdf
##prior <- list( alpha=1, nu0=1, tau1=0.1, tau2=0.1, tinv=diag(1), mub=0, Sb=diag(1) )  #try7.pdf

nburn <-50
nsave <- 1e4
nskip <-100
ndisplay <-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

state <- NULL
fit <- DPlmm( fixed=y~1, random=~1|f, mcmc=mcmc, state=state, status=TRUE, prior=prior )

# let me extract pairwise posterior co-cluster probability from saved stateo

u <- fit$save.state$randsave[,-(n+1)]
ok <- fit$save.state$thetasave[,5] == 5  ### 5 clusters


u <- fit$save.state$randsave[ok,-(n+1)]
B <- sum(ok)
bayes.clust <- matrix(NA, B, n)
for( b in 1:B )
 {
  tmp <- u[b,]
  bayes.clust[b,] <- match( tmp, unique(tmp) )
 }

ABayes <- matrix(0, n,n )
for( i in 1:B )
 {
  tmp <- outer( u[i,], u[i,], "==" )
  ABayes <- ABayes + tmp/B
 }

image(ABayes)

## so that seems to work; now try randomized clustering with Gamma(1/2,1) weights
## and maybe some distance shrinkage adjustment

#e.g
library(cluster)

dst <- dist(y)

dst.m <- (as.matrix(dst))^2

LOG_u = function(x){
    if(x > 0){
        return(log(x))
    }
    else{
        return(0)
    }
}

LOG = function(x){
    sapply(x,LOG_u)
}

LL = function(param, x){
    a0 = param[1]   #shape for prior
    d0 = param[2]   #rate for prior
    #a1 = param[3]   #shape for sampling model
    #d0 = 2 * sqrt(a0)
    a1 = 1
    
    n = length(x)
    
    nc = ncol(x)
    
    I = matrix(1,nc,nc)
    
    I = I - diag(nc)
    
    C = d0 + a1 * x
    
    res = (n - nc) * (lgamma(a0 + a1) - lgamma(a0) - lgamma(a1) + a0 * log(d0) + a1 * log(a1)) + sum((a1 - 1) * log(x + diag(nc))) - sum((a0 + a1) * log(C) * I)
    
    return(-res)
}

fit3 <- suppressMessages(nlminb( start=c(0.1,0.1), objective=LL, x=dst.m, lower=c(0,0) , upper = c(2,Inf)))

fit3 <- suppressMessages(nlminb( start=c(0.1,0.1), objective=LL, x=D_c, lower=c(0,0) , upper = c(1,Inf)))


a0 = fit3$par[1]
d0 = fit3$par[2]
#a1 = fit3$par[3]

a = a0 + 1

B <- 1000
Aboot <- matrix(0,n,n)
lam <- mean(dst)
lam2 <- mean(dst^2)
boot.clust <- matrix(NA, B, n )
for( b in 1:B )
 {
  print(b)
#   e <- rgamma(n,shape=(1/2))   # makes Gamma[1,] weights; try7.pdf and -b.pdf
#   e <- rgamma(n,shape=(1))   # makes Gamma[2,] weights;  -c.pdf
#   e <- rgamma(n,shape=(1/4), rate=(1/4) )   # makes Gamma[1/2,1/2] weights;  -d.pdf
#   e <- rgamma(n,shape=(1/5), rate=(1/5) )   # makes Gamma[2/5,2/5] weights;  -e.pdf
#e1 <- rgamma(n^2,shape=(a + 1), rate=(a) )   # makes Gamma[1/4,1/4] weights;  -f.pdf


  e2 = rgamma(n^2, shape =  (a0 + 1) , rate = a0)
  bar <- dst.m/matrix(e2, ncol = n)
  #bar = dst.m/(matrix(e1, ncol = n) + e3)
  dst.star <- as.dist(bar)
  #cstar <- pam(dst.star, k=5 ,diss =T)$clustering
  cstar = cutree(hclust(dst.star), k = 5)
  tmp <- outer(cstar,cstar,"==")
  boot.clust[b,] <- cstar
  Aboot <- Aboot + tmp/B
 }

## process

library(mcclust) ## looks at lots of clusterings and summarizes

diag(ABayes) <- 1; diag(Aboot) <- 1 ## removing some numerical inaccuracies
mode.bayes <- minbinder(ABayes)$cl

# having trouble with Aboot
tmp <- .99*Aboot + .01*diag(n)
mode.boot <- minbinder(tmp)$cl   ## slightly different modal clusters

prand.bayes <- numeric( nrow(bayes.clust))
prand.boot <- numeric( nrow(boot.clust))

for( b in 1:length(prand.bayes) )
 {
  prand.bayes[b] <- arandi( bayes.clust[b,] , true.clust )
 }

for( b in 1:length(prand.boot) )
 {
  prand.boot[b] <- arandi( boot.clust[b,] , true.clust )
 }

## distributions under random partitions, for comparison

B <- 1000
Arand <- matrix(0,n,n)
rand.clust <- matrix(NA, B, n )
prand.rand <- numeric( nrow(rand.clust))
for( b in 1:B )
 {
  print(b)
  e <- rgamma(n,shape=1)
  ##  bar <-  e %o% e
  bar <- outer( e, e, "+" )
  dst.star <- as.dist(bar)
  cstar <- pam(dst.star, k=5 )
  tmp <- outer(cstar$clustering,cstar$clustering,"==")
  rand.clust[b,] <- cstar$clustering
  Arand <- Arand + tmp/B
  prand.rand[b] <- arandi( rand.clust[b,] , true.clust )
 }

## some plots

#pdf( file="try7.pdf" )
#pdf( file="try7-b.pdf" )
#pdf( file="try7-c.pdf" )   ## Gamma[2,] weights; cor = 0.915
#pdf( file="try7-d.pdf" )   ## Gamma[1/2,] weights;  cor = 0.94
#pdf( file="try7-e.pdf" )   ## Gamma[2/5,2/5] weights;  cor = 0.945
pdf( file="try7-g.pdf" )   ## Gamma[1/4,1/4] weights;  cor = 0.947

par(mfrow=c(2,2), mgp=c(2.5,.5,0), mar=c(4,4,1,1) )

image(ABayes, axes=FALSE, main="Bayes")

image(Aboot, axes=FALSE, main="Random weight (divisors)")

plot( ABayes, Aboot )
abline(0,1,col="red", lwd=2 )

boxplot( rev(list( Bayes_MCMC=prand.bayes, Random_Dist=prand.boot, 
 Control=prand.rand )), horizontal=TRUE, las=1, cex.axis=.5, cex.label=.8,
		xlab="adjusted Rand index to true clustering" )


dev.off()

# it might be working

## also, cor( c(Aboot), c(ABayes) ) is high (

# check : plot( dst, Aboot ) and plot( dst, ABayes )
