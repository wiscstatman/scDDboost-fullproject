
## idea..like try1, but let's add a random matrix instead of deconvolving 

library(DPpackage)  ## for MCMC


# a toy data set

mu <- 2*c( rep(-3,15), rep(-1,20), rep(0,20), rep(1,30), rep(5,15) )
sig <- 1
true.clust <- match(mu, unique(mu))

n <- length(mu)
y <- rnorm(n,mean=mu, sd=sig)
f <- as.factor(1:n)


# try DPlmm

prior <- list( alpha=1, nu0=1, tau1=1, tau2=1, tinv=diag(1), mub=0, Sb=diag(1) )

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

dst.m <- as.matrix(dst)

B <- 1000
Aboot <- matrix(0,n,n)
lam <- mean(dst)
boot.clust <- matrix(NA, B, n )
for( b in 1:B )
 {
  print(b)
  tmp <- 1*sqrt(lam)*matrix(rgamma(n*n,shape=1),n,n)
  bar <- dst.m + (tmp+t(tmp))   ##symmtrize  ; compare to try1.pdf
  diag(bar) <- 0
  dst.star <- as.dist(bar)
  cstar <- pam(dst.star, k=5 )
  tmp <- outer(cstar$clustering,cstar$clustering,"==")
  boot.clust[b,] <- cstar$clustering
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

pdf( file="try2.pdf" )
par(mfrow=c(2,2), mgp=c(2.5,.5,0), mar=c(4,4,1,1) )

image(ABayes, axes=FALSE)

image(Aboot, axes=FALSE)

plot( ABayes, Aboot )
abline(0,1,col="red", lwd=2 )

boxplot( rev(list( Bayes_MCMC=prand.bayes, Random_Dist=prand.boot, 
 Control=prand.rand )), horizontal=TRUE, las=1, cex.axis=.5, cex.label=.8,
		xlab="adjusted Rand index to true clustering" )


dev.off()

## also do `plot(y)` to see the data and 5 groups
## or plot(y,mu)
## also, cor( c(Aboot), c(ABayes) ) is high (0.95 )


