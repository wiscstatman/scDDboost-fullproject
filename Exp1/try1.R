
library(DPpackage)  ## for MCMC


# a toy data set

mu <- 2*c( rep(-3,15), rep(-1,20), rep(0,20), rep(1,30), rep(5,15) )
sig <- 1
#sig <- 0.1
#sig <- 1/3

n <- length(mu)
y <- rnorm(n,mean=mu, sd=sig)
f <- as.factor(1:n)


# try DPlmm

prior <- list( alpha=1, nu0=1, tau1=1, tau2=1, tinv=diag(1), mub=0, Sb=diag(1) )

nburn <-50
nsave <-5000
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
for( b in 1:B )
 {
  print(b)
  e <- 2*sqrt(lam)*rgamma(n,shape=1/2)
  ##e <- 2*sqrt(lam)*rgamma(n,shape=1)
  bar <- dst.m + e %o% e
  dst.star <- as.dist(bar)
  cstar <- pam(dst.star, k=5 )
  tmp <- outer(cstar$clustering,cstar$clustering,"==")
  Aboot <- Aboot + tmp/B
 }
