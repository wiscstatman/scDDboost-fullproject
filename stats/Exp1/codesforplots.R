library(DPpackage)  ## for MCMC

#set.seed(75751)

# a toy data set

mu <- 2*c( rep(-3,15), rep(-1,20), rep(0,20), rep(1,30), rep(5,15) )
sig <- 1
true.clust <- match(mu, unique(mu))

n <- length(mu)
y <- rnorm(n,mean=mu, sd=sig)
f <- as.factor(1:n)


# try DPlmm

#prior <- list( alpha=1, nu0=1, tau1=1, tau2=1, tinv=diag(1), mub=0, Sb=diag(1) )
prior <- list( alpha=1, nu0=1, tau1=0.1, tau2=0.1, tinv=diag(1), mub=0, Sb=diag(1) )

nburn <-50
nsave <- 1e3
nskip <-100
ndisplay <-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)


DP_bayes = function(x, mcmc, prior, K){
  n = length(x)
  state = NULL
  fit = DPlmm( fixed=x~1, random=~1|f, mcmc=mcmc, state=state, status=TRUE, prior=prior )
  u <- fit$save.state$randsave[,-(n+1)]
  ok <- fit$save.state$thetasave[,5] == K
  u <- fit$save.state$randsave[ok,-(n+1)]
  B <- sum(ok)
  ABayes <- matrix(0, n,n )
  rand_bayes = rep(0,B)
  for( b in 1:B )
  {
    tmp <- outer( u[b,], u[b,], "==" )
    ABayes <- ABayes + tmp/B
    tmp = u[b,]
    tt <- match( tmp, unique(tmp) )
    print(length(tt))
    rand_bayes[b] = adjustedRandIndex(tt,true.clust)
  }
  res = list()
  res[[1]] = ABayes
  res[[2]] = rand_bayes
  return(res)
}

res_bayes = DP_bayes(y, mcmc, prior, 5)
ABayes = res_bayes[[1]]
rand_bayes =res_bayes[[2]]
image(ABayes)

boot = function(x, K, B){
  #mu = c(1,1)
  #sig = diag(2) * 0.1
  n = length(x)
  rand_boot = rep(0,B)
  Aboot <- matrix(0,n,n)
  D_ = as.matrix(dist(x))
  for( b_ in 1:B){
    #e <-rgamma(n,shape = 0.5 ,rate = 2)
    #e = rgamma(n * n, shape = 0.5, rate = 10)
    #j_ = sample(1:n,1)
    #e = rexp(n,0.5)
    e_ = runif(n,0,sum(D_) / (n * (n - 1)))
    #e_ = runifdisc(n,2)
    #e_ =cbind(e_[[3]],e_[[4]])
    e_[1] = 0
    #noise = as.matrix((dist(e_)))
    noise = matrix(0, n,n)
    noise[1,] = e_
    for(i in 2:(n - 1)){
      for(j in (i + 1):n){
        # print("I")
        # print(i - 1)
        a = 0
        b = 2000
        for(t_ in 1:(i - 1)){
          a = max(abs(noise[t_,i] - noise[t_,j]), a)
          #print(noise[t_,i] + noise[t_,j])
          b = min(noise[t_,i] + noise[t_,j], b)
        }
        # print("inter")
        # print(a)
        # print(b)
        if(j > i + 1){
          for(tt_ in (i + 1):(j - 1)){
            tmp_b = 10000
            for(t_ in 1:(i - 1)){
              tmp_b = min(noise[t_, tt_] + noise[t_, j], tmp_b)
            }
            a = max(noise[i,tt_] - tmp_b, a)
          }
        }
        # print("i")
        # print(i)
        # print("j")
        # print(j)
        # print("2222")
        # print(a)
        # print(b)
        # print("3333")
        tmp_e = runif(1,a,b)
        print("new")
        print(tmp_e)
        noise[i,j] = tmp_e
      }
    }
    noise = noise + t(noise)
    
    #mu_ = 0
    #p_ = rbinom(1,1,0.1)
    #e = rchisq(n * n, df = 1, ncp = mu_ * p_)
    #e = matrix(e,nrow = n)
    #bar <- outer( e, e, "+" ) + D_
    # bar = e + D_
    #bar = noise + D_
    bar = noise + D_
    dst.star <- as.dist( bar )
    hc = hclust(dst.star)
    clus = cutree(hc, k = K)
    #clus = pam(bar,K,diss = T)$clustering
    rand_boot[b_] = adjustedRandIndex(clus, true.clust)
    tmp = outer(clus,clus, "==")
    Aboot <- Aboot + tmp/B
  }
  res = list()
  res[[1]] = Aboot
  res[[2]] = rand_boot
  return(res)
}









y <- rnorm(n,mean=mu, sd=sig)
res = boot(y,5,100)
Aboot = res[[1]]
rand_boot = res[[2]]
res_bayes = DP_bayes(y, mcmc, prior, 5)
ABayes = res_bayes[[1]]
rand_bayes = res_bayes[[2]]
#res1 = boot1(y,5)
#Aboot = res1[[1]]
pdf("heat2.pdf")
par(mfrow = c(2,1))
#image(Aboot, axes = F, main = "2dDisc")
image(ABayes, axes=FALSE, main="Bayes")
image(Aboot,axes = F, main = "dist")
#plot( ABayes, Aboot )
#abline(0,1,col="red", lwd=2 )

#boxplot( rev(list( Bayes_MCMC=rand_bayes, Random_Dist=rand_boot 
              #     )), horizontal=TRUE, las=1, cex.axis=.5, cex.label=.8,
        # xlab="adjusted Rand index to true clustering" )
dev.off()
#boxplot(rand_boot, rand_bayes)



K = 5
res1 = c()
res2 = c()
N = 20
rand_bayes = list()
rand_boot = list()
for (i in 1:N){
  y <- rnorm(n,mean=mu, sd=sig)
  res_bayes = DP_bayes(y, mcmc, prior, K)
  res_boot = boot(y, K, 100)
  ABayes = res_bayes[[1]]
  rand_bayes[[i]] = res_bayes[[2]]
  Aboot = res_boot[[1]]
  rand_boot[[i]] = res_boot[[2]]
  tmp1 = as.vector(ABayes)
  tmp2 = as.vector(Aboot)
  print("spearman")
  print(cor(tmp1,tmp2,method="spearman"))
  print("pearson")
  print(cor(tmp1,tmp2))
  res1 = c(res1, tmp1)
  res2 = c(res2, tmp2)
}

pdf( file="try6.pdf" )
par(mfrow=c(2,1), mgp=c(2.5,.5,0), mar=c(4,4,1,1) )

hist(res1)

hist(res2)

dev.off()









