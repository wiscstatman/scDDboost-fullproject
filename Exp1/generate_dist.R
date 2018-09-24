

###generate dist based on dirichlet process
require(partitions)

get_max_seq = function(pt){
  #! param pt, a valid parition
  #return a mapped array m of pt. e.g. m[i] = max_{j <= i} pt[j]
  n = length(pt)
  m = rep(1L,n)
  tmp_max = 1L
  for(i in 1:n){
    tmp_max = max(tmp_max,pt[i])
    m[i] = tmp_max
  }
  return(m)
}


get_first_partition = function(n,K){
  #! param n, number of elements
  #! param K, number of groups
  #return, first partititon of n elements into K groups
  if(typeof(K) != "integer"){
    K = as.integer(K)
  }
  pt = rep(1L,n)
  pt[(n - K + 1):n] = 1L:K
  return(pt)
}

get_last_partition = function(n, K){
  #! param n, number of elements
  #! param K, number of groups
  #return, last partition of n elements into K groups
  if(typeof(K) != "integer"){
    K = as.integer(K)
  }
  pt = rep(K,n)
  pt[1:K] = 1L:K
  return(pt)
}

get_next_partition = function(n, K, current_pt){
  #! param n, number of elements
  #! param K, number of groups
  #! param current_pt, current partition of n elements into K groups
  #return, next partition of n elements into K groups
  m = get_max_seq(current_pt)
  if(typeof(K) != "integer"){
    K = as.integer(K)
  }
  for(i in n:2){
    if( (current_pt[i] < K) && (current_pt[i] <= m[i - 1]) ){
      current_pt[i] = current_pt[i] + 1L
      m[i] = max(m[i], current_pt[i])
      j = i + 1L
      while(j  <= (n - K + m[i])){
        current_pt[j] = 1L
        m[j] = m[i]
        j = j + 1L
      }
      j = (n - K + m[i] + 1L)
      while(j <= n){
        m[j] = K - n + j
        current_pt[j] = m[j]
        j = j + 1L
      }
      return(current_pt)
    }
  }
}



get_post = function(n,K){
  #! param n, number of elements
  #! param K, number of groups
  # return distance matrix of elements that two elements belong to same cluster 
  
  weight = c()
  D_ = matrix(0,nrow = n, ncol = n)
  start = get_first_partition(n,K)
  end = get_last_partition(n,K)
  first = paste(start,collapse = "")
  last = paste(end,collapse = "")
  ##index for weight
  i_ = 1
  while(first < last){
    D = 1 * outer(start,start,"==")
    weight = c(weight,exp(sum(lgamma(table(start)))))
    D_ = D_ + D * weight[i_]
    start = get_next_partition(n,K,start)
    first = paste(start,collapse = "")
    i_ = i_ + 1
  }
  D = 1 * outer(start,start,"==")
  weight = c(weight,exp(sum(lgamma(table(start)))))
  D_ = D_ + D * weight[i_]
  res = list()
  res[[1]] = D_
  res[[2]] = weight
  return(res)
}


count_partitions = function(n,k){
  #! param n, number of elements
  #! param k, number of groups
  # return number of partitions of n elements into k groups
  dp = matrix(0, nrow = n + 1, ncol = k + 1)
  
  for (i in 2:(n + 1)){
    for (j in 2:(k + 1)){
      if (j == 2 || i == j){
        dp[i,j] = 1
      }
      else{
        dp[i,j] = (j - 1) * dp[i-1,j] + dp[i-1,j-1]
      }
    }
  }
      
  return(dp[n + 1 , k + 1])
  
}


EB_prob = function(n,k){
  ##typically n is much larger than k
  #return the posterior prob that two elements belong to same cluster
  
  
  # base = count_partitions(n,k)
  ##
  # a = count_partitions(n - 2, k - 1)
  ##
  # b = count_partitions(n - 2, k)
  
  
  # return((b * k + a) / base)
  
  
  return( 1 / k)
}



mimic_dp = function(n, K){
  #! param n, number of elements
  #! param K, number of clusters
  alp = K - 1
  if(alp < 1){
    print("error, K should be greater than 1")
    return(0)
  }
  p_ = rep(0,n)
  #number of elements added, start with 1 element
  num = 1
  #number of elements per group
  num_per_cl = rep(0,n)
  num_per_cl[1] = 1
  #initial prob for assigning class under dirichlet process
  p_[1] = num_per_cl[1] / (num + alp)
  p_[2] = alp / (num + alp)
  
  #class label
  cl_label = rep(0, n)
  cl_label[0] = 1
  pool = 1:n
  for(i in 2:n){
    idx_ = sample(pool,1,F,p_)
    num_per_cl[idx_] = num_per_cl[idx_] + 1
    num = num + 1
    p_ = num_per_cl / (num + alp)
    # position for staring new class
    pos = which(p_ == 0)[1]
    p_[pos] = alp / (num + alp)
    cl_label[i] = idx_
  }
  noise = 1 * outer(cl_label,cl_label,"==")
  
  return(noise)
  
}



boot = function(n, K, x, B){
  # p_ = 1 / K
  D_ = as.matrix(dist(x))
  rand_boot = rep(0,B)
  Aboot <- matrix(0,n,n)
  
  for(b_ in 1:B){
    # noise = matrix(rbinom(n^2, prob=1-p_,size=1),nrow=n)
    # up_ = 1 * upper.tri(noise, diag = FALSE)
    # noise = noise * up_
    # noise = noise + t(noise)
    noise = mimic_dp(n,K)
    weight = 1
    weight = mean(D_) + mean(noise)
    noise = noise * weight
    bar = noise + D_
    dst.star <- as.dist( bar )
    hc = hclust(dst.star)
    clus = cutree(hc, k = K)
    # clus = pam(bar,K,diss = T)$clustering
    rand_boot[b_] = adjustedRandIndex(clus, true.clust)
    tmp = outer(clus,clus, "==")
    Aboot <- Aboot + tmp/B
    
  }
  res = list()
  res[[1]] = Aboot
  res[[2]] = rand_boot
  return(res)
  
}









####not working
boot3 = function(x, K, B, in_lam, out_lam, gm, theta){
  #mu = c(1,1)
  #sig = diag(2) * 0.1
  n = length(x)
  rand_boot = rep(0,B)
  Aboot <- matrix(0,n,n)
  D_ = as.matrix(dist(x))
  for(b_ in 1:B){
    noise = matrix(0, n, n)
    count = 1
    for(j in 2:n){
      flag = rbinom(1, 1, p = theta / (count + theta))
      if (flag  > 0){
        d = rexp(1,out_lam)
        noise[1,j] = d
        if(j < 3){
          next
        }
        for(i in 2:(j -1)){
          a = 0
          b = 2000
          for(t_ in 1:(i - 1)){
            if(abs(noise[t_,i] - noise[t_,j]) > a){
              a = abs(noise[t_,i] - noise[t_,j])
            }
            
            if(noise[t_,i] + noise[t_,j] < b){
              b = noise[t_,i] + noise[t_,j]
            }
          }
          tmp_e = runif(1,a,b)
          noise[i,j] = tmp_e
        }
      }
      else{
        d = rexp(1,in_lam)
        noise[1,j] = d
        if(j < 3){
          next
        }
        for(i in 2:(j -1)){
          a = 0
          b = 2000
          for(t_ in 1:(i - 1)){
            if(abs(noise[t_,i] - noise[t_,j]) > a){
              a = abs(noise[t_,i] - noise[t_,j])
            }
              
            if(noise[t_,i] + noise[t_,j] < b){
              b = noise[t_,i] + noise[t_,j]
            }
            
          }
        if(a > b){
          print(a)
          print(b)
        }
        tmp_e = runif(1,a,b)
        noise[i,j] = tmp_e
        }
      }
      count = count + 1
      }
    
    noise = noise + t(noise)
    noise = noise * gm #sum(D_) / (n * (n - 1)) 
    bar = noise + D_
    dst.star <- as.dist( bar )
    hc = hclust(dst.star)
    clus = cutree(hc, k = K)
    # clus = pam(bar,K,diss = T)$clustering
    rand_boot[b_] = adjustedRandIndex(clus, true.clust)
    tmp = outer(clus,clus, "==")
    Aboot <- Aboot + tmp/B
  }
  res = list()
  res[[1]] = Aboot
  res[[2]] = rand_boot
  return(res)
}


##simple uniform
boot4 = function(x, K, B, gm){
  #mu = c(1,1)
  #sig = diag(2) * 0.1
  n = length(x)
  rand_boot = rep(0,B)
  Aboot <- matrix(0,n,n)
  D_ = as.matrix(dist(x))
  for(b_ in 1:B){
    noise = matrix(0, n, n)
    noise[1,] = runif(n,0,1)
    noise[1,1] = 0
    for(j in 2:n){
        if(j < 3){
          next
        }
        for(i in 2:(j -1)){
          a = 0
          b = 2000
          for(t_ in 1:(i - 1)){
            if(abs(noise[t_,i] - noise[t_,j]) > a){
              a = abs(noise[t_,i] - noise[t_,j])
            }
            if(noise[t_,i] + noise[t_,j] < b){
              b = noise[t_,i] + noise[t_,j]
            }
          }
          tmp_e = runif(1,a,b)
          noise[i,j] = tmp_e
        }
      }
    
    
    noise = noise + t(noise)
    noise = noise * gm #sum(D_) / (n * (n - 1)) 
    bar = noise + D_
    dst.star <- as.dist( bar )
    hc = hclust(dst.star)
    clus = cutree(hc, k = K)
    # clus = pam(bar,K,diss = T)$clustering
    rand_boot[b_] = adjustedRandIndex(clus, true.clust)
    tmp = outer(clus,clus, "==")
    Aboot <- Aboot + tmp/B
  }
  res = list()
  res[[1]] = Aboot
  res[[2]] = rand_boot
  return(res)
}



