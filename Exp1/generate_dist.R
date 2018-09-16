

###generate dist based on dirichlet process
require(partitions)

get_invariant = function(pt){
  #! param pt, a valid parition
  #return a mapped array m of pt. e.g. m[i] = max_{j <= i} pt[j]
  n = length(pt)
  m = rep(1,n)
  tmp_max = 1
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
  pt = rep(1,n)
  pt[(n - K + 1):n] = 1:K
  return(pt)
}

get_last_partition = function(n, K){
  #! param n, number of elements
  #! param K, number of groups
  #return, last partition of n elements into K groups
  pt = rep(K,n)
  pt[1:K] = 1:K
  return(pt)
}

get_next_partition = function(n, K, current_pt){
  #! param n, number of elements
  #! param K, number of groups
  #! param current_pt, current partition of n elements into K groups
  #return, next partition of n elements into K groups
  
  
}



get_dist = function(n,K){
  partitions_space = setparts(K)
  point_mass = apply(partitions_space, 2, )
}





####not working
boot = function(x, K, B, in_lam, out_lam, gm, theta){
  #mu = c(1,1)
  #sig = diag(2) * 0.1
  n = length(x)
  rand_boot = rep(0,B)
  Aboot <- matrix(0,n,n)
  D_ = as.matrix(dist(x))
  for( b_ in 1:B){
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
          # print("I")
          # print(i - 1)
          a = 0
          b = 2000
          for(t_ in 1:(i - 1)){
            if(abs(noise[t_,i] - noise[t_,j]) > a){
              a = abs(noise[t_,i] - noise[t_,j])
              # print("noise_sub")
              # print("j")
              # print(j)
              # print("t_")
              # print(t_)
              # print("i")
              # print(i)
              # print(noise[t_,i])
              # print(noise[t_,j])
            }
            
            if(noise[t_,i] + noise[t_,j] < b){
              b = noise[t_,i] + noise[t_,j]
              # print("noise_add")
              # print("j")
              # print(j)
              # print("t_")
              # print(t_)
              # print("i")
              # print(i)
              # print(noise[t_,i])
              # print(noise[t_,j])
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
          # print("I")
          # print(i - 1)
          a = 0
          b = 2000
          for(t_ in 1:(i - 1)){
            if(abs(noise[t_,i] - noise[t_,j]) > a){
              a = abs(noise[t_,i] - noise[t_,j])
              # print("noise_sub")
              # print("j")
              # print(j)
              # print("t_")
              # print(t_)
              # print("i")
              # print(i)
              # print(noise[t_,i])
              # print(noise[t_,j])
            }
              
            
        
            #print(noise[t_,i] + noise[t_,j])
            if(noise[t_,i] + noise[t_,j] < b){
              b = noise[t_,i] + noise[t_,j]
              # print("noise_add")
              # print("j")
              # print(j)
              # print("t_")
              # print(t_)
              # print("i")
              # print(i)
              # print(noise[t_,i])
              # print(noise[t_,j])
            }
            
            # b = min(noise[t_,i] + noise[t_,j], b)
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
    #dst.star <- as.dist( bar )
    #hc = hclust(dst.star)
    #clus = cutree(hc, k = K)
    clus = pam(bar,K,diss = T)$clustering
    rand_boot[b_] = adjustedRandIndex(clus, true.clust)
    tmp = outer(clus,clus, "==")
    Aboot <- Aboot + tmp/B
  }
  res = list()
  res[[1]] = Aboot
  res[[2]] = rand_boot
  return(res)
}

  