source("./synData.R")

## empirical data used to estimate splatter internal parameters can be found via dropbox
params = est_param("./data",10)

## 4 pairs location and scale paramters
LOC = c(0.1,-0.1,-0.1,0.3)
SCALE = c(0.4,0.3,1,0.5)


####### K = 3
save_DIR = "user-specified-1"
## constraint for proportions of clusters
## phi_1 + phi_2 + phi_3 = psi_1 + psi_2 + psi_3
patP = c(1,1,1)
for(i in 1:4)
{
  synData(LOC[i],SCALE[i],params,patP,save_DIR,nrep = 10,
          Seed = 10, ncells = 400, de_prob = 0.1)
}

####### K = 7
save_DIR = "user-specified-2"
## constraint for proportions of clusters
## phi_1 + phi_2 = psi_1 + psi_2
## phi_3 + phi_4 + phi_5 = psi_3 + psi_4 + psi_5
## phi_6 + phi_7 = psi_6 + psi_7
patP = c(1,1,2,2,2,3,3)
for(i in 1:4)
{
  synData(LOC[i],SCALE[i],params,patP,save_DIR,nrep = 10,
          Seed = 10, ncells = 400, de_prob = 0.1)
}

####### K = 15
save_DIR = "user-specified-3"
## constraint for proportions of clusters
## phi_1 + phi_2 = psi_1 + psi_2
## phi_3 + phi_4 + phi_5 = psi_3 + psi_4 + psi_5
## phi_6 + ... + phi_10 = psi_6 + ... + psi_10
## phi_11 + ... + phi_15 = psi_11 + ... + psi_15
patP = c(1,1,2,2,2,3,3,3,3,3,4,4,4,4,4)
for(i in 1:4)
{
  synData(LOC[i],SCALE[i],params,patP,save_DIR,nrep = 10,
          Seed = 10, ncells = 400, de_prob = 0.1)
}
