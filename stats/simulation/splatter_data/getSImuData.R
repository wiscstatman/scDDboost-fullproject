###cases of 4 subtypes. set seed 10, pat phi_1 + phi_2 = psi_1 + psi_2
# at server
root_DIR = "~/data" 
params = est_param(root_DIR,10)
save_DIR = "~/simulation_data/K4/"

##location and scale paramters
loc = c(0.1,-0.1,-0.1,0.3)
scale = c(0.4,0.3,1,0.5)

patP = c(1,1,2,2)

for(i in 1:4)
{
  synData(loc[i],scale[i],params,patP,save_DIR,nrep = 2, 
          Seed = 10, ncells = 400, de_prob = 0.1)
}




###cases of 7 subtypes. pat phi_1 + phi_2 = psi_1 + psi_2, phi_3 + phi_4 + phi_5 = ....

save_DIR = "~/simulation_data/K7/"

##location and scale paramters
loc = c(0.1,-0.1,-0.1,0.3)
scale = c(0.4,0.3,1,0.5)

patP = c(1,1,2,2,2,3,3)

for(i in 1:4)
{
  synData(loc[i],scale[i],params,patP,save_DIR,nrep = 2, 
          Seed = 10, ncells = 400, de_prob = 0.1)
}


###cases of 15 subtypes. pat phi_1 + phi_2 = psi_1 + psi_2, phi_3 + phi_4 + phi_5 = ....

save_DIR = "~/simulation_data/K15/"

##location and scale paramters
loc = c(0.1,-0.1,-0.1,0.3)
scale = c(0.4,0.3,1,0.5)

patP = c(1,1,2,2,2,3,3,3,3,3,4,4,4,4,4)

for(i in 1:4)
{
  synData(loc[i],scale[i],params,patP,save_DIR,nrep = 2, 
          Seed = 10, ncells = 400, de_prob = 0.1)
}


