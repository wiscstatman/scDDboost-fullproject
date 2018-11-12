library(splatter)
library(Oscope)
library(scater)
library(BiocParallel)

##import loadDataset function for splatter simulation
source("/ua/xiuyu/simulation_data/loadDataset.r")

##load data for parameters estimation in splatter simulation
root <- "/ua/xiuyu/simulation_data/data"
datasets <- read_tsv(file.path(root, "datasets.txt"),
                     col_types = cols(.default = col_character(),
                                      NumCells = col_integer()
                     )
)
real <- loadDataset(datasets[3, ], root)

##estimation of parameters
set.seed(10)
real <- real[, sample(1:ncol(real), 200)]
real <- real[rowSums(real) > 0, ]
params <- splatEstimate(real)

#for each simulation setting, we have two replicates
bp <- BiocParallel::MulticoreParam(2)

synData <- function(Loc, Scale){
  sims <- bplapply(1:2, function(seed) {
    message("Simulating ", seed)
    sim <- splatSimulateGroups(params,
                               batchCells      = 400,
                               group.prob      = c(0.3,0.2,0.1,0.15,0.05,0.05,0.15),
                               de.prob         = 0.1,
                               de.facLoc       = Loc,
                               de.facScale     = Scale,
                               # dropout.present = FALSE,
                               seed            = seed)
    sim <- calculateQCMetrics(sim)
    return(sim)
  }, BPPARAM = bp)
  for(i in 1:2){
    if(Loc < 0){
      name_vec = c("sim","_","neg",abs(Loc),'_',Scale,'_',II,".rds")
      filename = paste(name_vec,collapse = '')
    }
    else{
      name_vec = c("sim","_",abs(Loc),'_',Scale,'_',II,".rds")
      filename = paste(name_vec,collapse = '')
    }
    saveDir = paste0("/ua/xiuyu/simulation_data/",filename)
    saveRDS(sims[[i]],  file = saveDir )
  }
}