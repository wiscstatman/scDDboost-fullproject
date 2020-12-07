library(splatter)
library(Oscope)
library(scater)
library(BiocParallel)
library(readr)
library(magrittr)
library(SingleCellExperiment)
library(MCMCpack)

loadDataset <- function(dataset, root) {
    
    #file <- file.path(root, dataset["Path"], dataset["CountsFile"])
    file <- file.path(root, dataset["CountsFile"])
    
    if (dataset["BpipePipeline"] == "Yes") {
        counts <- readr::read_tsv(file, comment = "#",
        col_types = readr::cols(
        .default = readr::col_integer(),
        Geneid = readr::col_character(),
        Chr = readr::col_character(),
        Start = readr::col_character(),
        End = readr::col_character(),
        Strand = readr::col_character()
        )
        )
        
        counts <- dplyr::select(counts, -Geneid, -Chr, -Start, -End, -Strand)
        counts <- as.matrix(counts)
        
    } else if (dataset["Conquer"] == "Yes") {
        counts <- readr::read_tsv(file,
        col_types = readr::cols(
        .default = readr::col_double(),
        Gene = readr::col_character()
        )
        )
        
        counts <- dplyr::select(counts, -Gene)
        counts <- as.matrix(counts)
        
    } else if (dataset["Dataset"] == "Grun") {
        counts <- readr::read_tsv(file,
        col_types = readr::cols(
        .default = readr::col_double(),
        GENENAME = readr::col_character()
        )
        )
        
        counts <- dplyr::select(counts, -GENENAME)
        counts <- as.matrix(counts)
    } else if (dataset["Dataset"] == "Klein") {
        counts <- readr::read_csv(file,
        col_types = readr::cols(
        .default = readr::col_integer(),
        Gene = readr::col_character()
        ),
        col_names = c("Gene", paste0("S", 1:239)),
        skip = 1)
        
        counts <- dplyr::select(counts, -Gene)
        counts <- as.matrix(counts)
    } else if (dataset["Dataset"] == "Tung") {
        counts <- readr::read_tsv(file,
        col_types = readr::cols(
        .default = readr::col_integer(),
        X1 = readr::col_character()
        ))
        
        counts <- dplyr::select(counts, -X1)
        counts <- as.matrix(counts)
    } else if (dataset["Dataset"] == "Zeisel") {
        counts <- read_tsv(file,
        col_types = readr::cols(
        .default = readr::col_character()
        ),
        col_names = FALSE,
        skip = 11)
        
        counts <- dplyr::select(counts, -X1)
        counts <- as.matrix(counts)
        counts <- apply(counts, 2, as.numeric)
    } else if (dataset["Dataset"] == "Zieg") {
        counts <- readr::read_tsv(file,
        col_types = readr::cols(
        .default = readr::col_integer(),
        Gene = readr::col_character()
        ),
        col_names = c("Gene", paste0("S", 1:482)),
        skip = 1)
        counts <- dplyr::select(counts, -Gene)
        counts <- as.matrix(counts)
    } else if (dataset["Dataset"] == "Velten") {
        counts <- readr::read_csv(file,
        col_types = readr::cols(
        .default = col_integer(),
        X1 = col_character()
        ))
        counts <- dplyr::select(counts, -X1)
        counts <- as.matrix(counts)
    } else {
        stop("Dataset not valid")
    }
    
    rownames(counts) <- paste0("G", seq_len(nrow(counts)))
    return(counts)
}

##estimation of parameters for simulation
est_param <- function(root_DIR, Seed){
    ##load Data for parameters estimation of splatter
    root <- root_DIR
    datasets <- read_tsv(file.path(root, "datasets.txt"),
    col_types = cols(.default = col_character(),
    NumCells = col_integer()
    )
    )
    real <- loadDataset(datasets[3, ], root)
    set.seed(Seed)
    real <- real[, sample(1:ncol(real), 200)]
    real <- real[rowSums(real) > 0, ]
    params <- splatEstimate(real)
    return(params)
}



#for each simulation setting, we have ten replicates
synData <- function(loc, scale, params, patP, save_DIR,
                    nrep = 2, Seed = 10,
                    ncells = 400, de_prob = 0.1){
  
    set.seed(Seed)
  
    #number of groups
    K = length(patP)
    
    b = table(patP)
    
    bb = as.numeric(b)
    
    subNames = as.numeric(names(b))
  
    BETA = 2 * as.numeric(b)  
    
    subK = length(b)
    
    PHI = rdirichlet(1, BETA)
    
    phi = list()
    
    psi = list()
    
    for(i in 1:subK)
    {
        alpha = rep(1,bb[i])
        
        phi[[i]] = rdirichlet(1,alpha) * PHI[i]
        
        psi[[i]] = rdirichlet(1,alpha) * PHI[i]
        
    }
  
    phi = as.numeric(do.call(cbind,phi))
    
    psi = as.numeric(do.call(cbind,psi))
    
    ##having zero prob that phi and psi have one component same unless in the patP specify some params to be the same
    
    ##prob for bernoulli indicator that a cell being sampled into condition 1
    
    PR = rep(0,K)
    
    for(i in 1:K)
    {
        PR[i] = phi[i] / (phi[i] + psi[i])
    }
    
    group_prob = (phi + psi) / 2
    
    bp <- BiocParallel::MulticoreParam(nrep)
    sims <- bplapply(1:nrep, function(seed) {
    message("Simulating ", seed)
    sim <- splatSimulateGroups(params,
                               batchCells      = ncells,
                               group.prob      = group_prob,
                               de.prob         = de_prob,
                               de.facLoc       = loc,
                               de.facScale     = scale,
                               # dropout.present = FALSE,
                               seed            = seed)
    sim <- calculateQCMetrics(sim)
    return(sim)
    }, BPPARAM = bp)
    
  
  for(i in 1:nrep){
      dat = sims[[i]]
      
      ##get simulated counts and group label
      data_counts_all = assays(dat)$counts
      
      group = colData(dat)$Group
      
      ##change Group1 to 1, convinient for further usage
      group = as.numeric(sapply(group,function(x) gsub("Group","",x)))
                                                               
      
      # indicator for a cell being conditon 1 or 2
      indic = rbinom(ncells, 1, PR[group])
      
      ### make first group belong to condition 1 and last group belong to condition2
      g1 = which(indic == 1)
      
      g2 = which(indic == 0)
      
      label1 = group[g1]
      
      label2 = group[g2]
      
      data_counts = data_counts_all[,c(g1,g2)]
      
      cd = c(rep(1,length(g1)),rep(2,length(g2)))
      
      trueLabel = c(label1,label2)
      
      ## set may have problem                          
      eq = c()
                                
      for(j in 1:K)
      {
          if(phi[j] == psi[j])
          {
              eq = c(eq, j)
          }
      }
      ###get index of DD genes
      tmp = c()
      for(I in 1:K){
          if(!(I %in% eq))
          {
              a = paste0("rowData(dat)$DEFacGroup",I)
              tmp = union(tmp,which(eval(parse(text = a)) != 1))
          }
      }
      
      DD = tmp
      ED = setdiff(1:nrow(data_counts),DD)
      
      
      if(loc < 0){
        name_vec = c("sim","_","neg",abs(loc),'_',scale,'_',i,".rds")
        filename = paste(name_vec,collapse = '')
      }
      else{
        name_vec = c("sim","_",abs(loc),'_',scale,'_',i,".rds")
        filename = paste(name_vec,collapse = '')
      }
      saveDir = paste0(save_DIR, filename)
      save(dat,data_counts,cd,g1,g2,group,DD,ED,label1,label2,trueLabel,phi,psi,file = saveDir)
  }
}
