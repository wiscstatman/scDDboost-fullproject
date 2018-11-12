library(splatter)
library(Oscope)
library(scater)
library(BiocParallel)


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



#for each simulation setting, we have two replicates
synData <- function(loc, scale, params, group_prob, save_DIR,
                    nrep = 2, Seed = 10,
                    ncells = 400, de_prob = 0.1){
  
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
  #number of groups
  K = length(group_prob)
  
  for(i in 1:nrep){
      ###get index of DD genes
      tmp = c()
      for(i in 1:K){
          a = paste0("rowData(dat)$DEFacGroup",i)
          tmp = union(tmp,which(eval(parse(text = a)) != 1))
      }
      
      DD = tmp
      ED = setdiff(1:nrow(data_counts),DD)
      
      ##get simulated counts and group label
      data_counts_all = assays(sim[[i]])$counts
      
      group = colData(sim[[i]])$Group
      
      ###make group1 belong to condition 1 and group 7 belong to condition2
      g11 = which(group == 1)
      
      g22 = which(group == 7)
      
      rest = setdiff(1:400,c(g11,g22))
      
      gp1 = sample(rest, 100 - length(g11))
      
      gp2 = setdiff(rest, gp1)
      
      g1 = c(g11,gp1)
      
      g2 = c(g22,gp2)
      
      label1 = c(group[g11],group[gp1])
      
      label2 = c(group[g22],group[gp2])
      
      data_counts = data_counts_all[,c(g1,g2)]
      cd = c(rep(1,200),rep(2,200))
      
      
      
      if(Loc < 0){
        name_vec = c("sim","_","neg",abs(loc),'_',scale,'_',i,".rds")
        filename = paste(name_vec,collapse = '')
      }
      else{
        name_vec = c("sim","_",abs(loc),'_',scale,'_',i,".rds")
        filename = paste(name_vec,collapse = '')
      }
      saveDir = paste0(save_DIR, filename)
      saveRDS(sim[[i]],data_counts,cd,g1,g2,group,DD,ED,label1,label2, file = saveDir)
  }
}
