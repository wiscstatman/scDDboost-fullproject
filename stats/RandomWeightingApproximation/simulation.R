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

root <- "/ua/xiuyu/simulation_data/data"

datasets <- read_tsv(file.path(root, "datasets.txt"),
                     col_types = cols(.default = col_character(),
                                      NumCells = col_integer()
                     )
)

real <- loadDataset(datasets[3, ], root)

set.seed(10)
real <- real[, sample(1:ncol(real), 200)]
real <- real[rowSums(real) > 0, ]

params <- splatEstimate(real)


bp <- BiocParallel::MulticoreParam(10)




sim.groups <- splatSimulate(group.prob = c(0.5, 0.5), method = "groups",
                            verbose = FALSE)
sim.groups <- normalise(sim.groups)

###updated splatter remove 
library("Oscope")
library("scater")
sims <- bplapply(1:20, function(seed) {
  message("Simulating ", seed)
  sim <- splatSimulateGroups(params,
                             batchCells      = 400,
                             group.prob      = c(0.3,0.2,0.1,0.15,0.05,0.05,0.15),
                             de.prob         = 0.1,
                             de.facLoc       = -0.3,
                             de.facScale     = 0.5,
                             # dropout.present = FALSE,
                             seed            = seed)
  sim <- calculateQCMetrics(sim)
  return(sim)
}, BPPARAM = bp)

###saving simulated data

for(i in 2:20){
  saveRDS(sims[[i]],  file = paste0(paste0("sim",i),".rds"))
}

library("MCMCpack")
##split simulated data into two conditions.
# group1_p = rdirichlet(1, rep(1,7))
# group2_p = rdirichlet(1, 1:7)
# split_p = group1_p / (group1_p + group2_p)

###as splatter can not control how genes are differentially expressed between groups

scd = readRDS(paste0(paste0("/ua/xiuyu/simulation_data/sim",i),".rds"))
#BatchCellMeans BaseCellMeans BCV CellMeans TrueCounts counts
ASSAY = assays(scd)
counts = ASSAY[[6]]
g1 = sample(1:400,200)
g2 = setdiff(1:400,g1)
clus = colData(scd)$Group
clusadj = sapply(clus,function(x) as.numeric(substr(x,6,6)))

c1 = clusadj[g1]
c2 = clusadj[g2]

dat1 = counts[,c1]
dat2 = counts[,c2]
dat = cbind(dat1,dat2)
rm(dat1)
rm(dat2)
cd = c(rep(1,200),rep(2,200))

D_c = cal_D(dat, 10)
K = 7

sz = MedianNorm(dat)
hp = c(1, rep(1,nrow(dat))) ##get hyper parameter

pDD7 = PDD(data = dat, cd = cd, ncores = 10, K = K, D = D_c,
           sz = sz, hp, pat(K)[[1]], 10, random = T, lambda = 0.5, nrandom = 50)

EDDb = which(pDD7 > 0.95)
length(EDDb)


####run deseq
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(monocle))

run_DESeq2census <- function(L) {
  message("DESeq2census")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      cds <- newCellDataSet(L$tpm, 
                            phenoData = new("AnnotatedDataFrame", 
                                            data = data.frame(condition = L$condt, 
                                                              row.names = colnames(L$tpm))))
      censuscounts <- relative2abs(cds)
      dds <- DESeqDataSetFromMatrix(countData = round(censuscounts), 
                                    colData = data.frame(condition = L$condt), 
                                    design = ~condition)
      dds <- DESeq(dds)
      res <- results(dds, contrast = c("condition", levels(factor(L$condt))[1], 
                                       levels(factor(L$condt))[2]), alpha = 0.05)
    })
    
    plotDispEsts(dds)
    plotMA(res)
    summary(res)
    
    list(session_info = session_info,
         timing = timing,
         res = res,
         df = data.frame(pval = res$pval,
                         padj = res$padj,
                         row.names = rownames(res)))
  }, error = function(e) {
    "DESeq2 results could not be calculated"
    list(session_info = session_info)
  })
}







