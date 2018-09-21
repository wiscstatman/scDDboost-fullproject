root <- "../data"

datasets <- read_tsv(file.path(root, "datasets.txt"),
                     col_types = cols(.default = col_character(),
                                      NumCells = col_integer()
                                      )
                     )

real <- loadDataset(datasets[3, ], root)

set.seed(1)
real <- real[, sample(1:ncol(real), 200)]
real <- real[rowSums(real) > 0, ]

