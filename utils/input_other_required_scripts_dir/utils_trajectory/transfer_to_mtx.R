###TopGO###
library(Matrix)

args <- commandArgs(T)

# arguments
input_sparse_fl <- as.character(args[1])

input_output_dir <- as.character(args[2])

prefix <- as.character(args[3])

activity <- read.table(input_sparse_fl,stringsAsFactors = T)

activity <- sparseMatrix(i=as.numeric(activity$V1),
                         j=as.numeric(activity$V2),
                         x=as.numeric(activity$V3),
                         dimnames=list(levels(activity$V1), levels(activity$V2)))
atac_mtx <- as(activity, "dgCMatrix")
saveRDS(atac_mtx,paste0(input_output_dir,'/',prefix,'.rds'))