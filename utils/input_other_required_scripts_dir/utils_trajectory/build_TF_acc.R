###TopGO###
library(Matrix)

args <- commandArgs(T)

# arguments
input_gene_acc_obj_rds_fl <- as.character(args[1])

input_TF_fl <- as.character(args[2])

input_output_dir <- as.character(args[3])


ipt_gene_acc_mtx <- readRDS(input_gene_acc_obj_rds_fl)$gene_acc

ipt_TF_dt <- read.delim(input_TF_fl)

ipt_TF_genes <- unique(ipt_TF_dt$Gene_ID)

shared_ids <- intersect(rownames(ipt_gene_acc_mtx),ipt_TF_genes)

opt_TF_acc_mtx <- ipt_gene_acc_mtx[shared_ids,]



saveRDS(opt_TF_acc_mtx,paste0(input_output_dir,'/TF.acc.rds'))