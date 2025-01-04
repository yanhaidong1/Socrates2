##this script is to store the all the prepared files to the previous soc object

library(Matrix)

# load rds files and pre-processing --------------------------------------
##create the Socrates obj
args <- commandArgs(T)

input_soc_object_fl <- as.character(args[1])

input_gene_tn5_sparse_fl <- as.character(args[2])

input_gene_500bp_sorted_fl <- as.character(args[3])

input_prefix <- as.character(args[4])

input_output_dir <- as.character(args[5])

##load the object
ipt_soc_object <- readRDS(input_soc_object_fl)

##load the tn5 gene sparse file
ipt_tn5_gene_sparse_dt <- read.table(input_gene_tn5_sparse_fl,stringsAsFactors = T)

##load the gene 500bp sorted fl
ipt_gene_500bp_bed_dt <- read.table(input_gene_500bp_sorted_fl)

# Adding a new item
final_obj <- append(ipt_soc_object, list(
  
  gene_Tn5 = ipt_tn5_gene_sparse_dt,
  gene_bed = ipt_gene_500bp_bed_dt
  
))

saveRDS(final_obj,file=paste0(input_output_dir,'/',input_prefix,'.atac.soc.rds'))

