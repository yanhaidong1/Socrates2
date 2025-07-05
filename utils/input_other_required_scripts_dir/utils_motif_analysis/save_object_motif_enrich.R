# run eFDR for genome
library(Matrix)

# load data
args <- commandArgs(T)
ipt_soc_obj_fl <- as.character(args[1])
ipt_enrich_res_fl <- as.character(args[2])
ipt_peak_tn5_fl <- as.character(args[3])
ipt_prefix <- as.character(args[4])
output_dir <- as.character(args[5])

ipt_soc_obj <- readRDS(ipt_soc_obj_fl)
a <- read.table(ipt_peak_tn5_fl,stringsAsFactors=T)
a <- sparseMatrix(i=as.numeric(a$V1),
                  j=as.numeric(a$V2),
                  x=as.numeric(a$V3),
                  dimnames=list(levels(a$V1), levels(a$V2)))

ipt_enrich_res_dt <- read.delim(ipt_enrich_res_fl)

final_obj <- append(ipt_soc_obj, list(
  peak_tn5_mtx = a,
  enrich_motif = ipt_enrich_res_dt
))

saveRDS(final_obj,file=paste0(output_dir,'/',ipt_prefix,'.atac.soc.rds'))