# run eFDR for genome


# load data
args <- commandArgs(T)
ipt_soc_obj_fl <- as.character(args[1])
output_dir <- as.character(args[2])

ipt_soc_obj <- readRDS(ipt_soc_obj_fl)

opt_meta_dt <- ipt_soc_obj$final_meta
write.table(opt_meta_dt,paste0(output_dir,'/','final_meta.txt'),sep = '\t',quote = F)

opt_peak_dt <- ipt_soc_obj$final_peak
write.table(opt_peak_dt,paste0(output_dir,'/final_peak.txt'),col.names = FALSE,row.names = FALSE,quote = F,sep = '\t')


