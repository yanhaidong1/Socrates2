# run eFDR for genome


# load data
args <- commandArgs(T)
ipt_soc_obj_fl <- as.character(args[1])
output_dir <- as.character(args[2])

ipt_soc_obj <- qs2::qs_read(ipt_soc_obj_fl)

opt_meta_dt <- ipt_soc_obj$meta
write.table(opt_meta_dt,paste0(output_dir,'/','temp_meta.txt'),sep = '\t',quote = F)



