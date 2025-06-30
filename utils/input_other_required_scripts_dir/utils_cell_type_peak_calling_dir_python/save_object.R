# run eFDR for genome


# load data
args <- commandArgs(T)
ipt_soc_obj_fl <- as.character(args[1])
ipt_peak_fl <- as.character(args[2])
output_dir <- as.character(args[3])
ipt_prefix <- as.character(args[4])

ipt_soc_obj <- readRDS(ipt_soc_obj_fl)
ipt_peak_dt <- read.table(ipt_peak_fl)

final_obj <- append(ipt_soc_obj, list(
  final_peak = ipt_peak_dt
))

saveRDS(final_obj,file=paste0(output_dir,'/',ipt_prefix,'.atac.soc.rds'))