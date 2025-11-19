# libraries
args <- commandArgs(T)
# vars
ipt_findcell_obj_fl  <- as.character(args[1]) 
output_dir <- as.character(args[2])
prefix <- as.character(args[3])

obj <- readRDS(ipt_findcell_obj_fl)

peak_fl = obj$acr

write.table(peak_fl,paste0(output_dir,'/temp_',prefix,'_acr_peak.txt'),quote = F,sep = '\t',row.names = FALSE,col.names = FALSE)

