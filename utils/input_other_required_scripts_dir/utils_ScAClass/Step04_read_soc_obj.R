##this script is to conduct the feature embeddings
##updating 063024

args = commandArgs(T)

input_soc_object_fl <- as.character(args[1])

input_output_dir <- as.character(args[2])

input_soc_object <- readRDS(input_soc_object_fl)

ipt_meta_dt <- input_soc_object$meta

write.table(ipt_meta_dt,paste0(input_output_dir,'/temp_meta_fl.txt'),quote = F,sep = '\t')



