# load libraries
library(Matrix)

# args
args <- commandArgs(T)

ipt_pmat_fl <- as.character(args[1])
ipt_gmat_fl <- as.character(args[2])
ipt_output_dir <- as.character(args[3])


ipt_pmat_dt <- read.table(ipt_pmat_fl)
ipt_gmat_dt <- read.table(ipt_gmat_fl)

colnames(ipt_pmat_dt) <- tolower(colnames(ipt_pmat_dt))
colnames(ipt_gmat_dt) <- tolower(colnames(ipt_gmat_dt))

shared_celltypes <- intersect(colnames(ipt_pmat_dt),colnames(ipt_gmat_dt))

# Remove all case-insensitive "unknown"
shared_celltypes <- shared_celltypes[tolower(shared_celltypes) != "unknown"]

ipt_pmat_dt <- ipt_pmat_dt[,shared_celltypes]
ipt_gmat_dt <- ipt_gmat_dt[,shared_celltypes]


write.table(ipt_pmat_dt, file=paste0(ipt_output_dir,'/opt_pmat.txt'), quote=F, row.names=T, col.names=T, sep="\t")
write.table(ipt_gmat_dt, file=paste0(ipt_output_dir,'/opt_gmat.txt'), quote=F, row.names=T, col.names=T, sep="\t")



