library(dplyr)
library(ggvenn)
library(vioplot)
library(sm)
library(RColorBrewer)
library(ggpubr)

###-----Update history------
#v4 20240105 check the quality of the links from v4 for nearby links only.
#v5 20240405 modified function "filter_distal_links2", about filter links with minor distal.
#            

#module load R/4.3.1-foss-2022a
#source("/scratch/xz24199/02Will82/scatac-seq/step09_ACR_gene_links/functions_check_links.R")
#source("./02_functions_check_links.R")

# load arguments ----------------------------------------------------------------------------------
args <- commandArgs(T)
#if(length(args)!=4){stop("Rscript Check_acr2gene_links.R <links.txt> <acr2gene.bed> <fdr> <prefix>")}

link <- as.character(args[1])
#gene_bed <- as.character(args[2])
acr2gene <- as.character(args[2])
#fdr <- as.numeric(args[3])
prefix <- as.character(args[3])
min.cor <- 0.25
output_dir <- as.character(args[4])
input_scripts_dir <- as.character(args[5])

ipt_fdr_cutoff <- as.numeric(args[6])
ipt_pval_cutoff <- as.numeric(args[7])


source(paste0(input_scripts_dir,'/',"02_functions_check_links.R"))

acr2gene.m <- read.table(acr2gene)
link.m.all <- read.table(link)
link.m.all <- left_join(link.m.all, acr2gene.m, by="peak_id")

#link.m <- link.m.all[link.m.all$FDR < fdr,]

##updating 080825
link.m <- link.m.all[link.m.all$p.value < ipt_pval_cutoff & link.m.all$FDR < ipt_fdr_cutoff & abs(link.m.all$cor.value) >= min.cor,]
link.m.ck <- link.m.all[link.m.all$p.value > 0.3 & abs(link.m.all$cor.value) < 0.05,]

#link.m <- add_distance(link.m, gene_bed)
#link.m <- filter_distal_links(link.m)
#link.m <- filter_distal_links1(link.m)
link.m$link.type <- ifelse(abs(link.m$pair_dis) < 260, "gLink", ifelse(abs(link.m$pair_dis) < 2000, "pLink", "dLink"))
link.m$distance.type <- link.m$link.type
link.m$distance <- abs(link.m$pair_dis)
link.m <- filter_distal_links2(link.m)
link.m <- link.m[link.m$Filter_distal == "keep",]
link.m <- check_cre_type(link.m)

write.table(link.m, file = paste0(output_dir,'/',prefix,".observed.gene-peak-links.annotated.txt"), sep = "\t", quote = F)

summary_link(output_dir,link.m=link.m, prefix = paste0(prefix))
