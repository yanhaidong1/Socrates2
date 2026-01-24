###################################################################################################
##               Predict ACR function by correlation ACR with nearby gene expression            ##
###################################################################################################
#
#20240402 v1. 


# load libraries ----------------------------------------------------------------------------------
#module load R/4.3.1-foss-2022a
library(Matrix)
library(parallel)
library(ggplot2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
#library(ArchR)
library(edgeR)
#library(GenomicRanges)
library(foreach)
library(doParallel)
library(dplyr)


#source("/scratch/xz24199/02Will82/scatac-seq/step09_ACR_gene_links/Liger_v4_MaxK80_cpm0.1/function_Liger_linkGenesAndPeak.R")
#source("./01_function_Liger_linkGenesAndPeak.R")

# #load arguments ----------------------------------------------------------------------------------
 args <- commandArgs(T)
#if(length(args)!=6){stop("Rscript ACR2gene_Liger_v3.2.R <scRNA_seurate.obj.rds> <scATAC.soc.processed.rds> <scatac.gene.sparse.rds> <peak_sparse.rds> <atac.meta> <prefix>")}
acr_pairs <- as.character(args[1])
cores <- as.numeric(args[2])
pmat <- as.character(args[3])
gmat <- as.character(args[4])
prefix <- as.character(args[5])
##updating 080725
opt_dir <- as.character(args[6])
ipt_required_scripts_dir <- as.character(args[7])

source(paste0(ipt_required_scripts_dir,"/01_function_Liger_linkGenesAndPeak.R"))

# #test input
# acr_pairs <- read.table("/scratch/xz24199/02Will82/scatac-seq/step06_peaks_sparse/03_annotate_peaks/Gm_atlas_all_tissues_cpm4.ACR.sorted_ACR2gene_pairs.bed")
# #cellNum <- 100 #The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses, eg number of k groups is TotalCellNum/cellNum
# #max_k <- 80
# cores=30
# pmat <- read.table("/scratch/xz24199/02Will82/scatac-seq/step09_ACR_gene_links/ACR_gene_cor/01_data/Gm_atlas_All_tissues_celltype_ACR_4cpm_celltype_cpm_matching.txt")
# gmat <- read.table("/scratch/xz24199/02Will82/scatac-seq/step09_ACR_gene_links/ACR_gene_cor/01_data/Gm_atlas_All_tissues_celltype_RNA_gene_celltype_cpm_matching.txt")
# prefix <- "Gm_atlas_All_tissues_celltype"

#test input
acr_pairs <- read.table(acr_pairs)
#cellNum <- 100 #The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses, eg number of k groups is TotalCellNum/cellNum
#max_k <- 80
#cores=30
pmat <- read.table(pmat)
gmat <- read.table(gmat)
#prefix <- "Gm_atlas_All_tissues_celltype"


#---process acr and gene pairs---
acr_pairs$peak_id <- paste(acr_pairs$V1,acr_pairs$V2,acr_pairs$V3, sep = "_")
acr_pairs$gene_id <- gsub(".Wm82.a4.v1","",x=acr_pairs$V8)
acr_pairs$pair_id <- paste(acr_pairs$peak_id, acr_pairs$gene_id, sep = ":")
acr_pairs$pair_dis <- acr_pairs$V11
acr_pairs$gene_pos <- paste(acr_pairs$V5,acr_pairs$V6,acr_pairs$V7, acr_pairs$V10, sep = "_")

#table overlap
overlap.g <- intersect(rownames(gmat), acr_pairs$gene_id)
gmat <- gmat[rownames(gmat) %in% overlap.g,]
acr_pairs <- acr_pairs[acr_pairs$gene_id %in% overlap.g,]
overlap.p <- intersect(rownames(pmat), acr_pairs$peak_id)
pmat <- pmat[rownames(pmat) %in% overlap.p,]
acr_pairs <- acr_pairs[acr_pairs$peak_id %in% overlap.p,]

acr_pairs <- acr_pairs[,c("gene_pos","pair_id","pair_dis")]

#gmat <- as(as.matrix(gmat), "sparseMatrix")
#pmat <- as(as.matrix(pmat), "sparseMatrix")

##save a object to store the gmat pmat and acr_pairs
final_obj <- list(
  gmat_mtx = gmat,
  pmat_mtx = pmat,
  acr_pairs_dt = acr_pairs
)
saveRDS(final_obj,file=paste0(opt_dir,'/','temp','.ipt.rds'))


# scan for gene-peak links
message(" - getting gene-peak linkages ...")
regnet <- linkGenesAndPeaks4(gene_counts = gmat,
                            peak_counts = pmat,
                            dist = 'spearman',
                            peak_genes_pairs = acr_pairs$pair_id,
							              cores=cores)
#merge all distance information.
regnet$pair_id <- paste(regnet$peak_id, regnet$gene_id, sep = ":")
regnet <- left_join(regnet, acr_pairs, by = "pair_id")

#regnet_p05 <- regnet[regnet$p.value < 0.05,]
#regnet_fdr05 <- regnet[regnet$FDR < 0.05,]
write.table(regnet, file=paste0(opt_dir,'/',prefix,"_ACR2gene_pairs",".observed.gene-peak-links.txt"), quote=F, sep="\t")
#write.table(regnet_p05, file=paste0(prefix,"_ACR2gene_pairs",".p05.observed.gene-peak-links.txt"), quote=F, sep="\t")
#write.table(regnet_fdr05, file=paste0(prefix,"_ACR2gene_pairs",".fdr05.observed.gene-peak-links.txt"), quote=F, sep="\t")



