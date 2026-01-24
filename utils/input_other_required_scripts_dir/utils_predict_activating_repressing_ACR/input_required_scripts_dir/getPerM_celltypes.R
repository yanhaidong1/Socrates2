# load libraries
library(Matrix)
library(gtools)
library(edgeR)
library(preprocessCore)
library(Seurat)

# args
args <- commandArgs(T)

if(length(args)==3){

  ipt_object_fl <- as.character(args[1])
  ipt_object <- readRDS(ipt_object_fl)
  a <- GetAssayData(ipt_object, slot = "counts")
  b <- ipt_object@meta.data

  #a <- read.table(ipt_peak_sparse_fl,stringsAsFactors = TRUE) 
  #a <- a_obj$gene_Tn5
  #a <- read.table(as.character(args[1]),stringsAsFactors = T)
  #b <- read.delim(as.character(args[2]),row.names = 1)
  #b <- read.table(as.character(args[2]))
  #parameter_fl <- as.character(args[3])
  #target_col <- as.character(args[3])
  #source(parameter_fl)
  target_col <- as.character(args[2])
  
  output_dir <- as.character(args[3])
  
  #prefix <- as.character(args[5]) ##rna or atac
}else{
  message(" ## wrong number of arguments ##")
  stop("Rscript differential_cluster_accessibility.R [meta] [sparse] [target_clust] [output_dir] [prefix]")
}

# load data
message(" - loading data ...")

# format
#a$V1 <- factor(a$V1)
#a$V2 <- factor(a$V2)
#a <- sparseMatrix(i=as.numeric(a$V1),j=as.numeric(a$V2),x=as.numeric(a$V3),dimnames=list(levels(a$V1),levels(a$V2)))

# iterate over clusters
clusts <- mixedsort(unique(b[[target_col]]))
mat <- matrix(0,nrow=nrow(a),ncol=length(clusts))

##get the intersect id
intersect_id <- intersect(colnames(a),rownames(b))
a <- a[,intersect_id]
b <- b[intersect_id,]

it <- 0
for (i in clusts){
  it <- it+1
  message(" - iterate over all peaks for cluster ", i)
  ids <- rownames(subset(b, b[[target_col]]==i))
  aa <- a[,colnames(a) %in% ids]
  mat[,it] <- Matrix::rowSums(aa)
}

colnames(mat) <- clusts
rownames(mat) <- rownames(a)

write.table(mat, file=paste0(output_dir,'/opt_rawcounts_peaks_accessibility_',target_col,".txt"), quote=F, row.names=T, col.names=T, sep="\t")

# adjust to per M
#options(mc.cores = 1)
#options("preprocessCore.nthreads" = 1)
mat <- as.matrix(mat)
mat <- cpm(mat, log=T, prior.count=5) #<- mat %*% diag(1e6/colSums(mat))
#mat <- normalize.quantiles(mat)
rownames(mat) <- rownames(a)
colnames(mat) <- clusts

write.table(mat, file=paste0(output_dir,'/opt_perM_peaks_accessibility_',target_col,".txt"), quote=F, row.names=T, col.names=T, sep="\t")







