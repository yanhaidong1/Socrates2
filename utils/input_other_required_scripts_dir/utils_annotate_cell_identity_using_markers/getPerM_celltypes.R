# load libraries
library(Matrix)
library(gtools)
library(edgeR)
library(preprocessCore)

# args
args <- commandArgs(T)

if(length(args)==5){
  a_obj <- readRDS(as.character(args[1]))
  a <- a_obj$gene_Tn5
  #a <- read.table(as.character(args[1]),stringsAsFactors = T)
  b <- read.delim(as.character(args[2]),row.names = 1)
  #b <- read.table(as.character(args[2]))
  target_col <- as.character(args[3])
  output_dir <- as.character(args[[4]])
  prefix <- as.character(args[5]) ##rna or atac
}else{
  message(" ## wrong number of arguments ##")
  stop("Rscript differential_cluster_accessibility.R [meta] [sparse] [target_clust] [output_dir] [prefix]")
}

# load data
message(" - loading data ...")

# format
a$V1 <- factor(a$V1)
a$V2 <- factor(a$V2)
a <- sparseMatrix(i=as.numeric(a$V1),j=as.numeric(a$V2),x=as.numeric(a$V3),dimnames=list(levels(a$V1),levels(a$V2)))

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
  message(" - iterate over all genes for cluster ", i)
  ids <- rownames(subset(b, b[[target_col]]==i))
  aa <- a[,colnames(a) %in% ids]
  mat[,it] <- Matrix::rowSums(aa)
}

colnames(mat) <- clusts
rownames(mat) <- rownames(a)

write.table(mat, file=paste0(output_dir,'/opt_',prefix,"_rawcounts_gene_accessibility_",target_col,".txt"), quote=F, row.names=T, col.names=T, sep="\t")

# adjust to per M
mat <- cpm(mat, log=T, prior.count=5) #<- mat %*% diag(1e6/colSums(mat))
mat <- normalize.quantiles(mat)
rownames(mat) <- rownames(a)
colnames(mat) <- clusts

write.table(mat, file=paste0(output_dir,'/opt_',prefix,"_perM_genes_accessibility_",target_col,".txt"), quote=F, row.names=T, col.names=T, sep="\t")







