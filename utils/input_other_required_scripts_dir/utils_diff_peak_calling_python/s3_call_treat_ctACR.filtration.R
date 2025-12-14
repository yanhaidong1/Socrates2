##this script is to perform the PerM for the treat and control peaks to perform the filtration

library(Matrix)
library(gtools)
library(edgeR)
library(preprocessCore)

# args
args <- commandArgs(T)

if(length(args)==3){
  
  #a <- readRDS(as.character(args[1]))
  ipt_peak_sparse_fl <- as.character(args[1])
  a <- read.table(ipt_peak_sparse_fl,stringsAsFactors = TRUE) 
  
  ##b is the group meta temp_add_group_meta.txt
  b <- read.delim(as.character(args[2]),row.names = 1)
  ##make a group cell type column
  b$group_celltype <- paste0(b$cell_identity,'.',b$group)
  target_cluster <- 'group_celltype'
  target_col <- target_cluster
  
  ##perform the filtration for the wrong cell type
  b$count <- 1
  aggregated_dt <- aggregate(count ~ cell_identity + group,data = b,sum)
  all_groups <- unique(aggregated_dt$group)
  cell_ids <- unique(aggregated_dt$cell_identity)
  missing_cell_types <- cell_ids[
    sapply(cell_ids, function(cid) {
      groups_present <- aggregated_dt$group[aggregated_dt$cell_identity == cid]
      length(setdiff(all_groups, groups_present)) > 0
    })
  ]
  b <- b[!(b$cell_identity %in% missing_cell_types),]

  
  
  output_dir <- as.character(args[3])
  
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
  message(" - iterate over all peaks for cluster ", i)
  ids <- rownames(subset(b, b[[target_col]]==i))
  aa <- a[,colnames(a) %in% ids]
  mat[,it] <- Matrix::rowSums(aa)
}

colnames(mat) <- clusts
rownames(mat) <- rownames(a)

write.table(mat, file=paste0(output_dir,'/opt_rawcounts_peaks_accessibility_',target_col,".txt"), quote=F, row.names=T, col.names=T, sep="\t")

##updating 121325
##obtain the final peaks
threshold <- 3
df <- as.data.frame(mat)
result <- do.call(
  rbind,
  lapply(colnames(df), function(col) {
    keep_rows <- rownames(df)[df[[col]] >= threshold]
    data.frame(
      column = col,
      row_name = keep_rows,
      stringsAsFactors = FALSE
    )
  })
)

write.table(result, file=paste0(output_dir,'/opt_',target_col,"filteration.peak.list.txt"), quote=F, row.names=F, col.names=F, sep="\t")




# adjust to per M
#options(mc.cores = 1)
#options("preprocessCore.nthreads" = 1)
mat <- as.matrix(mat)
mat <- cpm(mat, log=T, prior.count=5) #<- mat %*% diag(1e6/colSums(mat))
#mat <- normalize.quantiles(mat)
rownames(mat) <- rownames(a)
colnames(mat) <- clusts

write.table(mat, file=paste0(output_dir,'/opt_perM_peaks_accessibility_',target_col,".txt"), quote=F, row.names=T, col.names=T, sep="\t")









##Test
#df <- data.frame(
#  A = c(1, 4, 5, 2),
#  B = c(3, 1, 6, 4),
#  C = c(2, 7, 1, 5),
#  row.names = c("gene1", "gene2", "gene3", "gene4")
#)

#threshold <- 3

#result <- do.call(
#  rbind,
#  lapply(colnames(df), function(col) {
#    keep_rows <- rownames(df)[df[[col]] >= threshold]
#    data.frame(
#      column = col,
#      row_name = keep_rows,
#      stringsAsFactors = FALSE
#    )
#  })
#)












