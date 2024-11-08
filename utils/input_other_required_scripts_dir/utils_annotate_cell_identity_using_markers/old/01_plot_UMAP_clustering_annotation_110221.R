##updating 122521 do not filter cells
##updating 110221 we will change the b$unique to b$total
##udpating 110221 we will set the UMAP performance
##updating 081821 this is for the All.meta.txt from the annotation work
##updating 063021 plot the UMAP for the rice
##updating 021921 add the color information to the meta
##updation 121620 add intersection analysis
##for the mclust it does not matter how many of them will not be write out to a plot 
##since no filtration for the low number of cluster
##updation 113020 add the harmony filtration
##updation 112220 do some cluster size filtration for the plotting

################################################################
##develope a pipeline to check the clustering for the replicates
library(Matrix)
library(irlba)
#library(harmony)
library(uwot)
library(tcltk)
library(iterators)
library(itertools)
library(parallel)
library(Seurat)
library(FNN)
library(RColorBrewer)
library(scales)
library(viridis)
library(mclust)
library(gplots)
library(glmnet)
library(dplyr)
library(gtools)
library(densityClust)
library(dbscan)

suppressMessages(library(DESeq2))
suppressMessages(library(BiocParallel))
suppressMessages(library(caret))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(reshape2))

##updation we need to plot the replicate information
args <- commandArgs(T)

b <- read.delim(as.character(args[1])) ##ipt_all_meta_fl
prefix <- as.character(args[2]) ##cell_annotation_smooth
column <- as.character(args[3]) ##cell_annotation_smooth
target_cluster <- as.character(args[4]) ##cell_annotation_smooth

output_dir <- as.character(args[5])

#b <- read.table('All.LouvainCluster_FRiP03.txt',stringsAsFactors = T)
#head(b)

##check the replicate 
#column = 'LouvainClusters'
#cols <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(b[,column])))
#colv <- cols[factor(b[,column])]

#b <- read.table('PI30.LouvainCluster.txt')
#head(b)
#b <- read.table('All.LouvainCluster_aftrhm.txt',stringsAsFactors = T)

##this script will help to plot the UMAP cluster
##the mclust_size only not plot the cell number < 200 but do not conduct the filtration for the meta data
plotUMAP      <- function(b,newmeta_nm,output_dir, prefix="out", column="LouvainClusters_afthm", m1="umap1", m2="umap2",
                          target_cluster = 'LouvainClusters',mclust_size = 50,min.reads= 0.5e6){
  
  
  if (target_cluster == 'LouvainClusters_afthm') {
  
    ##important change the LouvainClusters to be factor
    #b$LouvainClusters <- factor(b$LouvainClusters)
    b$LouvainClusters_afthm <- factor(b$LouvainClusters_afthm)
    
    ##udpation 113020 chagne the LouvainClusters to the LouvainClusters_afthm
    #updation 112220 filter some clusters with size
    agg.reads <- aggregate(b$unique~b$LouvainClusters_afthm, FUN=sum)
    colnames(agg.reads) <- c("clusters","readDepth")
    clust.cnts <- table(b$LouvainClusters_afthm)
    agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    agg.pass <- subset(agg.reads, agg.reads$num_cells>=mclust_size & agg.reads$readDepth>=min.reads)
    b_filt <- b[b$LouvainClusters_afthm %in% agg.pass$clusters,]
  }
  if (target_cluster == 'LouvainClusters') {
    
    ##important change the LouvainClusters to be factor
    b$LouvainClusters <- factor(b$LouvainClusters)
    #b$LouvainClusters_afthm <- factor(b$LouvainClusters_afthm)
    
    ##udpation 113020 chagne the LouvainClusters to the LouvainClusters_afthm
    #updation 112220 filter some clusters with size
    agg.reads <- aggregate(b$unique~b$LouvainClusters, FUN=sum)
    colnames(agg.reads) <- c("clusters","readDepth")
    clust.cnts <- table(b$LouvainClusters)
    agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    agg.pass <- subset(agg.reads, agg.reads$num_cells>=mclust_size & agg.reads$readDepth>=min.reads)
    b_filt <- b[b$LouvainClusters %in% agg.pass$clusters,]
  }
  
  if (target_cluster == 'LouvainClusters_new') {
    
    ##important change the LouvainClusters_new to be factor
    b$LouvainClusters_new <- factor(b$LouvainClusters_new)
    #b$LouvainClusters_afthm <- factor(b$LouvainClusters_afthm)
    
    ##udpation 113020 chagne the LouvainClusters_new to the LouvainClusters_afthm
    #updation 112220 filter some clusters with size
    agg.reads <- aggregate(b$unique~b$LouvainClusters_new, FUN=sum)
    colnames(agg.reads) <- c("clusters","readDepth")
    clust.cnts <- table(b$LouvainClusters_new)
    agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    agg.pass <- subset(agg.reads, agg.reads$num_cells>=mclust_size & agg.reads$readDepth>=min.reads)
    b_filt <- b[b$LouvainClusters_new %in% agg.pass$clusters,]
  }
  
  if (target_cluster == 'LouvainClusters_comb') {
    
    ##important change the LouvainClusters_comb to be factor
    b$LouvainClusters_comb <- factor(b$LouvainClusters_comb)
    #b$LouvainClusters_afthm <- factor(b$LouvainClusters_afthm)
    
    ##udpation 113020 chagne the LouvainClusters_comb to the LouvainClusters_afthm
    #updation 112220 filter some clusters with size
    agg.reads <- aggregate(b$unique~b$LouvainClusters_comb, FUN=sum)
    colnames(agg.reads) <- c("clusters","readDepth")
    clust.cnts <- table(b$LouvainClusters_comb)
    agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    agg.pass <- subset(agg.reads, agg.reads$num_cells>=mclust_size & agg.reads$readDepth>=min.reads)
    b_filt <- b[b$LouvainClusters_comb %in% agg.pass$clusters,]
  }
  
  if (target_cluster == 'cell_annotation_smooth') {
    
    ##important change the cell_annotation_smooth to be factor
    b$cell_annotation_smooth <- factor(b$cell_annotation_smooth)
    #b$LouvainClusters_afthm <- factor(b$LouvainClusters_afthm)
    
    b_filt <- b
    
    ##updating 122521 no fitlration
    ##udpation 113020 chagne the cell_annotation_smooth to the LouvainClusters_afthm
    #updation 112220 filter some clusters with size
    
    #agg.reads <- aggregate(b$total~b$cell_annotation_smooth, FUN=sum)
    #colnames(agg.reads) <- c("clusters","readDepth")
    #clust.cnts <- table(b$cell_annotation_smooth)
    #agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    #agg.pass <- subset(agg.reads, agg.reads$num_cells>=mclust_size & agg.reads$readDepth>=min.reads)
    #b_filt <- b[b$cell_annotation_smooth %in% agg.pass$clusters,]
  }
  
  if (target_cluster == 'cell_annotation_knn') {
    
    ##important change the cell_annotation_knn to be factor
    b$cell_annotation_knn <- factor(b$cell_annotation_knn)
    #b$LouvainClusters_afthm <- factor(b$LouvainClusters_afthm)
    
    b_filt <- b
    
    ##udpation 113020 chagne the cell_annotation_knn to the LouvainClusters_afthm
    #updation 112220 filter some clusters with size
    
    
    #agg.reads <- aggregate(b$total~b$cell_annotation_knn, FUN=sum)
    #colnames(agg.reads) <- c("clusters","readDepth")
    #clust.cnts <- table(b$cell_annotation_knn)
    #agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    #agg.pass <- subset(agg.reads, agg.reads$num_cells>=mclust_size & agg.reads$readDepth>=min.reads)
    #b_filt <- b[b$cell_annotation_knn %in% agg.pass$clusters,]
  }
  
  if (target_cluster == 'cell_annotation_glmnet') {
    
    ##important change the cell_annotation_glmnet to be factor
    b$cell_annotation_glmnet <- factor(b$cell_annotation_glmnet)
    #b$LouvainClusters_afthm <- factor(b$LouvainClusters_afthm)
    
    b_filt <- b
    
    ##udpation 113020 chagne the cell_annotation_glmnet to the LouvainClusters_afthm
    #updation 112220 filter some clusters with size
    
    #agg.reads <- aggregate(b$total~b$cell_annotation_glmnet, FUN=sum)
    
    #colnames(agg.reads) <- c("clusters","readDepth")
    #clust.cnts <- table(b$cell_annotation_glmnet)
    #agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    #agg.pass <- subset(agg.reads, agg.reads$num_cells>=mclust_size & agg.reads$readDepth>=min.reads)
    #b_filt <- b[b$cell_annotation_glmnet %in% agg.pass$clusters,]
  }
  
  
  
  
  # plot space
  pdf(paste0(output_dir,'/',prefix,".UMAP_clusters.pdf"), width=6, height=6)
  
  # test if column is present
  if(!column %in% colnames(b_filt)){
    stop(" column header: ", column, ", is missing from input dataframe ...")
  }
  
  if(m1 != "umap1" | m2 != "umap2"){
    b_filt$umap1 <- b_filt[,m1]
    b_filt$umap2 <- b_filt[,m2]
  }
  
  # cols
  if (column == 'genotype') {
    b_filt <- b_filt[sample(nrow(b_filt)),] ##pre12 1:10
    cols <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(b_filt[,column])))
    colv <- cols[factor(b_filt[,column])]
  }else{
  
  if(is.factor(b_filt[,c(column)])){
    b_filt <- b_filt[sample(nrow(b_filt)),] ##pre22 1:12
    cols <- colorRampPalette(brewer.pal(12,"Paired")[1:12])(length(unique(b_filt[,column])))
    colv <- cols[factor(b_filt[,column])]
  }else if(is.character(b_filt[,column])){
    b_filt[,column] <- factor(b_filt[,column])
    b_filt <- b_filt[sample(nrow(b_filt)),]
    cols <- colorRampPalette(brewer.pal(24,"Paired")[1:10])(length(unique(b_filt[,column])))
    colv <- cols[factor(b_filt[,column])]
  }else if(is.numeric(b_filt[,column])){
    b_filt <- b_filt[order(b_filt[,column], decreasing=F),]
    cols <- viridis(100)
    colv <- cols[cut(b_filt[,column], breaks=101)]
  }
  }
  #unique(colv)

  
  # plot 
  ##pre cex is 0.5 pch=16 pre cex is 0.3
  p <- plot(b_filt$umap1, b_filt$umap2, pch=16, cex=0.3, col=colv,
       xlim=c(min(b_filt$umap1), max(b_filt$umap1)+(abs(max(b_filt$umap1))*0.5)))
  
  #print(p)
  
  if(is.factor(b_filt[,column])){
    legend("right", legend=sort(unique(b_filt[,column])),cex=0.5,
           fill=cols[sort(unique(b_filt[,column]))])
  }
  
  ##updating 021921 add the col information to the meta
  head(b_filt)
  b_filt$sub_cluster_color <- colv
  write.table(b_filt,paste0(output_dir,'/',newmeta_nm,'.txt'),sep = '\t',quote = F)
  
  
  #text(x='0', y='1', '1', adj=c(0,1))
  
  #?text
  #text()
  
  #?legend
  
  # device off
  dev.off()
  
}


##we need to set some parameters

##updating 081821
#b <- read.delim('All.meta.txt')

##smooth
#column='cluster_annotation_smooth'
#prefix='cluster_annotation_smooth'
m1="umap1"
m2="umap2"
mclust_size = 10
#target_cluster = 'cluster_annotation_smooth'
min.reads= 0.5e6


##smooth
#column='cluster_annotation_smooth'
#prefix='cluster_annotation_smooth'
#m1="umap1"
#m2="umap2"
#mclust_size = 10
#target_cluster = 'cluster_annotation_smooth'
#min.reads= 0.5e6

##knn
#column='cluster_annotation_knn'
#prefix='cluster_annotation_knn'
#m1="umap1"
#m2="umap2"
#mclust_size = 10
#target_cluster = 'cluster_annotation_knn'
#min.reads= 0.5e6

##knn
#column='cluster_annotation_glmnet'
#prefix='cluster_annotation_glmnet'
#m1="umap1"
#m2="umap2"
#mclust_size = 10
#target_cluster = 'cluster_annotation_glmnet'
#min.reads= 0.5e6

plotUMAP(b,'opt_final_meta_addcolr_cells',output_dir,column=column, prefix=prefix,m1="umap1", m2="umap2",mclust_size = 10,target_cluster = target_cluster)



