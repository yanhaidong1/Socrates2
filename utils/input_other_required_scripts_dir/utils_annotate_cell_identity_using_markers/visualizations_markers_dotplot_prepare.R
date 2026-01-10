###################################################################################################
##updating 011026 we will set the final_meta if we use the cell_identity column
##updating 010525 we will use the object
#' Visualizing marker genes using dot plot

#' Function to plot dot plot
#' @import Matrix
#' @import ggplot2

library(Matrix)
library(ggplot2)
library(gtools)

args <- commandArgs(T)

input_object_fl <- as.character(args[1])

#input_gene_body_acc_mtx_rds_fl <- as.character(args[1])

#input_gene_body_sparse_fl <- as.character(args[2])

#input_meta_fl <- as.character(args[3])

input_target_clust <- as.character(args[2])

input_output_dir <- as.character(args[3])

input_prefix <- as.character(args[4])


prepare_geneImputAcc_prop     <- function(input_object,
                                          input_output_dir,target_col){
  
  ####
  ##s1 load the files
  message(" -s1 load files ...")
  #acts <- readRDS(input_imputed_mtx_rds)
  acts <- input_object$gene_acc_smooth
  
  
  #meta_dt <- read.delim(input_meta_fl,row.names = 1)
  ##updating 011026
  if (target_col != 'LouvainClusters'){
    meta_dt <- input_object$final_meta
  }else{
    meta_dt <- input_object$meta
  }
  
  
  ####
  ##s2 calculat the avg mean values
  # verbose start-up
  message(" -s2 cal avg mean ...")
  if(is.null(meta_dt) | is.null(acts)){
    stop("ERROR: must provide metadata and activity matrix")
  }
  if(! "Cluster" %in% colnames(meta_dt)){
    meta_dt$Cluster <- meta_dt[[target_col]]
  }
  #rownames(markers) <- markers$geneID
  
  ##correspond cells to the meta
  meta_dt <- meta_dt[colnames(acts),]
  clust.o <- mixedsort(unique(meta_dt$Cluster))
  meta_dt$Cluster <- factor(meta_dt$Cluster, levels=clust.o)
  
  ##iterate over genes
  message(" - transversing gene activity ...")
  it <- 0
  df.acts <- as.data.frame(t(as.matrix(acts)))
  ##split cells based on the cluster
  df.splits <- split(df.acts, meta_dt$Cluster)
  ##calcualte the mean across all cells per gene
  difs.mean <- lapply(df.splits, function(z){colMeans(z)})
  difs.mean <- as.matrix(do.call(cbind, difs.mean))
  c.means <- difs.mean[,mixedorder(colnames(difs.mean))]
  
  ##return z score
  z <- t(as.matrix(scale(t(as.matrix(c.means)))))
  
  ####
  ##s3 calculate the prop per gene
  message(" -s3 cal prop per gene ...")
  #ipt_sparse_dt <- read.table(input_sparse_fl,stringsAsFactors = T)
  ipt_sparse_dt <- input_object$gene_Tn5
  
  ##correspond the gene name
  ipt_sparse_dt_t <- ipt_sparse_dt[ipt_sparse_dt$V1 %in% rownames(z),]
  ##correspond the cell name
  ipt_sparse_dt_t <- ipt_sparse_dt_t[ipt_sparse_dt_t$V2 %in% rownames(meta_dt),]
  
  if( "cellID.1" %in% colnames(meta_dt)){
    meta_dt$cellID <- meta_dt$cellID.1
  }
  
  ##updating 101524 debug add the cellID as some meta files do not have the cellID
  if (! 'cellID' %in% colnames(meta_dt)){
    meta_dt$cellID <- rownames(meta_dt)
    
  }
  
  
  ##merge meta and sparse file
  merged_dt <- merge(ipt_sparse_dt_t,meta_dt,by.x='V2',by.y = 'cellID') ##we will not use TRUE
  
  cell_type_list = as.character(unique(merged_dt[[target_col]]))
  
  outs <- lapply(cell_type_list,function(x){
    
    ##calculate total number of cells
    target_ct_cellnum <- length(meta_dt[meta_dt[[target_col]] == x,]$cellID)
    
    ##calculate the number of cells per genes
    merged_dt_target_ct <- merged_dt[merged_dt[[target_col]] == x,]
    
    merged_dt_target_ct$number <- 1
    
    gene_cellnum_dt <- aggregate(number ~ V1,data = merged_dt_target_ct, sum)
    
    gene_cellnum_dt$prop <- gene_cellnum_dt$number/target_ct_cellnum
    gene_cellnum_dt$organct <- x
    colnames(gene_cellnum_dt) <- c('Gene','AccNum','AccProp','CellType')
    
    return(gene_cellnum_dt)
  })
  
  combine_dt <- do.call(rbind,outs)
  
  ####
  ##s4 merge with z score
  message(" -s4 cal prop per gene ...")
  outs2 <- lapply(seq(1,nrow(combine_dt)), function(x){
    
    target_line <- combine_dt[x,]
    target_gene <- as.character(target_line$Gene)
    target_ctorgan <- as.character(target_line$CellType)
    
    #print(target_gene)
    #print(target_ctorgan)
    target_zscore <- z[target_gene,target_ctorgan]
    target_line$zscore <- target_zscore
    
    return(target_line)
    
  })
  
  ##return combine_dt_addzscore
  combine_dt_addzscore <- do.call(rbind,outs2)
  
  ##combine the information of current name
  ##return combine_dt_addzscore_addcomm
  #combine_dt_addzscore_addcomm <- merge(combine_dt_addzscore,target_gene_dt,by.x = 'Gene',by.y = 'geneID')
  
  
  return(list(zscore=z, dt=combine_dt_addzscore, cmeans = c.means))
  
}



##prepare the object 
#input_imputed_mtx_rds <- input_gene_body_acc_mtx_rds_fl
#input_sparse_fl <- input_gene_body_sparse_fl

input_object <- readRDS(input_object_fl)

opt_prepare_plot_obj <- prepare_geneImputAcc_prop(input_object,
                         input_output_dir,input_target_clust)


##udpating 010525 we will add the output to the oject
final_obj <- append(input_object, list(
  dotplot = opt_prepare_plot_obj
))

saveRDS(final_obj,file=paste0(input_output_dir,'/',input_prefix,'.atac.soc.rds'))


#saveRDS(opt_prepare_plot_obj,paste0(input_output_dir,'/opt_gene_zscore_accprop.rds'))







