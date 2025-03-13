##updating 031225 this script is to annotate the cell identity and plot the final UMAP


################################################################
##develope a pipeline to check the clustering for the replicates

library(RColorBrewer)
#library(viridis)


args <- commandArgs(TRUE)

ipt_output_soc_obj_fl <- as.character(args[1])

ipt_cell_annot_fl <- as.character(args[2])

ipt_output_dir <- as.character(args[3])

ipt_plot_cex_val <- as.numeric(args[4])

##read the file
ipt_soc_obj <- readRDS(ipt_output_soc_obj_fl)

ipt_meta_dt <- ipt_soc_obj$meta

ipt_cell_annot_dt <- read.delim(ipt_cell_annot_fl,head=FALSE)

colnames(ipt_cell_annot_dt) <- c('cluster','cell_identity')

##merge the meta file
merged_dt <- merge(ipt_meta_dt,ipt_cell_annot_dt,by.x='LouvainClusters',by.y='cluster')
head(merged_dt)

table(merged_dt$cell_identity)

rownames(merged_dt) <- merged_dt$cellID

table(merged_dt$LouvainClusters)

##this script will help to plot the UMAP cluster
##the mclust_size only not plot the cell number < 200 but do not conduct the filtration for the meta data
plotUMAP <- function(b, prefix="out", column="LouvainClusters_afthm", m1="umap1", m2="umap2",
                          target_cluster = 'LouvainClusters',knowncolor_TCP = 'no',knowncolor_lib = 'no',knowncolor_tissue = 'no',
                          mclust_size = 50,min.reads= 0.5e6,
                          newmeta_nm='addcol',target_unknown_name='no',openlegend='no',plot_bar = 'yes',
                          plot_cex = 0.15){
  
  b[[target_cluster]] <- factor(b[[target_cluster]])
  b_filt <- b
 
  # plot space
  pdf(paste0(ipt_output_dir,'/',prefix,".UMAP_clusters.pdf"), width=6, height=6)
  
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
  
  ##updating 011222
  if(is.factor(b_filt[,c(column)])){
    
    #if (target_unknown_name != 'Lseed_unknown'){
    
    if (knowncolor_TCP == 'yes'){
      
      colv <- as.character(b_filt[,'color_TCP'])
     
      
    }else{
      
      
      if (knowncolor_lib == 'yes'){
        
        cols <- as.character(unique(b_filt[,'color_lib']))
        colv <- as.character(b_filt[,'color_lib'])
        
        
      }else{
          
          if (knowncolor_tissue == 'yes'){
            
            
            cols <- as.character(unique(b_filt[,'color_tissue']))
            colv <- as.character(b_filt[,'color_tissue'])
            
          }else{
          
            if (target_unknown_name != 'no'){
              b_filt <- b_filt[sample(nrow(b_filt)),] ##pre22 1:12
              cols <- colorRampPalette(brewer.pal(12,"Paired")[1:12])(length(unique(b_filt[,column])))
              #colv <- cols[factor(b_filt[,column])]
              
              sorted_cellname <- sort(unique(factor(b_filt[,column])))
              class(sorted_cellname)
              
              #target_unknown_name <- 'seedling_unknown'
              
              target_name_order <- match(target_unknown_name, sorted_cellname)
              target_name_order_color <- cols[target_name_order]
              cols <- gsub(target_name_order_color,'grey',cols)
              colv <- cols[factor(b_filt[,column])]
            
            }else {
          
              print('the cluster is the number')
              b_filt <- b_filt[sample(nrow(b_filt)),] ##pre22 1:12
              cols <- colorRampPalette(brewer.pal(12,"Paired")[1:12])(length(unique(b_filt[,column])))
              colv <- cols[factor(b_filt[,column])]
            }
            
          }
      }
    }
    
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
  ##choose open or not
  #p <- plot(b_filt$umap1, b_filt$umap2, pch=16, cex=0.3, col=colv,
  #     xlim=c(min(b_filt$umap1), max(b_filt$umap1)+(abs(max(b_filt$umap1))*0.5)))
  
  ##previous is 0.2
  ##article is 0.06
  ##for the two lib is 0.3
  ##for the merged organs is 0.06
  ##for the new scRNA scATAC integration cex0.5
  
  if (plot_bar == 'yes'){
    
    p <- plot(b_filt$umap1, b_filt$umap2, pch=16, cex=plot_cex, col=colv,xaxt='n',yaxt='n',xlab="",ylab="",
              xlim=c(min(b_filt$umap1), max(b_filt$umap1)+(abs(max(b_filt$umap1))*0.5)))
    
    print('Only print the legend bar')
    cols <- as.character(unique(b_filt[,'color_TCP']))
    legend("right", legend=as.character(unique(b_filt[,column])),cex=0.5,
           fill=cols)
    
  }else {
    
    ##previous is 0.5
    ##default is 0.15
    p <- plot(b_filt$umap1, b_filt$umap2, pch=16, cex=plot_cex, col=colv,xaxt='n',yaxt='n',xlab="",ylab="",
              xlim=c(min(b_filt$umap1), max(b_filt$umap1)+(abs(max(b_filt$umap1))*0.5)))
    p
    
    ##choose open or not
    if (openlegend == 'yes'){
      if(is.factor(b_filt[,column])){
        legend("bottomright", legend=sort(unique(b_filt[,column])),cex=0.5,
               fill=cols[sort(unique(b_filt[,column]))])
      }
    }
    
    ##updating 021921 add the col information to the meta
    head(b_filt)
    b_filt$sub_cluster_color <- colv
    #write.table(b_filt,paste0(newmeta_nm,'.txt'),sep = '\t',quote = F)
    
    
    #text(x='0', y='1', '1', adj=c(0,1))
    
    #?text
    #text()
    
    #?legend
    
    # device off
  }
  dev.off()
  
  
}


##updating 031225
##plot the ScAClass prediction cell
b <- merged_dt
#target_unknown_name <- 'Unknown'
head(b)
##plot_cex = 0.2
plotUMAP(b,column='cell_identity', prefix='annot',m1="umap1", m2="umap2",mclust_size = 1,
         target_cluster = 'cell_identity',knowncolor_TCP = 'no',knowncolor_lib = 'no',knowncolor_tissue = 'no',
         min.reads= 1,newmeta_nm='addcol',target_unknown_name='no',openlegend = 'yes',plot_bar = 'no',
         plot_cex = ipt_plot_cex_val)

