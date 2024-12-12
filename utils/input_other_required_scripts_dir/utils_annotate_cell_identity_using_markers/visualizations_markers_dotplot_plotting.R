###################################################################################################
###################################################################################################
###################################################################################################

#' Visualizing marker genes using dot plot

#' Function to plot dot plot
#' @import Matrix
#' @import ggplot2

library(Matrix)
library(ggplot2)
library(gtools)


args <- commandArgs(T)

opt_prepare_plot_obj_fl <- as.character(args[1])

input_marker_gene_fl <- as.character(args[2])

input_output_dir <- as.character(args[3])

plot_width <- as.numeric(args[4])

plot_height <- as.numeric(args[5])


plot_markers <- function(prefix,input_marker_gene_fl,opt_prepare_plot_obj,
                         plot_width,plot_height,OpenaddCTtoName='no'){

  input_z_score_dt <- opt_prepare_plot_obj$zscore

  combine_dt_addzscore <- opt_prepare_plot_obj$dt
  
  markers <- read.delim(input_marker_gene_fl)
  
  ##merge the markers
  combine_dt_addzscore_addcomm <- merge(combine_dt_addzscore,markers,by.x = 'Gene',by.y = 'geneID')
  
  ##filter z score dt for the marker genes and top Var
  z.means <- input_z_score_dt[,mixedorder(colnames(input_z_score_dt))]
  z.means <- z.means[rownames(z.means) %in% rownames(markers),]
  markers <- markers[rownames(z.means),]
  
  ##do the clustering of z.means from markers
  v_clust <- hclust(dist(z.means %>% as.matrix() %>% t()))
  
  ##do the clustering of z.means_topVar
  #v_clust <- hclust(dist(z.means_topVar %>% as.matrix() %>% t()))
  cluster_ctorgan_list <- v_clust$labels[v_clust$order]
  z.means <- z.means[,cluster_ctorgan_list]
  dim(z.means)
  
  ##order the z.means
  row.o <- apply(z.means, 1, which.max)
  ##> head(row.o)
  #LOC_Os04g33860 LOC_Os01g41710 LOC_Os03g16430 LOC_Os05g37800 LOC_Os05g51510
  #1              2              2              2              2
  #LOC_Os09g17740
  #2
  ##after finding the order of gene with the max value, we will order the order value using decreasing and reorder the markers
  markers <-  markers[order(row.o, decreasing=F),]
  z.means <- z.means[order(row.o, decreasing=F),]
  
  ##reorder using the z.means file
  combine_dt_addzscore_addcomm$Gene <- factor(combine_dt_addzscore_addcomm$Gene, levels = as.character(rownames(z.means)))
  combine_dt_addzscore_addcomm <- combine_dt_addzscore_addcomm[order(combine_dt_addzscore_addcomm$Gene),]
  combine_dt_addzscore_addcomm$CellType <- factor(combine_dt_addzscore_addcomm$CellType, levels = as.character(colnames(z.means)))
  combine_dt_addzscore_addcomm$name <- factor(combine_dt_addzscore_addcomm$name, levels = as.character(unique(combine_dt_addzscore_addcomm$name)))
  
  ##do the plotting
  p <- ggplot(combine_dt_addzscore_addcomm, aes(y=CellType, x = name, 
                                                color = zscore, size = AccProp)) + 
    geom_point() + 
    theme(axis.text.x = element_text(angle = 90)) + 
    #scale_colour_gradient2(high = "gray0",  
    #                       mid = 'gray90',
    #                       low = "gray100",
    #                       midpoint = 0.5) +
    scale_colour_gradient2(midpoint = 0.5) +  ##original is 0.5
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.text.x = element_text(colour = "black", size=12,angle = 45, vjust = 1,hjust =1),
          axis.text.y = element_text(colour = "black", size=12),
          axis.line.x = element_line(color="black", size = 1), 
          axis.line.y = element_line(color="black", size = 1)) +
    ggtitle("Z Score of Markers - Across Annotations")
  p ##the middle point 
  pdf(paste0(input_output_dir,'/opt_dotplot_marker.pdf'),width = plot_width,height = plot_height)
  print(p)
  dev.off()
  
  ##updating 090722 we will set an option to check whether we will pot name and ct at the same time
  if (OpenaddCTtoName == 'yes'){
    
    combine_dt_addzscore_addcomm$CTnm <- paste0(combine_dt_addzscore_addcomm$name,'.',combine_dt_addzscore_addcomm$AAtype)
    combine_dt_addzscore_addcomm$CTnm <- factor(combine_dt_addzscore_addcomm$CTnm, levels = as.character(unique(combine_dt_addzscore_addcomm$CTnm)))
    
    
    p <- ggplot(combine_dt_addzscore_addcomm, aes(y=CellType, x = CTnm, 
                                                  color = zscore, size = AccProp)) + 
      geom_point() + 
      theme(axis.text.x = element_text(angle = 90)) + 
      #scale_colour_gradient2(high = "gray0",  
      #                       mid = 'gray90',
      #                       low = "gray100",
      #                       midpoint = 0.5) +
      scale_colour_gradient2(midpoint = 0.5) +  ##original is 0.5
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(), 
            axis.text.x = element_text(colour = "black", size=12,angle = 45, vjust = 1,hjust =1),
            axis.text.y = element_text(colour = "black", size=12),
            axis.line.x = element_line(color="black", size = 1), 
            axis.line.y = element_line(color="black", size = 1)) +
      ggtitle("Z Score of Markers - Across Annotations")
    p ##the middle point 
    pdf(paste0(input_output_dir,'/opt_',prefix,'_addCTnm_dotplot.pdf'),width = plot_width,height = plot_height)
    print(p)
    dev.off()
    
  }
  
}
#?scale_colour_gradient2


##plotting
prefix <- 'marker'

opt_prepare_plot_obj <- readRDS(opt_prepare_plot_obj_fl)
plot_markers(prefix,input_marker_gene_fl,opt_prepare_plot_obj,
             plot_width,plot_height,OpenaddCTtoName='no')









