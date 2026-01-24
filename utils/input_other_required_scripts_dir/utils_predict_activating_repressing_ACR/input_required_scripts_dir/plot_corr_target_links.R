# load libraries
library(Matrix)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# args
args <- commandArgs(T)

ipt_pmat_fl <- as.character(args[1])
ipt_gmat_fl <- as.character(args[2])
ipt_gene_link_annot_fl <- as.character(args[3])
ipt_top_link_num <- as.numeric(args[4])
ipt_output_dir <- as.character(args[5])
ipt_type_str <- as.character(args[6]) ##pLink
ipt_meta_fl <- as.character(args[7]) ##use for the heatmap cell type colors
ipt_target_pair_id_fl <- as.character(args[8]) ##na or have a file
target_col_in_meta <- 'cell_identity'
target_color_in_meta <- 'celltype_color'

####################################
##step01 plot the correlation panels

ipt_gmat_dt <- read.table(ipt_gmat_fl)
head(ipt_gmat_dt)
ipt_pmat_dt <- read.table(ipt_pmat_fl)
head(ipt_pmat_dt)

ipt_gene_link_annot_dt <- read.table(ipt_gene_link_annot_fl)
head(ipt_gene_link_annot_dt)
dim(ipt_gene_link_annot_dt)

ipt_gene_link_annot_dt <- ipt_gene_link_annot_dt %>% distinct(pair_id, .keep_all = TRUE)

if (ipt_target_pair_id_fl == 'na'){
  
  ipt_gene_link_annot_dt <- ipt_gene_link_annot_dt[ipt_gene_link_annot_dt$link.type == ipt_type_str,]
  
  ##sort the pvalue to get the top for the silencer and enhancer
  top_20 <- ipt_gene_link_annot_dt %>%
    group_by(cre.type) %>%
    arrange(p.value) %>%
    slice_head(n = ipt_top_link_num) %>%  
    ungroup()
  
  pair_id_list <- top_20$pair_id
  rownames(top_20) <- top_20$pair_id

}else{
  
  ipt_target_pair_id_dt <- read.delim('./ipt_target_pair_id_fl.txt',header = F)
  pair_id_list <- ipt_target_pair_id_dt$V1
  length(pair_id_list)
  
  top_20 <- ipt_gene_link_annot_dt[ipt_gene_link_annot_dt$pair_id %in% pair_id_list,]
  dim(top_20)
  
}

for (i in 1:length(pair_id_list)){
  
  ipt_pair_id <- pair_id_list[i]
  result <- strsplit(ipt_pair_id, ":")[[1]]
  
  peakID <- result[1]
  geneID <- result[2]
  
  ipt_pair_id_dt <- top_20[ipt_pair_id,]
  opt_cate <- ipt_pair_id_dt$cre.type
  opt_distance2gene <- ipt_pair_id_dt$dis2gene
  
  ipt_gmat_target_dt <- ipt_gmat_dt[geneID,]
  ipt_pmat_target_dt <- ipt_pmat_dt[peakID,]
  
  merged_dt <- rbind(ipt_gmat_target_dt,ipt_pmat_target_dt)
  
  mat <- t(merged_dt[c(peakID, geneID), ])
  
  
  df_plot <- data.frame(
    peak = mat[,1],
    gene = mat[,2]
  )
  
  p <- ggplot(df_plot, aes(x = peak, y = gene)) +
    geom_point(color = "blue",size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    ggtitle("Correlation between peak and gene") +
    annotate("text", x = min(df_plot$peak), y = max(df_plot$gene),
             label = paste("r =", round(cor(df_plot$peak, df_plot$gene), 2)),
             hjust = 0, vjust = 1)+
    theme(
      plot.title = element_text(face="bold.italic",size=20,hjust = 0.5),
      axis.title = element_text(size =20, face="bold"),
      #axis.text = element_text(angle=45, vjust=1,size=15,colour = 'black'),
      axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
      axis.text.y = element_text(size=20,colour = "black",face = 'bold'),
      axis.ticks = element_line(size = rel(2.5)),
      axis.ticks.length = unit(0.5, "cm"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      strip.text = element_text(size=20),
      #legend.title = element_blank(),
      #legend.position = "none",
      plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
    )
  
  pdf(file = paste0(ipt_output_dir,'/opt_',opt_cate,'_',peakID,'_',geneID,'_ds2gene',opt_distance2gene,'bp_',ipt_type_str,'.pdf'),width =7,height = 7)
  print(p)
  dev.off()
}

##############################################################
##step02 plot the heatmap for the accessibility and expression
ipt_gmat_dt <- read.table(ipt_gmat_fl)
head(ipt_gmat_dt)
ipt_pmat_dt <- read.table(ipt_pmat_fl)
head(ipt_pmat_dt)

ipt_gene_link_annot_dt <- read.table(ipt_gene_link_annot_fl)
head(ipt_gene_link_annot_dt)
dim(ipt_gene_link_annot_dt)

ipt_gene_link_annot_dt <- ipt_gene_link_annot_dt %>% distinct(pair_id, .keep_all = TRUE)

ipt_gene_link_annot_activate_ACR_dt <- ipt_gene_link_annot_dt[ipt_gene_link_annot_dt$cre.type == 'enhancer',]
ipt_gene_link_annot_repressing_ACR_dt <- ipt_gene_link_annot_dt[ipt_gene_link_annot_dt$cre.type == 'silencer',]

############################
##Selection#################
target_type <- 'activate'

target_type <- 'repressing'
############################


if (target_type == 'activate'){
  ipt_target_gene_annot_ACR_dt = ipt_gene_link_annot_activate_ACR_dt
  rownames(ipt_target_gene_annot_ACR_dt) <- ipt_target_gene_annot_ACR_dt$peak_id
}else{
  ipt_target_gene_annot_ACR_dt = ipt_gene_link_annot_repressing_ACR_dt
  rownames(ipt_target_gene_annot_ACR_dt) <- ipt_target_gene_annot_ACR_dt$peak_id
}

##plot the activate ACR gene accessiblity plot using cpm
target_peakID_list <- ipt_target_gene_annot_ACR_dt$peak_id
opt_peak_num <- length(unique(target_peakID_list)) ##2426

ipt_pmat_target_dt <- ipt_pmat_dt[target_peakID_list,]
mat <- ipt_pmat_target_dt

############
##plot peaks
z <- t(as.matrix(scale(t(as.matrix(mat)))))
dim(z)

o.order <- colnames(z)
col.clust <- hclust(as.dist(1-cor(z)))
col.o <- col.clust$order
col.dendro <- as.dendrogram(col.clust)
z <- z[,col.o]
z <- z[order(apply(z, 1, which.max), decreasing=F),]
z <- z[,o.order]
head(z)
n.range <- c(-4, -2, 0, 2, 4)

meta <- read.delim(ipt_meta_fl,row.names = 1)
head(meta)
meta[[target_col_in_meta]] <- tolower(meta[[target_col_in_meta]] )
shared_celltypes <- intersect(colnames(z),unique(meta[[target_col_in_meta]]))
meta <- meta[meta[[target_col_in_meta]]%in%shared_celltypes,]
meta_TCP <- meta[!duplicated(meta[[target_col_in_meta]]),]
m.cols <- as.character(meta_TCP[[target_color_in_meta]])
names(m.cols) <- as.character(meta_TCP[[target_col_in_meta]])
m.cols <- m.cols[colnames(z)]
ha.col = HeatmapAnnotation(clusterID = colnames(z), col = list(clusterID = m.cols),
                           annotation_legend_param=list(clusterID=list(ncol=3, 
                                                                       grid_height = unit(2, "mm"), grid_width = unit(2, "mm"),
                                                                       labels_gp = gpar(fontsize = 6))))

par(mar = c(7, 7, 7, 7))
pdf(paste0(ipt_output_dir,'/opt_heatmp_',target_type,'_ACRacc','_PeakNum',opt_peak_num,'.pdf'), width=7, height=12)
Heatmap(z, name = "ACR accessibility for each cell type", cluster_rows = F, cluster_columns=col.dendro, 
        top_annotation = ha.col,
        col=colorRamp2(as.numeric(n.range),c("dodgerblue4","deepskyblue","grey80","darkorange","firebrick3")),
        use_raster=T,
        show_row_names = F,
        show_column_names = T)
dev.off()

############
##plot genes
ordered_peaks <- rownames(z)
ordered_genes <- ipt_target_gene_annot_ACR_dt[ordered_peaks,]$gene_id

mat <- ipt_gmat_dt[ordered_genes,]

z_gene <- t(as.matrix(scale(t(as.matrix(mat)))))
dim(z_gene)

par(mar = c(7, 7, 7, 7))
pdf(paste0(ipt_output_dir,'/opt_heatmp_',target_type,'_Geneacc.pdf'), width=7, height=12)
Heatmap(z_gene, name = "ACR accessibility for each cell type", cluster_rows = F, cluster_columns=col.dendro, 
        top_annotation = ha.col,
        col=colorRamp2(as.numeric(n.range),c("dodgerblue4","deepskyblue","grey80","darkorange","firebrick3")),
        use_raster=T,
        show_row_names = F,
        show_column_names = T)
dev.off()











