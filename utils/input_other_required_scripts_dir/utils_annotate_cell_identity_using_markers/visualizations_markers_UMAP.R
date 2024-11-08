#Visualizing marker genes using UMAP
##load a function
library(Matrix)
library(ggplot2)

args <- commandArgs(T)

ipt_meta_fl <- as.character(args[1])

ipt_impute_mtx_fl <- as.character(args[2])

ipt_marker_fl <- as.character(args[3])

ipt_lim_val <- as.numeric(args[4])

input_output_dir <- as.character(args[5])

plot.act.scores    <- function(df,output_dir, 
                               acts=acts, 
                               info=NULL, 
                               top=NULL,
                               logT=F,
                               marker.dist=NULL,
                               outname="markerActivityScores.pdf", 
                               lim=0.95){
  
  # prep data
  df <- df[rownames(df) %in% colnames(acts),]
  acts <- acts[,which(rownames(df) %in% colnames(acts))]
  
  # reorder rows
  ##updating 102221 do not use the geneID to be the rownames
  ##we need to generaete a unique infor
  #unique(info[,c("geneID")])
  rownames(info) <- info$geneID
  
  info <- info[order(info$type),]
  info.genes <- rownames(info)
  #info.genes <- info$geneID ##the info.genes are not unique
  #info.genes <- unique(info.genes) ##make the genes to be unique
  
  act.genes <- rownames(acts)
  rd.cells <- rownames(df)
  
  # common genes
  common <- intersect(info.genes, act.genes)
  info <- info[which(rownames(info) %in% common),]
  #info <- info[which(info$geneID %in% common),] ##but the info contain not unique geneID
  
  info.ordered <- rownames(info)
  sub.scores <- acts[info.ordered,]
  gids <- info.ordered
  
  # setup plot size
  nrows <- ceiling(length(gids)/6)
  totals <- nrows*6
  ratio <- nrows/6
  
  # params
  png(file=paste0(output_dir,'/',outname), width=12, height=ratio*12, units="in", res=500, type="cairo")
  layout(matrix(c(1:totals), ncol=6, byrow=T))
  par(mar=c(2,2,1,1))
  
  # adjust cluster IDs
  message("begin plotting pre-defined markers...")
  for (i in 1:length(gids)){
    
    # copy meta data
    gene.index <- which(rownames(sub.scores) == gids[i])
    acv <- sub.scores[gene.index,]
    
    # set up plot cols/sizes
    orderRow <- order(acv, decreasing=F)
    #cols <- colorRampPalette(c("grey75","grey75","goldenrod2","firebrick3"), bias=1)(100)
    #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
    #cols <- inferno(100)
    #cols <- plasma(100)
    cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.5)(100)
    acv <- as.numeric(acv[orderRow])
    if(logT==T){
      acv <- log2(acv+1)
    }
    df2 <- df[orderRow,]
    acv[is.na(acv)] <- 0
    acv[is.infinite(acv)] <- 0
    upper.lim <- quantile(acv, lim)
    acv[acv > upper.lim] <- upper.lim
    if(!is.null(marker.dist)){
      message(" - # cells = ", length(acv), "| min: ", marker.dist[[gids[i]]][1], " | max: ",marker.dist[[gids[i]]][2])
      colvec <- cols[cut(acv, breaks=seq(from=marker.dist[[gids[i]]][1], to=marker.dist[[gids[i]]][2], length.out=101))]
    }else{
      min.acv <- min(acv) - (1e-6*min(acv))
      max.acv <- max(acv) + (1e-6*max(acv))
      message(" - # cells = ", length(acv), "| min: ", min.acv, " | max: ",max.acv)
      if(min.acv == max.acv){
        next
      }
      colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv, length.out=101))]
    }
    colvec[is.na(colvec) & acv > mean(acv)] <- cols[length(cols)]
    colvec[is.na(colvec) & acv == 0] <- cols[1]
    sizes <- rescale(acv, c(0.25, 0.3))
    
    # plot
    plot(df2$umap1, df2$umap2, col=colvec,
         main=paste(info$name[i],info$type[i],info$common[i],sep="-"),
         #main=info$name[i],
         xlab="", ylab="", bty="n",
         xaxt="n", yaxt="n", pch=16, cex=0.2) ##previous cex is 0.25 we set smaller
    
  }
  
  # turn device off
  dev.off()
  
}



UMAP_plot <- function(ipt_meta_fl,
                      ipt_impute_mtx_fl,
                      ipt_marker_fl,
                      input_output_dir,
                      lim = 0.98,
                      opt_name = 'output'){
    
  meta_dt <- read.table(ipt_meta_fl)
  #meta_dt <- obj_impute$meta
  impute.activity <- readRDS(ipt_impute_mtx_fl)
  #marker.info <- obj_impute$marker_infor
  marker.info <- read.delim(ipt_marker_fl,header=T)

  plot.act.scores(meta_dt, input_output_dir,acts=impute.activity,
                  info=marker.info,
                  logT=F,
                  lim=0.98,
                  marker.dist=NULL,
                  outname=opt_name)

}


UMAP_plot(ipt_meta_fl,
          ipt_impute_mtx_fl,
          ipt_marker_fl,
          input_output_dir,
          lim = ipt_lim_val,
          opt_name = 'opt_UMAP_marker.png')





