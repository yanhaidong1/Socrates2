###################################################################################################
###################################################################################################
##                                                                                               ##
##                  functions for plotting accessibility of markers from cicero                  ##
##                                                                                               ##
###################################################################################################
###################################################################################################

##updating 102221 we finally decide use the unique marker file since there will be causing problem to plot marker
##we next add the common name when plotting
##updating 102221 debug make sure smooth.data npcs should be same as previous ncol(pcs)
##updating 102121 we will filter cluster with only small number of cells
##updating 102021 we will set rownames that not unique by combining the gene name and cell type name
##updating 080921 we need to set another new cluster Louv
##updating 080921 save the all gene information
##updating 071521 we need to check if the act file exists or not 


# load libraries
library(viridis)
library(mclust)
library(irlba)
library(Matrix)
library(RANN)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(gplots)
library(scales)
library(varistran)
library(edgeR)
library(parallel)
library(png)



args <- commandArgs(T)

config <- as.character(args[1])

input_soc_obj_fl <- as.character(args[2])

input_output_dir <- as.character(args[3])


source(config)


###################################################################################################
# plot activity scores of de novo genes
plot.new.markers   <- function(input_soc_obj_fl,output_dir,
                               outname="ActivityScores.pdf",
                               row.o=NULL,
                               top=5, 
                               normT='prop.dif',
                               logT=F, 
                               threshold=0.01,
                               lim=0.95){
  
    
    ## read the soc obj
    ipt_soc_obj <- readRDS(input_soc_obj_fl)  
    
    acts <- ipt_soc_obj$gene_acc_smooth
    df <- ipt_soc_obj$meta
  
    # verbose start-up
    if(is.null(df) | is.null(acts)){
        stop("ERROR: must provide metadata and activity matrix")
    }
    if(! "Cluster" %in% colnames(df)){
        #df$Cluster <- df$umapclust
        df$Cluster <- df$LouvainClusters
    }

    # estimate significant differences
    df <- df[colnames(acts),]
    df$Cluster <- factor(df$Cluster)

    # iterate over genes
    message(" - transversing gene activity ...")
    it <- 0
    total.mean <- Matrix::rowMeans(acts)
    props <- acts
    props@x <- rep(1, length(props@x))
    total.prop <- Matrix::rowMeans(props)
    df.splits <- split(colnames(acts), df$Cluster)

    # estimate mean, proportions
    message(" - estimating difference in proportions and means ...")
    difs.prop <- lapply(df.splits, function(z){
        sub.act <- acts[,z]
        Matrix::rowMeans(sub.act >0)
    })
    difs.mean <- lapply(df.splits, function(z){
        sub.act <- acts[,z]
        rowMeans(sub.act)
    })
    difs.prop <- do.call(cbind, difs.prop)
    difs.mean <- do.call(cbind, difs.mean)
    colnames(difs.prop) <- names(df.splits)
    colnames(difs.mean) <- names(df.splits)

    # estimate pairwise log2 fold change
    message(" - estimating mean log2 fold changes")
    ave.difs.mean <- t(apply(difs.mean, 1, function(x){
        sapply(x, function(z){mean(log2((z+1)/(x+1)), na.rm=T)})
    }))

    # save output
    #write.table(ave.difs.mean, file=paste0(output_dir,'/',outname,".NormalizedAveLog2FC.txt"),
    #            quote=F, row.names=T, col.names=T, sep="\t")
    #write.table(difs.mean, file=paste0(output_dir,'/',outname,".NormalizedClusterMeans.txt"),
    #            quote=F, row.names=T, col.names=T, sep="\t")
    #write.table(difs.prop, file=paste0(output_dir,'/',outname,".NormalizedClusterProportions.txt"),
    #            quote=F, row.names=T, col.names=T, sep="\t")

    # output options
    if(normT=="mean.dif"){
        difs <- ave.difs.mean
    }else if(normT=="prop.dif"){
        difs <- difs.prop
    }else if(normT=="adj.dif"){
        difs <- ave.difs.mean * difs.prop
    }
    difs <- difs[,mixedorder(colnames(difs))]

    # # plot normalized heatmap
    # c.means <- difs.mean[,mixedorder(colnames(difs.mean))]
    # pdf(paste0(output,".normalized_heatmap.pdf"), width=5, height=6)
    # heatmap.2(c.means, trace='none', scale='row', Colv=F, Rowv=F, dendrogram='none',
    #           useRaster=T, col=colorRampPalette(c("dodgerblue4","white","firebrick4"))(100), labRow=F)
    # dev.off()

    # filter out genes by proportion in at least one cluster
    message(" - selecting genes passing minimum proportion threshold ...")
    pass.genes <- apply(difs.prop, 1, function(x){length(x[x>threshold])})
    difs <- difs[pass.genes > 0,]

    # melt
    reduced <- melt(difs)
    reduced <- reduced[order(reduced$value, decreasing=T),]
    top.genes <- Reduce(rbind, by(reduced, reduced$Var2, head, n=top))
    top.genes.out <- subset(reduced, reduced$value >= 1)
    #write.table(top.genes.out, file=paste0(output_dir,'/',outname,".log2FC1genes.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    topN <- top.genes$Var1
    message(" - # genes = ", nrow(top.genes), " | # per cluster: ", top)
    print(head(top.genes))
  
    ##updatin 080921
    message(' - save the rank difs')
    ##we need to write out all the top
    all_top.genes <- Reduce(rbind, by(reduced, reduced$Var2, head, n=1000000000))
    ##write out the all_top.genes information for each cluster
    #write.csv(all_top.genes,paste0(output_dir,'/',outname,".rank_difs.csv"),quote = F)
    
    ##save the rank file to the soc obj
    final_obj <- append(ipt_soc_obj, list(
      rank_diff_gene = all_top.genes
    ))
    
    saveRDS(final_obj,file=paste0(output_dir,'/',input_prefix,'.atac.soc.rds'))
    
    ##save the top genes
    top.defined.genes <- Reduce(rbind, by(reduced, reduced$Var2, head, n=top))
    write.table(top.defined.genes,paste0(output_dir,'/',input_prefix,'.top',top,"genes.rank_difs.txt"),sep='\t',quote = F)
    
    
    # plot top 50 -----------------------------------------------------------------------------
    nrows <- ceiling(length(topN)/top)
    totals <- nrows*top
    ratio <- nrows/top

    # params
    pdf(file=paste0(output_dir,'/',input_prefix,".denovo.pdf"), width=20, height=ratio*20)
    layout(matrix(c(1:totals), ncol=top, byrow=T))
    par(mar=c(2,2,1,1))

    # adjust cluster IDs
    message("begin plotting ...")
    for (i in 1:nrow(top.genes)){

        # copy meta data
        dff <- df
        geneID <- top.genes$Var1[i]
        acv <- as.numeric(acts[rownames(acts) %in% geneID,])
        if(logT==T){
            acv <- log2(acv+1)
        }
        message(" - gene: ", geneID)

        # set up plot cols/sizes
        orderRow <- order(acv, decreasing=F)
        #cols <- colorRampPalette(c("grey85",brewer.pal(9,"YlGnBu")), bias=1)(100)
        #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
        #cols <- inferno(100)
        #cols <- plasma(100)
        cols <- colorRampPalette(viridis(100), bias=0.5)(100)
        #cols <- colorRampPalette(c("grey75","darkorange","firebrick3"), bias=0.75)(100)
        acv <- acv[orderRow]
        df2 <- df[orderRow,]

        # set up upper limits
        upper.lim <- quantile(acv, lim)
        acv[acv > upper.lim] <- upper.lim

        min.acv <- -0.1
        max.acv <- max(acv) + (0.05*max(acv))
        colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv,length.out=101))]
        colvec[is.na(colvec)] <- cols[1]
        sizes <- rescale(acv, c(0.35, 0.4))

        # plot
        plot(df2$umap1, df2$umap2,
             col=colvec,
             main=paste(top.genes$Var1[i],top.genes$Var2[i],sep="-"),
             xlab="", ylab="", bty="n", xaxt="n", yaxt="n", pch=16, cex=0.5)

    }

    # turn device off
    dev.off()

}


plot.new.markers(input_soc_obj_fl,input_output_dir,
                 top=top_gene,
                 normT=normT_type,
                 threshold=threshold_prop,
                 logT=F,
                 lim=0.99,
                 outname=paste0(output,".imputed"))








