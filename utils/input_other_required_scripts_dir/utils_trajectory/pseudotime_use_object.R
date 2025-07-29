##updating 072825 we will set an option to plot umap if there is no color provided
##updating 071625 we will set an cutoff for the genes with FDR < 0.05 cutoff
##updating 051425 we will set an option to decide which matrix will be loaded 
##udpating 071223 we will change the color of each traj
##updating 071023 we will open the traject genes anlaysis add the beta for the lm
##updating 041623 we will set an option to open the equally cell number across different clusters
##updating 010623 do not add the gene item as it is too large
##updating 010523 we will add the other motif acr etc.
##updating 111722 sub_cluster_color we need to set an option to denote which column we will look at to draw the cluster
##updating 080822 we would check the pseudotime across different tissues
##updating 081121 decide which parts we will open
##updating 070521 open the plotting and change colors
##updating 070221 plot the target genes
##updating 063021 open all the tf and motif and have a checking
##updating 061421 open the motif to have a checking
##updating 061021 add a step to transfer peak sparse to matrix and save a temp file
## pseudotime ##
##this script will check the trajectory 
##1) regular checking
##2) motif
##3) tf 
##we will use the top_tissue_cluster 
##also change the LouvainClusters to LouvainClusters_afthm

##the mechainims of significant motifs
##is use one-way anova check whether trajetory time influence on the motif acc
##calculate the F value and check the p value


##we will first try the TFs

# load arguments
arg <- commandArgs(T)
print(arg)
#if(length(arg)!= 9){stop("pseudotime.R <sparse> <motif.deviations> <gene> <TFgene> <meta> <svd> <prefix> <config> <threads>")}

sparse <- as.character(arg[1]) ##peak matrix rds file
##sparse mtx rds


motif <- as.character(arg[2]) ## /scratch/hy17471/soybean_scATAC_100120/pipeline_analysis_110220/add_11_chromVAR_motif_dev_010221/02_run_chromVAR_041221/all.smoothed_motifs.txt
##mtx  should be smooth

gene <- as.character(arg[3]) ## /scratch/hy17471/soybean_scATAC_100120/pipeline_analysis_110220/add_06_geneAccessibility_UMAP_analysis_regmodel2_res01_shift_subc_comb_nbi_TPOnorm_newpeak_041221/all.normalizedActivity.sparse
##udpating 051425
##gene acc from the gene object

##three column file

tfgene <- as.character(arg[4])
##updating 051425
##this is the tfgene sparse mtx
##three column file
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/10_corssTissue_analysis_080222/05_predict_pesudoTime_080822/01_make_TF_acc_080822/output_dir/opt_TF_normalized_activity.sparse


ipt_object_fl <- as.character(arg[5])
#meta <- as.character(arg[5])
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/06_1_seperate_clustering_112421/vascular_080322/output_dir_win30000_PCs30_SVD/opt_tfidf_SVD_30PCs_win30000_res1_metadata.txt

#svd <- as.character(arg[6])
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/06_1_seperate_clustering_112421/vascular_080322/output_dir_win30000_PCs30_SVD/opt_tfidf_SVD_30PCs_win30000_res1_reduced_dimensions.txt

prefix <- as.character(arg[6])
##vascular

config <- as.character(arg[7]) ##this config records the target clusters: since we have different cell types we will use this information
##vascular_related_cells_config.txt
##it contains 
#traj <- c('Axis_Shoot_Apical_Meri&Root_Meristem_10','Cotyledon_Adaxial_Parenchyma_12','Cotyledon_Adaxial_Parenchyma_8','Cotyledon_Adaxial_Parenchyma_29')
#LC <- c(10,12,8,29)

threads <- as.numeric(arg[8])

output_dir <- as.character(arg[9])

target_cluster_col <- as.character(arg[10]) ##LouvainClusters
target_cluster_color <- as.character(arg[11]) ##na if no color

##updating loading gene 010623
#open_load_all_gene_acc <- as.character(arg[13]) ##yes or no
ipt_target_gene_list_fl <- as.character(arg[12]) ##if no we will use na to represent

input_number_top_plot_ACR_peudo <- as.numeric(arg[13])

##updating 041623
open_equal_cellnum_diff_clusters <- as.character(arg[14])

##updating 071023
whether_open_provided_reduced_meta_fl <- as.character(arg[15])

##updating 070221
#target_gene_set <- as.character(arg[11])

##updating 070621
#input_number_top_plot_ACR_peudo <- as.numeric(arg[12])

##updating 081121
#input_config_soure_fl <- as.character(arg[13])

##updating 081121
#collect_targetGN_set_opt_dir <- as.character(arg[8])

##updating 081121 
#source(input_config_soure_fl) ##this means whether we need to open or close some functions



##########
# defaults
#column <- "top_tissue_cluster"
column <- target_cluster_col
#openTF <- 'yes' ##this has been wrote in the config files
#LC <- 1
featureMin <- 0

##################
# load config file
source(config) ##there is a traj already in the config
##we will directly use the umap1 and umap2 to do the analysis
umap1 <- "umapPT_1"
umap2 <- "umapPT_2"

# load libraries
library(Matrix)
library(nabor)
library(dplyr)
library(viridis)
library(data.table)
library(scales)
library(mgcv)
library(gplots)
library(RColorBrewer)
library(parallel)
library(speedglm)
library(gtools)
library(uwot)
library(splines) ##provide the ns function
library(RANN)
library(phytools)



# load functions
smooth.data<- function(x, k=15, step=3, npcs=30, df=NULL, rds=NULL, verbose=F){
    
    # verbose
    if(verbose){message(" - imputing gene activity ...")}
    
    # input
    data.use <- x
    
    # verbose
    if(!is.null(rds)){
        
        if(!is.null(df)){
            if(verbose){message("   * using UMAP manifold for smoothing ...")}
            pcs <- df[,c("umap1","umap2")]
        }else{
            if(verbose){message("   * using prior PC space as manifold ...")}
            pcs <- rds[colnames(x),c(1:npcs)]
        }
    }else{
        
        # LSI
        if(verbose){message("   * PC manifold set to NULL, running LSI (TFIDF)...")}
        x[x>0] <- 1
        tf.idf <- tfidf(x)
        
        # get PCS
        if(verbose){message("   * PC manifold set to NULL, running LSI ...")}
        pc <- irlba(t(tf.idf), npcs)
        pcs <- pc$u 
        rownames(pcs) <- colnames(x)
        colnames(pcs) <- paste0("PC_", seq(1:ncol(pcs)))
        
        # do l2-norm
        pcs <- apply(pcs, 2, function(x){x/sqrt(sum(x^2))})
    }
    
    # get KNN
    if(verbose){message("   * finding knn graph ...")}
    knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
    j <- as.numeric(x = t(x = knn.graph))
    i <- ((1:length(x = j)) - 1) %/% k + 1
    edgeList = data.frame(i, j, 1)
    A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
    
    # Smooth graph
    if(verbose){message("   * smoothing graph ...")}
    A = A + t(A)
    A = A / Matrix::rowSums(A)
    step.size = step
    if(step.size > 1){
        for(i in 1:step.size){
            if(verbose){message("     ~ step ",i)}
            A = A %*% A
        }
    }
    
    # smooth data
    if(verbose){message("   * smoothing activity ...")}
    impute.activity <- t(A %*% t(data.use))
    colnames(impute.activity) <- colnames(x)
    rownames(impute.activity) <- rownames(x)
    
    # clean empty rows (if present) and round to two decimal places
    #if(verbose){message("   * remove empty rows/columns and scale to per million ...")}
    #impute.activity <- impute.activity[Matrix::rowSums(impute.activity)>0,]
    #impute.activity <- impute.activity[,Matrix::colSums(impute.activity)>0]
    #impute.activity <- impute.activity %*% Diagonal(x=1e6/Matrix::colSums(impute.activity))
    #impute.activity@x <- round(impute.activity@x, digits=2)
    
    
    # clean empty rows (if present) and round to two decimal places
    if(verbose){message("   * remove empty rows/columns and scale to per million ...")}
    ##this step will filter out many motifs that show large impute activity
    impute.activity <- impute.activity[Matrix::rowSums(impute.activity)>0,]
    impute.activity <- impute.activity[,Matrix::colSums(impute.activity)>0]
    
    ##updating 031924
    ##debug for add the name to the impute activity
    impute.activity_storecolnm <- impute.activity
    impute.activity <- impute.activity %*% Diagonal(x=1e6/Matrix::colSums(impute.activity))
    colnames(impute.activity) <- colnames(impute.activity_storecolnm)
    
    impute.activity@x <- round(impute.activity@x, digits=2)
    
    
    # return sparse Matrix
    return(impute.activity)
}



 
##updating 051425
loadData_option_version   <- function(ipt_object,sm_fl, mt_fl,gn_fl, tf_fl, doL2=0,
                                      traj,open_equal_cellnum_diff_clusters,column="top_tissue_cluster"){

  message(" - loading meta ...")
  #ma <- read.delim(ma_fl,row.names = 1)
  ma <- ipt_object$final_meta
  
  message(" - loading rd ...")
  rd <- ipt_object$svd
  #rd <- readRDS(rd_fl)$svd

  
  ##load the ACR
  if (sm_fl != 'na'){
    message(" - loading ACR ...")
    sm <- readRDS(sm_fl) ##sparse matrix of ACRs
    message(" - filter ACR mtx")
    sm <- sm[Matrix::rowSums(sm)>0,]
    sm <- sm[,Matrix::colSums(sm)>0]
  }
  
  ##load the motif
  if (mt_fl != 'na'){
    message(" - loading motif ...")
    
    mt <- t(ipt_object$smooth_motif_dev)
    ##allow the col is the cell
      
    #mt <- t(read.table(mt_fl))
  }else{
    mt <- 'na'
  }
  
  ##load the gene
  if (gn_fl != 'na'){
    message(" - loading gene ...")
    
    gn <- ipt_object$gene_acc
    ##gn must be the matrix
    
    #gn <- readRDS(gn_fl)$gene_acc
    
    #gn <- read.table(gn,stringsAsFactors = T)
    #gn <- sparseMatrix(i=as.numeric(gn$V1),
    #                   j=as.numeric(gn$V2),
    #                   x=as.numeric(gn$V3),
    #                   dimnames=list(levels(gn$V1),levels(gn$V2)))
  }
  
  ##load the tf
  if (tf_fl != 'na'){
    message(" - loading tf ...")
    tf <- readRDS(tf_fl)
    #tf <- read.table(tf,stringsAsFactors = T)
    #tf <- sparseMatrix(i=as.numeric(tf$V1),
    #                   j=as.numeric(tf$V2),
    #                   x=as.numeric(tf$V3),
    #                   dimnames=list(levels(tf$V1),levels(tf$V2)))
  }
    
  
  ##intersection
  message(" - intersect ids")
  shared_ids <- intersect(rownames(ma),rownames(rd))
  rd <- rd[shared_ids,]
  ma <- ma[shared_ids,]
  
  ##intersect between meta and sm
  if (sm_fl != 'na'){
    ids <- intersect(rownames(ma), colnames(sm))
    sm <- sm[,ids]
    ma <- ma[ids,]
    #ma <- ma[colnames(sm),]
  }else{
    sm <- 'na'
  }
  
  # subset cells by cluster
  message(" - subset cells by cluster")
  message (paste0(" - dimensino of ma is ",dim(ma)," dimension"))
  as.character(traj)
  message(" - filter cells")
  ma <- ma[as.character(ma[,column]) %in% as.character(traj),]
    
  ##updating 041623
  ##we will check if the ma should be downsample or not
  if (open_equal_cellnum_diff_clusters == 'yes'){
    
    if (whether_open_provided_reduced_meta_fl == 'na'){
      
      ct_num <- table(ma[[column]])
      min_ct_num <- min(ct_num)
      
      target_ct <- unique(ma[[column]])
      
      outs <- lapply(target_ct, function(x){
        t_ct = x
        t_ct_ma <- ma[as.character(ma[[column]]) %in% as.character(t_ct),]
        sample_t_ct_cells <- sample(rownames(t_ct_ma),min_ct_num,replace=FALSE)
        sample_t_ct_ma <- ma[sample_t_ct_cells,]
        return(sample_t_ct_ma)
      })
      combine_dt <- do.call(rbind,outs)
      ##we will define ma for a new ma
      ma <- combine_dt
      
    }else{
      
      ipt_provided_cell_meta_fl = whether_open_provided_reduced_meta_fl
      ipt_provided_cell_meta_dt <- read.delim(ipt_provided_cell_meta_fl,row.names = 1)
      
      shared_IDs <- rownames(ipt_provided_cell_meta_dt,rownames(ma))
      ##filter the ma
      ma <- ma[shared_IDs,]
      
    }
  }
  
  
  if (tf_fl != 'na'){  
    intersect_tf_ma_id <- intersect(colnames(tf),rownames(ma))
    tf <- tf[,intersect_tf_ma_id]
  }else{
    tf <- 'na'
  }
  
  if (gn_fl != 'na'){
    intersect_gn_ma_id <- intersect(colnames(gn),rownames(ma))
    gn <- gn[,intersect_gn_ma_id]
    gn <- gn[Matrix::rowSums(gn)>0,]
  }else{
    gn <- 'na'
  }
  
  if (sm_fl != 'na'){
    sm <- sm[,rownames(ma)]
  }else{
    sm <- 'na'
  }
  
  rd <- rd[rownames(ma),]
  #gn <- smooth.data(gn, npcs=ncol(rd), rds=rd)
  #tf <- smooth.data(tf, npcs=ncol(rd), rds=rd)
  
  # L2
  if(doL2==1){
    rd <- t(apply(rd, 1, function(x) x/sqrt(sum(x^2))))
  }else if(doL2==2){
    rd <- apply(rd, 2, function(x) x/sqrt(sum(x^2)))
  }
  
  # re-run UMAP
  new.umap <- umap(rd, min_dist=0.1, n_neighbors=15)
  colnames(new.umap) <- c("umapPT_1", "umapPT_2")
  rownames(new.umap) <- rownames(rd)
  ma <- cbind(ma, new.umap)
  
  # return list
  return(list(a=sm, m=mt, gn=gn, tf=tf, b=ma, d=rd))
}

calcPseudo <- function(b, d, column="top_tissue_cluster", traj=NULL, cell.dist1=0.9, cell.dist2=0.99,
                       dof=250, spar=1, new.column="trajectory"){
    
    
    #### the following code was repurposed from ArchR written by Jeff Granja & Ryan Corces             ####
    #### for the original source, see https://github.com/GreenleafLab/ArchR/blob/master/R/Trajectory.R ####
    
    # initiation checks
    if(is.null(traj)){
        stop(" - Argument (vector) is missing to object: traj ...")
    }
    if(!is.character(traj)){
        stop(" - traj argument must be a character vector ...")
    }
    if(!is.character(b[,c(column)])){
        stop(" - supplied columnID, ",column,", is not a character vector ...")
    }
    
    # hidden functions - taken from archR 
    .getQuantiles <- function(v = NULL, len = length(v)){
        if(length(v) < len){
            v2 <- rep(0, len)
            v2[seq_along(v)] <- v
        }else{
            v2 <- v
        }
        p <- trunc(rank(v2))/length(v2)
        if(length(v) < len){
            p <- p[seq_along(v)]
        }
        return(p)
    }
    
    # get average coordinates for specified clusters
    ##for each cluster
    ##select the target cells
    #x <- 'root_39'
    
    ########
    ##ideas:
    ##1. analyze on each cluster based on the clusters written in the traj config file
    ##2. extract all cells from that cluster
    ##3. collect the PC information (here is 25 PCs)
    ##4. calculate average for each PC across all the cells
    ##5. select the top cells with the largest PC differences, it means we only analyze cells with the largests difference among each other,
    ##since the PC is singular value that describes the differences among different cells
    ##it looks like d is singular values for all the patterns that describes the variations among different samples
    ##we therefore want to select cells that can have the largest PC variations within each cluster, which indicates describe the largest variations among the cells
    
    ##In summary, we just filter out 5% if we set cell.dist1 to be 95%,
    ##how to filter cells, we keep the the top 95% cells with the largest PC variations
    
    t.cells <- lapply(traj, function(x){
        b.sub <- b[b[,column]==as.character(x),]
        d.sub <- d[rownames(b.sub),]
        
        #> dim(d.sub)
        #[1] 193  25(PC dimention number)
        
        ##d.sub
        ##                                      PC_1      PC_2       PC_3
        #CB:Z:AAACTCGTCAGAACGG-Root2      0.2051595 -3.598139 -0.6020455
        #CB:Z:AAAGATGGTATCACAC-crownroot2 0.3815244 -3.880755 -0.7636675
        #CB:Z:AAAGGATAGGGTTCTT-Root2      1.0485984 -3.836079 -0.1842535
        
        ##get aves of each pc
        aves <- colMeans(d.sub)
        ##> aves
        #PC_1        PC_2        PC_3        PC_4        PC_5        PC_6
        #0.27657527 -3.67750115 -0.80183601 -0.02569824 -0.76595861  1.45755063
        #PC_7        PC_8        PC_9       PC_10       PC_11       PC_12
        #0.14540460  0.75898408 -0.19017462 -0.29891631 -0.31645131 -0.20448515
        
        # filter to top 5% of cells
        ##> t(d.sub)[1:3,1:3]
        #CB:Z:AAACTCGTCAGAACGG-Root2 CB:Z:AAAGATGGTATCACAC-crownroot2
        #PC_1                   0.2051595                        0.3815244
        #PC_2                  -3.5981393                       -3.8807548
        #PC_3                  -0.6020455                       -0.7636675
        #CB:Z:AAAGGATAGGGTTCTT-Root2
        #PC_1                   1.0485984
        #PC_2                  -3.8360794
        #PC_3                  -0.1842535
        
        #######
        ##this step means we will select the cells with the largest difference of PCs
        ##??why should we use this step
        
        ##or say the 25 PCs have the largest differences
        
        ##calculate variance
        per.traj <- sqrt(colSums((t(d.sub) - aves)^2))
        ##> per.traj[1:3]
        #CB:Z:AAACTCGTCAGAACGG-Root2 CB:Z:AAAGATGGTATCACAC-crownroot2
        #1.164982                         2.695434
        #CB:Z:AAAGGATAGGGTTCTT-Root2
        #1.766252
        
        ##set the top 95% if cell.dist1 is 0.95
        idx <- which(per.traj <= quantile(per.traj, cell.dist1))
        ##> length(idx)
        #[1] 183
        
        ids.keep <- rownames(d.sub)[idx]
        #> length(ids.keep)
        #[1] 183
        ##[179] "CB:Z:TTGCGGGGTCCAGACC-Root2" "CB:Z:TTGGTCCAGCGGACAT-Root2"
        #[181] "CB:Z:TTGTTCAAGACAGCTG-Root1" "CB:Z:TTGTTGTAGCATACCT-Root2"
        #[183] "CB:Z:TTTGGTTCACGCTCAG-Root2"
        
        return(ids.keep)
    })
    names(t.cells) <- traj
    ##t.cells is a list 
    ##cluster name and then it contains a vector of cells
    
    ########
    ##ideas:
    ##1. analyze each cluster
    ##2. extract PCs for all cells from each cluster
    ##3. 
    
    # get initial trajectory
    #> seq_along(traj)
    #[1] 1 2 3 4 
    
    #x <- 1
    
    ##analyze one each cluster in the traj starting from 1
    ##> seq_along(traj)
    ##[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    
    ##Basic ideas
    ##the basic idea is to calculate the PCs differences among each nearby clusters by using dif <- sqrt(colSums((t(c.i) - c.j)^2))
    ##and transfer this difference to be the psedo time by using pt <- (1 - .getQuantiles(dif)) + x
    
    raw.pseudo <- unlist(lapply(seq_along(traj), function(x){
        
        # current cluster
        ##when the x is 1
        c.i <- d[t.cells[[x]],]
        
        ##check all the PC information from this cluster
        ##> c.i[1:3,1:3]
        #PC_1      PC_2       PC_3
        #CB:Z:AAACTCGTCAGAACGG-Root2 0.2051595 -3.598139 -0.6020455
        #CB:Z:AAAGGATAGGGTTCTT-Root2 1.0485984 -3.836079 -0.1842535
        #CB:Z:AAATGAGCAGAACAGC-Root2 0.2509564 -3.798713 -0.6890806
        
        # compute trajectory 
        ##if x != 4 which means x == 1 or 2 or 3 
        if(x != length(traj)){
          
            ##c.j is the colMean of the second one
            c.j <- colMeans(d[t.cells[[(x+1)]],])
            
            ##> c.j
            #PC_1        PC_2        PC_3        PC_4        PC_5        PC_6
            #-0.36138596 -3.38457210 -0.81722269 -0.47880210 -0.06121096  1.04534735
            #PC_7        PC_8        PC_9       PC_10       PC_11       PC_12
            #0.23318704  0.73971569 -0.59670760 -0.20875787 -0.62498055  0.58308108
            #PC_13       PC_14       PC_15       PC_16       PC_17       PC_18
            #0.35614442 -0.89672287  1.57776326 -0.23365715 -0.02372611  0.88390081
            #PC_19       PC_20       PC_21       PC_22       PC_23       PC_24
            #0.31831858  0.01045418  0.57493199  0.09552533 -0.55618413  0.93099445
            #PC_25
            #0.89456591
        
            ##> t(c.i)[1:3,1:3]
            #CB:Z:AAACTCGTCAGAACGG-Root2 CB:Z:AAAGGATAGGGTTCTT-Root2
            #PC_1                   0.2051595                   1.0485984
            #PC_2                  -3.5981393                  -3.8360794
            #PC_3                  -0.6020455                  -0.1842535
            #CB:Z:AAATGAGCAGAACAGC-Root2
            #PC_1                   0.2509564
            #PC_2                  -3.7987126
            #PC_3                  -0.6890806
            
            ##calculate the difference between the first one to the second value
            ##calculate the euclidean distance
            dif <- sqrt(colSums((t(c.i) - c.j)^2))
            
            pt <- (1 - .getQuantiles(dif)) + x
        }else{
            ##when the x is the last one
            ##c.j is the one before the last
            c.j <- colMeans(d[t.cells[[(x-1)]],])
            dif <- sqrt(colSums((t(c.i) - c.j)^2))
            pt <- .getQuantiles(dif) + x
        }
        
        return(pt)
        
    }))
    
    ##> raw.pseudo[1:10]
    #CB:Z:AAACTCGTCAGAACGG-Root2 CB:Z:AAAGGATAGGGTTCTT-Root2
    #1.765027                    1.120219
    #CB:Z:AAATGAGCAGAACAGC-Root2 CB:Z:AAATGAGGTGTGACCC-Root2
    #1.508197                    1.322404
    #CB:Z:AAATGAGTCCTATCCG-Root2 CB:Z:AACAGTCTCGGTGATT-Root2
    #1.672131                    1.568306
    #CB:Z:AACCGATGTTTAGGAA-Root2 CB:Z:AACCGATTCTGGAAGG-Root1
    #1.273224                    1.846995
    #CB:Z:AACCTGATCGTTGTAG-Root2 CB:Z:AATGGAATCACTGATG-Root2
    #1.464481                    1.863388
    ##raw.pseudo is the combination of all the cellswith the value of pt
    
    ####################
    ##step 2 smooth data
    
    # fit spline
    ##names(raw.pseudo) are the cellnames
    d.filt <- d[names(raw.pseudo),]
    
    ##extract the original PC value
    #> d.filt[1:3,1:3]
    #PC_1      PC_2       PC_3
    #CB:Z:AAACTCGTCAGAACGG-Root2 0.2051595 -3.598139 -0.6020455
    #CB:Z:AAAGGATAGGGTTCTT-Root2 1.0485984 -3.836079 -0.1842535
    #CB:Z:AAATGAGCAGAACAGC-Root2 0.2509564 -3.798713 -0.6890806
    
    d.spline <- lapply(seq_len(ncol(d.filt)), function(x){
        
        ##analyze on each PC
      
        # fit split
        ##it will fit a model that shows the relaiton between raw traject and the raw PCA 
        stats::smooth.spline(x = raw.pseudo, 
                             y = d.filt[,x], 
                             df = dof, 
                             spar = spar)[[2]]
        
        
    }) %>% Reduce("cbind",.) %>% data.frame()
    
    ##> d.spline[1:3,1:3]
    #init        V2         V3
    #1 0.5643304 -3.798877 -0.6492334
    #2 0.5604523 -3.797623 -0.6514531
    #3 0.5565742 -3.796370 -0.6536729
    
    
    # KNN fit versus actual fit
    knnObj <- nabor::knn(data = d.spline, query = d.filt, k = 3)
    
    # place along trajectory
    knnIdx <- knnObj[[1]]
    knnDist <- knnObj[[2]]
    knnDiff <- ifelse(knnIdx[,2] > knnIdx[,3], 1, -1)
    knnDistQ <- .getQuantiles(knnDist[,1])
    
    #Filter Outlier Cells to Trajectory for High Resolution
    idxKeep <- which(knnDist[,1] <= quantile(knnDist[,1], cell.dist2))
    dfTrajectory <- data.frame(row.names = rownames(d.filt),
                               distance = knnDist[, 1],
                               distanceIdx = knnIdx[, 1] + knnDiff * knnDist)[idxKeep, , drop = FALSE]
    dfTrajectory$trajectory <- 100 * .getQuantiles(dfTrajectory[,2])
    dfTrajectory <- dfTrajectory[rownames(b),]
    b[,c(new.column)] <- dfTrajectory$trajectory
    
    # return
    return(b)
}

plotPT     <- function(meta,output_dir, prefix,subsetCluster=NULL, t.id="trajectory", umap1="umapsub_1", umap2="umapsub_2",
                       smoother=5, addArrow=T, cex=0.5, xlab="umap1",ylab="umap2", bty='o'){
    
    # functions
    .centerRollMean <- function(v = NULL, k = NULL){
        o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
        if(k%%2==0){
            o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
        }else if(k%%2==1){
            o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
        }else{
            stop("Error!")
        }
        o2
    }
    
    # plot test
    if(!is.null(subsetCluster)){
        meta <- meta[meta[[target_cluster_col]] %in% subsetCluster,]
    }
    
    # define pseudotime coloring
    cols <- viridis(100)
    
    plot.new()
    add.color.bar(0.1, cols,lwd=10, title=NULL, lims=c(-1,1),prompt=F,x=0.3,y=0.5,outline=TRUE,digits = 1)
    dev.off()
    
    
    # separate by NA
    test.1 <- meta[is.na(meta[,c(t.id)]),]
    test.2 <- meta[!is.na(meta[,c(t.id)]),]
    test.3 <- test.2
    
    # layout
    pdf(paste0(output_dir,'/',prefix,'_PTplot.pdf'),width = 10,height = 10)
    
    layout(matrix(c(1:2), nrow=1))
    
    # plot grey first
    plot(test.1[,c(umap1)], test.1[,c(umap2)], xlim=range(meta[,umap1]), ylim=range(meta[,c(umap2)]),
         col=NA, pch=16, cex=cex, xlab=xlab, ylab=ylab, bty=bty)
    
    # plot trajectory values
    points(test.2[,c(umap1)], test.2[,c(umap2)], col=cols[cut(test.2[,c(t.id)], breaks=101)], pch=16, cex=cex)
    
    ##updating 072925 add the legend bar
    usr <- par("usr")
    xleft  <- usr[2] - 0.3 * diff(usr[1:2])  # 30% from right edge
    xright <- usr[2] - 0.05 * diff(usr[1:2]) # 5% from right edge
    ybottom <- usr[3] + 0.05 * diff(usr[3:4])
    ytop    <- ybottom + 0.02 * diff(usr[3:4])
    
    # Draw the gradient bar using rect()
    n <- length(cols)
    x_seq <- seq(xleft, xright, length.out = n + 1)
    for (i in 1:n) {
      rect(xleft = x_seq[i], xright = x_seq[i + 1],
           ybottom = ybottom, ytop = ytop,
           col = cols[i], border = NA)
    }
    
    # Add axis ticks
    axis_labels <- c(0, 25, 50, 75, 100)
    axis_positions <- xleft + (axis_labels / 100) * (xright - xleft)
    text(x = axis_positions,
         y = ybottom - 0.01 * diff(usr[3:4]),
         labels = axis_labels,
         cex = 0.8, xpd = NA)
    
    # Optional: label for the color bar
    text(x = (xleft + xright) / 2, y = ytop + 0.01 * diff(usr[3:4]),
         labels = "Pseudotime", cex = 0.9, font = 2, xpd = NA)
    ##################
    
    
    
    
    # add arrow
    if(addArrow){
        test.2 <- test.2[c(umap1, umap2, t.id)]
        dfArrow <-  split(test.2, floor(test.2[,c(t.id)] / 1.01)) %>% 
            lapply(colMeans) %>% Reduce("rbind",.) %>% data.frame(row.names = NULL)
        dfArrow[,c(umap1)] <- .centerRollMean(dfArrow[,c(umap1)], smoother)
        dfArrow[,c(umap2)] <- .centerRollMean(dfArrow[,c(umap2)], smoother)
        dfArrow.smooth <- smooth.spline(dfArrow[,umap1],dfArrow[,umap2], spar=1)
        
        # plot
        lines(dfArrow[,umap1], dfArrow[,umap2], lwd=3, col=alpha("black", 0.8))
        dfArrow <- dfArrow[!duplicated(dfArrow[,c(umap1,umap2)]),]
        len.df <- nrow(dfArrow)
        a.x1 <- dfArrow[,umap1][len.df]
        a.y1 <- dfArrow[,umap2][len.df]
        dir1 <- a.x1 - mean(dfArrow[,umap1][(len.df-1):(len.df-1)])
        dir2 <- a.y1 - mean(dfArrow[,umap2][(len.df-1):(len.df-1)])
        a.x2 <- a.x1+dir1
        a.y2 <- a.y1+dir2
        arrows(a.x1, a.y1, a.x2, a.y2, lwd=3, length=0.1, col=alpha("black",0.8))
    }
    
    # plot celltypes
    test.2 <- test.3
    print(head(test.2))
    #plot(test.2[,c(umap1)], test.2[,c(umap2)], col=as.character(test.2$sub_cluster_color), pch=16,
    #     xlim=range(meta[,umap1]), ylim=range(meta[,c(umap2)]),
    #     cex=cex, xlab=xlab, ylab=ylab, bty=bty)
    ##updating 111722 change the sub_cluster_color to target_cluster_col
    
    ##updating 072825
    ##check if target_cluster_color in the test.2
    if (length(intersect(colnames(test.2),target_cluster_color)) == 0){
      cols <- colorRampPalette(brewer.pal(12,"Paired")[1:10])(length(unique(test.2[,target_cluster_col])))
      
      cluster_colors <- setNames(cols, LC)  # named vector: cluster → color
      
      #colv <- cols[factor(test.2[,target_cluster_col])]
      colv <- cluster_colors[as.character(test.2[[target_cluster_col]])]
      
      plot(test.2[,c(umap1)], test.2[,c(umap2)], col=colv, pch=16,
           xlim=range(meta[,umap1]), ylim=range(meta[,c(umap2)]),
           cex=cex, xlab=xlab, ylab=ylab, bty=bty)
      
      ##updating 111722 we will close it since it overlaps with the figure
      ##updating 072825 we will open it 
      #legend("topleft", legend=sort(unique(as.character(test.2[[target_cluster_col]]))),
      #       fill=cols)
      
      legend("topleft",
             legend = LC,
             fill = cluster_colors[LC],
             border = NA,
             bty = "n")
      
      
      
    }else{
      
      #########
      test.2[[target_cluster_col]] <- as.character(test.2[[target_cluster_col]])
      test.2[[target_cluster_color]] <- as.character(test.2[[target_cluster_color]])
      #########
      
      plot(test.2[,c(umap1)], test.2[,c(umap2)], col=as.character(test.2[[target_cluster_color]]), pch=16,
           xlim=range(meta[,umap1]), ylim=range(meta[,c(umap2)]),
           cex=cex, xlab=xlab, ylab=ylab, bty=bty)
      
      ##updating 111722 we will close it since it overlaps with the figure
      ##updating 072825 we will open it 
      #cluster_labels <- sort(unique(as.character(test.2[[target_cluster_col]])))
      
      
      cluster_colors <- sapply(LC, function(cl) {
        test.2[test.2[[target_cluster_col]] == cl, target_cluster_color][1]
      })
      
      #legend("topleft", legend=sort(unique(as.character(test.2[[target_cluster_col]]))), 
      #       fill=cluster_colors)
      legend("topleft",
             legend = LC,
             fill = cluster_colors,
             border = NA,
             bty = "n")
      
    
    }
    
    #legend("topright", legend=sort(unique(as.character(test.2$celltypeID))), 
    #       col=as.character(test.2$subcluster_color)[factor(sort(unique(as.character(test.2$celltypeID))),levels=sort(unique(as.character(test.2$celltypeID))))])
    
    # add arrow
    if(addArrow){
        test.2 <- test.2[c(umap1, umap2, t.id)]
        dfArrow <-  split(test.2, floor(test.2[,c(t.id)] / 1.01)) %>% 
            lapply(colMeans) %>% Reduce("rbind",.) %>% data.frame(row.names = NULL)
        dfArrow[,c(umap1)] <- .centerRollMean(dfArrow[,c(umap1)], smoother)
        dfArrow[,c(umap2)] <- .centerRollMean(dfArrow[,c(umap2)], smoother)
        dfArrow.smooth <- smooth.spline(dfArrow[,umap1],dfArrow[,umap2], spar=1)
        
        # plot
        lines(dfArrow[,umap1], dfArrow[,umap2], lwd=3, col=alpha("black", 0.8))
        dfArrow <- dfArrow[!duplicated(dfArrow[,c(umap1,umap2)]),]
        len.df <- nrow(dfArrow)
        a.x1 <- dfArrow[,umap1][len.df]
        a.y1 <- dfArrow[,umap2][len.df]
        dir1 <- a.x1 - mean(dfArrow[,umap1][(len.df-1):(len.df-1)])
        dir2 <- a.y1 - mean(dfArrow[,umap2][(len.df-1):(len.df-1)])
        a.x2 <- a.x1+dir1
        a.y2 <- a.y1+dir2
        arrows(a.x1, a.y1, a.x2, a.y2, lwd=3, length=0.1, col=alpha("black",0.8))
    }
    
    dev.off()
}

sigPseudo  <- function(obj, meta, n.pseudo.cells=NULL, num.bins=NULL, threads=1, type="ACRs", test="binomial"){
    
    # hidden functions
    .estimate_sf_sparse <- function(counts,
                                    round_exprs=TRUE,
                                    method="mean-geometric-mean-total"){
        if (round_exprs)
            counts <- round(counts)
        
        if(method == 'mean-geometric-mean-total') {
            cell_total <- Matrix::colSums(counts)
            sfs <- cell_total / exp(mean(log(cell_total)))
        }else if(method == 'mean-geometric-mean-log-total') {
            cell_total <- Matrix::colSums(counts)
            sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
        }
        
        sfs[is.na(sfs)] <- 1
        sfs
    }
    
    # use
    if(type=="ACRs"){
        use <- obj$a
    }else if(type=="TFs"){
        use <- obj$tf
    }else if(type=="MTs"){
        use <- obj$m - min(obj$m, na.rm=T)
        print(head(use[,1:5]))
    }else if(type=="Genes"){
        use <- obj$gn
    }
    
    # order cells
    meta <- meta[!is.na(meta$trajectory),]
    meta <- meta[order(meta$trajectory, decreasing=F),]
    if(is.null(n.pseudo.cells) & !is.null(num.bins)){
        num.per <- ceiling(nrow(meta)/num.bins)
        meta$pseudoCells <- ceiling(1:nrow(meta)/num.per)
    }else if(!is.null(n.pseudo.cells)){
        meta$pseudoCells <- ceiling(1:nrow(meta)/n.pseudo.cells)
        if(max(meta$pseudoCells) < 10){
            num.per <- ceiling(nrow(meta)/10)
            meta$pseudoCells <- ceiling(1:nrow(meta)/num.per)
        }
    }else{
        meta$pseudoCells <- ceiling(1:nrow(meta)/200)
    }
    pseudocell <- split(meta, factor(meta$pseudoCells))
    
    # collapse into pseudocells
    mat <- matrix(nrow=nrow(use), ncol=length(unique(meta$pseudoCells)), 
                  dimnames=list(rownames(use), paste0("pseudo_",names(pseudocell))))
    p.cells <- lapply(names(pseudocell), function(x){
        df <- pseudocell[[x]]
        ids <- rownames(df)
        a.sub <- use[,colnames(use) %in% ids]
        if(type == "ACRs"){
            a.sum <- Matrix::rowSums(a.sub)
        }else if(type=="TFs"){
            a.sum <- Matrix::rowMeans(a.sub)
        }else if(type == "MTs"){
            a.sum <- Matrix::rowMeans(a.sub)
        }else if(type== "Genes"){
            a.sum <- Matrix::rowMeans(a.sub)
        }
        p.id <- paste0("pseudo_",x)
        mat[,p.id] <<- a.sum
        mean.site <- sum(a.sum)
        mean.traj <- mean(df$trajectory)
        outs <- c(mean.site, mean.traj, ncol(a.sub))
        names(outs) <- c("nSites", "trajectory", "num_cells")
        return(outs)
    })
    p.cells <- as.data.frame(do.call(rbind, p.cells))
    rownames(p.cells) <- paste0("pseudo_", names(pseudocell))
    mat <- Matrix(mat, sparse=T)
    mat <- mat[Matrix::rowSums(mat)>0,]
    p.cells$size_factors <- .estimate_sf_sparse(mat)

    # run tests for each site
    outs <- mclapply(seq(1:nrow(mat)), function(x){
        if((x %% 1000)==0){message(" - iterated over ",x, " sites ...")}
        df.sub <- mat[x,rownames(p.cells)]
        df.sub <- cbind(p.cells, df.sub)
        colnames(df.sub) <- c("nSites","trajectory", "num_cells","size_factors","accessibility")
        if(test=="binomial"){
            df.dat <- as.matrix(data.frame(succ=df.sub$accessibility, failure=df.sub$num_cells-df.sub$accessibility))
            #df.sub$trajectory <- as.factor(df.sub$trajectory)
            mod <- glm(df.dat ~ trajectory + nSites, 
                       data=df.sub, 
                       family=stats::binomial())
        }else if(type != "ACRs"){
            df.sub$norm_access <- df.sub$accessibility / df.sub$size_factors
            df.sub$norm_access <- round(df.sub$norm_access)
            mod <- glm(norm_access ~ trajectory + nSites, 
                       data=df.sub, 
                       family=stats::quasipoisson())
        }
        mod.sum <- summary(mod)
        traj <- mod.sum$coefficients["trajectory",]
        names(traj) <- c("Estimate","se","Tval","pval")
        return(traj)
    }, mc.cores=threads)
    outs <- as.data.frame(do.call(rbind, outs))
    rownames(outs) <- rownames(mat)
    
    # estimate q-values
    outs1 <- subset(outs, outs$pval == 1)
    outs2 <- subset(outs, outs$pval != 1)
    outs1$qval <- rep(1, nrow(outs1))
    outs2$qval <- p.adjust(outs2$pval, method="fdr")
    outs <- rbind(outs1, outs2)
    outs <- outs[mixedorder(rownames(outs), decreasing=F),]
    # return
    return(list(pseudotests=outs, aggmat=mat, pseudometa=p.cells))
    
}
sigPseudo2 <- function(obj, meta, type="ACRs", threads=1){
    
    # use
    if(type=="ACRs"){
        use <- obj$a
    }else if(type=="TFs"){
        use <- obj$tf
        use <- use %*% Diagonal(x=1e5/Matrix::colSums(use))
        colnames(use) <- colnames(obj$tf)
    }else if(type=="MTs"){
        use <- obj$m
    }else if(type=="Genes"){
        use <- obj$gn
        use <- use %*% Diagonal(x=1e5/Matrix::colSums(use))
        ##updating 081824
        ##debug as this function will deplete the colname
        colnames(use) <- colnames(obj$gn)
    }
    
    ##use is a matrix
    ##row is motif name and col is cell names
  
    # align cells
    use <- as.matrix(use)
  
    ids <- intersect(colnames(use), rownames(meta))
    use <- use[,ids]
    meta <- meta[ids,]
    meta[[target_cluster_col]] <- factor(meta[[target_cluster_col]])
    meta[[target_cluster_col]] <- droplevels(meta[[target_cluster_col]])
    
    lmp <- function (modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
    }
    
    meta$lib_size <- Matrix::colSums(use)
    
    # run anova per gene
    #x <- 1
    ##use the multiple threads to have analyses
    ##note the use would be the 0 or 1 to match the binomial test
    outs <- mclapply(seq(1:nrow(use)), function(x){
        if((x%%1000)==0){message("   ~ iterated over ",x," gene/peak/motifs ...")}
        if(type=="ACRs"){
            peak1 <- as.numeric(residuals(glm(as.numeric(use[x,])~meta$log10nSites, family=binomial()), type="response"))
        }else{
            peak1 <- as.numeric(use[x,])
        }
        ##use the y~ns(x,df) function yields a smooth curve given by a linear combination
        ##it is kind of splines
        mod <- lm(peak1~ns(meta$trajectory, df=6))
        res <- lmp(mod) ##permutation tests for the linear model
        names(res) <- rownames(use)[x]
        return(res)
    }, mc.cores=threads)
    pval <- do.call(c, outs)
    ids <- names(pval)
    qval <- p.adjust(pval, method="fdr")
    
    ##updating 071023 we will add the coef
    outs <- mclapply(seq(1:nrow(use)), function(x){
      if((x%%1000)==0){message("   ~ iterated over ",x," gene/peak/motifs ...")}
      if(type=="ACRs"){
        peak1 <- as.numeric(residuals(glm(as.numeric(use[x,])~meta$log10nSites, family=binomial()), type="response"))
      }else{
        peak1 <- as.numeric(use[x,])
      }
      ##use the y~ns(x,df) function yields a smooth curve given by a linear combination
      ##it is kind of splines
      mod <- lm(peak1~ns(meta$trajectory, df=6))
      res <- summary(mod)$coefficients[2]
      names(res) <- rownames(use)[x]
      return(res)
    }, mc.cores=threads)
    coef_val <- do.call(c, outs)
    
    df <- data.frame(pval=pval,qval=qval,beta=coef_val,row.names=ids)
    return(df)
}
plotTrajHM <- function(obj, pt, cluster=1, prefix="temp", threads=1, top=30000, featureMin=0, tests=NULL, FDR=0.05){
    
    # subset traj
    message(" - cleaning input ...")
    pt <- subset(pt, pt[[target_cluster_col]] %in% cluster)
    pt <- pt[!is.na(pt$trajectory),]
    pt <- pt[order(pt$trajectory, decreasing=F),]
    binary <- obj$a[,rownames(pt)]
    
    # filter by tests
    if(!is.null(tests)){
        message(" - filter by differential testing ...")
        tests <- subset(tests, tests$qval < FDR)
        binary <- binary[rownames(binary) %in% rownames(tests),]
    }
    binary <- binary[Matrix::rowSums(binary)>featureMin,]
    binary <- binary[,Matrix::colSums(binary)>0]
    pt <- pt[colnames(binary),]
    message(" - filtered zero-sum columns/rows of binary cell x ACR matrix : ",ncol(binary)," | ", nrow(binary))
    
    # generalized additive model for logistic regression
    message(" - running generalized additive model ...")
    newdat <- data.frame(p.time=seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500))
    fit <- mclapply(seq(1:nrow(binary)),function(x){
        df <- data.frame(acc=as.numeric(binary[x,]), p.time=pt$trajectory, lib=pt$log10nSites)
        mod.scores <- glm(acc~lib, data=df)
        df$res.acc <- residuals(mod.scores)
        mod <- gam(res.acc~s(p.time, bs="cr"), data=df)
        pred <- predict(mod, newdat, type="response")
        zscore <- (pred - mean(pred, na.rm=T))/sd(pred, na.rm=T)
        return(zscore)
    }, mc.cores=threads)
    names(fit) <- rownames(binary)
    fit <- do.call(rbind, fit)
    row.o <- apply(fit, 1, which.max)
    fit <- fit[order(row.o, decreasing=F),]
    write.table(fit, file=paste0(output_dir,'/',prefix,".ACR_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    ##updating use this scaled would work
    # plot top 50% by variance
    fit <- t(apply(fit, 1, function(x){scales::rescale(x, c(-1,1))}))
    
    ##updating use this scaled would work
    #fit <- t(apply(fit, 1, function(x) {
    #  if (all(is.na(x)) || diff(range(x, na.rm = TRUE)) == 0) {
    #    rep(0, length(x))
    #  } else {
    #    scales::rescale(x, to = c(-1, 1))
    #  }
    #}))
    
    
    # reformat output
    row.o <- apply(fit, 1, which.max)
    fit <- fit[order(row.o, decreasing=F),]
    write.table(fit, file=paste0(output_dir,'/',prefix,".ACR_pt_sorted.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    
    # plot
    message(" - plotting cell trajectory ...")
    #cols <- colorRampPalette(c(rev(brewer.pal(8,'RdGy'))))(100)
    cols <- colorRampPalette(c("paleturquoise4", "white","palevioletred3"))(100)
    pdf(paste0(output_dir,'/',prefix,".trajectoryACR.pdf"), width=10, height=10)
    heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
              scale="none", labRow = NA, labCol=NA, useRaster=T, 
              ylab=paste("ACRs", paste0("(n=",nrow(fit),")"), sep=" "),
              key = TRUE,                # show color key (legend bar)
              key.title = '',      # title for the legend
              key.xlab = "ACR chromatin accessibility",      # x-axis label for legend
              keysize = 1,            # size of the legend bar
              density.info = "none"
    )
    dev.off()
    
    ##add the legend bar plot
    #pdf(paste0(output_dir,'/',prefix,".trajectoryACR_lb.pdf"), width=5, height=5)
    #plot.new()
    #add.color.bar(0.1, cols,lwd=10, title=NULL, lims=c(-1,1),prompt=F,x=0.3,y=0.5,outline=TRUE,digits = 1)
    #dev.off()
    
    
    # return
    return(fit)
    
}
plotTrajMT <- function(obj, pt, cluster=1, prefix="temp", threads=1, tests=NULL, FDR=0.05){
    
    # subset traj
    message(" - cleaning input ...")
    pt <- subset(pt, pt[[target_cluster_col]] %in% cluster)
    pt <- pt[!is.na(pt$trajectory),]
    pt <- pt[order(pt$trajectory, decreasing=F),]
    shared_cells <- intersect(rownames(pt),colnames(obj$m))
    
    motif <- obj$m[,shared_cells]
    
    # filter by tests
    if(!is.null(tests)){
        message(" - filter by differential testing ...")
        tests <- subset(tests, tests$qval < FDR)
        motif <- motif[rownames(motif) %in% rownames(tests),]
    }
    pt <- pt[colnames(motif),]
    message(" - filtered zero-sum columns/rows of cell x motif matrix : ",ncol(motif)," | ", nrow(motif))
    
    # generalized additive model for logistic regression
    message(" - running generalized additive model ...")
    newdat <- data.frame(p.time=seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500))
    fit <- mclapply(seq(1:nrow(motif)),function(x){
        df <- data.frame(acc=as.numeric(motif[x,]), p.time=pt$trajectory)
        mod <- gam(acc~s(p.time, bs="cr"), data=df)
        pred <- predict(mod, newdat, type="response")
        zscore <- (pred - mean(pred, na.rm=T))/sd(pred, na.rm=T)
        return(zscore)
    }, mc.cores=threads)
    names(fit) <- rownames(motif)
    fit <- do.call(rbind, fit)
    
    fit <- t(apply(fit, 1, function(x){scales::rescale(x, c(-1,1))}))
    
   
    
    
    
    write.table(fit, file=paste0(output_dir,'/',prefix,".Mt_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    # reformat output
    row.o <- apply(fit, 1, which.max)
    fit <- fit[order(row.o, decreasing=F),]
    write.table(fit, file=paste0(output_dir,'/',prefix,".Mt_pt_sorted.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    fit <- as.matrix(fit)
    
    # plot
    message(" - plotting cell trajectory ...")
    cols <- colorRampPalette(c(rev(brewer.pal(8,'PiYG'))))(100)
    #cols <- colorRampPalette(c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3"))(100)
    pdf(paste0(output_dir,'/',prefix,".trajectoryMT.pdf"), width=10, height=10)
    heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
              scale="none", labRow = NA, labCol=NA, useRaster=T,
              ylab=paste("Motifs", paste0("(n=",nrow(fit),")"), sep=" "),
              key = TRUE,                # show color key (legend bar)
              key.title = '',      # title for the legend
              key.xlab = "Motif deviation scores",      # x-axis label for legend
              keysize = 1,            # size of the legend bar
              density.info = "none")
    dev.off()
    
    ##add the legend bar plot
    ##do not add the bar as we already have in the heatmap.2
    #pdf(paste0(output_dir,'/',prefix,".trajectoryMT_lb.pdf"), width=5, height=5)
    #plot.new()
    #add.color.bar(0.1, cols,lwd=10, title=NULL, lims=c(-1,1),prompt=F,x=0.3,y=0.5,outline=TRUE,digits = 1)
    #dev.off()
    
    
    # return
    return(fit)
    
}
plotTrajTF <- function(obj, pt, cluster=1, prefix="temp", threads=1, tests=NULL, FDR=0.05){
  
    ##basic ideas
    ##1. select the target clusters and order them from small to large across all cells
    ##2. smooth the the matrix of tf (row is gene ID and col is the cell name)
    ##3. only select the TFs less than 0.05 FDR test (this is the reason why we see ony a part of TFs could be plotted)
    ##4. create a new data frame to generate 500 points from min(pt$trajectory, na.rm=T) to max(pt$trajectory, na.rm=T)
    ##for example, > seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500)
  #[1]   0.1272265   0.3273723   0.5275181   0.7276640   0.9278098   1.1279557
  #[7]   1.3281015   1.5282473   1.7283932   1.9285390   2.1286849   2.3288307
  #[13]   2.5289765   2.7291224   2.9292682   3.1294140   3.3295599   3.5297057
  # ....100
  #.....[100]
    ##we totally have 500 points
    
    ##5. For each TF, it has a vector of acc and a vector of pseudo time
    ##use mod <- gam(acc ~ s(p.time,bs='cr'),data=df)
    ##next, we will predict the acc using our previous created 500 points (it has same value difference among different point)
    ##pred <- predict(mod, newdat, type="response")
    ##we transfer the pred acc to be zcore
    ##6 plot
  
    ##In summary, we will use the known relation of acc and pseduo time to predict the acc using the simulated pseduo tiem points
    ##Using this analysis, we will can predict acc change according to the pseudo time change.
  


   
  
    ##fit a logistic regression to plot TFs  
    
    # subset traj
    message(" - cleaning input ...")
    ##only consider the target cluster based on specifying the target_cluster_col
    pt <- subset(pt, pt[[target_cluster_col]] %in% cluster)
    pt <- pt[!is.na(pt$trajectory),]
    ##now the trajectory is from small to large
    pt <- pt[order(pt$trajectory, decreasing=F),]
    ##smooth the data using the same coordinate and pcs
    obj$tf <- smooth.data(obj$tf, npcs=ncol(obj$d), rds=obj$d)
    
    shared_IDs <- intersect(colnames(obj$tf),rownames(pt))
    pt <- pt[shared_IDs,]
    binary <- obj$tf[,shared_IDs]
    
    #binary <- obj$tf[,rownames(pt)]
    ##> binary[1:3,1:3]
    #3 x 3 sparse Matrix of class "dgCMatrix"
    #CB:Z:TGCTATTCATGGTTTG-PI93_seed_rep1
    #Glyma.01G000600                                80.63
    #Glyma.01G002100                                51.52
    #Glyma.01G003000                                12.86
    #CB:Z:CCACGTTGTAAAGCTA-PI93_seed_rep2
    #Glyma.01G000600                                87.79
    #Glyma.01G002100                                54.68
    #Glyma.01G003000                                13.45
    #CB:Z:CCGAAGCGTCCGTCGA-PI93_seed_rep2
    #Glyma.01G000600                                61.81
    #Glyma.01G002100                                50.79
    #Glyma.01G003000                                12.98
    
    # filter by tests
    if(!is.null(tests)){
        message(" - filter by differential testing ...")
        tests <- subset(tests, tests$qval < FDR)
        binary <- binary[rownames(binary) %in% rownames(tests),]
    }
    ##filter by the colSums and rowSums
    binary <- binary[,Matrix::colSums(binary)>0]
    binary <- binary[Matrix::rowSums(binary)>0,]
    pt <- pt[colnames(binary),]
    message(" - filtered zero-sum columns/rows of binary cell x gene matrix : ",ncol(binary)," | ", nrow(binary))
    
    # generalized additive model for logistic regression
    message(" - running generalized additive model ...")
    ##select 500 from the small to the max of the p.time
    newdat <- data.frame(p.time=seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500))
    #> head(newdat)
    #p.time
    #1 0.08025682
    #2 0.28049679
    #3 0.48073675
    #4 0.68097672
    
    ##analyze one each gene
    fit <- mclapply(seq(1:nrow(binary)),function(x){
        ##generate a df that includes the smooth data
        ##because pt is ordered by from small to max and binary is ordered exactly followed by the pt$trajectory
        df <- data.frame(acc=as.numeric(binary[x,]), p.time=pt$trajectory)
        ##it looks like if we directly plot there might be no obvious trend
        ##so in this case, we will fit a model 
        ##     acc     p.time
        #1  80.63 0.08025682
        #2  87.79 0.16051364
        #3  61.81 0.24077047
        #4  70.78 0.32102729
        #5  76.23 0.40128411
        #6 128.49 0.48154093
        
        mod <- gam(acc~s(p.time, bs="cr"), data=df)
        pred <- predict(mod, newdat, type="response")
        zscore <- (pred-mean(pred, na.rm=T))/sd(pred, na.rm=T)
    }, mc.cores=threads)
    names(fit) <- rownames(binary)
    fit <- do.call(rbind, fit)
    write.table(fit, file=paste0(output_dir,'/',prefix,".TF_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    # filter by variances
    ##rescale the x to be 0,1
    fit <- t(apply(fit, 1, function(x){
      scales::rescale(x, c(0,1))
    }))
    
    ##> fit[1:3,1:3]
    #1         2         3
    #Glyma.01G015900 1 0.9924176 0.9848399
    #Glyma.01G016600 1 0.9884767 0.9769618
    #Glyma.01G043300 1 0.9877446 0.9754974
    ##> dim(fit)
    #[1] 3490  500
    ##3490 genes and 500 cells
  
    # reformat output
    ##we will reorder 
    ##find the max order of each gene or say each row
    row.o <- apply(fit, 1, which.max)
    #length(row.o) is 3490 cells
    ##the rationale here:
    ##1. we want to find the order of max value of each gene
    ##2. re order the gene to allow the smallest ranking value to be the most left and the max ranking value to the right
    ##3. this is a way to plot
    
    #order(c(1,2,3,4,2,5),decreasing = F)
    ##[1] 1 2 5 3 4 6
    ##get the order of this list from small value to large value 
    
    fit <- fit[order(row.o, decreasing=F),]
    
    write.table(fit, file=paste0(output_dir,'/',prefix,".TF_pt_sorted.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    # plot
    message(" - plotting cell trajectory ...")
    ##we can choose other color plot
    cols <- colorRampPalette(c('grey80','grey75',brewer.pal(8,'BuPu')[2:8]))(100)
    #cols <- colorRampPalette(c("grey80", "grey75",brewer.pal(7, "YlGnBu")[2:7]))(100)
    pdf(paste0(output_dir,'/',prefix,".trajectoryTF.pdf"), width=10, height=10)
    heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
              scale="none", labRow = NA, labCol=NA, useRaster=T,
              ylab=paste("TFs", paste0("(n=",nrow(fit),")"), sep=" "),
              key = TRUE,                # show color key (legend bar)
              key.title = '',      # title for the legend
              key.xlab = "TF gene accessibility",      # x-axis label for legend
              keysize = 1,            # size of the legend bar
              density.info = "none"
    )
    
    ##create a small plot to create a legend bar
    dev.off()
    
    #pdf(paste0(output_dir,'/',prefix,".trajectoryTF_lb.pdf"), width=5, height=5)
    #plot.new()
    #add.color.bar(0.1, cols,lwd=10, title=NULL, lims=c(0,1),prompt=F,x=0.3,y=0.5,outline=TRUE,digits = 1)
    #dev.off()
    
    
    # return
    return(fit)
    
}
plotTrajGN <- function(obj, pt, cluster=1, prefix="temp", threads=1, tests=NULL, FDR=0.05){
    
    # subset traj
    message(" - cleaning input ...")
    pt <- subset(pt, pt[[target_cluster_col]] %in% cluster)
    pt <- pt[!is.na(pt$trajectory),]
    pt <- pt[order(pt$trajectory, decreasing=F),]
    obj$gn <- smooth.data(obj$gn, npcs=ncol(obj$d), rds=obj$d)
    binary <- obj$gn[,rownames(pt)]
    
    # filter by tests
    if(!is.null(tests)){
        message(" - filter by differential testing ...")
        tests <- subset(tests, tests$qval < FDR)
        binary <- binary[rownames(binary) %in% rownames(tests),]
    }
    binary <- binary[,Matrix::colSums(binary)>0]
    binary <- binary[Matrix::rowSums(binary)>0,]
    pt <- pt[colnames(binary),]
    message(" - filtered zero-sum columns/rows of binary cell x gene matrix : ",ncol(binary)," | ", nrow(binary))
    
    # generalized additive model for logistic regression
    message(" - running generalized additive model ...")
    newdat <- data.frame(p.time=seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500))
    fit <- mclapply(seq(1:nrow(binary)),function(x){
        df <- data.frame(acc=as.numeric(binary[x,]), p.time=pt$trajectory)
        mod <- gam(acc~s(p.time, bs="cr"), data=df)
        pred <- predict(mod, newdat, type="response")
        zscore <- (pred-mean(pred, na.rm=T))/sd(pred, na.rm=T)
    }, mc.cores=threads)
    names(fit) <- rownames(binary)
    fit <- do.call(rbind, fit)
    write.table(fit, file=paste0(output_dir,'/',prefix,".Gene_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    # filter by variances
    fit <- t(apply(fit, 1, function(x){
      scales::rescale(x, c(0,1))
    }))
    
    # reformat output
    row.o <- apply(fit, 1, which.max)
    fit <- fit[order(row.o, decreasing=F),]
    
    write.table(fit, file=paste0(output_dir,'/',prefix,".Gene_pt_sorted.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    
    # plot
    message(" - plotting cell trajectory ...")
    cols <- viridis(100)
    pdf(paste0(output_dir,'/',prefix,".trajectoryGene.pdf"), width=10, height=10)
    heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
              scale="none", labRow = NA, labCol=NA, useRaster=T,
              ylab=paste("Genes", paste0("(n=",nrow(fit),")"), sep=" "))
    dev.off()
    
    # return
    return(fit)
    
}

##updating 070221
##also keep the FDR < 0.05
plotTargetGN <- function(obj, pt, collect_targetGN_set_opt_dir,cluster=1, prefix="temp", threads=1, tests=NULL,target_gene_set=NULL, FDR=0.05){
  
  # subset traj
  message(" - cleaning input ...")
  pt <- subset(pt, pt[[target_cluster_col]] %in% cluster)
  pt <- pt[!is.na(pt$trajectory),]
  pt <- pt[order(pt$trajectory, decreasing=F),]
  obj$gn <- smooth.data(obj$gn, npcs=ncol(obj$d), rds=obj$d)
  
  ##updating 101224
  ##debug allow the same cells between obj$gn to the obj$d
  shared_cells <- intersect(colnames(obj$gn),rownames(pt))
  obj$gn <- obj$gn[,shared_cells]
  pt <- pt[shared_cells,]
  
  binary <- obj$gn[,rownames(pt)]
  
  # filter by tests
  if(!is.null(tests)){
    message(" - filter by differential testing ...")
    tests <- subset(tests, tests$qval < FDR)
    
    ##updating 071625
    if (nrow(tests) > 10000) {
      tests <- tests[order(tests$qval), ][1:10000, ]
    }
    
    binary <- binary[rownames(binary) %in% rownames(tests),]
  }
  
  ##filter by target gene sets
  if (!is.null(target_gene_set)){
    message(" - filter by target genes ...")
    binary <- binary[rownames(binary) %in% rownames(target_gene_set),]
  }
  
  binary <- binary[,Matrix::colSums(binary)>0]
  binary <- binary[Matrix::rowSums(binary)>0,]
  pt <- pt[colnames(binary),]
  message(" - filtered zero-sum columns/rows of binary cell x gene matrix : ",ncol(binary)," | ", nrow(binary))
  
  # generalized additive model for logistic regression
  message(" - running generalized additive model ...")
  newdat <- data.frame(p.time=seq(from=min(pt$trajectory, na.rm=T), to=max(pt$trajectory, na.rm=T), length.out=500))
  fit <- mclapply(seq(1:nrow(binary)),function(x){
    df <- data.frame(acc=as.numeric(binary[x,]), p.time=pt$trajectory)
    mod <- gam(acc~s(p.time, bs="cr"), data=df)
    pred <- predict(mod, newdat, type="response")
    zscore <- (pred-mean(pred, na.rm=T))/sd(pred, na.rm=T)
  }, mc.cores=threads)
  names(fit) <- rownames(binary)
  fit <- do.call(rbind, fit)
  write.table(fit, file=paste0(collect_targetGN_set_opt_dir,'/',prefix,".Gene_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  
  # filter by variances
  fit <- t(apply(fit, 1, function(x){
    scales::rescale(x, c(0,1))
  }))
  
  # reformat output
  row.o <- apply(fit, 1, which.max)
  fit <- fit[order(row.o, decreasing=F),]
  write.table(fit, file=paste0(collect_targetGN_set_opt_dir,'/',prefix,".Gene_pt_sorted.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  
  
  #library(scales)
  #tb <- read.table('Cotyledon_Adaxial_Parenchyma.Gene_pt.txt')
  #fit <- tb
  #prefix <- 'Cotyledon_Adaxial_Parenchyma'
  #fit <- t(apply(fit, 1, function(x){
  #  rescale(x, c(0,1))
  #}))
  #row.o <- apply(fit, 1, which.max)
  #fit <- fit[order(row.o, decreasing=F),]
  #write.table(fit, file=paste0(prefix,".Gene_pt_sorted.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  
  
  
  # plot
  message(" - plotting cell trajectory ...")
  cols <- viridis(100)
  cols <- colorRampPalette(c('grey80','grey75',brewer.pal(8,'YlOrRd')[2:8]))(100)
  pdf(paste0(collect_targetGN_set_opt_dir,'/',prefix,".trajectoryGene.pdf"), width=10, height=10)
  heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
            scale="none", labRow = NA, labCol=NA, useRaster=T,
            ylab=paste("Genes", paste0("(n=",nrow(fit),")"), sep=" "),
            key = TRUE,                # show color key (legend bar)
            key.title = '',      # title for the legend
            key.xlab = "Gene accessibility",      # x-axis label for legend
            keysize = 1,            # size of the legend bar
            density.info = "none"
  )
  dev.off()
  
  #pdf(paste0(collect_targetGN_set_opt_dir,'/',prefix,".trajectoryGene_lb.pdf"), width=5, height=5)
  #plot.new()
  #add.color.bar(0.1, cols,lwd=10, title=NULL, lims=c(0,1),prompt=F,x=0.3,y=0.5,outline=TRUE,digits = 1)
  #dev.off()
  
  # return
  return(fit)
  
}

plotSites  <- function(acr, mt, tf, x){
    
    # get sites
    acr.1 <- as.numeric(acr[which(rownames(acr) %in% x[1]),])
    mt.1 <- as.numeric(mt[which(rownames(mt) %in% x[2]),])
    tf.1 <- as.numeric(tf[which(rownames(tf) %in% x[3]),])
    xvals <- seq(from=0, to=100, length.out=length(acr.1))
    
    # plot
    plot(xvals, acr.1, type="l", col="darkorchid", lwd=2)
    lines(xvals, mt.1, type="l", col="darkorange", lwd=2)
    lines(xvals, tf.1, type="l", col="dodgerblue", lwd=2)
    legend("topright", legend=x, pch=16, col=c("darkorchid", "darkorange", "dodgerblue"))
}
findCor    <- function(acr, mt, tf, x){
    acr.use <- as.numeric(acr[x,])
    best.mt <- cor(t(mt), acr.use, method="pearson")
    best.tf <- cor(t(tf), acr.use, method="pearson")
    top.mt <- rownames(best.mt)[which.max(best.mt[,1])]
    top.tf <- rownames(best.tf)[which.max(best.tf[,1])]
    out <- c(x, top.mt, top.tf)
    return(out)
}

###################################################################################################
# load and process raw data -----------------------------------------------------------------------
###################################################################################################
message(paste0('analyze ',prefix))


ipt_object <- readRDS(ipt_object_fl)

##updating 062725 we will reload everything if users want to focuse on one more specific items 
#if (file.exists(paste0(output_dir,'/',prefix,".pseudotime_allgene.RDS"))){
#  obj <- readRDS(paste0(output_dir,'/',prefix,".pseudotime_allgene.RDS"))
#}else{
obj <- loadData_option_version(ipt_object,sparse, motif, gene, tfgene, doL2=0, traj,open_equal_cellnum_diff_clusters,column=target_cluster_col)
saveRDS(obj, file=paste0(output_dir,'/',prefix,".pseudotime_allgene.RDS"))
#}

#obj <- loadData_c(sparse,meta,svd,motif,doL2=0, traj)
#saveRDS(obj, file=paste0(prefix,".pseudotime.rds"))

if (file.exists(paste0(output_dir,'/',prefix,".trajectory.txt"))){
  out <- read.delim(paste0(output_dir,'/',prefix,".trajectory.txt"))
  plotPT(out, output_dir,prefix,subsetCluster=LC, smoother=5, t.id="trajectory", umap1=umap1, umap2=umap2)
}else{
  # iterate over types of trajectories
  obj$b[,column] <- as.character(obj$b[,column])
  out <- calcPseudo(obj$b, obj$d, column=column, traj=traj, cell.dist1=0.95, cell.dist2=0.95)
  # plot UMAP traj
  #pdf(paste0(output_dir,'/',prefix,".trajectory.pdf"), width=10, height=5)
  plotPT(out, output_dir,prefix,subsetCluster=LC, smoother=5, t.id="trajectory", umap1=umap1, umap2=umap2)
  #dev.off()
  write.table(out, file=paste0(output_dir,'/',prefix,".trajectory.txt"), quote=F, row.names=T, col.names=T, sep="\t")
}


# find significant peaks across pseudotime
if (sparse != 'na'){
  
  message(" - finding significant ACRs across pseudotime ...")
  
  if (file.exists(paste0(output_dir,'/',prefix,".ACR.diffTestsQVALs.txt"))) {
    diff.peaks <- read.delim(paste0(output_dir,'/',prefix,".ACR.diffTestsQVALs.txt"))
    
    ##uddating 060721 add the top ranking diff ACRs
    diff.peaks <- diff.peaks[order(diff.peaks$qval,decreasing=F),]
    diff.peaks <- diff.peaks[1:input_number_top_plot_ACR_peudo,]
    
  }else{
    diff.peaks <- sigPseudo2(obj, out, threads=threads, type="ACRs")
    write.table(diff.peaks, file=paste0(output_dir,'/',prefix,".ACR.diffTestsQVALs.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    
    ##uddating 060721 add the top ranking diff ACRs
    diff.peaks <- diff.peaks[order(diff.peaks$qval,decreasing=F),]
    diff.peaks <- diff.peaks[1:input_number_top_plot_ACR_peudo,]
    
  }
}

# find significant TFs across pseudotime
if (tfgene != 'na'){
  message(" - finding significant TFs across pseudotime ...")
  if (file.exists(paste0(output_dir,'/',prefix,".TFs.diffTestsQVALs.txt"))) {
    diff.TFs <- read.delim(paste0(output_dir,'/',prefix,".TFs.diffTestsQVALs.txt"))
  }else{
    diff.TFs <- sigPseudo2(obj, out, threads=threads, type="TFs")
    write.table(diff.TFs, file=paste0(output_dir,'/',prefix,".TFs.diffTestsQVALs.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  }
}

# find significant Genes across pseudotime
if (gene != 'na'){
  message(" - finding significant Genes across pseudotime ...")
  if (file.exists(paste0(output_dir,'/',prefix,".Genes.diffTestsQVALs.txt"))) {
      diff.Gns <- read.delim(paste0(output_dir,'/',prefix,".Genes.diffTestsQVALs.txt"))
  }else{
    diff.Gns <- sigPseudo2(obj, out, threads=threads, type="Genes")
    write.table(diff.Gns, file=paste0(output_dir,'/',prefix,".Genes.diffTestsQVALs.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  }
}

# find significant motifs across pseudotime
if (motif != 'na'){
  message(" - finding significant MTs across pseudotime ...")
  if (file.exists(paste0(output_dir,'/',prefix,".Motifs.diffTestsQVALs.txt"))) {
    diff.MTs <- read.delim(paste0(output_dir,'/',prefix,".Motifs.diffTestsQVALs.txt"))
  }else{
    diff.MTs <- sigPseudo2(obj, out, threads=threads, type="MTs")
    write.table(diff.MTs, file=paste0(output_dir,'/',prefix,".Motifs.diffTestsQVALs.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  }
}

##open them when necessary
##updating 081121 we will decide which parts we will open
# plot heatmap traj
if (openACR == 'yes'){
  acr <- plotTrajHM(obj, out, cluster=LC, prefix=prefix, threads=threads, featureMin=featureMin, tests=diff.peaks, FDR=0.05)
}
if (openMT == 'yes'){
  mt <- plotTrajMT(obj, out, cluster=LC, prefix=prefix, threads=threads, tests=diff.MTs, FDR=0.05)
}
if (openTF == 'yes'){
  tf <- plotTrajTF(obj, out, cluster=LC, prefix=prefix, threads=threads, tests=diff.TFs, FDR=0.05)
}
#if (openGN == 'yes'){
#  gn <- plotTrajGN(obj, out, cluster=LC, prefix=prefix, threads=threads, tests=diff.Gns, FDR=0.05)
#}
if (openGN == 'yes'){
  #tg_set <- read.delim(target_gene_set,row.names=1)
  collect_targetGN_set_opt_dir <- output_dir
  tgn <- plotTargetGN(obj, out, collect_targetGN_set_opt_dir,cluster=LC, prefix=prefix, threads=threads, tests=diff.Gns, target_gene_set = NULL,FDR=0.05)
}
if (openTGN == 'yes'){
  collect_targetGN_set_opt_dir <- output_dir
  tg_set <- read.delim(target_gene_set,row.names=1)
  tgn <- plotTargetGN(obj, out, collect_targetGN_set_opt_dir,cluster=LC, prefix=prefix, threads=threads, tests=diff.Gns, target_gene_set = tg_set,FDR=0.05)
}



