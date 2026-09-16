# run Socrates on merged socrates object #

##updating 013826 we will set an option to rename library
##updating 031225 we will not build the subclustering as it will complex the processing, high res will solve the issue
##updating 031125 we will open the win sparse in order to perform the subclustering
##updating 020225 we will add the parameter of harmony function
##updating 010525 we will return a prefix.atac.soc.qs2 file
##updating 011124 we will set an option to remove the temp file
##updating 041722 we need to modify npc number detect doublet
##updating 121221 we set an argument to add a harmony function
##updating 112521 change the cex of plotting UMAP
##updating 112021 add the res as the out
##updating 111621 we will set the normalized way
##updating 111121 we will check whether the rds file is existing
##updating 110921 we will set the parameters as arguments
##udpating 102221 we will save a win matrix since later we will use this to do the subclustereing
##updating 101821 we will filter windows first \
##updating 101521 no harmony information here
##updating 101421 add the detect doublet function
##updating 101021 set the NMF and SVD
##udpating 100821 we will update the the clustering
##updating 100621 we will add the harmony 
##updating we will directly use the function within in the Socrates

# libraries
#library(Socrates)
#library(harmony)
#library(RcppML)
##ml rgdal/1.4-8-foss-2019b-R-4.0.0
#.libPaths(c("/apps/eb/rgdal/1.4-8-foss-2019b-R-4.0.0", .libPaths()))
#.libPaths(c("/home/hy17471/R/x86_64-pc-linux-gnu-library/4.0", .lƒibPaths()))


library("here")
library(devtools)
library(Seurat)
library(harmony)
library(RcppML)

#load_all('/home/jpm73279/Socrates')




# load rds files and pre-processing --------------------------------------
##create the Socrates obj
args <- commandArgs(T)
# vars
#input_scobj_dir  <- as.character(args[1]) ##opt_combine_fltmeta.sparse
##analyze in the current dir cp the filtered.soc.qs2 to the current dir

input_obj_dir <- as.character(args[1])
##allow the filtered obj into the input object directory

##updating 111124
path_to_preload_R_script <- as.character(args[2])

##updating 111124
config <- as.character(args[3])

input_output_dir <- as.character(args[4])

source(config)

##set parameters
##set the number of var of the selected windows
#num_var <- as.numeric(args[2])

if (num_var_final == 0){
  num_var = NULL
}else{
  num_var = num_var_final
}

##set the rdType
#rdType <- as.character(args[3]) ##NMF or SVD
##set the numPCs
#numPCs <- as.numeric(args[4])

##updating 111621
#normalize_way <- as.character(args[5]) ##tfidf or regmodel

##updating 112021
#res_clust <- as.numeric(args[6])

##updating 121221
#do_harmony <- as.character(args[7])

##updating 121221
#do_clust_saving_step <- as.character(args[8])

##updating 111124
#remove_temp <- as.character(args[8])


##load all the required R script under Socrates2
load_all(path_to_preload_R_script)

##load the configure file
source(config)


##updating 091326 remove the user-specified clusters and re-cluster
##the labels come from the final object of a previous run (written by clustering.py into the config file)
if (exists('remove_clusters_final') && as.character(remove_clusters_final) != 'no'){
  remove_clusters_vec <- trimws(strsplit(as.character(remove_clusters_final), ',')[[1]])
  remove_clusters_vec <- remove_clusters_vec[remove_clusters_vec != '']
  do_remove_clusters <- length(remove_clusters_vec) > 0
}else{
  remove_clusters_vec <- NULL
  do_remove_clusters <- FALSE
}

##defensive guard: this flag is meaningless together with only_cluster == 'yes'
if (do_remove_clusters && only_cluster == 'yes'){
  stop(' -remove_clusters cannot be combined with -only_cluster yes: the merged object has to be re-processed.')
}

##namespace the intermediate files, so that a run with -remove_clusters never overwrites
##the caches produced by a run without it (this keeps the removal reversible)
rmc_suffix <- ''
if (do_remove_clusters && exists('remove_clusters_suffix')){
  rmc_suffix <- as.character(remove_clusters_suffix)
}


##set out name
out <- paste0('opt_',normalize_way,'_',rdType,'_',as.character(numPCs),'PCs_','win',as.character(num_var),rmc_suffix)
aft_clust_out <- paste0('opt_',normalize_way,'_',rdType,'_',as.character(numPCs),'PCs_','win',as.character(num_var),'_res',as.character(res_clust),rmc_suffix)

##merget the object
mergeSocObjects <- function(obj.list){
  
  # functions
  .merge.sparse <- function(cnt.list) {
    
    cnnew <- character()
    rnnew <- character()
    x <- vector()
    i <- numeric()
    j <- numeric()
    
    for (M in cnt.list) {
      
      cnold <- colnames(M)
      rnold <- rownames(M)
      
      cnnew <- union(cnnew,cnold)
      rnnew <- union(rnnew,rnold)
      
      cindnew <- match(cnold,cnnew)
      rindnew <- match(rnold,rnnew)
      ind <- summary(M)
      i <- c(i,rindnew[ind[,1]])
      j <- c(j,cindnew[ind[,2]])
      x <- c(x,ind[,3])
    }
    
    sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
  }
  
  # separate counts and meta
  counts <- lapply(obj.list, function(x){
    x$counts
  })
  counts <- .merge.sparse(counts)
  
  # meta
  metas <- lapply(obj.list, function(x){
    x$meta
  })
  metas <- do.call(rbind, metas)
  rownames(metas) <- metas$cellID
  
  
  ##updating 012826
  ##if rename_library_final is no we will open 
  if (rename_library_final == 'no'){
  
    ##updating 061925
    ##we will still open this one
    metas$library <- data.frame(do.call(rbind, strsplit(rownames(metas), "-")))[,2]
    
    ##udpating 021925 we will mute this library as we have developed this library before
    #metas$library <- data.frame(do.call(rbind, strsplit(rownames(metas), "-")))[,2]
  }
  
  metas <- metas[colnames(counts),]
  
  # new object
  new.obj <- list(counts=counts, meta=metas)
  return(new.obj)
  
}
plotComponents <- function(obj, comp=1, embedding="NMF", umapslot="Clusters_NMF"){
  
  cols <- colorRampPalette(c("grey80", "grey80", brewer.pal(9, "PuRd")[2:9]))(100)
  df <- obj[[umapslot]]
  comps <- obj[[embedding]]
  ids <- intersect(rownames(df), rownames(comps))
  df <- df[ids,]
  comps <- comps[ids,]
  row.o <- order(comps[,comp], decreasing=F)
  df <- df[row.o,]
  comps <- comps[row.o,]
  plot(df$umap1, df$umap2, col=cols[cut(comps[,comp], breaks=101)], 
       pch=16, 
       cex=0.2,
       main=paste0(embedding,"_component_",comp))
  
  
}


##decide whether we will only perform the clustering 
if (only_cluster == 'no'){
  
  ##updating 111121
  if (file.exists(paste0(input_output_dir,'/',out,".aft_rD_UMAP_rDb.qs2"))){
    message (' - aft_rD_UMAP_rDb object is existed, and we will directly read it')
    soc.obj <- qs2::qs_read(paste0(input_output_dir,'/',out,".aft_rD_UMAP_rDb.qs2"))
  }else{
  
    if (file.exists(paste0(input_output_dir,'/',out,".aft_rD_UMAP.qs2"))){
      message (' - project_UMAP object is existed, and we will directly read it')
      soc.obj <- qs2::qs_read(paste0(input_output_dir,'/',out,".aft_rD_UMAP.qs2"))
      
      #############################
      message(" - detect doublets")
      soc.obj <- detectDoublets(
        obj = soc.obj,
        nTrials = 5,
        nSample = 1000,
        k = 10,
        #n.pcs = 50,
        threads = 1
      )
      
      if (remove_temp == 'no'){
        qs2::qs_save(soc.obj,file=paste0(input_output_dir,'/',out,".aft_rD_UMAP_rDb.qs2"))
      }
      
    }else{
    
      if (file.exists(paste0(input_output_dir,'/',out,".aft_rD.qs2"))){
        message (' - reduceDim object is existed, and we will directly read it')
        soc.obj <- qs2::qs_read(paste0(input_output_dir,'/',out,".aft_rD.qs2"))
        
        ##########################
        message(' - project_UMAP')
        soc.obj <- projectUMAP(soc.obj, 
                               verbose=T,
                               k.near=50,
                               m.dist=0.01,
                               #svd_slotName=rdType,  ##do not need this one
                               umap_slotName=paste0('UMAP_',rdType))
        
        if (remove_temp == 'no'){
          qs2::qs_save(soc.obj, file=paste0(input_output_dir,'/',out,".aft_rD_UMAP.qs2"))
        }
        
        #############################
        message(" - detect doublets")
        soc.obj <- detectDoublets(
          obj = soc.obj,
          nTrials = 5,
          nSample = 1000,
          k = 10,
          #n.pcs = 50,
          threads = 1
        )
        
        if (remove_temp == 'no'){
          qs2::qs_save(soc.obj,file=paste0(input_output_dir,'/',out,".aft_rD_UMAP_rDb.qs2"))
        }
        
      }else{
        
        message(' - no previous analyses have been done. We need start from the merged obj')
        # load rds files and pre-processing --------------------------------------
        message(' - merge all the objects')
        if(file.exists(paste0(input_output_dir,'/',"opt_merged_obj.qs2"))){
          message (' - merged object is existed, and we will directly read it')
          soc.obj <- qs2::qs_read(paste0(input_output_dir,'/',"opt_merged_obj.qs2"))
        }else{
          message ('- merged object is not existed, and we will create it')
          
          dat <- list.files(path= input_obj_dir,pattern = '.filtered.soc.qs2')
          #dat <- list.files(pattern=paste0(input_obj_dir,'/','*.filtered.soc.qs2'))
          dat <- lapply(dat, function(x){
            obj <- qs2::qs_read(paste0(input_obj_dir,'/',x))
            
            ##udpating 010426
            obj$counts <- obj$counts[Matrix::rowMeans(obj$counts > 0)>minimum_matrix_rowmean_final,]
            
            #obj$counts <- obj$counts[Matrix::rowMeans(obj$counts > 0)>0.01,]
            obj$counts <- obj$counts[,Matrix::colSums(obj$counts)>0]
            obj$meta <- obj$meta[colnames(obj$counts),]
            
            
            ##updating 012826 open when rename_library_final is yes
            if (rename_library_final == 'yes'){
            
              ##updating 061925 close the following
              ##updating 021925
              ##modify the meta file add the lib to the cell barcode
              libnm <- gsub('.filtered.soc.qs2','',x)
              obj_meta <- obj$meta
              obj_counts <- obj$counts
              
              obj_meta$library <- libnm
              obj_meta$cellID <- gsub('-.+','',obj_meta$cellID)
              obj_meta$cellID <- paste0(obj_meta$cellID,'-',libnm)
              rownames(obj_meta) <- obj_meta$cellID
              
              colnames(obj_counts) <- gsub('-.+','',colnames(obj_counts))
              colnames(obj_counts) <- paste0(colnames(obj_counts),'-',libnm)
              
              obj$counts <- obj_counts
              obj$meta <- obj_meta
              ####################
            }
              
              
            return(obj)
          })
          
          soc.obj <- mergeSocObjects(dat)
          qs2::qs_save(soc.obj, file=paste0(input_output_dir,'/',"opt_merged_obj.qs2"))

        }

        ##updating 091326 remove the user-specified clusters before the dimension reduction
        if (do_remove_clusters){

          ref_fl <- if (exists('remove_clusters_ref_final') && as.character(remove_clusters_ref_final) != ''){
            as.character(remove_clusters_ref_final)
          }else{
            paste0(input_output_dir,'/',input_prefix,'.atac.soc.qs2')
          }
          if (!file.exists(ref_fl)){
            stop(' -remove_clusters: the reference object ', ref_fl, ' does not exist. Please run the clustering ',
                 'once with the same -o, -prefix and parameters before using -remove_clusters.')
          }

          ref_obj  <- qs2::qs_read(ref_fl)
          ref_meta <- ref_obj$meta

          if (!('LouvainClusters' %in% colnames(ref_meta))){
            stop(' -remove_clusters: no LouvainClusters column found in ', ref_fl, '.')
          }
          if (!('library' %in% colnames(ref_meta))){
            stop(' -remove_clusters: no library column found in ', ref_fl, '.')
          }

          ref_clust <- as.character(ref_meta$LouvainClusters)
          names(ref_clust) <- rownames(ref_meta)
          if (anyDuplicated(names(ref_clust))){
            warning(' -remove_clusters: duplicated cell IDs in the reference metadata; only the first occurrence is used.')
            ref_clust <- ref_clust[!duplicated(names(ref_clust))]
          }

          ##guard 1: the requested labels must exist in the reference
          requested <- intersect(remove_clusters_vec, unique(ref_clust))
          if (length(requested) == 0){
            stop(' -remove_clusters: none of the requested clusters (',
                 paste(remove_clusters_vec, collapse=','), ') exist in ', ref_fl,
                 '. Available clusters: ', paste(sort(unique(ref_clust)), collapse=','))
          }
          if (length(requested) < length(remove_clusters_vec)){
            warning(' -remove_clusters: cluster(s) ',
                    paste(setdiff(remove_clusters_vec, requested), collapse=','),
                    ' were not found in the reference and are ignored.')
          }

          ##map the reference labels onto the cells of the merged object
          all_cells  <- colnames(soc.obj$counts)
          bad_cells  <- intersect(names(ref_clust)[ref_clust %in% requested], all_cells)
          n_total    <- length(all_cells)
          n_labelled <- sum(all_cells %in% names(ref_clust))

          message(' - removing cluster(s) ', paste(requested, collapse=','),
                  ': ', length(bad_cells), ' cell(s) matched in the merged object')
          message('   * the reference covers ', n_labelled, ' of ', n_total,
                  ' cell(s); the ', n_total - n_labelled,
                  ' unlabelled cell(s) have no cluster assignment and are kept')

          ##guard 2: a total mismatch means the reference does not belong to this merged object
          if (length(bad_cells) == 0){
            stop(' -remove_clusters: the requested clusters matched 0 cells in the merged object. ',
                 'The reference was probably produced from a different merged object ',
                 '(check -ipt_obj_dir, -rename_lib, -minimum_matrix_rowmean and -prefix).')
          }

          ##audit trail: which cells were removed and how they distribute across the libraries
          audit <- ref_meta[match(bad_cells, rownames(ref_meta)), , drop=FALSE]
          write.table(audit, file=paste0(input_output_dir,'/',out,'_removed_cells.txt'),
                      quote=F, row.names=T, col.names=T, sep="\t")
          write.table(as.data.frame.matrix(table(audit$LouvainClusters, audit$library)),
                      file=paste0(input_output_dir,'/',out,'_removed_cluster_library_composition.txt'),
                      quote=F, row.names=T, col.names=T, sep="\t")

          ##guard 3: never remove every cell
          if (length(bad_cells) >= length(all_cells)){
            stop(' -remove_clusters: removing cluster(s) ', paste(requested, collapse=','),
                 ' would remove every cell in the merged object.')
          }

          keep_cells <- all_cells[!(all_cells %in% bad_cells)]
          soc.obj$counts <- soc.obj$counts[, keep_cells, drop=FALSE]
          ##re-index by column name so that counts and meta can never drift apart
          soc.obj$meta <- soc.obj$meta[colnames(soc.obj$counts), , drop=FALSE]
          if (anyNA(rownames(soc.obj$meta))){
            stop(' -remove_clusters: counts and meta could not be aligned after filtering.')
          }

          ##guard 4: harmony needs at least two libraries
          if (do_harmony == 'yes' && length(unique(soc.obj$meta$library)) < 2){
            stop(' -remove_clusters: only one library remains after filtering, so harmony cannot be run. ',
                 'Please re-run with -harmony no.')
          }
        }

        write.table(soc.obj$meta, file=paste0(input_output_dir,'/',out,'_combine_metadata.txt'), quote=F, row.names=T, col.names=T, sep="\t")
        
        # get per cell feature counts --------------------------------------------
        cell.counts <- log10(Matrix::colSums(soc.obj$counts))  # count number of features with Tn5 insertions per cell
        cell.counts.z <- as.numeric(scale(cell.counts)) # convert features counts into Z-scores
        cell.counts.threshold <- max(c((10^cell.counts[cell.counts.z < -1]), 1000)) # minimum feature counts (greater of 1 std or 1000)  #filter cells
        
        # clean sparse counts matrix ---------------------------------------------
        
        ##updating 111224
        ##here we will directly set the number of minimum cell count other than consideration of std
        
        soc.ori <- soc.obj
        soc.obj <- cleanData(soc.obj, 
                             min.c=min_cell_final,
                             #min.c=cell.counts.threshold,  # minimum number of accessible features per cell   ##filter cells 
                             min.t=min_ft_freq_final,
                             max.t=max_ft_freq_final,
                             #min.t=0.005, # minimum feature frequency across cells  ##filter windows
                             #max.t=0.005,  # maximum feature frequency across cells ##filter windows
                             verbose=T)
        
        soc.ori$meta$anchors <- as.factor(ifelse(rownames(soc.ori$meta) %in% colnames(soc.obj$counts), "anchor", "projection"))
        
        
        if (normalize_way == 'tfidf') { 
          # normalize with TFIDF ---------------------------------------------------
          soc.obj <- tfidf(soc.obj)
        }else{
          if (normalize_way == 'regmodel'){
            soc.obj <- regModel(soc.obj)
          }
        }
        
        
        # project with NMF -------------------------------------------------------
        #######################
        message(' - reduceDim')
        soc.obj <- reduceDims(soc.obj,
                              method=rdType, 
                              n.pcs=numPCs,  ##50
                              num.var=num_var, ##50000
                              verbose=T,
                              scaleVar=T,
                              doSTD=T,
                              doL1=F,
                              doL2=F,
                              refit_residuals=F)
        
        if (remove_temp == 'no'){
          qs2::qs_save(soc.obj, file=paste0(input_output_dir,'/',out,".aft_rD.qs2"))
        }
          
          
        # L2 norm loadings
        soc.obj$PCA <- t(apply(soc.obj$PCA, 1, function(x) x/sqrt(sum(x^2))))
        
        # reduce to 2-dimensions with UMAP ---------------------------------------
        ##########################
        message(' - project_UMAP')
        soc.obj <- projectUMAP(soc.obj, 
                               verbose=T,
                               k.near=50,
                               m.dist=0.01,
                               #svd_slotName=rdType,  ##do not need this one
                               umap_slotName=paste0('UMAP_',rdType))
        
        if (remove_temp == 'no'){
          qs2::qs_save(soc.obj, file=paste0(input_output_dir,'/',out,".aft_rD_UMAP.qs2"))
        }
          
        #############################
        message(" - detect doublets")
        soc.obj <- detectDoublets(
                                obj = soc.obj,
                                nTrials = 5,
                                nSample = 1000,
                                k = 10,
                                #n.pcs = 50,
                                threads = core_doublet_final
        )
        
        #if (remove_temp == 'no'){
          qs2::qs_save(soc.obj,file=paste0(input_output_dir,'/',out,".aft_rD_UMAP_rDb.qs2"))
        #}
      }
    }
  }

}
  
   
##updating 121221
##check whether we will do the harmony
if (do_harmony == 'yes'){
  
  if (only_cluster == 'no'){
    
    nmf.rd <- soc.obj$PCA
    meta_dt <- soc.obj$meta
    
    harmony_embeddings <- HarmonyMatrix(
      nmf.rd, meta_dt, 'library',
      lambda=lambda_val,
      theta=theta_val,
      sigma=sigma_val,
      do_pca = FALSE
    )
    ##add the harmony_embeddings to the obj
    soc.obj[['HM_EMBS']] <- harmony_embeddings
    
  
    message(' - project_UMAP')
    soc.obj <- projectUMAP(soc.obj, 
                           verbose=T,
                           k.near=50,
                           m.dist=0.01,
                           svd_slotName='HM_EMBS',
                           #svd_slotName=rdType,  ##do not need this one
                           umap_slotName=paste0('UMAP_',rdType))
    #if (remove_temp == 'no'){
      qs2::qs_save(soc.obj, file=paste0(input_output_dir,'/',out,".aft_rD_rHM_UMAP.qs2"))
    #}
    
  }
    
  soc.obj <- qs2::qs_read(paste0(input_output_dir,'/',out,".aft_rD_rHM_UMAP.qs2"))
  # identify clusters using neighborhood graph -----------------------------
  message(' - call clustering')
  soc.obj <- callClusters(soc.obj, 
                          res=res_clust, ##default is 0.3
                          verbose=T,
                          #svd_slotName=rdType, ##do not need this one
                          umap_slotName=paste0('UMAP_',rdType),
                          cluster_slotName=paste0('Clusters_',rdType),
                          cleanCluster=F,
                          cl.method=3,
                          e.thresh=5)
  
  #plot cluster membership on UMAP embedding ------------------------------
  pdf(paste0(input_output_dir,'/',aft_clust_out,'.aftHM.UMAP.clusters.pdf'), width=10, height=10)
  plotUMAP(soc.obj, cluster_slotName=paste0('Clusters_',rdType), cex=0.2)
  dev.off()
  pdf(paste0(input_output_dir,'/',aft_clust_out,".aftHM.UMAP.library.pdf"), width=10, height=10)
  plotUMAP(soc.obj, cluster_slotName=paste0('Clusters_',rdType), column="library", cex=0.2)
  dev.off()
  pdf(paste0(input_output_dir,'/',aft_clust_out,".aftHM.UMAP.nSites.pdf"), width=10, height=10)
  plotUMAP(soc.obj, cluster_slotName=paste0('Clusters_',rdType), column="log10nSites", cex=0.2)
  dev.off()
  
  if (remove_temp == 'no'){
    qs2::qs_save(soc.obj, file=paste0(input_output_dir,'/',aft_clust_out,".rHM.processed_aftplotting.qs2"))
  }
    
  # # output text files
  if (rdType == 'NMF'){
    nmf.meta <- soc.obj$Clusters_NMF
  }
  if (rdType == 'SVD') {
    nmf.meta <- soc.obj$Clusters_SVD
  }
  
  ##since we did not specify the slot name so we need to use the PCA
  nmf.rd <- soc.obj$HM_EMBS
  
  ##updating 010425 only store the meta and svd
  final_obj <- list(
    meta = nmf.meta,
    svd = nmf.rd
  )
  qs2::qs_save(final_obj,file=paste0(input_output_dir,'/',input_prefix,'.atac.soc.qs2'))
  
  #write.table(nmf.meta, file=paste0(input_output_dir,'/',aft_clust_out,'_rHM_metadata.txt'), quote=F, row.names=T, col.names=T, sep="\t")
  #write.table(nmf.rd, file=paste0(input_output_dir,'/',aft_clust_out,'_rHM_reduced_dimensions.txt'), quote=F, row.names=T, col.names=T, sep="\t")
  # # save data --------------------------------------------------------------
  #qs2::qs_save(soc.obj, file=paste0(input_output_dir,'/',aft_clust_out,".rHM.processed_final.qs2"))
  
  
  ##return a win matrix
  #mtx_m <- soc.obj$counts
  #ia <- as.data.frame(summary(mtx_m))
  #ia$i <- rownames(mtx_m)[as.numeric(ia$i)]
  #ia$j <- colnames(mtx_m)[as.numeric(ia$j)]
  #write.table(ia, file=paste0(input_output_dir,'/',aft_clust_out,'.rHM.win.sparse'), quote=F, row.names=F, col.names=F, sep="\t")
  
  
}else{

  
  soc.obj <- qs2::qs_read(paste0(input_output_dir,'/',out,".aft_rD_UMAP_rDb.qs2"))
  # identify clusters using neighborhood graph -----------------------------
  soc.obj <- callClusters(soc.obj, 
                          res=res_clust, ##default is 0.3
                          verbose=T,
                          #svd_slotName=rdType, ##do not need this one
                          umap_slotName=paste0('UMAP_',rdType),
                          cluster_slotName=paste0('Clusters_',rdType),
                          cleanCluster=F,
                          cl.method=3,
                          e.thresh=5)
  
  #plot cluster membership on UMAP embedding ------------------------------
  pdf(paste0(input_output_dir,'/',aft_clust_out,'.UMAP.clusters.pdf'), width=10, height=10)
  plotUMAP(soc.obj, cluster_slotName=paste0('Clusters_',rdType), cex=0.2)
  dev.off()
  pdf(paste0(input_output_dir,'/',aft_clust_out,".UMAP.library.pdf"), width=10, height=10)
  plotUMAP(soc.obj, cluster_slotName=paste0('Clusters_',rdType), column="library", cex=0.2)
  dev.off()
  pdf(paste0(input_output_dir,'/',aft_clust_out,".UMAP.nSites.pdf"), width=10, height=10)
  plotUMAP(soc.obj, cluster_slotName=paste0('Clusters_',rdType), column="log10nSites", cex=0.2)
  dev.off()
  
  if (remove_temp == 'no'){
    qs2::qs_save(soc.obj, file=paste0(input_output_dir,'/',aft_clust_out,".processed_aftplotting.qs2"))
  }
  # plot components
  #png(paste0(input_output_dir,'/',aft_clust_out,".components.png"), width=25, height=25, type="cairo", units="in", res=300)
  #layout(matrix(c(1:25), nrow=5, byrow=T))
  #for(i in 1:25){
    ##updating 101221
  #  plotComponents(soc.obj, embedding=rdType, comp=i)
  #}
  #dev.off()
  # # output text files
  if (rdType == 'NMF'){
    nmf.meta <- soc.obj$Clusters_NMF
  }
  if (rdType == 'SVD') {
    nmf.meta <- soc.obj$Clusters_SVD
  }
  
  ##since we did not specify the slot name so we need to use the PCA
  nmf.rd <- soc.obj$PCA
  
  #if (rdType == 'NMF'){
  #  nmf.rd <- soc.obj$NMF
  #}
  #if (rdType == 'SVD') {
  #  nmf.rd <- soc.obj$SVD
  #}
  
  ##updating 010425 only store meta and rd
  final_obj <- list(
    meta = nmf.meta,
    svd = nmf.rd
  )
  qs2::qs_save(final_obj,file=paste0(input_output_dir,'/',input_prefix,'.atac.soc.qs2'))
  

  #write.table(nmf.meta, file=paste0(input_output_dir,'/',aft_clust_out,'.metadata.txt'), quote=F, row.names=T, col.names=T, sep="\t")
  #write.table(nmf.rd, file=paste0(input_output_dir,'/',aft_clust_out,'.reduced_dimensions.txt'), quote=F, row.names=T, col.names=T, sep="\t")
  # # save data --------------------------------------------------------------
  #qs2::qs_save(soc.obj, file=paste0(input_output_dir,'/',aft_clust_out,".processed_final.qs2"))
  
  ##return a win matrix
  #mtx_m <- soc.obj$counts
  #ia <- as.data.frame(summary(mtx_m))
  #ia$i <- rownames(mtx_m)[as.numeric(ia$i)]
  #ia$j <- colnames(mtx_m)[as.numeric(ia$j)]
  #write.table(ia, file=paste0(input_output_dir,'/',aft_clust_out,'.win.sparse'), quote=F, row.names=F, col.names=F, sep="\t")

}


