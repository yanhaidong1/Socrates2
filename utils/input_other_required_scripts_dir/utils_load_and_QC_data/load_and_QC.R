# run Socrates on merged socrates object #

##updating 020125 add the max reads
##updating 100821 we will do Socrates from very beginning
##updating 100621 we will add the harmony 
##updating we will directly use the function within in the Socrates

# libraries
library("here")
library(devtools)
library(Seurat)
library(harmony)
library(RcppML)

#load_all('/scratch/hy17471/software/Socrates/R/')

args <- commandArgs(T)
# vars
bed  <- as.character(args[1]) 
chr   <- as.character(args[2])
ann <- as.character(args[3])
config <- as.character(args[4])
output_dir <- as.character(args[5])

open_build_object_final <- as.character(args[6])
open_find_cells_step_final <- as.character(args[7])
open_is_cells_step_final <- as.character(args[8])

path_to_preload_R_script <- as.character(args[9])

load_all(path_to_preload_R_script)

##load all the parameters
source(config)

##set the condition to run each step
if (open_build_object_final == 'yes'){

  ##load the obj
  obj <- loadBEDandGenomeData(bed, ann, chr)
  ##count organellar reads
  obj <- countRemoveOrganelle(obj, org_scaffolds=c(Pt, Mt), remove_reads=T)  ##ChrPt and ChrMt for the rice
  ##call peaks 
  obj <- callACRs(obj, genomesize=genomesize, 
                  shift= macs2_shift_final, 
                  extsize= macs2_extsize_final,
                  fdr= macs2_fdr_final,
                  output=paste0(out,"_peaks"), 
                  tempdir=paste0(output_dir,'/',out, '_macs2_temp'), 
                  verbose=T)
  ##build meta data
  obj <- buildMetaData(obj,
                       tss.window= tss_window_size_final, 
                       verbose=TRUE)
  saveRDS(obj,paste0(output_dir,'/',out,'.raw_obj.rds'))
  meta_fl <-obj$meta
  write.table(meta_fl, file=paste0(output_dir,'/',out, ".raw_obj.metadata.txt"), quote=F, row.names=T, col.names=T, sep="\t")
}


if (open_find_cells_step_final == 'yes'){
  
  if (open_build_object_final == 'yes'){

    message(" - finding cells ...")
    
    ##udpating 031925
    if (open_knee_plot_filter_cell == 'yes'){
      obj <- findCells(obj, 
                       output_dir,
                       doplot=T,
                       min.cells= min_cell_val_final,
                       max.cells= max_cell_val_final,
                       set.tn5.cutoff=NULL,
                       min.tn5=min_tn5_val,
                       max.tn5=max_tn5_val,
                       org.filter.thresh=organelle_filter_cutoff_final,
                       tss.min.freq=tss_min_freq_val_final,
                       tss.z.thresh=tss_z_thresh_final,
                       frip.min.freq=frip_min_freq_final,
                       frip.z.thresh=frip_z_thresh_final,
                       prefix=out)
    }else{
      obj <- findCells_no_knee_plot_filter(obj, 
                       output_dir,
                       doplot=T,
                       min.cells= min_cell_val_final,
                       max.cells= max_cell_val_final,
                       set.tn5.cutoff=NULL,
                       min.tn5=min_tn5_val,
                       max.tn5=max_tn5_val,
                       org.filter.thresh=organelle_filter_cutoff_final,
                       tss.min.freq=tss_min_freq_val_final,
                       tss.z.thresh=tss_z_thresh_final,
                       frip.min.freq=frip_min_freq_final,
                       frip.z.thresh=frip_z_thresh_final,
                       prefix=out)
      
    }
    
    
    
    obj <- generateMatrix(obj,
                          filtered=T, 
                          windows=win_size, 
                          peaks=F, 
                          verbose=T)
  }else{
    
    message(" - finding cells ...")
    message(" - load the raw object")
    obj <- readRDS(paste0(output_dir,'/',out,'.raw_obj.rds'))
    #obj <- findCells(obj, 
    #                 output_dir,
    #                 doplot=T,
    #                 max.cells=16000,
    #                 set.tn5.cutoff=1000,
    #                 frip.min.freq=0.2,
    #                 frip.z.thresh=1,
    #                 prefix=out)
    
    ##udpating 031925
    if (open_knee_plot_filter_cell == 'yes'){
      
      obj <- findCells(obj, 
                       output_dir,
                       doplot=T,
                       min.cells= min_cell_val_final,
                       max.cells= max_cell_val_final,
                       set.tn5.cutoff=NULL,
                       min.tn5=min_tn5_val,
                       max.tn5=max_tn5_val,
                       org.filter.thresh=organelle_filter_cutoff_final,
                       tss.min.freq=tss_min_freq_val_final,
                       tss.z.thresh=tss_z_thresh_final,
                       frip.min.freq=frip_min_freq_final,
                       frip.z.thresh=frip_z_thresh_final,
                       prefix=out)
    }else{
      
      obj <- findCells_no_knee_plot_filter(obj, 
                       output_dir,
                       doplot=T,
                       min.cells= min_cell_val_final,
                       max.cells= max_cell_val_final,
                       set.tn5.cutoff=NULL,
                       min.tn5=min_tn5_val,
                       max.tn5=max_tn5_val,
                       org.filter.thresh=organelle_filter_cutoff_final,
                       tss.min.freq=tss_min_freq_val_final,
                       tss.z.thresh=tss_z_thresh_final,
                       frip.min.freq=frip_min_freq_final,
                       frip.z.thresh=frip_z_thresh_final,
                       prefix=out)
      
    }
    
    
    
    obj <- generateMatrix(obj,
                          filtered=T, 
                          windows=win_size, 
                          peaks=F, 
                          verbose=T)
  }
  
  saveRDS(obj,paste0(output_dir,'/',out,'.findcell_obj.rds'))
  
}

if (open_is_cells_step_final == 'yes'){
  
  if (open_find_cells_step_final == 'yes'){
    
    message(" - filter cell by comparing background ...")
    obj <- isCell(obj,
                  num.test=num_test_val_final,
                  num.tn5=num_tn5_val_final,
                  num.ref=num_ref_val_final,
                  background.cutoff=min_tn5_val+1000,
                  min.pTSS=0.2,
                  min.FRiP=0.2,
                  min.pTSS.z= -2,
                  min.FRiP.z= -2,
                  verbose=F)

  }else{
    
    message(" - filter cell by comparing background ...")
    message(" - laod the findcell object ...")
    obj <- readRDS(paste0(output_dir,'/',out,'.findcell_obj.rds'))
    obj <- isCell(obj,
                  num.test=num_test_val_final,
                  num.tn5=num_tn5_val_final,
                  num.ref=num_ref_val_final,
                  background.cutoff=min_tn5_val+1000,
                  min.pTSS=0.2,
                  min.FRiP=0.2,
                  min.pTSS.z= -2,
                  min.FRiP.z= -2,
                  verbose=F)
  }
}else{
  
  message(" - close the isCell")
  obj <- readRDS(paste0(output_dir,'/',out,'.findcell_obj.rds'))
  
  
}



# save QC data
saveRDS(obj, file=paste0(output_dir,'/',out, ".QC.rds"))
# convert to Socrates format for downstream analysis. 
soc.obj <- convertSparseData(obj, verbose=T)
# save QC object
saveRDS(soc.obj, file=paste0(output_dir,'/',out,".filtered.soc.rds"))




