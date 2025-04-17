###################################################################################################
###                             chromVAR analysis of scATAC data                                ###
###################################################################################################
##updating 040925 we will load the species which has been built the BSgenome 
##updating 071323 add the ara BS
##updating 042923 we will debug to keep the peak order same as the matrix
##updating 011923 add a step to check whether the input peak file contains four columns
##updating 072422 debug make sure meta$depth <- meta$total depth is total or unique
##updating 062121 add the output_dir rm the wd
##updating 021721 generate a new function to getJasparMotifs

##add note information

#.libPaths(c('/home/hy17471/.conda/envs/r4.2lib_071922/lib/R/library',.libPaths()))


# arguments
args <- commandArgs(TRUE)
threads <- as.numeric(args[1])
input.sp <- as.character(args[2])
peaks.in <- as.character(args[3])


#species <- as.character(args[4]) ##eg. soybean;rice_MSU;rice_323

meta <- read.delim(as.character(args[4]),row.names = 1)
#rownames(meta) <- meta$cellID
output_dir <- as.character(args[5])

##parameter
#min_fragments_per_peak_val <- as.numeric(args[6]) ##for all organ version it is 2500
min_fragments_per_peak_val <- 2000

input_BSgenome_config_fl <- as.character(args[6])
input_GCbias_config_fl <- as.character(args[7])
input_motifmatch_config_fl <- as.character(args[8])



# setwd
#setwd(wd)

# set variables


# load libraries
library(chromVAR)
library(motifmatchr)
library(BiocParallel)

##updating 040925 we will source the configure file
##updating 062121 decide the library we will load for the species
source(input_BSgenome_config_fl)
#if (species == 'soybean') {
# library(BSgenome.Gmax.a4.v1)
#}
#if (species == 'rice_MSU') {
#  library(BSgenome.Osa.MSU.v7)
#}
#if (species == 'rice_323') {
#  library(BSgenome.Osa.323.v7)
#}
#if (species == 'maize_v4') {
#  library(BSgenome.maize.v4)
#}
#if (species == 'Ath.TAIR10'){
#  library(BSgenome.Ara.TAIR10)
#}



library(Matrix)
library(SummarizedExperiment)
library(GenomicAlignments)
library(dplyr)
library(TFBSTools)
library(JASPAR2024)
library(pheatmap)
library(circlize)

# functions
loadPeaks <- function(filename, extra_cols=4){

	# load bed
	bed <- read.delim(file = filename, header = FALSE, sep = "\t",
            stringsAsFactors = FALSE)[, c(1:3, extra_cols)]

	# convert 2 GR
	colnames(bed) <- c("chr", "start", "end", names(extra_cols))
	bed[, "start"] <- bed[, "start"] + 1
        bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

	# sort
  sorted_bed <- sortSeqlevels(bed)
  sorted_bed <- sort(sorted_bed, ignore.strand = TRUE)
  return(sorted_bed)
}

# set number of cores
register(MulticoreParam(threads))

# verbose
message("########################################")
message("########################################")
message("")
message("============================")
message("     running chromVAR       ")
message("============================")
message("")


###################################################################################################
### load and process data									   
###################################################################################################

# input files


if(file.exists(paste0(output_dir,'/temp_peak_sparse_mtx.rds'))){
  
  message("peak sparse exists, directly load")
  
  a <- readRDS(paste0(output_dir,'/temp_peak_sparse_mtx.rds'))
  
}else{
  # build counts matrix
  message("Loading count matrix ...")
  a <- read.table(input.sp,stringsAsFactors=T)
  a <- sparseMatrix(i=as.numeric(a$V1),
                    j=as.numeric(a$V2),
                    x=as.numeric(a$V3),
                    dimnames=list(levels(a$V1), levels(a$V2)))
  #a <- a[as.character(peaks2$V4),] 
  
  
  saveRDS(a,paste0(output_dir,'/temp_peak_sparse_mtx.rds'))
  
}



message("Loading peak information ...")

##updating 011923
message('check whether peaks have four columns or not')
temp_peak_dt <- read.delim(peaks.in,header = F)
temp_peak_dt$peaknm <- paste0(temp_peak_dt$V1,'_',temp_peak_dt$V2,'_',temp_peak_dt$V3)
temp_peak_dt <- temp_peak_dt[,c('V1','V2','V3','peaknm')]
colnames(temp_peak_dt) <- c('V1','V2','V3','V4')

if (ncol(temp_peak_dt) == 3){
  
  temp_peak_dt$V4 <- paste0(temp_peak_dt$V1,'_',temp_peak_dt$V2,'_',temp_peak_dt$V3)
  
  ##updating 031123 we will make sure there is no chr
  rownames(a) <- gsub('chr','',rownames(a))

  ##updating 030823
  ##we need to make sure the peaks match between two sets
  ##there are peaks missing as we downsample cells so in the peak sparse files, there are some peaks are lost just like 100 peaks
  intersect_ids <- intersect(rownames(a),temp_peak_dt$V4)
  temp_peak_dt <- temp_peak_dt[temp_peak_dt$V4 %in% intersect_ids,]
  
  write.table(temp_peak_dt,paste0(output_dir,'/temp_modi_peak_fl.txt'),row.names=F,col.names = F,quote = F,sep = '\t')
  peaks.in = paste0(output_dir,'/temp_modi_peak_fl.txt')
}else{
  
  ##updating 030823
  ##if the temp_peak_dt has 4 col
  if (ncol(temp_peak_dt) == 4){
    
    rownames(a) <- gsub('chr','',rownames(a))
    
    intersect_ids <- intersect(rownames(a),temp_peak_dt$V4)
    temp_peak_dt <- temp_peak_dt[temp_peak_dt$V4 %in% intersect_ids,]
    write.table(temp_peak_dt,paste0(output_dir,'/temp_modi_peak_fl.txt'),row.names=F,col.names = F,quote = F,sep = '\t')
    peaks.in = paste0(output_dir,'/temp_modi_peak_fl.txt')
    
  
  }else{
    print('Wrong: the temp peak dt does not have four column')
  }
}



##filter a
a <- a[intersect_ids,]

peaks <- loadPeaks(peaks.in, extra_cols=c(4))

##order the matrix to be same as the peaks
peaks_bed <- as.data.frame(peaks)
s.ids <- paste(peaks_bed$seqnames,peaks_bed$start-1,peaks_bed$end,sep="_")
a <- a[s.ids,]

peaks2 <- read.table(peaks.in)

saveRDS(peaks,paste0(output_dir,'/peaks.rds'))
saveRDS(a, paste0(output_dir,'/a.rds'))

# load meta.data
message("Loading meta data ...")
#meta <- read.table(as.character(args[5]))
#meta$depth <- meta$unique
a <- a[,colnames(a) %in% rownames(meta)]
meta <- meta[colnames(a),]
print(head(a[,1:5]))
print(head(meta))
#meta$depth <- meta$total  ##here is total not the unique
meta$depth <- Matrix::colSums(a)
message("cells = ",ncol(a), " | peaks = ", nrow(a))
message('peaks number = ',length(peaks))

# create frag counts object
message("Creating experiment object ...")
##first enerate fragment information from the function SummarizedExperiment
fragment_counts <- SummarizedExperiment(assays = list(counts = a),
                                        rowRanges = peaks,
                                        colData = meta)

saveRDS(fragment_counts,paste0(output_dir,'/fragment_counts.rds'))

# clean-up memory
rm(a)
rm(peaks2)

# add GC data
message("Estimating GC bias ...")

source(input_GCbias_config_fl)

#if (species == 'soybean') {
#  fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.Gmax.a4.v1)
#}
#if (species == 'rice_MSU') {
#  fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.Osa.MSU.v7)
#}
#if (species == 'rice_323') {
#  fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.Osa.323.v7)
#}
#if (species == 'maize_v4') {
#  fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.maize.v4)
#}
#if (species == 'Ath.TAIR10') {
#  fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.Ara.TAIR10)
#}


#fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.Gmax.a4.v1)

# filter cells
message("Filtering samples ...")
#this is from the ChromVAR
##filter peaks
filtered_counts <- filterSamples(fragment_counts, min_depth=100, min_in_peaks=0.1, shiny=F)

saveRDS(filtered_counts,paste0(output_dir,'/filtered_counts.rds'))


# filter peaks
message("Filtering peaks ...")
#filtered_counts <- filterPeaks(filtered_counts, non_overlapping=T, min_fragments_per_peak=100)
##100 does not work for the filtration
filtered_counts <- filterPeaks(filtered_counts, non_overlapping=T, min_fragments_per_peak= min_fragments_per_peak_val)
message('dimension of filtered_counts is')
print(dim(filtered_counts))
##> dim(filtered_counts_fltpeaks)
#[1] 97563 69485

filtered_peaks <- rownames(filtered_counts)
write.csv(as.data.frame(filtered_peaks),paste0(output_dir,'/temp_filtered_peaks.csv'),quote = F)
##save to another form
#filtered_peaks_dt <- as.data.frame(filtered_peaks)



###############################################################################
## kmer deviation
###############################################################################
#message("Running Kmer analysis...")
#kmers            <- matchKmers(7, filtered_counts, genome = BSgenome.Gmax.a4.v1)
#dev.kmers        <- computeDeviations(object = filtered_counts, annotations = kmers)
#dev.kmers.scores <- deviationScores(dev.kmers)
#kmer.devs        <- deviations(dev.kmers)
#write.table(dev.kmers.scores, file="kmer.scores.JBv3.txt", quote=F, row.names=T, col.names=T, sep="\t")
#write.table(kmer.devs, file="kmer.deviations.JBv3.txt", quote=F, row.names=T, col.names=T, sep="\t")

# kmer plots
#variability <- computeVariability(dev.kmers)
#pdf("kmer.variability.JBv3.pdf", width=6, height=4)
#plotVariability(variability, use_plotly = FALSE, n=10)
#dev.off()

# kmer cor
#sample_cor <- getSampleCorrelation(dev.kmers)
#pdf(file="kmer.correlation.pdf", width=6, height=6)
#pheatmap(as.dist(sample_cor), 
#         clustering_distance_rows = as.dist(1-sample_cor), 
#         clustering_distance_cols = as.dist(1-sample_cor))
#dev.off()

# kmer tSNE
#tsne_results <- deviationsTsne(dev.kmers, threshold = 0.5, perplexity = 50, 
#                               shiny = FALSE)
#tsne_plots   <- plotDeviationsTsne(dev.kmers, tsne_results, 
#                                   sample_column = "library", shiny = FALSE)
#pdf(file="kmer.tSNE.JBv3.pdf", width=6, height=6)
#tsne_plots[[1]]
#dev.off()

# Z-score heatmap
#mat <- dev.kmers.scores
#mat[is.na(mat)] <- 0
#ha=HeatmapAnnotation(bar=data.frame(do.call(rbind, strsplit(colnames(dev.kmers.scores), "\\.")))[,4])
#cols <- colorRamp2(c(-2,0,2),c("darkorchid4", "grey85", "forestgreen"))
#pdf(file="kmer.heatmap.Zscores.pdf", width=6, height=6)
#Heatmap(mat, col=cols, use_raster=T, show_row_names=F, show_column_names=F, top_annotation = ha)
#dev.off()


###############################################################################
## motif deviation
###############################################################################

# estimate deviations
message("Running motif analysis ...")
##the chromVAR has function to get the JasparMotifs

##updating 021721
getJasparMotifs_new <- function (species = "Homo sapiens", collection = "CORE", ...)
{
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
    names(out) <- paste(names(out), TFBSTools::name(out),
                        sep = "_")
  return(out)
}

message("get the motif devs score")
jaspmotifs       <- getJasparMotifs_new(species = "Arabidopsis thaliana")

source(input_motifmatch_config_fl)

#if (species == 'soybean') {
#  motif            <- matchMotifs(jaspmotifs, filtered_counts, genome = BSgenome.Gmax.a4.v1)
  
#  saveRDS(motif,paste0(output_dir,'/motif.rds'))
#}
#if (species == 'rice_MSU') {
#  motif            <- matchMotifs(jaspmotifs, filtered_counts, genome = BSgenome.Osa.MSU.v7)
#}
#if (species == 'rice_323') {
#  motif            <- matchMotifs(jaspmotifs, filtered_counts, genome = BSgenome.Osa.323.v7)
#}
#if (species == 'maize_v4') {
#  motif            <- matchMotifs(jaspmotifs, filtered_counts, genome = BSgenome.maize.v4)
#}
#if (species == 'Ath.TAIR10') {
#  motif            <- matchMotifs(jaspmotifs, filtered_counts, genome = BSgenome.Ara.TAIR10)
  
#}


#motif            <- matchMotifs(jaspmotifs, filtered_counts, genome = BSgenome.Gmax.a4.v1)
saveRDS(motif, file=paste0(output_dir,'/opt_motif_matches.rds'))
##save the motif to the matrix
motif_mtx <- as.matrix(assays(motif)$motifMatches)
motif_mtx[motif_mtx = FALSE] <- 0
saveRDS(motif_mtx,file=paste0(output_dir,'/opt_motif_matches_mtx.rds'))

dev.motif        <- computeDeviations(object = filtered_counts, annotations = motif)
dev.motif.scores <- deviationScores(dev.motif)
motif.devs       <- deviations(dev.motif)
write.table(t(dev.motif.scores), file=paste0(output_dir,'/',"motif.scores.txt"), quote=F, row.names=T, col.names=T, sep="\t")
write.table(t(motif.devs), file=paste0(output_dir,'/',"motif.deviations.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# plot motif
#variability <- computeVariability(dev.motif)
#pdf(paste0(output_dir,"/motif.variability.pdf"), width=6, height=4)
#plotVariability(variability, use_plotly = FALSE, n=10)
#dev.off()

# motif tSNE
#tsne_results <- deviationsTsne(dev.motif, threshold = 0.5, perplexity = 50, 
#                               shiny = FALSE)
#tsne_plots   <- plotDeviationsTsne(dev.motif, tsne_results, 
#                                   sample_column = "library", shiny = FALSE)
#pdf(file=paste0(output_dir,"/motif.tSNE.pdf"), width=6, height=6)
#tsne_plots[[1]]
#dev.off()

## background peaks
#bbpeaks <- getBackgroundPeaks(filtered_counts)
#write.table(bbpeaks, file=paste0(output_dir,'/',"backgroundPeaks.mat.txt"), quote=F, row.names=T, col.names=T, sep="\t")
message("--Finished--")
