## sample ##

##updating 032622 use the NMF way to do the classification
##

##this script is to do the embeddings of cells and features between training and independent testing dataset 

##updating 032522 this script is to check the scPred

##1. we need to obtain the reference embedding first
#1) we first use all the train set (90%) the other 10% is the independent testing dataset to do the training
#2) now we already have train set (PCs by cells) and train cell meta

##redo the processing
##1) all data (gene by cells) divided into independent and training
##2) transfer the train set to be Seurat object without any normalization
##3) use the Seurat to do the reducation to obtain the gene and cell embeddings respectively
##4) use the cell embeddings to have a training
##5) use the gene embedings to make a new embedings for the independent testing datasets.
##6) do the prediction

##load library
library("scPred")
library("Seurat")
library("magrittr")
library('matrixStats')
library('Matrix')
library('harmony')
library(RcppML)

args = commandArgs(T)

input_training_exp_fl <- as.character(args[1])
##'/scratch/hy17471/rice_altlas_scATAC_seq_042021/08_training_classifiers_031822/02_SPmarker_032322/working_dir_topCTver/Step1_prepare_data_dir/Step1_3_split_data_dir/output_dir/step1_split_train_indetest_opt_dir/opt_exp_train.csv'

input_testing_exp_fl <- as.character(args[2])
##'/scratch/hy17471/rice_altlas_scATAC_seq_042021/08_training_classifiers_031822/02_SPmarker_032322/working_dir_topCTver/Step1_prepare_data_dir/Step1_3_split_data_dir/output_dir/step1_split_train_indetest_opt_dir/opt_exp_test.csv'

input_training_meta_fl <- as.character(args[3])

input_output_dir <- as.character(args[4])

npcs <- as.numeric(args[5])

reduction_major_type <- as.character(args[6]) ##SVD or NMF

top_variant_ft_num <- as.numeric(args[7]) ##default 5000


##set some parameters
reduction_type <- 'pca'
#npcs = 100
max.iter.harmony = 20
seed=123

##read data
training_dt <- readRDS(input_training_exp_fl)
#rownames(training_dt) <- training_dt[,1]  ## set rownames
#training_dt <- training_dt[, -1]           ## remove the first variable

testing_obj <- readRDS(input_testing_exp_fl)
##updating 050725 we will use the soc to be the test
testing_dt <- testing_obj$gene_acc

#rownames(testing_dt) <- testing_dt[,1]  ## set rownames
#testing_dt <- testing_dt[, -1]           ## remove the first variable



##########################
##process training dataset
##create seurat project
reference_obj <- CreateSeuratObject(counts = training_dt, project = "training", min.cells = 0, min.features = 0)
#> reference_obj
#An object of class Seurat
#41935 features across 10201 samples within 1 assay
#Active assay: RNA (41935 features, 0 variable features)

##scale matrix

reference_obj <- NormalizeData(reference_obj)
reference_obj <- ScaleData(reference_obj, features = rownames(reference_obj))


##############
##do reduction 
##SVD ways
if (reduction_major_type == 'SVD'){
  reference_obj <- RunPCA(reference_obj, features = rownames(reference_obj),npcs=npcs)
  reduction_ref_data <- Reductions(reference_obj, slot = reduction_type)
  loadings_ref <- Loadings(reduction_ref_data) ##loadings_ref rowname is gene feature name and PCs is the column
  cellEmbeddings_ref <- Embeddings(reduction_ref_data) ##cellembeddings rowname is cell name and PCs is the column
}

if (reduction_major_type == 'NMF'){
  ##############################################
  ##use the nmf to calculate the loadings of ref
  ##obtain the matrix
  assay_mtx_ref <- as.matrix(GetAssay(reference_obj[,])[,])
  ##do the nmf reduction
  ft_pcs <- RcppML::nmf(t(assay_mtx_ref), npcs, verbose=T)
  dim(ft_pcs$h) ##20 32838
  ##scale
  ft_pcs$u <- t(ft_pcs$h)
  ft_pcs$v <- ft_pcs$w
  ft_pc <- t(ft_pcs$h) %*% Diagonal(x=1/ft_pcs$d)
  dim(ft_pc) ##32838genes 20PCs
  ft_pc[is.na(ft_pc)] <- 0
  ##rename
  ##get the reverse of ft_pc since the pcs output order of column is from small to large
  ft_pc <- as.matrix(rev(as.data.frame(ft_pc)))
  rownames(ft_pc) <- rownames(assay_mtx_ref)
  colnames(ft_pc) <- paste0("PC_", seq(1:ncol(ft_pc)))
  loadings_ref <- as.matrix(ft_pc)
  
  
  ##do the nmf reduction cell embeddings
  cells_pcs <- RcppML::nmf(assay_mtx_ref, npcs, verbose=T)
  dim(cells_pcs$h) ##20 3500
  ##scale
  cells_pcs$u <- t(cells_pcs$h)
  cells_pcs$v <- cells_pcs$w
  cells_pc <- t(cells_pcs$h) %*% Diagonal(x=1/cells_pcs$d)
  dim(cells_pc) ##3500cells 20PCs
  cells_pc[is.na(cells_pc)] <- 0
  ##rename
  ##get the reverse
  cells_pc <- as.matrix(rev(as.data.frame(cells_pc)))
  rownames(cells_pc) <- colnames(assay_mtx_ref)
  colnames(cells_pc) <- paste0("PC_", seq(1:ncol(cells_pc)))
  cellEmbeddings_ref <- as.matrix(cells_pc)
}


#########################
##process testing dataset
##obtain the shared ID feature

##note: some features have zero variance when we calculate the pca so we will filter it
reference_PCA_obj <- RunPCA(reference_obj, features = rownames(reference_obj),npcs=npcs)
reduction_ref_data <- Reductions(reference_PCA_obj, slot = reduction_type)
loadings_ref_PCA <- Loadings(reduction_ref_data) ##loadings_ref rowname is gene feature name and PCs is the column
shared_ref_features_PCA <- intersect(rownames(loadings_ref_PCA),rownames(loadings_ref))
loadings_ref <- loadings_ref[shared_ref_features_PCA,]


rownames(testing_dt) <- gsub('_','-',rownames(testing_dt))
shared_ID_test_ref <- intersect(rownames(testing_dt),rownames(loadings_ref))

assay_ref <- GetAssay(reference_obj)
#means_ref <- rowMeans(assay_ref)

assay_mtx_ref <- as.matrix(LayerData(object = reference_obj, layer = "counts"))
#assay_mtx_ref <- as.matrix(GetAssay(reference_obj[,])[,]) 

means_ref <- rowMeans(assay_mtx_ref)

std_ref <- rowSds(assay_mtx_ref)

names(std_ref) <- rownames(assay_ref) -> names(std_ref)

##obtain the shared feature ID information
means_ref <- means_ref[shared_ID_test_ref]
std_ref <- std_ref[shared_ID_test_ref]

loadings_ref <- loadings_ref[shared_ID_test_ref, ]

testing_data <- as.matrix(testing_dt[shared_ID_test_ref,])

##scale testing data
testing_data <- Matrix::t(testing_data)
scaled_testing_data <- scale(testing_data, means_ref, std_ref) ##row is cell and col is gene

##no need for this command
#scaled_testing_data[is.na(scaled_testing_data)] <- 0

##obtain the test embeddings
test_embeddings <- scaled_testing_data %*% loadings_ref

##create a dataset of reference and testing
dataset <- factor(c(rep("reference", nrow(cellEmbeddings_ref)), rep("test", nrow(test_embeddings))), 
                  levels = c("reference", "test"))

rownames(cellEmbeddings_ref) <- paste0("ref_", rownames(cellEmbeddings_ref))
rownames(test_embeddings) <- paste0("test_", rownames(test_embeddings))

eigenspace <- as.data.frame(rbind(cellEmbeddings_ref, test_embeddings))
#> head(eigenspace)
#PC_1        PC_2       PC_3       PC_4
#ref_CCACGGACATGCTGGC-1   3.5462928  4.24607540  0.7342539  0.4024824
#ref_GATTCAGGTCTCCACT-1   5.1154612 -1.93469009 -2.5112122 -0.1026473
#ref_TGTATTCTCTATGTGG-1 -32.3288611 -0.03741952 -4.1121667  7.3990149
#ref_GTCTCGTGTTCACGGC-1  -0.9768477 -6.43608814 14.1769009  0.2691508
#ref_TATCAGGCACCAACCG-1   0.7108476  8.61825733  0.8995970 -0.1301161
#ref_TTTGCGCAGTAAGTAC-1   5.6391219 -2.81575179 -2.1088155 -0.2297144

meta_data <- data.frame(rownames(eigenspace), dataset = dataset)
#> head(meta_data)
#rownames.eigenspace.   dataset
#1 ref_CCACGGACATGCTGGC-1 reference
#2 ref_GATTCAGGTCTCCACT-1 reference
#3 ref_TGTATTCTCTATGTGG-1 reference
#4 ref_GTCTCGTGTTCACGGC-1 reference
#5 ref_TATCAGGCACCAACCG-1 reference
#6 ref_TTTGCGCAGTAAGTAC-1 reference

##do the harmony embeddings for the testing dataset
set.seed(seed)
harmony_embeddings <- HarmonyMatrix(eigenspace, 
                                    meta_data, 
                                    'dataset', 
                                    do_pca = FALSE, 
                                    reference_values = "reference",
                                    max.iter.harmony = max.iter.harmony)

test_embeddings_aligned <- harmony_embeddings[dataset == "test", , drop = FALSE]

##save the training and testing embeddings
##remove the 'ref_' and 'test_'
test_embeddings_aligned_t <- t(test_embeddings_aligned)
colnames(test_embeddings_aligned_t) <- gsub('test_','',colnames(test_embeddings_aligned_t))
write.csv(test_embeddings_aligned_t,paste0(input_output_dir,'/opt_indep_testing_embeddings.csv'),quote = F)

cellEmbeddings_ref_t <- t(cellEmbeddings_ref)
colnames(cellEmbeddings_ref_t) <- gsub('ref_','',colnames(cellEmbeddings_ref_t))


##udpating 063024 
##we will intersect with the meta file
##updating 100425 we will update shared IDs
dt <- read.delim(input_training_meta_fl,header=F)
dt$V1 <- gsub(':','.',dt$V1)
dt$V1 <- gsub('-','.',dt$V1)
shared_IDs <- intersect(dt$V1,colnames(cellEmbeddings_ref_t))

if ((length(shared_IDs)) != 0){

  cellEmbeddings_ref_t <- cellEmbeddings_ref_t[,shared_IDs]
  
  write.csv(cellEmbeddings_ref_t,paste0(input_output_dir,'/opt_training_embeddings.csv'),quote = F)
}else{
  dt <- read.delim(input_training_meta_fl,header=F)
  shared_IDs <- intersect(dt$V1,colnames(cellEmbeddings_ref_t))
  cellEmbeddings_ref_t <- cellEmbeddings_ref_t[,shared_IDs]
  
  write.csv(cellEmbeddings_ref_t,paste0(input_output_dir,'/opt_training_embeddings.csv'),quote = F)
  
}






