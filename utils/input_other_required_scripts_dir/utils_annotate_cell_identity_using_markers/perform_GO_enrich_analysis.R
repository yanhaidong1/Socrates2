###TopGO###
library("topGO")

args <- commandArgs(T)

# arguments
input_soc_obj_fl <- as.character(args[1])

input_annot_fl <- as.character(args[2])

input_log2fc_cutoff <- as.numeric(args[3])

prefix <- as.character(args[4])

output_dir <- as.character(args[5])



#########
##step 01 read the soc and build the list of target gene files
obj_dt <- readRDS(input_soc_obj_fl)

rank_diff_gene_dt <- obj_dt$rank_diff_gene

colnames(rank_diff_gene_dt) <- c('gene','cluster','log2fc')

rank_diff_gene_flt_dt <- rank_diff_gene_dt[rank_diff_gene_dt$log2fc > input_log2fc_cutoff,]

all_cluster_list <- unique(rank_diff_gene_flt_dt$cluster)

#########
##step 02 perform the Top GO

outs <- lapply(all_cluster_list, function(x){
  
  rank_diff_gene_flt_cluster_dt <- rank_diff_gene_flt_dt[rank_diff_gene_flt_dt$cluster == x,]
  
  interest_gene_list <- as.character(rank_diff_gene_flt_cluster_dt$gene)
  
  geneID2GO <- readMappings(file = input_annot_fl)
  
  ######################
  ##do the GO enrichment
  geneUniverse <- names(geneID2GO) 
  ##input the the genes that we will analyze, change if we want to change the list of interesting genes
  genesOfInterest <- as.character(interest_gene_list) 
  ##Then we need to tell TopGO where these interesting genes appear in the 'geneUniverse' vector:
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  ##Putting the data together into an R object
  ##We need to put our data in an object of type 'topGOdata'. 
  ##This will contain the list of genes of interest, the GO annotations, and the GO hierarchy itself.
  myGOdata <- new("topGOdata", description="My project", ontology="BP", 
                  allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
  #########Performing enrichment tests
  ##One type of test that topGO does is a Fisher's exact test based on gene counts:
  resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
  ##We can list the top ten significant results found:
  #allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", 
  #                   ranksOf = "classicFisher", topNodes = 50,rm.one=TRUE)
  allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", 
                     ranksOf = "classicFisher", topNodes = 50)
  ##rm.one if it exists, it will get wrong report
  
  ##add qval
  allRes$qval <- p.adjust(allRes$classicFisher, method="fdr")
  
  ##add the cluster
  allRes$cluster <- paste0('cluster_',x)
  
  return(allRes)
  
  #write.table(allRes,paste0(output_dir,'/',prefix,"_annot.txt"),sep = '\t',quote = F)
})

combine_dt <- do.call(rbind,outs)
combine_dt <- as.data.frame(combine_dt)

write.table(combine_dt,paste0(output_dir,'/',prefix,"_denovo_genes_GO_enrich.txt"),sep = '\t', quote = F)







