##updating 050221
##change the provider version is TAIR10 and provider is NCBI 

##updation 010221 
##this script is trying to forge BSgenome

##soybean BSgenom
##modify the seed_file.txt first
##the current dir contains the seed_file.txt information
library(BSgenome)

args <- commandArgs(TRUE)
ipt_seed_file <- as.character(args[1])
output_dir <- as.character(args[2])
ipt_genome_name <- as.character(args[3])
#peaks.in <- as.character(args[3])

my_file <- read.dcf(ipt_seed_file, fields = NULL, all = FALSE, keep.white = NULL) 
my_file
write.dcf(my_file , file = paste0(output_dir,'/',"seed.dcf"), append = FALSE, useBytes = FALSE, indent = 0.1 * getOption("width"), width = 0.9 * getOption("width"), keep.white = NULL)
unlink(c(ipt_genome_name), recursive = TRUE, force = TRUE)
forgeBSgenomeDataPkg(paste0(output_dir,'/',"seed.dcf"),destdir=output_dir)

##in the terminal
##R CMD build BSgenome.Osa.323.v7
##R CMD check BSgenome.Osa.323.v7_1.0.0.tar.gz
##ignore checking the pdf was wrong
##R CMD INSTALL BSgenome.Osa.323.v7_1.0.0.tar.gz
##library(BSgenome.Osa.323.v7)
##Osa (type Gmax give the information)


##updation 010221
##how to prepare the seqnames
## paste("chr", c(1:20, "X", "M", "Un", paste(c(1:20, "X", "Un"), "_random", sep="")), sep="")
## seqnames: paste("chr", c(1:2), sep="")
##generate the real genome chromosome information for the soybean
#c('Gm01','Gm02','Gm03','Gm04','Gm05','Gm06','Gm07','Gm08','Gm09','Gm10','Gm11','Gm12','Gm13','Gm14','Gm15','Gm16','Gm17','Gm18','Gm19','Gm20','ChrCp','ChrMt')









##Test
#my_file <- read.dcf("seed_file.txt", fields = NULL, all = FALSE, keep.white = NULL) 
#my_file
#write.dcf(my_file , file = "seed.dcf", append = FALSE, useBytes = FALSE, indent = 0.1 * getOption("width"), width = 0.9 * getOption("width"), keep.white = NULL)

#library(BSgenome)
#unlink(c("OrfTxdb"), recursive = TRUE, force = TRUE)
#forgeBSgenomeDataPkg("seed.dcf")


