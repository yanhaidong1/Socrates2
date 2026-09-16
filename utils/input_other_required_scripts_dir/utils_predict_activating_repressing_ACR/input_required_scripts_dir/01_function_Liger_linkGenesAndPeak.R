getFDR <- function(observed, expected, fdr=0.01, grid=1000, verbose=F){
    
    # convert to data.frames
    obs <- as.data.frame(summary(observed))
    obs$i <- rownames(observed)[obs$i]
    obs$j <- colnames(observed)[obs$j]
    exp <- as.data.frame(summary(expected))
    exp$i <- rownames(expected)[exp$i]
    exp$j <- colnames(expected)[exp$j]
    
    # split into +/-
    n.exp <- subset(exp, exp$x < 0)
    p.exp <- subset(exp, exp$x > 0)
    
    # get counts
    pos.nexp <- nrow(p.exp)
    neg.nexp <- nrow(n.exp)
    pos.nobs <- nrow(subset(obs, obs$x > 0))
    neg.nobs <- nrow(subset(obs, obs$x < 0))
    if(verbose){message(" - number expected links = (+) ",pos.nexp, " | (-) ",neg.nexp)}
    if(verbose){message(" - number observed links = (+) ",pos.nobs, " | (-) ",neg.nobs)}
    
    # generate range of thresholds
    p.vals <- seq(from=0, to=1, length.out=grid)
    n.vals <- seq(from=0, to= -1, length.out=grid)
    
    # iterate over grid
    if(verbose){message(" - scanning positive thresholds ...")}
    p.thresh <- c()
    for(i in p.vals){
        num.exp <- sum(p.exp$x > as.numeric(i))
        c.fdr <- num.exp/pos.nexp
        if(is.na(c.fdr)){
            c.fdr <- 0
        }
        p.thresh <- c(p.thresh, c.fdr)
        message(" - (+) correlation threshold = ", i, " | FDR = ", c.fdr)
    }
    if(verbose){message(" - scanning negative thresholds ...")}
    n.thresh <- c()
    for(i in n.vals){
        num.exp <- sum(n.exp$x < as.numeric(i))
        c.fdr <- num.exp/(neg.nexp)
        if(is.na(c.fdr)){
            c.fdr <- 0
        }
        n.thresh <- c(n.thresh, c.fdr)
        message(" - (-) correlation threshold = ", i, " | FDR = ", c.fdr)
    }
    
    # select cut-offs
    p.threshold <- min(p.vals[which(p.thresh <= fdr)])
    n.threshold <- max(n.vals[which(n.thresh <= fdr)])
    
    # filter
    obs <- subset(obs, obs$x > p.threshold | obs$x < n.threshold)
    
    # verbose number of +/- linkages
    pos.links <- nrow(subset(obs, obs$x > 0))
    neg.links <- nrow(subset(obs, obs$x < 0))
    message(" - found ",pos.links, " + and ", neg.links," - gene-peak links ...")
    
    # return
    return(obs)
}
linkGenesAndPeaks2 <- function (gene_counts, peak_counts, genes.list = NULL, dist = "spearman",
                                alpha = 0.05, max_dist=200000, path_to_coords){
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Package \"GenomicRanges\" needed for this function to work. Please install it by command:\n",
             "BiocManager::install('GenomicRanges')", call. = FALSE)
    }
    if (!requireNamespace("IRanges", quietly = TRUE)) {
        stop("Package \"IRanges\" needed for this function to work. Please install it by command:\n",
             "BiocManager::install('IRanges')", call. = FALSE)
    }
    #peak.names <- strsplit(rownames(peak_counts), "[:-_]")
	peak.names <- strsplit(rownames(peak_counts), "_")
    chrs <- Reduce(append, lapply(peak.names, function(peak) {
        peak[1]
    }))
    chrs.start <- Reduce(append, lapply(peak.names, function(peak) {
        peak[2]
    }))
    chrs.end <- Reduce(append, lapply(peak.names, function(peak) {
        peak[3]
    }))
    peaks.pos <- GenomicRanges::GRanges(seqnames = chrs, ranges = IRanges::IRanges(as.numeric(chrs.start),
                                                                                   end = as.numeric(chrs.end)))
    gene.names <- read.csv2(path_to_coords, sep = "\t", header = FALSE,
                            stringsAsFactors = F)
    gene.names <- gene.names[complete.cases(gene.names), ]
    genes.coords <- GenomicRanges::GRanges(seqnames = gene.names$V1,
                                           ranges = IRanges::IRanges(as.numeric(gene.names$V2),
                                                                     end = as.numeric(gene.names$V3)))
    names(genes.coords) <- gene.names$V4
    gene_counts <- t(gene_counts)
    peak_counts <- t(peak_counts)
    if (missing(genes.list)) {
        genes.list <- colnames(gene_counts)
    }
    missing_genes <- !genes.list %in% names(genes.coords)
    if (sum(missing_genes) != 0) {
        print(paste0("Removing ", sum(missing_genes), " genes not found in given gene coordinates..."))
    }
    genes.list <- genes.list[!missing_genes]
    genes.coords <- genes.coords[genes.list]
    print("Calculating correlation for gene-peak pairs...")
    each.len <- 0
    elements <- lapply(seq(length(genes.list)), function(pos){
        gene.use <- genes.list[pos]
        gene.loci <- GenomicRanges::trim(suppressWarnings(GenomicRanges::promoters(GenomicRanges::resize(genes.coords[gene.use],
                                                                                                         width = 1, fix = "start"), upstream = max_dist, downstream = max_dist)))
        peaks.use <- S4Vectors::queryHits(GenomicRanges::findOverlaps(peaks.pos,gene.loci))
        if ((x <- length(peaks.use)) == 0L) {
            return(list(NULL, as.numeric(each.len), NULL))
        }
        res <- suppressWarnings(psych::corr.test(x = gene_counts[,gene.use], 
                                                 y = as.matrix(peak_counts[, peaks.use]),
                                                 method = dist, 
                                                 adjust = "holm", 
                                                 ci = FALSE, 
                                                 use = "complete"))
        pick <- res[["p"]] < alpha
        pick[is.na(pick)] <- FALSE
        if(sum(pick) == 0){
            return(list(NULL, as.numeric(each.len), NULL))
        }
        else{
            res.corr <- as.numeric(res[["r"]][pick])
            peaks.use <- peaks.use[pick]
        }
        assign("each.len", each.len + length(peaks.use), envir = parent.frame(2))
        return(list(as.numeric(peaks.use), as.numeric(each.len),res.corr))
    })
    i_index <- Reduce(append, lapply(elements, function(ele) {
        ele[[1]]
    }))
    p_index <- c(0, Reduce(append, lapply(elements, function(ele) {
        ele[[2]]
    })))
    value_list <- Reduce(append, lapply(elements, function(ele) {
        ele[[3]]
    }))
    regnet <- sparseMatrix(i = i_index, p = p_index, x = value_list,
                           dims = c(ncol(peak_counts), length(genes.list)), dimnames = list(colnames(peak_counts),
                                                                                            genes.list))
    return(regnet)
}

#----------linkGenesAndPeaks3-----------
#out put correlation score and p value.
# using pvalue from permutation test.
linkGenesAndPeaks3 <- function (gene_counts, peak_counts, genes.list = NULL, dist = "spearman",
                                alpha = 0.05, max_dist=200000, path_to_coords, num.sim=10000){
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Package \"GenomicRanges\" needed for this function to work. Please install it by command:\n",
             "BiocManager::install('GenomicRanges')", call. = FALSE)
    }
    if (!requireNamespace("IRanges", quietly = TRUE)) {
        stop("Package \"IRanges\" needed for this function to work. Please install it by command:\n",
             "BiocManager::install('IRanges')", call. = FALSE)
    }
    #peak.names <- strsplit(rownames(peak_counts), "[:-_]")
	peak.names <- strsplit(rownames(peak_counts), "_")
    chrs <- Reduce(append, lapply(peak.names, function(peak) {
        peak[1]
    }))
    chrs.start <- Reduce(append, lapply(peak.names, function(peak) {
        peak[2]
    }))
    chrs.end <- Reduce(append, lapply(peak.names, function(peak) {
        peak[3]
    }))
    peaks.pos <- GenomicRanges::GRanges(seqnames = chrs, ranges = IRanges::IRanges(as.numeric(chrs.start),
                                                                                   end = as.numeric(chrs.end)))
    gene.names <- read.csv2(path_to_coords, sep = "\t", header = FALSE,
                            stringsAsFactors = F)
    gene.names <- gene.names[complete.cases(gene.names), ]
    genes.coords <- GenomicRanges::GRanges(seqnames = gene.names$V1,
                                           ranges = IRanges::IRanges(as.numeric(gene.names$V2),
                                                                     end = as.numeric(gene.names$V3)))
    names(genes.coords) <- gene.names$V4
    gene_counts <- t(gene_counts)
    peak_counts <- t(peak_counts)
    if (missing(genes.list)) {
        genes.list <- colnames(gene_counts)
    }
    missing_genes <- !genes.list %in% names(genes.coords)
    if (sum(missing_genes) != 0) {
        print(paste0("Removing ", sum(missing_genes), " genes not found in given gene coordinates..."))
    }
    genes.list <- genes.list[!missing_genes]
    genes.coords <- genes.coords[genes.list]
    print("Calculating correlation for gene-peak pairs...")
    elements <- lapply(seq(length(genes.list)), function(pos){
        gene.use <- genes.list[pos]
        gene.loci <- GenomicRanges::trim(suppressWarnings(GenomicRanges::promoters(GenomicRanges::resize(genes.coords[gene.use],
                                                                                                         width = 1, fix = "start"), upstream = max_dist, downstream = max_dist)))
        peaks.use <- S4Vectors::queryHits(GenomicRanges::findOverlaps(peaks.pos,gene.loci))
        if ((y <- length(peaks.use)) == 0L) {
		res.list <- c("NA", gene.use, "NA", "NA")
		res.list <- as.data.frame(res.list)
		res.list <- t(res.list)
        }else{
        res <- lapply(peaks.use, function(x){
		a = gene_counts[,gene.use]
		b = peak_counts[, x]
		cor.value <- cor(a,b, method=dist)
		p.value <- NULL
		res.p <- NULL
		if (cor.value >0)
		{
		 res.p <- perm.cor.test(a, b, method=dist, alternative="greater", num.sim=num.sim)
		 p.value <- res.p$p.value
		}else{
		 res.p <- perm.cor.test(a, b, method=dist, alternative="less", num.sim=num.sim)
		 p.value <- res.p$p.value
		}
		res.list <- c(colnames(peak_counts)[x], gene.use, cor.value, p.value)
		res.list <- as.data.frame(res.list)
		res.list <- t(res.list)
		})
		res <- do.call(rbind, res)
	}
    })
     elements <- do.call(rbind, elements)
     colnames(elements) <- c("peak_id","gene_id","cor.value","p.value")
	 elements <- as.data.frame(elements)
	 elements <- elements[elements$p.value != "NA",]
	 #elements <- elements[elements$p.value < alpha,]
	 elements$FDR <- p.adjust(elements$p.value, method = "BH")
	return(elements)
}

#----------linkGenesAndPeaks3.2-----------
#out put correlation score and p value.
# using pvalue from permutation test.
# make a parallel version.
linkGenesAndPeaks3.2 <- function (gene_counts, peak_counts, genes.list = NULL, dist = "spearman",
                                alpha = 0.05, max_dist=200000, path_to_coords, num.sim=10000, cores=24){
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Package \"GenomicRanges\" needed for this function to work. Please install it by command:\n",
             "BiocManager::install('GenomicRanges')", call. = FALSE)
    }
    if (!requireNamespace("IRanges", quietly = TRUE)) {
        stop("Package \"IRanges\" needed for this function to work. Please install it by command:\n",
             "BiocManager::install('IRanges')", call. = FALSE)
    }
    #peak.names <- strsplit(rownames(peak_counts), "[:-_]")
	peak.names <- strsplit(rownames(peak_counts), "_")
    chrs <- Reduce(append, lapply(peak.names, function(peak) {
        peak[1]
    }))
    chrs.start <- Reduce(append, lapply(peak.names, function(peak) {
        peak[2]
    }))
    chrs.end <- Reduce(append, lapply(peak.names, function(peak) {
        peak[3]
    }))
    peaks.pos <- GenomicRanges::GRanges(seqnames = chrs, ranges = IRanges::IRanges(as.numeric(chrs.start),
                                                                                   end = as.numeric(chrs.end)))
    gene.names <- read.csv2(path_to_coords, sep = "\t", header = FALSE,
                            stringsAsFactors = F)
    gene.names <- gene.names[complete.cases(gene.names), ]
    genes.coords <- GenomicRanges::GRanges(seqnames = gene.names$V1,
                                           ranges = IRanges::IRanges(as.numeric(gene.names$V2),
                                                                     end = as.numeric(gene.names$V3)))
    names(genes.coords) <- gene.names$V4
    gene_counts <- t(gene_counts)
    peak_counts <- t(peak_counts)
    if (missing(genes.list)) {
        genes.list <- colnames(gene_counts)
    }
    missing_genes <- !genes.list %in% names(genes.coords)
    if (sum(missing_genes) != 0) {
        print(paste0("Removing ", sum(missing_genes), " genes not found in given gene coordinates..."))
    }
    genes.list <- genes.list[!missing_genes]
    genes.coords <- genes.coords[genes.list]
    print("Calculating correlation for gene-peak pairs...")
	cl <- makeCluster(cores)
	registerDoParallel(cl)
#    elements <- lapply(seq(length(genes.list)), function(pos){
	elements <- foreach( pos = seq(length(genes.list)), .combine=rbind) %dopar% {
        library(Matrix)
        library(GenomicRanges)
	    #permutation funciton is from jmuOutlier packages.
       .perm.cor.test <- function(x, y=NULL, alternative=c("two.sided", "less", "greater"), 
                          method=c("pearson", "spearman"), num.sim=20000 ) {
           # A permutation correlation test is performed.
           # P-values are approximated based on randomly sampling the permutations,
           #   where \code{x} and \code{y} are the paired vectors of observations.
           # 'x': Numeric vector of design variable.
           # 'y': Numeric vector of response variable.
           # 'alternative': A character string specifying the alternative hypothesis, and
           #   must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.  
           #   Only the initial letter needs to be specified.
           # 'num.sim': The upper limit on the number of permutations generated.
           # `method' is a character string specifying the type of correlation, and
           #   must be one of \code{"pearson"} (default) or \code{"spearman"}.  
           #   Only the initial letter needs to be specified.
           # `num.sim' is the number of permutations generated.
           # Examples:   x = c( 4, 6, 8, 11 ) ;   y = c( 19, 44, 15, 13 )
           #             perm.cor.test( x, y, "less", "pearson" ) 
           #             perm.cor.test( x, y, "less", "spearman" ) 
           if ( !is.numeric(x) )   stop("'x' must be numeric.")
           if ( !is.numeric(y) & !is.null(y) )   stop("'y' must be numeric or NULL.")
         #  alternative <- abbreviation(alternative, c("two.sided", "less", "greater"))
           if ( !(alternative %in% c("two.sided", "less", "greater")) ) 
                stop("'alternative' must be 'two.sided', 'less', or 'greater'.")
         #  method <- abbreviation(method, c("pearson", "spearman"))
           if ( !(method %in% c("pearson", "spearman")) ) 
                stop("'method' must be 'pearson' or 'spearman'.")
           if ( !is.numeric(num.sim) )  stop("'num.sim' must be numeric.")
           num.sim = floor(num.sim);   if ( num.sim < 1 )  stop("'num.sim' must be at least 1.")
           if (is.null(y))   {y <- x[,2] ;  x <- as.numeric(x[,1])}
           if (length(x)!=length(y)) stop("'x' and 'y' must have the same length.")
           if (method=="spearman")  {x <- rank(x);  y <- rank(y)}
           if (method %in% c("pearson", "spearman")) {
             test.stat0 <- cor(x,y);  test.stat <- NULL
             for (i in 1:num.sim)  {test.stat <- c(test.stat, cor(x, sample(y)))}
             if (alternative=="two.sided") p.value <- mean(abs(test.stat)>=abs(test.stat0))
             if (alternative=="less") p.value <- mean(test.stat<=test.stat0)
             if (alternative=="greater") p.value <- mean(test.stat>=test.stat0)   }
           output1 <- paste("Permutation correlation test.  Method is", method)
           output2 <- paste("p-value was estimated based on", num.sim, "simulations.")
           structure(list(output1, output2, alternative=alternative, p.value=p.value, cor.value=test.stat0))
           }

        gene.use <- genes.list[pos]
        gene.loci <- GenomicRanges::trim(suppressWarnings(GenomicRanges::promoters(GenomicRanges::resize(genes.coords[gene.use],
                                                                                                         width = 1, fix = "start"), upstream = max_dist, downstream = max_dist)))
        peaks.use <- S4Vectors::queryHits(GenomicRanges::findOverlaps(peaks.pos,gene.loci))
        if ((y <- length(peaks.use)) == 0L) {
		res.list <- c("NA", gene.use, "NA", "NA")
		res.list <- as.data.frame(res.list)
		res.list <- t(res.list)
		return(res.list)
        }else{
        res <- lapply(peaks.use, function(x){
		a = gene_counts[,gene.use]
		b = peak_counts[, x]
		cor.value <- cor(a,b, method=dist)
		p.value <- NULL
		res.p <- NULL
		if (cor.value >0)
		{
		 res.p <- .perm.cor.test(a, b, method=dist, alternative="greater", num.sim=num.sim)
		 p.value <- res.p$p.value
		}else{
		 res.p <- .perm.cor.test(a, b, method=dist, alternative="less", num.sim=num.sim)
		 p.value <- res.p$p.value
		}
		res.list <- c(colnames(peak_counts)[x], gene.use, cor.value, p.value)
		res.list <- as.data.frame(res.list)
		res.list <- t(res.list)
		return(res.list)
		})
		res <- do.call(rbind, res)
		return(res)
	}
    }
#     elements <- do.call(rbind, elements)
     stopCluster(cl)
     colnames(elements) <- c("peak_id","gene_id","cor.value","p.value")
	 elements <- as.data.frame(elements)
	 elements <- elements[elements$p.value != "NA",]
	 #elements <- elements[elements$p.value < alpha,]
	 elements$FDR <- p.adjust(elements$p.value, method = "BH")
	return(elements)
}
	
#----------linkGenesAndPeaks4-----------
# out put correlation score and p value.
# using pvalue from permutation test.
# make a parallel version.
# v5 test the given ACR and gene pairs
linkGenesAndPeaks4 <- function (gene_counts, peak_counts, peak_genes_pairs = NULL, dist = "spearman",
                                num.sim=10000, cores=24){
	#take overlaps
    # overlap.g <- intersect(rownames(gmat), acr_pairs$gene_id)
    # gmat <- gmat[rownames(gmat) %in% overlap.g,]
    # acr_pairs <- acr_pairs[acr_pairs$gene_id %in% overlap.g,]
    # overlap.p <- intersect(rownames(pmat), acr_pairs$peak_id)
    # pmat <- pmat[rownames(pmat) %in% overlap.p,]
    # acr_pairs <- acr_pairs[acr_pairs$peak_id %in% overlap.p,]
	gene_counts <- t(gene_counts)
    peak_counts <- t(peak_counts)
    print("Calculating correlation for gene-peak pairs...")
	cl <- makeCluster(cores)
	registerDoParallel(cl)
#    elements <- lapply(seq(length(genes.list)), function(pos){
	elements <- foreach( x = peak_genes_pairs, .combine=rbind) %dopar% {
        #library(Matrix)
	    #permutation funciton is from jmuOutlier packages.
       .perm.cor.test <- function(x, y=NULL, alternative=c("two.sided", "less", "greater"), 
                          method=c("pearson", "spearman"), num.sim=20000 ) {
           # A permutation correlation test is performed.
           # P-values are approximated based on randomly sampling the permutations,
           #   where \code{x} and \code{y} are the paired vectors of observations.
           # 'x': Numeric vector of design variable.
           # 'y': Numeric vector of response variable.
           # 'alternative': A character string specifying the alternative hypothesis, and
           #   must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.  
           #   Only the initial letter needs to be specified.
           # 'num.sim': The upper limit on the number of permutations generated.
           # `method' is a character string specifying the type of correlation, and
           #   must be one of \code{"pearson"} (default) or \code{"spearman"}.  
           #   Only the initial letter needs to be specified.
           # `num.sim' is the number of permutations generated.
           # Examples:   x = c( 4, 6, 8, 11 ) ;   y = c( 19, 44, 15, 13 )
           #             perm.cor.test( x, y, "less", "pearson" ) 
           #             perm.cor.test( x, y, "less", "spearman" ) 
           if ( !is.numeric(x) )   stop("'x' must be numeric.")
           if ( !is.numeric(y) & !is.null(y) )   stop("'y' must be numeric or NULL.")
         #  alternative <- abbreviation(alternative, c("two.sided", "less", "greater"))
           if ( !(alternative %in% c("two.sided", "less", "greater")) ) 
                stop("'alternative' must be 'two.sided', 'less', or 'greater'.")
         #  method <- abbreviation(method, c("pearson", "spearman"))
           if ( !(method %in% c("pearson", "spearman")) ) 
                stop("'method' must be 'pearson' or 'spearman'.")
           if ( !is.numeric(num.sim) )  stop("'num.sim' must be numeric.")
           num.sim = floor(num.sim);   if ( num.sim < 1 )  stop("'num.sim' must be at least 1.")
           if (is.null(y))   {y <- x[,2] ;  x <- as.numeric(x[,1])}
           if (length(x)!=length(y)) stop("'x' and 'y' must have the same length.")
           if (method=="spearman")  {x <- rank(x);  y <- rank(y)}
           if (method %in% c("pearson", "spearman")) {
             test.stat0 <- cor(x,y);  test.stat <- NULL
             for (i in 1:num.sim)  {test.stat <- c(test.stat, cor(x, sample(y)))}
             if (alternative=="two.sided") p.value <- mean(abs(test.stat)>=abs(test.stat0))
             if (alternative=="less") p.value <- mean(test.stat<=test.stat0)
             if (alternative=="greater") p.value <- mean(test.stat>=test.stat0)   }
           output1 <- paste("Permutation correlation test.  Method is", method)
           output2 <- paste("p-value was estimated based on", num.sim, "simulations.")
           structure(list(output1, output2, alternative=alternative, p.value=p.value, cor.value=test.stat0))
           }
		x <- unlist(strsplit(x,":"))
        peaks.use <- x[1]
        gene.use <- x[2]
		
		a = gene_counts[,gene.use]
		b = peak_counts[, peaks.use]
		cor.value <- cor(a,b, method=dist)
		p.value <- NULL
		res.p <- NULL
		if (cor.value >0)
		{
		 res.p <- .perm.cor.test(a, b, method=dist, alternative="greater", num.sim=num.sim)
		 p.value <- res.p$p.value
		}else{
		 res.p <- .perm.cor.test(a, b, method=dist, alternative="less", num.sim=num.sim)
		 p.value <- res.p$p.value
		}
		res.list <- c(peaks.use, gene.use, cor.value, p.value)
		res.list <- as.data.frame(res.list)
		res.list <- t(res.list)
		return(res.list)
	}
     stopCluster(cl)

     colnames(elements) <- c("peak_id","gene_id","cor.value","p.value")
	 elements <- as.data.frame(elements)
	 #elements <- elements[elements$p.value != "NA",]
	 #elements <- elements[elements$p.value < alpha,]
	 elements$FDR <- p.adjust(elements$p.value, method = "BH")
	return(elements)
}

#-------calculate gene expression specificity-----
fTau <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  } 
  return(res)
}



