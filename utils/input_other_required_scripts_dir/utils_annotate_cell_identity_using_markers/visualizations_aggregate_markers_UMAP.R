#########################
## per cell annotation ##
#########################

##updating 040325 we will make a version that allow the duplicate gene
##updating 010525 make the umap and store the meta into obejct and also store the meta directly in the output
##updating 110424 
##updating 110221 add the output dir
##updating 081821 replace the knn and log using the smooth data
##updating 080221
##also use the nflt markers for all the other methods
##updating 080321 modify smooth data to increase the performance of showing the annotation

# load arguments ----------------------------------------------------------------------------------

# enable command line arguments
args <- commandArgs(T)

# usage
#if(length(args) != 7){stop("Rscript annotate_cells.R <gene_accessibility> <meta> <markers> <svd> <prefix> <output_dir>")}

# arguments
input_soc_object_fl <- as.character(args[1])

#ga <- as.character(args[1])
#meta <- as.character(args[2])
markers <- as.character(args[2])
#svd <- as.character(args[4])
prefix <- as.character(args[3])
openAnnot <- as.character(args[4])
output_dir <- as.character(args[5])
input_prefix <- as.character(args[6])


# load libraries ----------------------------------------------------------------------------------
library(parallel)
library(doMC)
library(Matrix)
library(irlba)
library(glmnet)
library(FNN)
library(RANN)
library(kknn)
library(gtools)
library(RColorBrewer)
library(png)
library(scales)


# load functions ----------------------------------------------------------------------------------

# loading data into obj
loadData <- function(input_soc_object,
                     markers){
    
    #input_soc_object <- readRDS(input_soc_object_fl)
  
    #message(gene_access)
    #message(metadata)
    message(markers)
    #message(pcs)
    
    # files
    message(" - loading data ...")
    
    #x <- read.table(gene_access)
    #y <- read.table(metadata)
    y <- input_soc_object$meta
    
    ##updating 061923
    #y <- read.delim(metadata,row.names = 1)
    yy <- y
    #z <- read.table(markers, header=T)
    z <- read.delim(markers,header=T)
    
    ##updating 040325 we will not allow them to be the row name that will avoid the duplicate
    #rownames(z) <- z$geneID
    
    #pcs <- read.table(pcs)
    pcs <- input_soc_object$svd
    
    # make sure columns 1-2 are factors and 3 is numeric
    #message(" - reformating input data ...")
    #x[,1] <- as.factor(x[,1])
    #x[,2] <- as.factor(x[,2])
    
    # convert 2 sparse matrix
    #sm <- sparseMatrix(i=as.numeric(x[,1]), 
    #                   j=as.numeric(x[,2]), 
    #                   x=as.numeric(x[,3]), 
    #                   dimnames=list(levels(x[,1]), levels(x[,2])))
    #sm <- sm %*% Diagonal(x=10000/Matrix::colSums(sm))
    
    #sm <- readRDS(gene_access)
    sm <- input_soc_object$gene_acc
    
    ids <- colnames(sm)
    
    # verbose
    message("     * loaded gene accessibility matrix | n.rows = ", nrow(sm), " | n.cols = ", ncol(sm))
    
    # get shared cells with meta data
    message(" - cleaning input data by meta.data...")
    shared.cells <- intersect(colnames(sm), rownames(y))
    sm <- sm[,shared.cells]
    sm <- sm[Matrix::rowSums(sm)>0,]
    y <- y[shared.cells,]
    yy <- yy[shared.cells,]
    message("     * filtered gene accessibility matrix | n.rows = ", nrow(sm), " | n.cols = ", ncol(sm))
    
    # filter markers
    z <- z[as.character(z$geneID) %in% rownames(sm),]
    
    # check that LouvainClusters is a column
    if(! c("LouvainClusters") %in% colnames(y)){
        y$LouvainClusters <- y$iNMF_clusters
        yy$LouvainClusters <- yy$iNMF_clusters
    }
    
    # return
    return(list(sm=sm, meta=y, markers=z, original.meta=yy, pcs=pcs[rownames(yy),]))
    
}

# filter data
filterMarkers <- function(obj, 
                          enrich.threshold=0.5, 
                          permutations=100, 
                          k=50, 
                          threads=1){
    
    # entropy eq
    .entropy <- function(p){
            if (min(p) < 0 || sum(p) <= 0) return(NA)
            p.norm <- p[p>0]/sum(p)
            -sum(log2(p.norm)*p.norm)
    }
    
    # get centroids
    clusts <- unique(obj$meta$LouvainClusters)
    cents <- lapply(clusts, function(z){
        cl <- subset(obj$meta, obj$meta$LouvainClusters==z)
        x.val <- mean(cl$umap1)
        y.val <- mean(cl$umap2)
        dists <- lapply(seq(1:nrow(cl)), function(x){
            sqrt(((cl[x,]$umap1 - x.val)^2) + ((cl[x,]$umap2 - y.val)^2))
        })
        dists <- do.call(c, dists)
        rownames(cl)[which.min(dists)]
        
    })
    cents <- do.call(c, cents)
    
    # get matrix of UMAP coords
    umap.mat <- as.matrix(obj$meta[,c("umap1", "umap2")])
    rownames(umap.mat) <- rownames(obj$meta)
    
    # bin cells into groups by KNN
    knn.out <- get.knn(umap.mat, k=50)$nn.index
    rownames(knn.out) <- rownames(umap.mat)
    colnames(knn.out) <- paste0("KNN", seq(1:ncol(knn.out)))
    knn.out <- apply(knn.out, 2, function(z){
        rownames(knn.out)[z]
    })
    rownames(knn.out) <- rownames(umap.mat)
    cents.knn <- knn.out[cents,]
    cents.ga <- lapply(cents, function(z){
        Matrix::rowMeans(obj$sm[,cents.knn[z,]])
    })
    cents.ga <- do.call(cbind, cents.ga)
    colnames(cents.ga) <- clusts
    # cents.ga <- cents.ga * 100
    
    # estimate entropy across cluster
    se.score <- lapply(1:nrow(cents.ga), function(z){
        .entropy(cents.ga[z,])
    })
    se.score <- do.call(c, se.score)
    names(se.score) <- rownames(cents.ga)
    se.score[is.na(se.score)] <- mean(se.score, na.rm=T)
    
    # get entropy of x random genes
    #perms <- mclapply(seq(1:permutations), function(z){
    #    rand.genes <- sample(se.score, nrow(obj$markers))
    #    mean(rand.genes)
    #}, mc.cores=threads)
    #perms <- do.call(c, perms)
    ave.perm <- mean(se.score, na.rm=T)
    std.perm <- sd(se.score, na.rm=T)
    
    # estimate enrichment
    marker <- obj$markers
    marker$se.score.enrich <- (se.score[as.character(marker$geneID)] - ave.perm)/std.perm
    marker$se.score <- se.score[as.character(marker$geneID)]
    
    # filter markers
    marker <- marker[which(marker$se.score.enrich < enrich.threshold),]
    obj$sm <- obj$sm[unique(as.character(marker$geneID)),]
    obj$sm <- obj$sm[Matrix::rowSums(obj$sm)>0,]
    obj$sm <- obj$sm[,Matrix::colSums(obj$sm)>0]
    obj$markers <- marker[as.character(marker$geneID) %in% rownames(obj$sm),]
    obj$meta <- obj$meta[colnames(obj$sm),]
    obj$se.scores <- pnorm((se.score - ave.perm)/(std.perm))
    
    # return
    return(obj)
    
}

# cell identity enrichment
cellIdentityEnrichment <- function(obj, 
                                   permutations=1000, 
                                   threads=1, 
                                   weight=T){
    
    # transform to z-score (dense matrix)
    r.mat <- as.matrix(t(scale(t(as.matrix(obj$sm)))))
    #r.mat <- t(apply(as.matrix(obj$sm), 1, function(z){
        #ranks <- rank(z)
        #ranks/sum(ranks)
    #}))
    
    # filter markers
    if(weight){
        
        # check marker correlation
        weights <- lapply(unique(obj$markers$type), function(y){
            ct.markers <- subset(obj$markers, obj$markers$type==y)
            if(nrow(ct.markers) < 3){
                cor.vals <- rep(1, nrow(ct.markers))
                names(cor.vals) <- rownames(ct.markers)
                return(cor.vals)
            }else{
                cors <- cor(t(z.mat[as.character(ct.markers$geneID),]), method="spearman")
                diag(cors) <- rep(0, length(diag(cors)))
                cors[cors > 0] <- 1
                cors[cors < 0] <- 0
                cor.vals <- colSums(cors, na.rm=T)
                names(cor.vals) <- colnames(cors)
                return(cor.vals)
            }
        })
        weights <- do.call(c, weights)
        weights <- weights[weights > 0]
        shared <- intersect(rownames(z.mat), names(weights))
        z.mat <- z.mat[shared,]
        
        # filter marker genes
        obj$markers <- obj$markers[as.character(obj$markers$geneID) %in% shared,]
        z.mat <- z.mat[,colSums(obj$sm[shared,])>0]
    }
    
    # iterate over each cell type
    message(" - begin per cell annotation...")
    types <- unique(obj$markers$type)
    print(types)
    outs <- lapply(types, function(z){
        
        # verbose
        message("      * current cell type: ", z)
        
        # get gene for/against cell type
        m.genes <- as.character(subset(obj$markers, obj$markers$type == z)$geneID)
        nm.genes <- as.character(subset(obj$markers, obj$markers$type != z)$geneID)
        
        # num.genes each
        num.m.genes <- length(m.genes)
        
        # estimate enrichment
        if(num.m.genes > 1){
            enrich <- colMeans(r.mat[m.genes,])
        }else{
            enrich <- r.mat[m.genes,]
        }
        
        # permute enrichment 
        perms <- mclapply(seq(1:permutations), function(y){

            #rand.genes <- sample(nm.genes, num.m.genes)
 	    if ((length(nm.genes)-num.m.genes) < 0){
              rand.genes <- sample(nm.genes, num.m.genes,replace = T)
            }else{
              rand.genes <- sample(nm.genes, num.m.genes)
            }

            if(num.m.genes > 1){
                colMeans(r.mat[rand.genes,])
            }else{
                r.mat[rand.genes,]
            }
        }, mc.cores=threads)
        perms <- as.matrix(do.call(rbind, perms))
        
        # get Z-score over observed
        ave.perms <- colMeans(perms)
        std.perms <- apply(perms, 2, sd)
        
        # cell type enrichment
        ct.enrich <- (enrich - ave.perms)/(std.perms)
        
        # return
        return(ct.enrich)
    })
    outs <- do.call(rbind, outs)
    rownames(outs) <- types
    outs <- t(outs)
    return(outs)
    
}

# smooth identities using diffusion graph
smoothIdentity <- function(obj, 
                           idents, 
                           normIdents=F, 
                           embedding="PCs",
                           step=3,
                           k=15,
                           npcs=20){
    
    # functions
    .smooth.data       <- function(x, obj, k=25, step=3, npcs=30, embedding="PCs"){
        
        # verbose
        message(" - formating input data ...")
        
        # input
        x <- as.data.frame(x)
        x <- x[colnames(obj$sm),]
        rownames(x) <- colnames(obj$sm)
        x[is.na(x)] <- 0
        x <- as.matrix(x)
        data.use <- t(Matrix(x, sparse=T))
        
        # check if pcs are given
        message(" - collecting embedding ...")
        if(embedding=="UMAP"){
            
                # get UMAP coords
                pcs <- obj$meta[,c("umap1","umap2")]
                colnames(pcs) <- c("PC_1", "PC_2")
            
         }else if(embedding=="PCs"){
            
             # check if PCs are in object
             if(npcs > ncol(obj$pcs)){
                 npcs <- ncol(obj$pcs)
             }
             pcs <- obj$pcs[colnames(data.use),1:npcs]
             
         }else{
             
                
             # verbose
             message(" - generating new PCs ...")
            
             # get PCS
             pc <- irlba(t(obj$sm), nv=npcs)
             pcs <- pc$u %*% diag(pc$d)
             rownames(pcs) <- colnames(data.use)
             colnames(pcs) <- paste0("PC_", seq(1:ncol(pcs)))
                
             # do l2-norm
             pcs <- t(apply(pcs, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
        }
        
        # get KNN
        message(" - building graph ...")
        knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
        j <- as.numeric(x = t(x = knn.graph))
        i <- ((1:length(x = j)) - 1) %/% k + 1
        edgeList = data.frame(i, j, 1)
        A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
        
        
        ##updating 112625 debug
        if (nrow(A) != ncol(A)){
          n <- nrow(pcs)
          A <- sparseMatrix(
            i = edgeList[,1],
            j = edgeList[,2],
            x = edgeList[,3],
            dims = c(n, n)
          )
        }
        
        
        # Smooth graph
        message(" - smoothing identities ...")
        A = A + t(A)
        A = A / Matrix::rowSums(A)
        step.size = step
        if(step.size > 1){
            for(i in 1:step.size){
                A = A %*% A
            }
        }
        
        # smooth data
        impute.activity <- t(A %*% t(data.use))
        colnames(impute.activity) <- colnames(data.use)
        rownames(impute.activity) <- rownames(data.use)
        
        # normalize scores
        norm.ann <- as.matrix(t(impute.activity))
        norm.ann <- t(apply(norm.ann, 1, function(z){
            z / sum(z)
        }))
        
        # return sparse Matrix
        return(norm.ann)
    }
    
    # normalize enriched scores to positive values
    idents <- t(apply(idents, 1, function(z){
        if(normIdents){
            z <- z - min(z, na.rm=T)        
            z <- z/sum(z)
        }else{
            z <- z - min(z, na.rm=T)
        }
        return(z)
    }))
    
    # smooth
    smooth.celltype <- .smooth.data(idents, obj, k=k, step=step, npcs=npcs, embedding=embedding)
    
    # return
    return(smooth.celltype)
    
    
}

# train classifier on whole-genome patterns of the top X% of cell annotations (by cell type)
trainIdentity <- function(obj, 
                          idents, 
                          perc=0.05, 
                          threads=1, 
                          smoothdata=F){
    
    # functions
    .normalizeIdentity <- function(idt){
        
        # normalize identities
        outs <- apply(idt, 1, function(z){
            z <- z - (min(z))
            z/sum(z)
        })
        
        # return
        return(t(outs))
    }
    .smooth.data       <- function(obj, k=15, step=3, npcs=30, embedding="PCs"){
        
        # verbose
        message("   * formating input data ...")
        
        # input
        data.use <- obj$sm
        message("   * ncol = ", ncol(data.use), " | nrow = ", nrow(data.use))
        
        # check if pcs are given
        if(embedding=="UMAP"){
            
            # get UMAP coords
            message("   * UMAP collecting embedding ...")
            pcs <- obj$meta[,c("umap1","umap2")]
            colnames(pcs) <- c("PC_1", "PC_2")
            
        }else if(embedding=="PCs"){
            
            # check if PCs are in object
            if(! c("pcs") %in% names(obj)){
                stop("   ! slot 'pcs' missing from input object ...")
            }
            if(npcs > ncol(obj$pcs)){
                npcs <- ncol(obj$pcs)
            }
            message("   * collecting SVD/PCA embedding ...")
            pcs <- obj$pcs[colnames(data.use),1:npcs]
            
        }else{
            
            # verbose
            message("   * generating new PCs ...")
            
            # get PCS
            pc <- irlba(t(obj$sm), nv=npcs)
            pcs <- pc$u %*% diag(pc$d)
            rownames(pcs) <- colnames(data.use)
            colnames(pcs) <- paste0("PC_", seq(1:ncol(pcs)))
            
            # do l2-norm
            pcs <- t(apply(pcs, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
        }
        
        # get KNN
        message("   * building graph ...")
        knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
        j <- as.numeric(x = t(x = knn.graph))
        i <- ((1:length(x = j)) - 1) %/% k + 1
        edgeList = data.frame(i, j, 1)
        A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
        
        # Smooth graph
        message("   * smoothing profiles ...")
        A = A + t(A)
        A = A / Matrix::rowSums(A)
        step.size = step
        if(step.size > 1){
            for(i in 1:step.size){
                A = A %*% A
            }
        }
        
        # smooth data
        impute.activity <- t(A %*% t(data.use))
        colnames(impute.activity) <- colnames(data.use)
        rownames(impute.activity) <- rownames(data.use)
        
        # return sparse Matrix
        return(impute.activity)
    }
    
    # if parallel
    if(threads > 1){
        doPar <- T
        require(doMC)
        if(10 < threads){
            threads <- 10
        }
        registerDoMC(cores=threads)
    }else{
        doPar <- F
    }
    
    # normalize identity scores
    norm.idents <- as.data.frame(.normalizeIdentity(idents))
    
    # get top cells
    message(" - finding cells to train classifier ...")
    outs1 <- lapply(seq(1:nrow(norm.idents)), function(zz){
        z <- norm.idents[zz,]
        ct <- names(z)[which.max(z)]
        score <- max(z)
        data.frame(celltype=ct, score=score)
    })
    outs1 <- do.call(rbind, outs1)
    rownames(outs1) <- rownames(norm.idents)
    
    # use raw scores for picking cells
    outs <- lapply(seq(1:nrow(idents)), function(zz){
        z <- idents[zz,]
        ct <- names(z)[which.max(z)]
        score <- max(z)
        data.frame(celltype=ct, score=score)
    })
    outs <- do.call(rbind, outs)
    rownames(outs) <- rownames(idents)
    
    # get top X% of each cell type
    types <- unique(outs$celltype)
    top.cells <- lapply(types, function(z){
        df <- subset(outs, outs$celltype==z)
        cut.off <- quantile(df$score, c(1-perc))
        pass <- subset(df, df$score > cut.off)
        if(nrow(pass) < 20 & nrow(df) >= 20){
            df <- df[order(df$score, decreasing=T),]
            pass <- head(df, n=20)
        }else if(nrow(pass) < 20 & nrow(df) < 20){
            pass <- df
        }
        return(pass)
    })
    top.cells <- do.call(rbind, top.cells)
    top.cell.cnts <- table(top.cells$celltype)
    top.cell.cnts <- top.cell.cnts[top.cell.cnts > 20]
    keep.celltypes <- names(top.cell.cnts)
    top.cells <- top.cells[as.character(top.cells$celltype) %in% keep.celltypes,]
    top.cells$celltype <- factor(as.character(top.cells$celltype))
    
    # smooth profiles?
    if(smoothdata){
        ga.use <- .smooth.data(obj)
    }else{
        ga.use <- obj$sm
    }
    
    # get gene accessibility for top cells
    ga.top <- ga.use[,rownames(top.cells)]
    ga.top <- t(ga.top[Matrix::rowSums(ga.top)>0,])
    ct.top <- factor(top.cells$celltype)
    ct.top <- droplevels(ct.top)
    
    # cross-validated model
    message(" - training logistic classifier with cross-validation ...")
    cvfit <- cv.glmnet(ga.top, ct.top, 
                       family = "multinomial", 
                       alpha = 1,
                       parallel = doPar,
                       trace.it = TRUE)
    
    # predict
    message(" - predicting celltypes from trained classifier ...")
    predicted <- as.data.frame(predict(cvfit, newx = t(ga.use[colnames(ga.top),]), s = "lambda.min", type = "response"))
    predicted <- as.matrix(predicted)
    predicted[is.na(predicted)] <- 0
    predicted[is.infinite(predicted) & predicted > 0] <- 1
    predicted[is.infinite(predicted) & predicted < 0] <- 0
    colnames(predicted) <- gsub("\\.1","",colnames(predicted))
    
    # get top cell identities
    outs2 <- lapply(seq(1:nrow(predicted)), function(zz){
        z <- predicted[zz,]
        ct <- names(z)[which.max(z)]
        score <- max(z)
        data.frame(celltype=ct, score=score)
    })
    outs2 <- do.call(rbind, outs2)
    rownames(outs2) <- rownames(predicted)
    
    # join
    colnames(outs2) <- c("celltype_glmnet", "probability_glmnet") 
    colnames(outs1) <- c("celltype_enrich", "probability_enrich")
    shared <- intersect(rownames(outs1), rownames(outs2))
    outs1 <- outs1[shared,]
    outs2 <- outs2[shared,]
    cts <- cbind(outs1, outs2)
    
    # report accuracy between enrich and predicted classifications
    message(" - accuracy = ", mean(as.character(cts$celltype_enrich)==as.character(cts$celltype_glmnet)))
    
    # return
    return(predicted)
    
}

# weighted knn on filtered annotations to predicted missing cell types
knnCellTypes <- function(obj,
                         idents,
                         train.top = 0.05,
                         k=25){
    
    # select cells to train on
    t.cells <- lapply(seq(1:nrow(idents)), function(zz){
        z <- idents[zz,]
        ct <- names(z)[which.max(z)]
        score <- max(z)
        data.frame(celltype=ct, score=score)
    })
    t.cells <- do.call(rbind, t.cells)
    rownames(t.cells) <- rownames(idents)
    
    # get top X% of each cell type
    types <- unique(t.cells$celltype)
    top.cells <- lapply(types, function(z){
        df <- subset(t.cells, t.cells$celltype==z)
        cut.off <- quantile(df$score, c(1-train.top))
        pass <- subset(df, df$score > cut.off)
        if(nrow(pass) < 20 & nrow(df) >= 20){
            df <- df[order(df$score, decreasing=T),]
            pass <- head(df, n=20)
        }else if(nrow(pass) < 20 & nrow(df) < 20){
            pass <- df
        }
        return(pass)
    })
    top.cells <- do.call(rbind, top.cells)
    not.top.cells <- rownames(idents)[!rownames(idents) %in% rownames(top.cells)]
    
    # get PCs
    t.pcs <- as.data.frame(obj$pcs[rownames(top.cells),])
    n.pcs <- as.data.frame(obj$pcs[not.top.cells,])
    t.pcs$celltype <- factor(top.cells$celltype)
    
    # run weighted knn
    knn.pred <- kknn(celltype~., train=t.pcs, test=n.pcs, k=k)

    # annotations
    pred.types <- knn.pred$prob
    rownames(pred.types) <- not.top.cells
    
    # shared cell types
    known.cells <- idents[rownames(top.cells),]
    known.cells <- t(apply(known.cells, 1, function(x){
        ifelse(x==max(x, na.rm=T), 1, 0)
    }))
    known.cells <- known.cells[,colnames(pred.types)]
    
    # output
    out <- rbind(known.cells, pred.types)
    
    # return
    return(out)
    
}

# plot cell identities by cell type
plotIdentities <- function(idents, 
                           obj, 
                           verbose=T, 
                           prefix="cell_identities.png", 
                           u.lim=0.99){
    
    # number of rows
    ncols <- 5
    nrows <- ceiling(ncol(idents)/5)
    
    # set up plot
    png(file=paste0(output_dir,'/',prefix), width=ncols*3, height=nrows*3, units="in", res=600, type="cairo")
    layout(matrix(c(1:(ncols*nrows)), ncol=ncols, nrow=nrows, byrow=T))
    par(mar=c(0.5,0.5,0.5,0.5))
    
    # color palette
    cols <- colorRampPalette(c("grey80","grey75", brewer.pal(9, "RdPu")[4:9]), bias=0.5)(100)
    
    # plot each idents
    for (i in mixedsort(colnames(idents))){
        
        # verbose
        if(verbose){
            message(" - plotting annotation scores for: ",i)
        }
        
        # extract coordinates and activity
        shared <- intersect(rownames(obj$meta), rownames(idents))
        if(length(shared) < nrow(obj$meta)){
            idents <- as.data.frame(idents)
            idents[rownames(obj$meta),]
            idents[is.na(idents)] <- 0
        }
        umap.coords <- obj$meta[rownames(idents),c("umap1","umap2")]
        umap.coords$celltype <- as.numeric(idents[,i])
        ct.range <- quantile(umap.coords$celltype, u.lim)
        umap.coords$celltype[umap.coords$celltype > ct.range] <- ct.range
        umap.coords <- umap.coords[order(umap.coords$celltype, decreasing=F),]
        
        # set up range
        ct.std <- sd(umap.coords$celltype, na.rm=T)
        min.ct <- min(umap.coords$celltype, na.rm=T) - (0.01*ct.std)
        max.ct <- max(umap.coords$celltype, na.rm=T) + (0.01*ct.std)
        
        # plot
        plot(umap.coords$umap1, umap.coords$umap2, 
             pch=16, 
             cex=0.25, 
             col=cols[cut(umap.coords$celltype, breaks=seq(from=min.ct, to=max.ct, length.out=101))],
             main=i,
             axes=F,
             xlab="",
             ylab="",
             bty="n")
        
    }
    
    # device off
    dev.off()
    
}

# average scores
aveScores <- function(scores.list,
                      norm=T,
                      add=T){
    
    # get all ids
    columns <- lapply(scores.list, function(x){
        colnames(x)
    })
    columns <- do.call('c', columns)
    columns <- unique(columns)
    rows <- lapply(scores.list, function(x){
        rownames(x)
    })
    rows <- do.call('c', rows)
    rows <- unique(rows)
    
    # reformat
    scores <- lapply(scores.list, function(df){
        if(sum(as.matrix(df)[1,]) != 1){
            df <- as.matrix(df)
            df <- t(apply(df, 1, function(x){
                if(any(x < 0)){
                    x <- x - min(x, na.rm=T)
                }
                x/sum(x)
            }))
            df <- data.frame(df)
        }
        dif.cols <- columns[!columns %in% colnames(df)]
        dif.rows <- rows[!rows %in% rownames(df)]
        if(length(dif.cols) > 0){
            for(i in dif.cols){
                df[i] <- 0
            }
        }
        if(length(dif.rows) > 0){
            for(i in dif.rows){
                vec <- rep(0, ncol(df))
                df <- rbind(df, vec)
                rownames(df)[nrow(df)] <- i
            }
        }
        df[is.na(df)] <- 0
        rownames(df) <- rows
        colnames(df) <- columns
        return(df)
    })
    total <- Reduce('+', scores)
    ave <- total/length(scores)
    ave <- t(apply(ave, 1, function(z){
        z/sum(z)
    }))
    
    # return
    return(ave)
    
}

# return results 
updateMeta <- function(obj, 
                       smooth.idents, 
                       knn.idents,
                       pred.idents, 
                       ext.cluster=T, 
                       cutoff=c(0.5,0.5,0.5,0.5)){
    
    # assign identities (glmnet)
    calls <- apply(pred.idents,1, function(x){
        val <- max(x)
        if(val > cutoff[4]){
            return(names(x)[which.max(x)])
        }else{
            return("unknown")
        }
    })
    obj$original.meta$cell_annotation_glmnet <- calls[rownames(obj$original.meta)]
    obj$original.meta$cell_annotation_glmnet[is.na(obj$original.meta$cell_annotation_glmnet)] <- "unknown"
    
    # assign identities (knn)
    calls <- apply(knn.idents,1, function(x){
        val <- max(x)
        if(val > cutoff[4]){
            return(names(x)[which.max(x)])
        }else{
            return("unknown")
        }
    })
    obj$original.meta$cell_annotation_knn <- calls[rownames(obj$original.meta)]
    obj$original.meta$cell_annotation_knn[is.na(obj$original.meta$cell_annotation_knn)] <- "unknown"
    
    ##updating 010525 do not assign the enrich identity
    # assign identities (enrich)
    #calls <- apply(enrich.idents,1, function(x){
    #    val <- max(x)
    #    if(val > cutoff[4]){
    #        return(names(x)[which.max(x)])
    #    }else{
    #        return("unknown")
    #    }
    #})
    #obj$original.meta$cell_annotation_enrich <- calls[rownames(obj$original.meta)]
    #obj$original.meta$cell_annotation_enrich[is.na(obj$original.meta$cell_annotation_enrich)] <- "unknown"
    
    # assign identities (smooth)
    calls <- apply(smooth.idents,1, function(x){
        val <- max(x)
        if(val > cutoff[4]){
            return(names(x)[which.max(x)])
        }else{
            return("unknown")
        }
    })
    obj$original.meta$cell_annotation_smooth <- calls[rownames(obj$original.meta)]
    obj$original.meta$cell_annotation_smooth[is.na(obj$original.meta$cell_annotation_smooth)] <- "unknown"
    
    # extend to cluster
    if(ext.cluster){
        obj$original.meta$LouvainClusters <- as.character(obj$original.meta$LouvainClusters)
        
        dat <- subset(obj$original.meta, obj$original.meta$cell_annotation_glmnet!="unknown")
        props <- prop.table(table(dat$LouvainClusters, dat$cell_annotation_glmnet), 1)
        calls <- apply(props, 1, function(x) names(x)[which.max(x)])
        p.cells <- apply(props, 1, max)
        obj$original.meta$cluster_annotation_glmnet <- calls[as.character(obj$original.meta$LouvainClusters)]
        
        dat <- subset(obj$original.meta, obj$original.meta$cell_annotation_knn!="unknown")
        props <- prop.table(table(dat$LouvainClusters, dat$cell_annotation_knn), 1)
        calls <- apply(props, 1, function(x) names(x)[which.max(x)])
        p.cells <- apply(props, 1, max)
        obj$original.meta$cluster_annotation_knn <- calls[as.character(obj$original.meta$LouvainClusters)]
        
        dat <- subset(obj$original.meta, obj$original.meta$cell_annotation_smooth!="unknown")
        props <- prop.table(table(dat$LouvainClusters, dat$cell_annotation_smooth), 1)
        calls <- apply(props, 1, function(x) names(x)[which.max(x)])
        p.cells <- apply(props, 1, max)
        obj$original.meta$cluster_annotation_smooth <- calls[as.character(obj$original.meta$LouvainClusters)]
        
        dat <- subset(obj$original.meta, obj$original.meta$cell_annotation_enrich!="unknown")
        props <- prop.table(table(dat$LouvainClusters, dat$cell_annotation_enrich), 1)
        calls <- apply(props, 1, function(x) names(x)[which.max(x)])
        p.cells <- apply(props, 1, max)
        obj$original.meta$cluster_annotation_enrich <- calls[as.character(obj$original.meta$LouvainClusters)]
    }
    
    # return
    return(obj$original.meta)
    
}


# ANALYSIS ----------------------------------------------------------------------------------------

# load data #####################################
#ga.obj <- loadData(ga, meta, markers, svd)

input_soc_object <- readRDS(input_soc_object_fl)
ga.obj <- loadData(input_soc_object,markers)


#saveRDS(ga.obj,paste0(output_dir,'/temp_store_obj.rds'))

# filter markers ################################
#ga.obj.filt <- filterMarkers(ga.obj, enrich.threshold=1, threads=30)

if (openAnnot == 'yes'){

  ##open when necessary
  # estimate cell identity enrichment #############
  ga.identities <- cellIdentityEnrichment(ga.obj, permutations=1000, threads=30, weight=F)
  #write.table(ga.identities, file=paste0(output_dir,'/',prefix,".enriched_celltypes.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  ##do not plot the enriched cell types
  #plotIdentities(ga.identities, ga.obj, prefix=paste0(prefix,".cell_identities.enriched.png"), u.lim=0.99)
  
  # smooth cell identity by neighborhood ##########
  ga.smooth.identities <- smoothIdentity(ga.obj, ga.identities, k=15, step=3, npcs=30)
  #write.table(ga.smooth.identities, file=paste0(output_dir,'/',prefix,".smoothed_celltypes.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  plotIdentities(ga.smooth.identities, ga.obj, prefix=paste0(prefix,".cell_identities.smooth.png"), u.lim=0.999)
  
  ##updating 080321 add a version to allow the value/max
  #smooth_dt <- read.table(paste0(output_dir,'/',prefix,".smoothed_celltypes.txt"))
  #ga.smooth.identities.norm = smooth_dt/apply(smooth_dt,1,max)
  #message('plot norm smooth data of ori')
  #plotIdentities(ga.smooth.identities.norm, ga.obj, prefix=paste0(prefix,".cell_identities_norm.smooth.png"), u.lim=0.9999)
  
  # knn identities ################################
  ga.knn.identities <- knnCellTypes(ga.obj, ga.identities, train.top=0.1, k=30)
  ga.knn.smooth <- smoothIdentity(ga.obj, ga.knn.identities, k=15, step=3, npcs=30)
  #write.table(ga.knn.smooth, file=paste0(output_dir,'/',prefix,".knn_celltypes.smooth.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  plotIdentities(ga.knn.smooth, ga.obj, prefix=paste0(prefix,".cell_identities.wknn.png"), u.lim=0.999)
  
  ##updating 080321 add a version to do the normalization
  #smooth_dt <- read.table(paste0(output_dir,'/',prefix,".knn_celltypes.smooth.txt"))
  #ga.knn.smooth.norm = smooth_dt/apply(smooth_dt,1,max)
  #message('plot norm smooth data of knn method')
  #plotIdentities(ga.knn.smooth.norm, ga.obj, prefix=paste0(prefix,".cell_identities_norm.wknn.png"), u.lim=0.9999)
  
  # train multinomial logistic regression cell type classifier on genome-wide accessibility
  ga.pred.identities <- trainIdentity(ga.obj, ga.identities, perc=0.1, threads=10, smoothdata=F)
  ga.pred.smooth <- smoothIdentity(ga.obj, ga.pred.identities, k=15, step=3, npcs=30)
  #write.table(ga.pred.smooth, file=paste0(output_dir,'/',prefix,".predicted_celltypes.smooth.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  plotIdentities(ga.pred.smooth, ga.obj, prefix=paste0(prefix,".cell_identities.predicted.png"), u.lim=0.999)
  
  ##updating 080321 add a version to doe the normalization
  #smooth_dt <- read.table(paste0(output_dir,'/',prefix,".predicted_celltypes.smooth.txt"))
  #ga.pred.smooth.norm = smooth_dt/apply(smooth_dt,1,max)
  #message('plot norm smooth data of pred method')
  #plotIdentities(ga.pred.smooth.norm, ga.obj, prefix=paste0(prefix,".cell_identities_norm.predicted.png"), u.lim=0.9999)
  
  ##directly perform the annot other than make the temp files
  cell.data <- updateMeta(ga.obj, 
                          ga.smooth.identities,
                          ga.knn.smooth, 
                          ga.pred.smooth, 
                          cutoff=c(0,0,0,0))
  
}

##we need to assign identities
##do not need the enrich ones
#ga.identities <- read.table(paste0(output_dir,'/',prefix,".enriched_celltypes.txt"))

#ga.smooth.identities <- read.table(paste0(output_dir,'/',prefix,'.smoothed_celltypes.txt'))
#ga.knn.smooth <- read.table(paste0(output_dir,'/',prefix,'.knn_celltypes.smooth.txt'))
#ga.pred.smooth <- read.table(paste0(output_dir,'/',prefix,'.predicted_celltypes.smooth.txt'))


# assign identities #############################
#cell.data <- updateMeta(ga.obj, 
#                        ga.smooth.identities,
#                        ga.knn.smooth, 
#                        ga.pred.smooth, 
#                        cutoff=c(0,0,0,0))


# return results ----------------------------------------------------------------------------------
write.table(cell.data, file=paste0(output_dir,'/',prefix,".meta.txt"), quote=F, row.names=T, col.names=T, sep="\t")


##updating 112125
##we do not create the object

##save the annotate to the meta of the final object
#final_obj <- append(input_soc_object, list(
#  meta_annot = cell.data
#))


#saveRDS(final_obj,file=paste0(output_dir,'/',input_prefix,'.atac.soc.rds'))















