# This script trains a model on HQ cells using lineage specific genes identified via FindAllMarkers 
# The model is then used to predict lineage annotations on LQ dataset
# parallelization require base R to run efficiently. Radian or Jupyter are not optimal and may flood memory
# splitCommon raises the following error when a cell type has < 3 cells: "Error in sample.int(length(x), size, replace, prob) : invalid 'size' argument"

# packages
library(Seurat)
library(qs)
library(singleCellNet)
library(data.table)
library(stringr)
library(parallel)
library(pROC)

setwd('/Users/ibishara/Desktop/FELINE_C1/')
numCores <- detectCores()
numCores 

# data
lineage.markers <- read.table('Annotation.lineage.markers.txt', sep = '\t' )
celltype.markers <- read.table('Annotation.celltype.markers.txt', sep = '\t' )

# High quality FELINE C1 data
seu_HQ <- qread(file = "seu_HQ.qs", nthreads = numCores)
seu_HQ <- subset(x = seu_HQ, subset = Celltype != "Normal epithelial cells")   ## Removes normal epithelial cells. Genes unique to normal epi cells are removed from analysis downstream
meta <- seu_HQ@meta.data
seu.HQ.counts <- GetAssayData(seu_HQ, assay = "RNA")

# Functions 
# This function is to reduce the number of count per cell above threshold
# function iterates over columns (cells). caluclulated total counts per cell. multiply gene counts by a fraction to yield total counts per cell above threshold
# x = column in filtered dataframe
# y = max total counts per cell threshold
fl00r <- function(x, y){
    total <- sum(x) # total counts for current cell 
    x <- floor((y/total)*x) # mutliply all gene counts by a fraction. genes with 0 counts are floored back to 0. other gene counts are floored to nearest integer
    return(x)
}

# function works on each column/cell to convert read counts to binary, remove n genes above threshold then export a binary matrix
bin = function(x, y){
    x[x > 0] <- 1  # convert reads to binary
    total = sum(x) # number of expressed genes 
    if(total > y){ 
        pre.index <- which(x == 1) # index of expressed genes 
        x[sample(pre.index , total - y)] <- 0 # random convert a number of genes over threshold from 1 -> 0
         }
    return(x)
}

# function works on each column/cell to convert read counts to binary, remove n genes above threshold then convert back to non-binary counts
nonbin = function(x, y){
    orig <- x # maitain count matrix 
    x[x > 0] <- 1  # convert reads to binary
    total = sum(x) # number of expressed genes 
    if(total > y){ 
        pre.index <- which(x == 1) # index of expressed genes 
        x[sample(pre.index , total - y)] <- 0 # random convert a number of genes over threshold from 1 -> 0
    }
    x <- ifelse(x == 0, 0, orig)   # convert expressed genes back to their counts
    return(x)
}



# This function run trains a classifier based off method, then loop over different thresholds to produce AUC values 
# Arguments: 
# class <- 'Celltype' , 'Lineage'
# method <- 'floor', 'binary', 'non-binary'

RF_run <- function (class, method) {

    # class <- 'Lineage' # diagnostic 
    # method <- 'floor' # diagnostic
    # i <- 2000 # diagnostic 
    # Create export directories 
    output.dir <- paste('SCN', method, class, sep='/')
    dir.create(output.dir, recursive = TRUE)
    sub.dir.perf <- paste(output.dir, '/model_performance', sep='')
    dir.create(sub.dir.perf)
    sub.dir.down <- paste(output.dir, '/downsample', sep='')
    dir.create(sub.dir.down)

    # Select for top genes for class determination 
    if ( class == 'Lineage'){ 
        common.genes <- intersect(rownames(seu.HQ.counts), lineage.markers$gene)
    } else if ( class == 'Celltype') {
        common.genes <- intersect(rownames(seu.HQ.counts), celltype.markers$gene)
        }

    
    set.seed(100)
    # split training set
    stList = splitCommon(sampTab = meta, ncells = Tncells, dLevel = class) # At certain thresholds, there's not enough remaining cells for training 
    stTrain = stList[[1]]
    expTrain = seu.HQ.counts[, rownames(stTrain)]
    expTrain <- expTrain[common.genes, ]

    if (method == 'binary'){
        expTrain[expTrain > 0] <- 1 # transform training counts to binary
    }

    # model training
    system.time(class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, 
                            nTopGenes = nGenes, nRand = 50, nTrees = 1000, nTopGenePairs = nGenes*2, 
                            dLevel = class, colName_samp = "Cell"))
    # Save model 
    qsave(class_info, file= paste(output.dir, '/Trained_model_for_', class, '_', method,'.qs', sep='' ), nthreads= numCores)

    # Cells available for validation set
    counts <- as.data.frame(seu.HQ.counts[, rownames(stList[[2]])])   # counts minus training set 
    counts.binary <- as.data.frame(mclapply(counts, FUN = function(x) {ifelse( x > 0, 1, 0)}, mc.cores= numCores))  # convert reads to binary
    pdf('SCN/Untransformed_reads_genes_plot_val.pdf')
    plot(log10(colSums(counts)), log10(colSums(counts.binary)), pch = 20, cex = 0.2,  main = 'Untransformed matrix', xlab = "number of reads", ylab = "number of unique genes" )
    dev.off()
    genes <- rownames(counts)

    # create empty summary and distribution dataframes 
    summ.out <- data.frame()
    dist.reads <- data.frame()
    dist.genes <- data.frame()

    # Transform count matrix by loop over thresholds 
    for (i in threshold_list){
        threshold <- format(i/1000, nsmall=1)
        print(noquote(paste('Processing threshold', threshold)))
        if (method == 'floor'){ 
            table_type <- 'reads'

            # filters for cells with total reads above threshold (to be transformed)
            x <- counts[, which(colSums(counts) > i)]
            transformed <- as.data.frame(mclapply(X = x, FUN = fl00r, y = i, mc.cores= numCores )) # run 2nd function to reduce the number of count per cell above threshold. iterates over columns (cells)
            total.reads <- colSums(transformed)
            t.binary <- transformed
            t.binary <- as.data.frame(mclapply(t.binary, FUN = function(x) {ifelse( x > 0, 1, 0)}, mc.cores= numCores))  # convert reads to binary
            total.genes <- colSums(t.binary)
            rownames(transformed) <- genes
            total <- total.reads

  
        } else if (method == 'binary'){
            table_type <- 'genes'
            # Transform genes tables 
            transformed <- as.data.frame(mclapply(counts, FUN = bin, i, mc.cores= numCores)) # binary output
            total.genes <- colSums(transformed)
            rownames(transformed) <- genes
            total.reads <- rep(0, ncol(transformed)) 
            total <- total.genes

        } else if (method == 'non-binary'){ 
            table_type <- 'genes'
            # Transform genes tables 
            transformed <- as.data.frame(mclapply(counts, FUN = nonbin, i, mc.cores= numCores)) # normal output
            total.reads <- colSums(transformed)
            t.binary <- transformed
            t.binary <- as.data.frame(mclapply(t.binary, FUN = function(x) {ifelse( x > 0, 1, 0)}, mc.cores= numCores))# convert reads to binary
            total.genes <- colSums(t.binary)
            rownames(transformed) <- genes
            total <- total.genes
        }
    
    # export transformed data
    path <- paste(sub.dir.down, '/', table_type, '_down_', threshold, '_', class, sep='')
    fwrite(transformed, paste(path, '.txt', sep=''), sep='\t', nThread = numCores, row.names = TRUE)

    # export hist and stats 
    sdat <- summary(total)   
    summStr <- paste(names(sdat), format(sdat, digits = 2), collapse = "; ")
    pdf(paste(path,'.pdf', sep=''), onefile =FALSE)
    plot(hist(total), xlab = paste('n', table_type, '/cell', sep='') , main = paste('cut-off', i, table_type), sub = summStr, col="#1e72d2") 
    dev.off()

    # Validation on transformed data
    stTestList = splitCommon(sampTab = stList[[2]][ stList[[2]]$Cell %in% colnames(transformed),], ncells = Vncells, dLevel = class)  # filter for transformed cells from validation set
    stTest = stTestList[[1]]
    expTest = transformed[ , rownames(stTest)]
    expTest <- expTest[common.genes, ] 

    # SCN prediction
    classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat = expTest, nrand = 50)  # number of training and validation cells must be equal. genes in model must be in validation set. | issue: some dropped genes lead to error
    # SCN model assessment 
    tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, 
                                    dLevelSID = "Cell", classTrain = class, classQuery = class, nRand = 50)
    AUC.SCN <- tm_heldoutassessment$AUPRC_w # get AUC value

    # Alternative model assessment (pROC package)
    classRes_val_labels <- apply(classRes_val_all[,1:(ncol(classRes_val_all)-50) ], MARGIN = 2, function(x) { x <- names(which(x == max(x))) }) # generate labels based off highest probabilities (excluding Random)
    
    classRes_val_labels <- classRes_val_labels[order(names(classRes_val_labels))] # order predicted labels alphabetically by cell 
    stTest <- stTest[ order(rownames(stTest)), ]    # order true labels alphabetically by cell
    if ( is.list(classRes_val_labels) ){ # at some thresholds, classRes_val_labels is a list instead of a vector
        classRes_val_labels <- as.data.frame(classRes_val_labels)
        AUC.pROC <- multiclass.roc(stTest[, class], factor(classRes_val_labels[1,], ordered = TRUE))$auc[1]
    } else {
        AUC.pROC <- multiclass.roc(stTest[, class], factor(classRes_val_labels, ordered = TRUE))$auc[1]}

    # changed from : 
    #     AUC.pROC <- multiclass.roc(stTest[, class], factor(classRes_val_labels[1,], ordered = TRUE))$auc[1]
    #     }
    # AUC.pROC <- multiclass.roc(stTest[, class], factor(classRes_val_labels, ordered = TRUE))$auc[1]

    ## Plot performance metrics 
    print(noquote('Generating plots'))
    # plots 
    pdf(paste(sub.dir.perf, '/', method, '_', class, '_', threshold, '.pdf', sep = ''))
            hist(total.reads, main = paste(table_type, '_', method, '_', class, '_', threshold, sep = ''))         
            hist(total.genes, main = paste(table_type, '_', method, '_', class, '_', threshold, sep = ''))
            plot(plot_PRs(tm_heldoutassessment))
            plot(plot_attr(classRes = classRes_val_all, sampTab=stTest, nrand=50, dLevel=class, sid="Cell"))
            plot(plot_metrics(tm_heldoutassessment))
            if (method == 'floor' | method == 'non-binary') {
                coeff <- round(cor(total.reads, total.genes), 2)
                plot(log10(total.reads), log10(total.genes), pch = 20, cex = 0.2,  
                main = paste('Transformed using', method, 'at threshold', i, table_type), sub = paste("Pearson's coefficient =", coeff), 
                xlab = "log10 number of reads", ylab = "log10 number of unique genes" )}
    dev.off()

    avg.reads <- mean(total.reads)
    avg.genes <- mean(total.genes)

    summ <- c(class, table_type, threshold, method, AUC.SCN, AUC.pROC, Vncells, nGenes, round(avg.reads), round(avg.genes)) 
    summ.out <- rbind(summ.out, summ)

    names(total.reads) <- NULL
    total.reads <- c(threshold, total.reads)
    dist.reads <- rbind(dist.reads, total.reads)

    names(total.genes) <- NULL
    total.genes <- c(threshold, total.genes)
    dist.genes <- rbind(dist.genes, total.genes)
    }

    print(noquote('Generating summary table'))
    colnames(summ.out) <- c( 'class', 'source', 'threshold','method', 'AUC_SCN', 'AUC_pROC', 'VnCells', 'nTopGenes', 'Avg.Reads', 'Avg.Genes')
    write.table(summ.out, paste(getwd() , '/SCN_Performance_summary_', method, '_', class, '.txt' , sep=''), col.names = TRUE, sep = '\t') 
    # export the reads and genes ditribution at each threshold 
    print(noquote('Generating distribution tables'))
    write.table(dist.reads, paste(getwd() , '/SCN_reads_distribution_', method, '.txt' , sep=''), col.names = TRUE, sep = '\t') # Lineage and celltype outputs are identical. They're re-exported for validation only. 
    write.table(dist.genes, paste(getwd() , '/SCN_genes_distribution_', method, '.txt' , sep=''), col.names = TRUE, sep = '\t') 
}

# parameters 
nGenes <- 25 # nTopGenes for model training 
Tncells <- 200 # nCells/class for training dataset
Vncells <- 400 # nCells/class for testing dataset
threshold_list <- c(200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 4000, 5000) # thresholds to be tested


RF_run('Lineage', 'floor')
RF_run('Celltype', 'floor')

RF_run('Lineage', 'non-binary')
RF_run('Celltype', 'non-binary')

RF_run('Lineage', 'binary')
RF_run('Celltype', 'binary')


