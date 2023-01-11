# Author: "Isaac Bishara"

# This script trains SingleCellNet classifier on high quality (HQ) cells using marker genes identified using 001_marker_genes.r
# The model is then used to predict lineage & cell type annotations on the low quality (LQ) cells
# parallelization require base R to run efficiently. Radian or Jupyter are not optimal and may lead to memory error

# packages
library(Seurat)
library(qs)
library(singleCellNet)
library(data.table)
library(stringr)
library(parallel)
library(pROC)

setwd('/Users/ibishara/Desktop/FELINE_C1/')
numCores <- detectCores()-2
numCores 

# load data 
# dataset <- c( 'feline', 'combes', 'pbmc' )
dataset <- 'feline' # select dataset to be analyzed 

if (dataset == 'feline'){
    #########High quality FELINE C1 dataset###############
    seu_HQ <- qread(file = "seu_HQ_no_id.qs", nthreads = numCores)
    seu_HQ <- subset(seu_HQ, subset = Celltype != "Normal epithelial cells")   ## Removes normal epithelial cells. Genes unique to normal epi cells are removed from analysis downstream
    seu_HQ <- subset(seu_HQ, subset = nCount_RNA < 15000 ) # filter out cells with > 15k UMI 
    meta <- seu_HQ@meta.data
    counts <- GetAssayData(seu_HQ, assay = "RNA")
    meta["Cell"] <- rownames(meta)
    lineage.markers <- read.table(paste(dataset,'Annotation.lineage.markers.txt', sep='/'), sep = '\t' )
    celltype.markers <- read.table(paste(dataset, 'Annotation.celltype.markers.txt', sep='/'), sep = '\t' )

}else if(dataset == 'combes'){
    #########Combes dataset###############
    combes <- qread('/Users/ibishara/Desktop/FELINE_C1/Additional_datasets/Combes_whole_blood/sobj_SCG1-2-3.qs', nthreads = numCores )
    combes <- subset(combes, subset = Celltype != 'RBC' ) # removing RBCs from downstream analysis
    meta <- combes@meta.data
    counts <- GetAssayData(combes, assay = "RNA")
    meta["Cell"] <- rownames(meta)
    celltype.markers <- read.table(paste(dataset, 'Annotation.celltype.markers.txt', sep='/'), sep = '\t' )

}else if(dataset == 'pbmc'){
    #########pbmc3k dataset###############
    # After running pbmc3k_tutorial.rmd #
    pbmc <- readRDS('/Users/ibishara/Desktop/FELINE_C1/Additional_datasets/PBMC3k/pbmc3k_final.rds')
    meta <- pbmc@meta.data
    counts <- GetAssayData(pbmc, assay = "RNA")
    meta["Cell"] <- rownames(meta)
    celltype.markers <- read.table(paste(dataset, 'Annotation.celltype.markers.txt', sep='/'), sep = '\t' )

}



# Functions 
# Reduce UMI by multiplying gene counts by a fraction. Uses a poisson distribution where the tolerance level is determined by selected theshold
# to randomize the fractions. This yields counts with mean around threshold
# x = vector/column/cell in filtered dataframe (test set)
# y = threshold
transform.pois <- function(x, y){
    x <- round(x * rpois(n=length(x), lambda=y/sum(x)))
    return(x)
    }

# function converts each column/cell to  read counts to binary, remove n genes above threshold then export a binary matrix
# x = vector/column/cell in filtered dataframe
# y = threshold
bin = function(x, y){
    x[x > 0] <- 1  # convert UMI to binary
    total = sum(x) # number of expressed genes 
    if(total > y){
        pre.index <- which(x == 1) # index of expressed genes 
        x[sample(pre.index , total - y)] <- 0 # randomly convert a number of genes over threshold from 1 -> 0
         }
    return(x)
}

# function works on each column/cell to convert read counts to binary, remove n genes above threshold then convert back to non-binary counts
# x = vector/column/cell in filtered dataframe
# y = threshold
nonbin = function(x, y){
    orig <- x # maitain count matrix 
    x[x > 0] <- 1  # convert UMI to binary
    total = sum(x) # number of expressed genes 
    if(total > y){ 
        pre.index <- which(x == 1) # index of expressed genes 
        x[sample(pre.index , total - y)] <- 0 # random convert a number of genes over threshold from 1 -> 0
    }
    x <- ifelse(x == 0, 0, orig)   # convert expressed genes back to their counts
    return(x)
}


# This function trains a classifier based off method, then loop over different thresholds to produce AUC values 
# Arguments: 
# level <- 'Celltype' , 'Lineage'
# method <- 'binary', 'non-binary', 'poisson'
SCN_run <- function (level, method) {

    # Create export directories 
    experiment <- 'SCN'
    output.dir <- paste(dataset, experiment, method, level, sep='/')
    dir.create(output.dir, recursive = TRUE)
    # sub.dir.perf <- paste(output.dir, '/model_performance', sep='')
    # dir.create(sub.dir.perf)
    sub.dir.down <- paste(output.dir, '/downsample', sep='')
    dir.create(sub.dir.down)

    # Select marker genes identified by Seurat::FindAllMarkers() for level determination 
    if ( level == 'Lineage'){ 
        common.genes <- intersect(rownames(counts), lineage.markers$gene)
    } else if ( level == 'Celltype') {
        common.genes <- intersect(rownames(counts), celltype.markers$gene)
    }
    
    set.seed(100)
    # Split 50 / 50 
    stList = splitCommon(sampTab = meta, ncells = ncells, dLevel = level) # At certain thresholds, there's not enough remaining cells for training 
    stSub = stList[[1]]
    stTrain = stSub[sample(nrow(stSub), round(nrow(stSub)/2)), ] # splits training set
    stTest = stSub[! rownames(stSub) %in% rownames(stTrain) ,]# splits test set (make sure not in training)


    expTrain.full = counts[, rownames(stTrain)] # training count matrix - all genes 
    expTrain <- expTrain.full[common.genes, ] # filtered training count matrix

    expTest = counts[ , rownames(stTest)] # test count matrix
    control.test <- expTest # untransformed UMI as a test control\


    if (method == 'binary'){
        expTrain[expTrain > 0] <- 1 # transform training counts to binary
       # expTest[expTest > 0] <- 1 # unnesessary since counts are transformed to binary via bin function 
    }

    # model training
    cat('Training model ..', "\n" )
    level_info <- scn_train(stTrain = stTrain, expTrain = expTrain, 
                    nTopGenes = nGenes, nRand = 50, nTrees = 1000, nTopGenePairs = nGenes*2, 
                    dLevel = level, colName_samp = "Cell")
    # Save model 
    qsave(level_info, file= paste(output.dir, '/Trained_model_for_', level, '_', method,'.qs', sep='' ), nthreads= numCores)

    genes <- rownames(expTest)

    # create empty summary and distribution dataframes 
    summ.out <- data.frame()
    dist.UMI <- data.frame() # for UMI distribution plots 
    dist.genes <- data.frame() # for gene distribution plots 
    counts <- expTest  # variable 'counts' now contains test count matrix - all genes

    # Transform count matrix by loop over thresholds 
    for (i in threshold_list){
        threshold <- format(i/1000, nsmall=1)
        cat(paste('Processing threshold', threshold), "\n")
       if (method == 'poisson'){ 
            table_type <- 'UMI'
            transformed <- apply(X = counts, MARGIN = 2, FUN = transform.pois, y = i) # run 2nd function to reduce the number of count per cell above threshold. iterates over columns (cells)
            total.UMI <- colSums(transformed)
            total.genes <- apply(transformed, MARGIN = 2, function(x) sum(x > 0))  # convert UMI to binary to pull n genes 
            rownames(transformed) <- genes
            total <- total.UMI # 

  
        } else if (method == 'binary'){
            table_type <- 'genes'
            # Transform genes tables 
            transformed <- apply(X = counts, MARGIN = 2, FUN = bin, y = i) # binary output
            total.genes <- colSums(transformed)
            rownames(transformed) <- genes
            transformed.nonbin <- ifelse(transformed == 0, 0, counts)   # convert expressed genes back to their counts
            total.UMI <- colSums(transformed.nonbin)
            total.genes <- colSums(transformed)
            total <- total.genes

        } else if (method == 'non-binary'){ 
            table_type <- 'genes'
            # Transform genes tables 
            transformed <- apply(X = counts, MARGIN = 2, FUN = nonbin, y = i) # non-binary output
            total.UMI <- colSums(transformed)
            total.genes <- apply(transformed, MARGIN = 2, function(x) sum(x > 0))  # convert UMI to binary to pull n genes          
            rownames(transformed) <- genes
            total <- total.genes
        }
        
        # export transformed data
        path <- paste(sub.dir.down, '/', table_type, '_down_', threshold, '_', level, sep='')
        fwrite(transformed, paste(path, '.txt', sep=''), sep='\t', nThread = numCores, row.names = TRUE)

        # export hist and stats 
        sdat <- summary(total)   
        summStr <- paste(names(sdat), format(sdat, digits = 2), collapse = "; ")
        pdf(paste(path,'.pdf', sep=''), onefile =FALSE)
            plot(hist(total), xlab = paste('n', table_type, '/cell', sep='') , main = paste('threshold', i, table_type), sub = summStr, col="#1e72d2") 
        dev.off()


        expTest <- as.data.frame(transformed[common.genes, ]) # filter test set for common genes 
        # SCN prediction
        levelRes_val_all = scn_predict(cnProc=level_info[['cnProc']], expDat = expTest, nrand = 0)  # Removed rand # number of training and validation cells must be equal. genes in model must be in validation set. 

        #  model assessment (pROC package)
        ## Remove Rand 
        levelRes_val_all <- as.data.frame(levelRes_val_all)
        levelRes_val_all <- levelRes_val_all[!rownames(levelRes_val_all) %in% 'rand' ,] # remove 'rand' category 

        # generate labels based off highest probabilities (excluding Random)
        levelRes_val_labels <- unlist( apply(levelRes_val_all, MARGIN = 2, function(x) { x <- names(which(x == max(x))[1]) })  ) # if two labels assigned equal probabilities, 2 prediction per cell generated, to avoid issue, it choose the top prediction
        levelRes_val_labels <- levelRes_val_labels[order(names(levelRes_val_labels))] # order predicted labels alphabetically by cell |         
       
        stTest <- stTest[order(rownames(stTest)),  ]    # order true labels alphabetically by cell
        test <- stTest[, level]
        AUC.pROC <- multiclass.roc(as.numeric(factor(test)), as.numeric(factor(levelRes_val_labels)))$auc[1] #  AUC value from pROC package

        cat(paste('SCN-AUC =', AUC.SCN), "\n") # diagnostic
        cat(paste('pROC-AUC =', AUC.pROC), "\n") # diagnostic

        avg.UMI <- mean(total.UMI)
        avg.genes <- mean(total.genes)

        summ <- c(level, table_type, threshold, method, AUC.SCN, AUC.pROC, ncells, nGenes, round(avg.UMI), round(avg.genes)) 
        summ.out <- rbind(summ.out, summ)

        names(total.UMI) <- NULL
        total.UMI <- c(threshold, total.UMI)
        dist.UMI <- rbind(dist.UMI, total.UMI)

        names(total.genes) <- NULL
        total.genes <- c(threshold, total.genes)
        dist.genes <- rbind(dist.genes, total.genes)
    } # end of loop 


    ## add AUC for untransformed control ##
    threshold <- 'untransformed'
    if (method == 'poisson'){ table_type <- 'UMI' }
    total.UMI <- colSums(control.test)
    total.genes <- apply(control.test, MARGIN = 2, function(x) sum(x > 0))  
    avg.UMI <- mean(total.UMI)
    avg.genes <- mean(total.genes)
    rownames(control.test) <- genes
    total <- total.UMI 
    # SCN prediction
    control.test <- control.test[common.genes, ] # filter for common genes 
    levelRes_val_all = scn_predict(cnProc=level_info[['cnProc']], expDat = control.test, nrand = 0) 
    # SCN model assessment | remove for deployment
    tm_heldoutassessment = assess_comm(ct_scores = levelRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "Cell", classTrain = level, classQuery = level, nRand = 0)
    AUC.SCN <- tm_heldoutassessment$AUPRC_w # get AUC value

    #  model assessment for control (pROC package)
    levelRes_val_all <- as.data.frame(levelRes_val_all)
    levelRes_val_all <- levelRes_val_all[!rownames(levelRes_val_all) %in% 'rand' ,] # remove 'rand' category 
    levelRes_val_labels <- unlist(apply(levelRes_val_all, MARGIN = 2, function(x) { x <- names(which(x == max(x)))[1] })) # generate labels based off highest probabilities (excluding Random)
    levelRes_val_labels <- levelRes_val_labels[order(names(levelRes_val_labels))] # order predicted labels alphabetically by cell | 
    
    stTest <- stTest[order(rownames(stTest)),  ]   # order true labels alphabetically by cell
    test <- stTest[, level]
    AUC.pROC <- multiclass.roc(as.numeric(factor(test)), as.numeric(factor(levelRes_val_labels)))$auc[1]

    summ <- c(level, table_type, threshold, method, AUC.SCN, AUC.pROC, ncells, nGenes, round(avg.UMI), round(avg.genes)) 
    summ.out <- rbind(summ.out, summ)

    names(total.UMI) <- NULL
    total.UMI <- c(threshold, total.UMI)
    dist.UMI <- rbind(dist.UMI, total.UMI)

    names(total.genes) <- NULL
    total.genes <- c(threshold, total.genes)
    dist.genes <- rbind(dist.genes, total.genes)
    #######


    cat('Generating summary table', "\n")
    colnames(summ.out) <- c( 'level', 'source', 'threshold','method', 'AUC_SCN', 'AUC_pROC', 'VnCells', 'nTopGenes', 'Avg.Reads', 'Avg.Genes')
    write.table(summ.out, paste0(dataset , '/', paste(dataset, experiment, 'performance_summary', method,  level, sep = '_'), '.txt'), col.names = TRUE, sep = '\t') 
    # export the UMI and genes ditribution at each threshold 
    cat('Generating distribution tables', "\n")
    write.table(dist.UMI, paste0(dataset, '/', paste(dataset, 'UMI_distribution', method, sep = '_'), '.txt' ), col.names = TRUE, sep = '\t') # Level is irrelevant since outputs are identical per dataset. 
    write.table(dist.genes, paste0(dataset , '/', paste(dataset, 'gene_distribution', method, sep = '_'), '.txt' ), col.names = TRUE, sep = '\t') 
}

# parameters 
nGenes <- 100 # nTopGenes for model training 
ncells <- 400 # nCells/level for training & testing dataset
threshold_list <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000) # thresholds to be tested

# Note: Lineage annotations are only available for the FELINE dataset
SCN_run('Lineage', 'poisson')
SCN_run('Celltype', 'poisson')

SCN_run('Lineage', 'non-binary')
SCN_run('Celltype', 'non-binary')

SCN_run('Lineage', 'binary')
SCN_run('Celltype', 'binary')
