# Author: "Isaac Bishara"

# This script trains SingleR classifier on high quality (HQ) cells using all available genes
# The model is then used to predict lineage & cell type annotations on the low quality (LQ) cells
# The output is AUROCC as a performance metric for different transformations e.g., 'binary', 'poisson'
# parallelization require base R to run efficiently. Radian or Jupyter are not optimal and may lead to memory error

# packages
library(Seurat)
library(qs)
library(data.table)
library(parallel)
library(SingleR)
library(singleCellNet)
library(pROC)
library(BiocParallel)

setwd('/Users/ibishara/Desktop/FELINE_C1/')
numCores <- detectCores()-2
numCores 

# dataset <- c( 'feline', 'combes', 'pbmc' )
dataset <- 'pbmc' # select dataset to be analyzed 

if (dataset == 'feline'){
    #########High quality FELINE C1 dataset###############
    seu_HQ <- qread(file = "seu_HQ_no_id.qs", nthreads = numCores)
    seu_HQ <- subset(seu_HQ, subset = Celltype != "Normal epithelial cells")   ## Removes normal epithelial cells. Genes unique to normal epi cells are removed from analysis downstream
    seu_HQ <- subset(seu_HQ, subset = nCount_RNA < 15000 ) # filter out cells with > 15k UMI 
    meta <- seu_HQ@meta.data
    counts <- GetAssayData(seu_HQ, assay = "RNA")

}else if(dataset == 'combes'){
    #########Combes dataset###############
    combes <- qread('/Users/ibishara/Desktop/FELINE_C1/Additional_datasets/Combes_whole_blood/sobj_SCG1-2-3.qs', nthreads = numCores )
    combes <- subset(combes, subset = Celltype != 'RBC' ) # removing RBCs from downstream analysis
    meta <- combes@meta.data
    counts <- GetAssayData(combes, assay = "RNA")

}else if(dataset == 'pbmc'){
    #########pbmc3k dataset###############
    # After running pbmc3k_tutorial.rmd #
    pbmc <- readRDS('/Users/ibishara/Desktop/FELINE_C1/Additional_datasets/PBMC3k/pbmc3k_final.rds')
    meta <- pbmc@meta.data
    counts <- GetAssayData(pbmc, assay = "RNA")
}


# Functions 
# Reduce UMI by multiplying gene counts by a fraction. Uses a poisson distribution where the tolerance level is determined by selected theshold
# to randomize the fractions. This yields counts with mean around threshold
# x = vector/column/cell in filtered dataframe
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
SR_run <- function (level, method) {

    ############### Diagnostic only ############
    # level <- 'Celltype' 
    # method <- 'poisson' 
    # i <- 4000 
    ############################################
    
    # Create export directories 
    experiment <- 'SR'
    output.dir <- paste(dataset, experiment, method, level, sep='/')
    dir.create(output.dir, recursive = TRUE)
    # sub.dir.perf <- paste(output.dir, '/model_performance', sep='')
    # dir.create(sub.dir.perf)
    sub.dir.down <- paste(output.dir, '/downsample', sep='')
    dir.create(sub.dir.down)

    
    set.seed(100)
    
    # Split 50 / 50 
    stList = splitCommon(sampTab = meta, ncells = ncells, dLevel = level) # At certain thresholds, there's not enough remaining cells for training 
    stSub = stList[[1]]
    stTrain = stSub[sample(nrow(stSub), round(nrow(stSub)/2)), ]
    stTest = stSub[! rownames(stSub) %in% rownames(stTrain) ,]

    expTrain = counts[, rownames(stTrain)] # train on all genes 

    expTest = counts[ , rownames(stTest)]
    control.test <- expTest # untransformed UMIs as a test control\

    if (method == 'binary'){
        expTrain[expTrain > 0] <- 1 # transform training counts to binary
       # expTest[expTest > 0] <- 1 # unnesessary since counts are transformed to binary via bin function 
    }

    # model training
    if ( level == 'Lineage'){ 
        level_info <- trainSingleR(expTrain, stTrain$Lineage)
    } else if ( level == 'Celltype') {
        level_info <- trainSingleR(expTrain, stTrain$Celltype)
    }

    # Save model 
    qsave(level_info, file= paste(output.dir, '/Trained_model_for_', level, '_', method,'.qs', sep='' ), nthreads= numCores)

    # # pre-transformation Diagnostics
    # tot_counts_train <- unlist(mclapply(as.data.frame(expTrain), function (x) sum(x), mc.cores= numCores ))
    # tot_genes_train <- unlist(mclapply(as.data.frame(expTrain), function (x) sum(x > 0), mc.cores= numCores ))
    # tot_counts_test <- unlist(mclapply(as.data.frame(expTest), function (x) sum(x), mc.cores= numCores ))
    # tot_genes_test <- unlist(mclapply(as.data.frame(expTest), function (x) sum(x > 0), mc.cores= numCores ))
    # pdf(paste0(experiment,'/Train_Test_sets_corr_plot.pdf'))
    #     plot(tot_counts_train, tot_genes_train, main= paste('Training set, n =', length(tot_counts_train))); abline(lm(tot_genes_train ~ tot_counts_train))
    #     plot(tot_counts_test, tot_genes_test, main= paste('Testing set, n =', length(tot_counts_test))); abline(lm(tot_genes_test ~ tot_counts_test))
    # dev.off()

    genes <- rownames(expTest)

    # create empty summary and distribution dataframes 
    summ.out <- data.frame()
    dist.UMI <- data.frame()
    dist.genes <- data.frame()
    counts <- expTest

    # Transform count matrix by loop over thresholds 
    for (i in threshold_list){
        threshold <- format(i/1000, nsmall=1)
        cat(paste('Processing threshold', threshold), "\n")
       if (method == 'poisson'){ 
            table_type <- 'UMI'
            transformed <- apply(X = counts, MARGIN = 2, FUN = transform.pois, y = i) # run 2nd function to reduce the number of count per cell above threshold. iterates over columns (cells)
            total.UMI <- colSums(transformed)
            total.genes <- apply(transformed, MARGIN = 2, function(x) sum(x > 0))  # convert UMI count to binary to pull n genes 
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
            total.genes <- apply(transformed, MARGIN = 2, function(x) sum(x > 0))  # convert UMIs to binary to pull n genes          
            rownames(transformed) <- genes
            total <- total.genes
        }
    
    # export transformed data
    path <- paste(sub.dir.down, '/', table_type, '_down_', threshold, '_', level, sep='')
    fwrite(transformed, paste(path, '.txt', sep=''), sep='\t', nThread = numCores, row.names = TRUE)

    # ## Plot performance metrics - Diagnostic only
    # cat('Generating plots', "\n")
    # # plots 
    # pdf(paste(sub.dir.perf, '/', method, '_', level, '_', threshold, '.pdf', sep = ''))
    #         hist(total.UMI, main = paste(table_type, '_', method, '_', level, '_', threshold, sep = ''))         
    #         hist(total.genes, main = paste(table_type, '_', method, '_', level, '_', threshold, sep = ''))
    #         if (max(total.UMI) != 0) {  
    #             coeff <- round(cor(total.UMI, total.genes), 2)
    #             plot(log10(total.UMI), log10(total.genes), pch = 20, cex = 0.2,  
    #             main = paste('Transformed using', method, 'at threshold', i, table_type), sub = paste("Pearson's coefficient =", coeff), 
    #             xlab = "log10 number of UMI", ylab = "log10 number of unique genes" )}
    # dev.off()

    # export hist and stats 
    sdat <- summary(total)   
    summStr <- paste(names(sdat), format(sdat, digits = 2), collapse = "; ")
    pdf(paste(path,'.pdf', sep=''), onefile =FALSE)
    plot(hist(total), xlab = paste('n', table_type, '/cell', sep='') , main = paste('threshold', i, table_type), sub = summStr, col="#1e72d2") 
    dev.off()


    expTest <- as.data.frame(transformed) # run on all genes 

    # SingleR prediction
    levelRes_val_all = classifySingleR(expTest, level_info, fine.tune = FALSE, prune = FALSE, BPPARAM = MulticoreParam(numCores)) 

    # SingleR model evaluation 
    stTest <- stTest[ order(rownames(stTest)), ]    # order true labels alphabetically by cell
    levelRes_val_all <- levelRes_val_all[order(rownames(levelRes_val_all)),] # order predicted labels alphabetically by cell 
    AUC.pROC <- multiclass.roc(as.numeric(factor(stTest[, level])), as.numeric(factor(levelRes_val_all$labels)))$auc[1]

    cat(paste('pROC-AUC =', AUC.pROC), "\n")

    avg.UMI <- mean(total.UMI)
    avg.genes <- mean(total.genes)

    summ <- c(level, table_type, threshold, method, AUC.pROC, ncells, round(avg.UMI), round(avg.genes)) 
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

    # SingleR prediction
    levelRes_val_all = classifySingleR(control.test, level_info, fine.tune = FALSE, prune = FALSE, BPPARAM = MulticoreParam(numCores)) 

    # model assessment (pROC package)
    stTest <- stTest[ order(rownames(stTest)), ]    # order true labels alphabetically by cell
    levelRes_val_all <- levelRes_val_all[order(rownames(levelRes_val_all)),] # order predicted labels alphabetically by cell 
    AUC.pROC <- multiclass.roc(as.numeric(factor(stTest[, level])), as.numeric(factor(levelRes_val_all$labels)))$auc[1]
    cat(paste('pROC-AUC =', AUC.pROC), "\n")

    summ <- c(level, table_type, threshold, method, AUC.pROC, ncells, round(avg.UMI), round(avg.genes)) 
    summ.out <- rbind(summ.out, summ)

    names(total.UMI) <- NULL
    total.UMI <- c(threshold, total.UMI)
    dist.UMI <- rbind(dist.UMI, total.UMI)

    names(total.genes) <- NULL
    total.genes <- c(threshold, total.genes)
    dist.genes <- rbind(dist.genes, total.genes)
    #######


    cat('Generating summary table', "\n")
    colnames(summ.out) <- c( 'level', 'source', 'threshold','method', 'AUC_pROC', 'nCells', 'Avg.Reads', 'Avg.Genes')
    write.table(summ.out, paste0(dataset , '/', paste(dataset, experiment, 'performance_summary', method,  level, sep = '_'), '.txt'), col.names = TRUE, sep = '\t') 
    # export the UMI and genes ditribution at each threshold 
    cat('Generating distribution tables', "\n")
    write.table(dist.UMI, paste0(dataset, '/', paste(dataset, 'UMI_distribution', method, sep = '_'), '.txt' ), col.names = TRUE, sep = '\t') # Level is irrelevant since outputs are identical per dataset. 
    write.table(dist.genes, paste0(dataset , '/', paste(dataset, 'gene_distribution', method, sep = '_'), '.txt' ), col.names = TRUE, sep = '\t') 
}

# parameters 
ncells <- 400 # nCells/level for training & testing dataset
threshold_list <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000) # thresholds to be tested


# Note: Lineage annotations are only available for the FELINE dataset
SR_run('Lineage', 'poisson')
SR_run('Celltype', 'poisson') 

SR_run('Lineage', 'non-binary')
SR_run('Celltype', 'non-binary')

SR_run('Lineage', 'binary')
SR_run('Celltype', 'binary')


