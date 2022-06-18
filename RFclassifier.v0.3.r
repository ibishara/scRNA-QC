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
seu_HQ <- qread(file = "seu_HQ_no_id3.qs", nthreads = numCores)
# seu_HQ_rename <- qread(file = "seu_HQ_no_id2.qs", nthreads = numCores)

# # ###########################
# meta <- seu_HQ@meta.data
# seu.HQ.counts <- GetAssayData(seu_HQ, assay = "RNA")

# meta_rename <- seu_HQ_rename@meta.data
# seu.HQ.counts_rename <- GetAssayData(seu_HQ_rename, assay = "RNA")

# dim(meta)
# dim(meta_rename)
# all(meta$Lineage == meta_rename$Lineage)
# all(seu.HQ.counts == seu.HQ.counts_rename)
# all(rownames(meta_rename) == colnames(seu.HQ.counts_rename))
# all(rownames(meta) == colnames(seu.HQ.counts))

# # ############################


seu_HQ <- subset(x = seu_HQ, subset = Celltype != "Normal epithelial cells")   ## Removes normal epithelial cells. Genes unique to normal epi cells are removed from analysis downstream
seu_HQ <- subset(seu_HQ, subset = nCount_RNA < 15000 ) # filter out cells with > 15k reads 
meta <- seu_HQ@meta.data
seu.HQ.counts <- GetAssayData(seu_HQ, assay = "RNA")

# Functions 
# Reduce reads by multiplying gene counts by a fraction. Uses a poisson distribution where the tolerance level is determined by selected theshold
# to randomize the fractions. This yields counts with mean around threshold
# x = vector/column/cell in filtered dataframe
# y = threshold
red.reads <- function(x, y){
    x <- round(x * rpois(n=length(x), lambda=y/sum(x)))
    return(x)
    }

# function converts each column/cell to  read counts to binary, remove n genes above threshold then export a binary matrix
# x = vector/column/cell in filtered dataframe
# y = threshold
bin = function(x, y){
    x[x > 0] <- 1  # convert reads to binary
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
    x[x > 0] <- 1  # convert reads to binary
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
# class <- 'Celltype' , 'Lineage'
# method <- 'binary', 'non-binary', 'poisson'
RF_run <- function (class, method) {

    class <- 'Lineage' # diagnostic 
    method <- 'poisson' # diagnostic
    i <- 2000 # diagnostic 

    # Create export directories 
    experiment <- 'SCN'
    output.dir <- paste(experiment, method, class, sep='/')
    dir.create(output.dir, recursive = TRUE)
    sub.dir.perf <- paste(output.dir, '/model_performance', sep='')
    dir.create(sub.dir.perf)
    sub.dir.down <- paste(output.dir, '/downsample', sep='')
    dir.create(sub.dir.down)

    # Select marker genes identified by Seurat::FindAllMarkers() for class determination 
    if ( class == 'Lineage'){ 
        common.genes <- intersect(rownames(seu.HQ.counts), lineage.markers$gene)
    } else if ( class == 'Celltype') {
        common.genes <- intersect(rownames(seu.HQ.counts), celltype.markers$gene)
        }
    
    set.seed(100)
    # Split 50 / 50 
    stList = splitCommon(sampTab = meta, ncells = ncells, dLevel = class) # At certain thresholds, there's not enough remaining cells for training 
    stSub = stList[[1]]
    stTrain = stSub[sample(nrow(stSub), round(nrow(stSub)/2)), ]
    stTest = stSub[! rownames(stSub) %in% rownames(stTrain) ,]

    expTrain.full = seu.HQ.counts[, rownames(stTrain)]
    expTrain <- expTrain.full[common.genes, ]

    expTest = seu.HQ.counts[ , rownames(stTest)]
    control.test <- expTest # untransformed reads as a test control\


    if (method == 'binary'){
        expTrain[expTrain > 0] <- 1 # transform training counts to binary
       # expTest[expTest > 0] <- 1 # unnesessary since counts are transformed to binary via bin function 
    }

    # model training
    print(noquote('Training model ..'))
    class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, 
                    nTopGenes = nGenes, nRand = 50, nTrees = 1000, nTopGenePairs = nGenes*2, 
                    dLevel = class, colName_samp = "Cell")
    # Save model 
    qsave(class_info, file= paste(output.dir, '/Trained_model_for_', class, '_', method,'.qs', sep='' ), nthreads= numCores)

    # pre-transformation Diagnostics
    tot_counts_train <- unlist(mclapply(as.data.frame(expTrain), function (x) sum(x), mc.cores= numCores ))
    tot_genes_train <- unlist(mclapply(as.data.frame(expTrain), function (x) sum(x > 0), mc.cores= numCores ))
    tot_counts_test <- unlist(mclapply(as.data.frame(expTest), function (x) sum(x), mc.cores= numCores ))
    tot_genes_test <- unlist(mclapply(as.data.frame(expTest), function (x) sum(x > 0), mc.cores= numCores ))
    pdf('SCN/Train_Test_sets_corr_plot.pdf')
        plot(tot_counts_train, tot_genes_train, main= paste('Training set, n =', length(tot_counts_train))); abline(lm(tot_genes_train ~ tot_counts_train))
        plot(tot_counts_test, tot_genes_test, main= paste('Testing set, n =', length(tot_counts_test))); abline(lm(tot_genes_test ~ tot_counts_test))
    dev.off()

    genes <- rownames(expTest)

    # create empty summary and distribution dataframes 
    summ.out <- data.frame()
    dist.reads <- data.frame()
    dist.genes <- data.frame()
    counts <- expTest


    # Transform count matrix by loop over thresholds 
    for (i in threshold_list){
        threshold <- format(i/1000, nsmall=1)
        print(noquote(paste('Processing threshold', threshold)))
       if (method == 'poisson'){ 
            table_type <- 'reads'
            transformed <- apply(X = counts, MARGIN = 2, FUN = red.reads, y = i) # run 2nd function to reduce the number of count per cell above threshold. iterates over columns (cells)
            total.reads <- colSums(transformed)
            total.genes <- apply(transformed, MARGIN = 2, function(x) sum(x > 0))  # convert reads to binary to pull n genes 
            rownames(transformed) <- genes
            total <- total.reads # 

  
        } else if (method == 'binary'){
            table_type <- 'genes'
            # Transform genes tables 
            transformed <- apply(X = counts, MARGIN = 2, FUN = bin, y = i) # binary output
            total.genes <- colSums(transformed)
            rownames(transformed) <- genes
            transformed.nonbin <- ifelse(transformed == 0, 0, counts)   # convert expressed genes back to their counts
            total.reads <- colSums(transformed.nonbin)
            total.genes <- colSums(transformed)
            total <- total.genes

        } else if (method == 'non-binary'){ 
            table_type <- 'genes'
            # Transform genes tables 
            transformed <- apply(X = counts, MARGIN = 2, FUN = nonbin, y = i) # non-binary output
            total.reads <- colSums(transformed)
            total.genes <- apply(transformed, MARGIN = 2, function(x) sum(x > 0))  # convert reads to binary to pull n genes          
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
            plot(hist(total), xlab = paste('n', table_type, '/cell', sep='') , main = paste('threshold', i, table_type), sub = summStr, col="#1e72d2") 
        dev.off()
####

        expTest <- as.data.frame(transformed[common.genes, ]) # filter test set for common genes 
        # SCN prediction
        classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat = expTest, nrand = 0)  # Removed rand # number of training and validation cells must be equal. genes in model must be in validation set. | issue: some dropped genes lead to error
        # SCN model assessment | remove for deployment
        tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "Cell", classTrain = class, classQuery = class, nRand = 0)
        AUC.SCN <- tm_heldoutassessment$AUPRC_w # get AUC value

        #  model assessment (pROC package)
        ## Remove Rand 
        classRes_val_all <- as.data.frame(classRes_val_all)
        classRes_val_all <- classRes_val_all[!rownames(classRes_val_all) %in% 'rand' ,] # remove 'rand' category 

        # generate labels based off highest probabilities (excluding Random)
        classRes_val_labels <- unlist( apply(classRes_val_all, MARGIN = 2, function(x) { x <- names(which(x == max(x))[1]) })  ) # if two labels assigned equal probabilities, 2 prediction per cell generated, to avoid issue, it choose the top prediction
        classRes_val_labels <- classRes_val_labels[order(names(classRes_val_labels))] # order predicted labels alphabetically by cell |         
       
        stTest <- stTest[order(rownames(stTest)),  ]    # order true labels alphabetically by cell
        test <- stTest[, class]
        AUC.pROC <- multiclass.roc(as.numeric(factor(test)), as.numeric(factor(classRes_val_labels)))$auc[1]

        print(paste('SCN-AUC =', AUC.SCN)) # diagnostic
        print(paste('pROC-AUC =', AUC.pROC)) # diagnostic
####
        ## Plot performance metrics 
        print(noquote('Generating plots'))
        # plots 
        pdf(paste(sub.dir.perf, '/', method, '_', class, '_', threshold, '.pdf', sep = ''))
                hist(total.reads, main = paste(table_type, '_', method, '_', class, '_', threshold, sep = ''))         
                hist(total.genes, main = paste(table_type, '_', method, '_', class, '_', threshold, sep = ''))
                plot(plot_PRs(tm_heldoutassessment))
                plot(plot_metrics(tm_heldoutassessment))
                if (max(total.reads) != 0) {  
                    coeff <- round(cor(total.reads, total.genes), 2)
                    plot(log10(total.reads), log10(total.genes), pch = 20, cex = 0.2,  
                    main = paste('Transformed using', method, 'at threshold', i, table_type), sub = paste("Pearson's coefficient =", coeff), 
                    xlab = "log10 number of reads", ylab = "log10 number of unique genes" )}
        dev.off()

        avg.reads <- mean(total.reads)
        avg.genes <- mean(total.genes)

        summ <- c(class, table_type, threshold, method, AUC.SCN, AUC.pROC, ncells, nGenes, round(avg.reads), round(avg.genes)) 
        summ.out <- rbind(summ.out, summ)

        names(total.reads) <- NULL
        total.reads <- c(threshold, total.reads)
        dist.reads <- rbind(dist.reads, total.reads)

        names(total.genes) <- NULL
        total.genes <- c(threshold, total.genes)
        dist.genes <- rbind(dist.genes, total.genes)
    } # end of loop 


    ## add AUC for untransformed control ##
    threshold <- 'untransformed'
    if (method == 'poisson'){ table_type <- 'reads' }
    total.reads <- colSums(control.test)
    total.genes <- apply(control.test, MARGIN = 2, function(x) sum(x > 0))  
    avg.reads <- mean(total.reads)
    avg.genes <- mean(total.genes)
    rownames(control.test) <- genes
    total <- total.reads 
    # SCN prediction
    control.test <- control.test[common.genes, ] # filter for common genes 
    classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat = control.test, nrand = 0) 
    # SCN model assessment | remove for deployment
    tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "Cell", classTrain = class, classQuery = class, nRand = 0)
    AUC.SCN <- tm_heldoutassessment$AUPRC_w # get AUC value

    #  model assessment for control (pROC package)
    ## Remove Rand 
    classRes_val_all <- as.data.frame(classRes_val_all)
    classRes_val_all <- classRes_val_all[!rownames(classRes_val_all) %in% 'rand' ,] # remove 'rand' category 
    classRes_val_labels <- unlist(apply(classRes_val_all, MARGIN = 2, function(x) { x <- names(which(x == max(x))) })) # generate labels based off highest probabilities (excluding Random)
    classRes_val_labels <- classRes_val_labels[order(names(classRes_val_labels))] # order predicted labels alphabetically by cell | 
    # classRes_val_labels <- classRes_val_labels[rownames(stTest)] ## Almost randomly, a cell or two are added with a digit after bar code, this step is to remove these extra cells until debugged
    
    stTest <- stTest[order(rownames(stTest)),  ]   # order true labels alphabetically by cell
    test <- stTest[, class]
    AUC.pROC <- multiclass.roc(as.numeric(factor(test)), as.numeric(factor(classRes_val_labels)))$auc[1]

    summ <- c(class, table_type, threshold, method, AUC.SCN, AUC.pROC, ncells, nGenes, round(avg.reads), round(avg.genes)) 
    summ.out <- rbind(summ.out, summ)

    names(total.reads) <- NULL
    total.reads <- c(threshold, total.reads)
    dist.reads <- rbind(dist.reads, total.reads)

    names(total.genes) <- NULL
    total.genes <- c(threshold, total.genes)
    dist.genes <- rbind(dist.genes, total.genes)
    #######


    print(noquote('Generating summary table'))
    colnames(summ.out) <- c( 'class', 'source', 'threshold','method', 'AUC_SCN', 'AUC_pROC', 'VnCells', 'nTopGenes', 'Avg.Reads', 'Avg.Genes')
    write.table(summ.out, paste(getwd() , '/', experiment, '_Performance_summary_', method, '_', class, '.txt' , sep=''), col.names = TRUE, sep = '\t') 
    # export the reads and genes ditribution at each threshold 
    print(noquote('Generating distribution tables'))
    write.table(dist.reads, paste(getwd() , '/', experiment, '_reads_distribution_', method, '.txt' , sep=''), col.names = TRUE, sep = '\t') # Lineage and celltype outputs are identical. They're re-exported for validation only. 
    write.table(dist.genes, paste(getwd() , '/', experiment, '_genes_distribution_', method, '.txt' , sep=''), col.names = TRUE, sep = '\t') 
}

# parameters 
nGenes <- 100 # nTopGenes for model training 
ncells <- 400 # nCells/class for training & testing dataset
threshold_list <- c(2000, 4000) # thresholds to be tested
# threshold_list <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500) # thresholds to be tested
# threshold_list <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000) # thresholds to be tested


RF_run('Lineage', 'poisson')
RF_run('Celltype', 'poisson')

RF_run('Lineage', 'non-binary')
RF_run('Celltype', 'non-binary')

RF_run('Lineage', 'binary')
RF_run('Celltype', 'binary')


