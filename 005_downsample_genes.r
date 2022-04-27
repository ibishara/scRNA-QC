
library(data.table)
library(dplyr)
library(qs)
library(Seurat)
setwd('/Users/ibishara/Desktop/FELINE_C1/')

# data
seu_HQ <- qread(file="seu_HQ.qs", nthreads=16)
seu_HQ <- subset(x = seu_HQ, subset = Celltype != "Normal epithelial cells")   ## Removes normal epithelial cells. Genes unique to normal epi cells are removed from analysis downstream

counts <- as.data.frame(GetAssayData(seu_HQ, assay = "RNA"))
set.seed(100)
numCores <- detectCores()
numCores

# functions 
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


# export
threshold <- c(200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 4000)

for (i in threshold){

    y <- i
    b <- as.data.frame(mclapply(counts, FUN = bin, y, mc.cores= numCores)) # binary output
    rownames(b) <- rownames(counts)
    n <- as.data.frame(mclapply(counts, FUN = nonbin, y, mc.cores= numCores)) # normal output
    rownames(n) <- rownames(counts)

    bpath <- paste('downsample/genes_downsample/binary/genes_down_', format(y/1000, nsmall=1), sep='')
    npath <- paste('downsample/genes_downsample/non-binary/genes_down_', format(y/1000, nsmall=1), sep='')

    fwrite(b, paste(bpath, '.txt', sep=''), sep='\t', nThread = 16, row.names = TRUE)
    fwrite(n, paste(npath, '.txt', sep=''), sep='\t', nThread = 16, row.names = TRUE) 

    # export binary hist and stats 
    btotal <- colSums(b)
    sdat <- summary(btotal)   
    summStr <- paste(names(sdat), format(sdat, digits = 2), collapse = "; ")
    pdf(paste(bpath,'.pdf', sep=''), onefile =FALSE)
    plot(hist(btotal), xlab = 'nGenes/cell', main = paste('cut-off', y, '# genes'), sub = summStr, col="#53b478") 
    dev.off()


    # export non-binary hist and stats 
    n[n > 0] <- 1  # convert reads to binary to calculate n genes
    ntotal <- colSums(n)
    sdat <- summary(ntotal)   
    summStr <- paste(names(sdat), format(sdat, digits = 2), collapse = "; ")
    pdf(paste(npath,'.pdf', sep=''), onefile =FALSE)
    plot(hist(ntotal), xlab = 'nGenes/cell', main = paste('cut-off', y, '# genes'), sub = summStr, col="#539cb4") 
    dev.off()

}

