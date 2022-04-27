
library(data.table)
library(dplyr)
library(qs)
library(Seurat)
setwd('/Users/ibishara/Desktop/FELINE_C1/')

# data
seu_HQ <- qread(file="seu_HQ.qs", nthreads=16)
seu_HQ <- subset(x = seu_HQ, subset = Celltype != "Normal epithelial cells")   ## Removes normal epithelial cells. Genes unique to normal epi cells are removed from downstream analysis 

counts <- as.data.frame(GetAssayData(seu_HQ, assay = "RNA"))
set.seed(100)
numCores <- detectCores()
numCores

# functions
# function extract cells to be transformed ( total counts above threshold ) then hand filtered dataframe to transformation function
# x = dataframe to be transformed 
# y = max total counts per cell threshold
parent <- function(x, y, fun){
    x <- x[, which(colSums(x) > y)] # filters for cells with total counts above threshold (to be transformed)
    x <- as.data.frame(mclapply(X = x, FUN= fun, y = y, mc.cores= numCores )) # run 2nd function to reduce the number of count per cell above threshold. iterates over columns (cells)

    return(x)
}

# Function to reduce the number of count per cell above threshold
# function iterates over columns (cells). caluclulated total counts per cell. multiply gene counts by a fraction to yield total counts per cell above threshold
# x = column in filtered dataframe
# y = max total counts per cell threshold
fl00r <- function(x, y){
    total <- sum(x) # total counts for current cell 
    x <- floor((y/total)*x) # mutliply all gene counts by a fraction. genes with 0 counts are floored back to 0. other gene counts are floored to nearest integer

    return(x)
}

nofl00r <- function(x, y){
    total <- sum(x) # total counts for current cell 
    x <- round((y/total)*x, 1) # mutliply all gene counts by a fraction with rounding to 2 significant figures (reduce size)

    return(x)
}

r0und <- function(x, y){
    total <- sum(x) # total counts for current cell 
    x <- round((y/total)*x, 0) # mutliply all gene counts by a fraction rounding to integers

    return(x)
}


# export
threshold <- c(200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 4000, 5000)

for (i in threshold){
    y <- i
    f <- parent(x = counts, y = y, fun = fl00r) # transformed dataframe (floored)
    # n <- parent(x = counts, y = y, fun = nofl00r) # transformed dataframe (no-floor)
    # r <- parent(x = counts, y = y, fun = r0und) # transformed dataframe (round to integers)

    fpath <- paste('downsample/reads_downsample/floor/reads_down_', format(y/1000, nsmall=1), sep='')
    # npath <- paste('downsample/reads_downsample/nofloor/reads_down_', format(y/1000, nsmall=1), sep='')
    # rpath <- paste('downsample/reads_downsample/round/reads_down_', format(y/1000, nsmall=1), sep='')

    fwrite(f, paste(fpath, '.txt', sep=''), sep='\t', nThread = 16, row.names = TRUE)
    # fwrite(n, paste(npath, '.txt', sep=''), sep='\t', nThread = 16, row.names = TRUE) 
    # fwrite(r, paste(rpath, '.txt', sep=''), sep='\t', nThread = 16, row.names = TRUE) 


    # export floor hist and stats 
    ftotal <- colSums(f)
    sdat <- summary(ftotal)   
    summStr <- paste(names(sdat), format(sdat, digits = 2), collapse = "; ")
    op <- par(mar = c(7,4,4,2) + 0.1, cex = 0.5)
    pdf(paste(fpath,'.pdf', sep=''), onefile =FALSE)
    plot(hist(ftotal), xlab = 'nReads/cell', main = paste('cut-off', y, 'UMI'), sub = summStr, col="#1e72d2") # ;     par(op)
    dev.off()


    # # export nofloor hist and stats 
    # ntotal <- colSums(n)
    # sdat <- summary(ntotal)   
    # summStr <- paste(names(sdat), format(sdat, digits = 2), collapse = "; ")
    # op <- par(mar = c(7,4,4,2) + 0.1, cex = 0.5)
    # pdf(paste(npath,'.pdf', sep=''), onefile =FALSE)
    # plot(hist(ntotal), xlab = 'nReads/cell', main = paste('cut-off', y, 'UMI'), sub = summStr, col="#d74a4a") # ;     par(op)
    # dev.off()

    # # export round hist and stats 
    # rtotal <- colSums(r)
    # sdat <- summary(rtotal)   
    # summStr <- paste(names(sdat), format(sdat, digits = 2), collapse = "; ")
    # op <- par(mar = c(7,4,4,2) + 0.1, cex = 0.5)
    # pdf(paste(rpath,'.pdf', sep=''), onefile =FALSE)
    # plot(hist(rtotal), xlab = 'nReads/cell', main = paste('cut-off', y, 'UMI'), sub = summStr, col="#56b453") # ;     par(op)
    # dev.off()

}

 