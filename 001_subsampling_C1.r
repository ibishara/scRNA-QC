## This script randomly subsamples 20% cells from each batch pre- and post-filtering from metadata and raw counts (FELINE cohort 1)


# Author: "Isaac Bishara"

setwd('/Users/ibishara/Desktop/FELINE_C1/')

# library
library(data.table)
library(stringr)
library(tidyverse)


# data
ref1 <- as.data.frame(fread('post-filter/Feline_metadata_101421.txt', sep='auto')) # post-filtered metadata as a reference for "high-quality" cells | 17 celltype classes
ref1 <- ref1[ref1$Cohort == 'C1', c('Cell', 'Patient', 'nCount_RNA', 'nFeature_RNA', 'Percent.Mitochondria', 'Batch', 'Annotation')] # filter for FELINE C1 
rownames(ref1) <- ref1$Cell

ref2 <- as.data.frame(fread('post-filter/FEL001046_scRNA.metadata_JF.txt', sep='auto')) # post-filtered metadata as a reference for "high-quality" cells | 11 celltype classes
colnames(ref2)[which(names(ref2) == "Cell.ID")] <- "Cell"
rownames(ref2) <- ref2$Cell
common.cells <- intersect(ref1$Cell, ref2$Cell)
ref2 <- ref2[common.cells, ]
ref1 <- ref1[common.cells, ]

# filtered_HQ_cellID <- ref1[ref1$Platform == '10x','Cell.ID'] # cell ids filtered for 10X technology 
meta_HQ <- left_join(ref1, ref2[, c('Cell','Celltype')], by = 'Cell' ) # combine both cell annotations


set.seed(100)
## metadata subsampling
meta_HQ_sub <- meta_HQ[sample(nrow(meta_HQ), round(0.2 * nrow(meta_HQ))), ] # subsample 10% cells 

## Subsample count table for matching cell ids 
## Takes forever to run in the notebook. Run in standard R or Radian
## creates and index of subsamples cells in metadata for each raw count table then subsample each table only for matching cells
all.counts <- list.files(path = "/Users/ibishara/Desktop/FELINE_C1/raw/FELINE_cellranger_premRNA/", pattern = "*counts.txt", recursive = TRUE) # create a list of raw count tables from all batches
count.filter.list <- lapply(paste('raw/FELINE_cellranger_premRNA/', all.counts, sep=''), FUN = function(x) {
count_sub <- fread(x, sep='auto', select= c("Gene Symbol",  meta_HQ_sub$Cell  )) # cell ids subsampled from the metadata
return(count_sub)  # returns warnings for every cell id NOT in each count table
})

count_batches_merged <- as.data.frame(do.call(cbind, count.filter.list))  #join list of tables 
count_batches_merged <- count_batches_merged[,!duplicated(colnames(count_batches_merged))] # remove duplicate cell id
count_batches_merged <- count_batches_merged[!duplicated(count_batches_merged$"Gene Symbol"),] # remove duplicate genes
count_batches_merged <- count_batches_merged[count_batches_merged$'Gene Symbol' != '' ,] # remove null genes

# ordering counts and meta cells alphabetically 
temp <- count_batches_merged[,-1]
sorted_counts <- temp[, order(colnames(temp))] 
sorted_meta <- meta_HQ_sub[order(meta_HQ_sub$Cell) ,]

meta_HQ_sub$nCount_RNA = colSums(sorted_counts)  # corrected nCount_RNA
sorted_counts[sorted_counts > 0] <- 1
meta_HQ_sub$nFeature_RNA = colSums(sorted_counts)  # corrected nFeature_RNA

fwrite(meta_HQ_sub, 'metadata_subsample_HQ.txt', sep='\t') # Export subsampled metadata
fwrite(count_batches_merged, 'raw_counts_subsample_HQ.txt', sep='\t') # Export subsampled raw counts
