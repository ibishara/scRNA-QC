# attach hpca, create 'Lineage' annotations and create seurat objects for HQ and LQ datasets 


setwd('/Users/ibishara/Desktop/FELINE_C1/')

# packages
library(data.table)
library(tidyverse)
library(Seurat)
library(qs)

# data 
counts <- as.data.frame(fread('raw_counts_subsample_HQ.txt', sep='\t',  check.names =FALSE)) # raw counts (subsampled)
counts <- counts[counts$"Gene Symbol" != '' ,] 
rownames(counts) <- counts$"Gene Symbol"
counts$"Gene Symbol" <- NULL

meta_sub <- as.data.frame(fread('metadata_subsample_HQ.txt', sep='auto')) # Jinfeng's high quality metadata 
rownames(meta_sub) <- meta_sub$Cell
meta_sub$V1 <- NULL

# export metadata
fwrite(meta_sub,'metadata_subsample_anno.txt', sep='\t', nThread = 16, row.names = TRUE) # metadata with singleR annotations (subsampled) | contains hpca and blueprint annotations

# construct seurat object
seu_HQ <- CreateSeuratObject(counts= counts, min.features= 0, min.cells = 0, names.delim= "_", meta.data= meta_sub) 
qsave(seu_HQ, file="seu_HQ.qs")


LQ_count_batches_merged <- as.data.frame(do.call(cbind, count.filter.list))  #join list of tables 
LQ_count_batches_merged <- LQ_count_batches_merged[,!duplicated(colnames(LQ_count_batches_merged))] # remove duplicate cell id
LQ_count_batches_merged <- LQ_count_batches_merged[!duplicated(LQ_count_batches_merged$"Gene Symbol"),] # remove duplicate genes
LQ_count_batches_merged <- LQ_count_batches_merged[LQ_count_batches_merged$'Gene Symbol' != '' ,] # remove null genes

# ordering counts and meta cells alphabetically 
temp <- LQ_count_batches_merged[,-1]