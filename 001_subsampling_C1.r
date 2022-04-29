## This script consists of 3 sections :
## 1. Subsample n cells from High Quality cells metadata including cell type annotaions (cite). Cell type annotations are used to construct lineage annotations.
## 2. Subsample read count matrix for corresponding high quality cells. HQ counts and metadata are used to construct HQ_seurat obj.
## 3. Subsample n cells from non-HQ cells (low quality cells) metadata.
## Note: nCounts and nFeature columns are recalculated from corresponding count matrix 

## HQ cells use JF's annotations 'Celltype', 'Annotation'
## LQ cells use hpca SingleR annotations
# Author: "Isaac Bishara"

setwd('/Users/ibishara/Desktop/FELINE_C1/')

# library
library(data.table)
library(stringr)
library(tidyverse)
library(parallel)
library(Seurat)
library(qs)

numCores <- detectCores()
numCores 

# data
meta.data  <- as.data.frame(fread('raw/FELINE_C1_raw_singler_metadata.txt', sep='\t')) # Full FELINE C1 metadata + SingleR annotations
ref1 <- as.data.frame(fread('post-filter/Feline_metadata_101421.txt', sep='auto')) # post-filtered metadata as a reference for "high-quality" cells | 17 celltype classes
ref1 <- ref1[ref1$Cohort == 'C1', c('Cell', 'Patient', 'nCount_RNA', 'nFeature_RNA', 'Percent.Mitochondria', 'Batch', 'Annotation')] # filter for FELINE C1 
rownames(ref1) <- ref1$Cell

ref2 <- as.data.frame(fread('post-filter/FEL001046_scRNA.metadata_JF.txt', sep='auto')) # post-filtered metadata as a reference for "high-quality" cells | 11 celltype classes
colnames(ref2)[which(names(ref2) == "Cell.ID")] <- "Cell"
rownames(ref2) <- ref2$Cell
common.cells <- intersect(ref1$Cell, ref2$Cell)
ref2 <- ref2[common.cells, ]
ref1 <- ref1[common.cells, ]

# FELINE C1 high quality metadata with cell type annotations 
meta_HQ <- left_join(ref1, ref2[, c('Cell','Celltype')], by = 'Cell' ) # combine both cell type annotations

set.seed(100)
sample.ncell <- 35000 # number of cells to sample for HQ and LQ

## HQ metadata subsampling ##
meta_HQ_sub <- meta_HQ[sample(nrow(meta_HQ), sample.ncell), ] # subsample 20% cells 

# Attach hpca annotation for HQ 
meta_HQ_sub <- left_join(meta_HQ_sub, meta.data[, c('Cell', 'hpca')])
rownames(meta_HQ_sub) <- meta_HQ_sub$Cell


# # initially, hpca labels are used for lineage annotations 
# meta_HQ_sub[["Lineage"]] <- meta_HQ_sub$Celltype
# meta_HQ_sub <- meta_HQ_sub %>% mutate( Lineage = case_when(
#     meta_HQ_sub[["hpca"]] == "Epithelial_cells"  ~ 'Epithelial_cells',
#     meta_HQ_sub[["hpca"]] == "Fibroblasts"  ~ 'Mesenchymal_cells',
#     meta_HQ_sub[["hpca"]] == "Smooth_muscle_cells"  ~ 'Mesenchymal_cells',
#     meta_HQ_sub[["hpca"]] == "Endothelial_cells"  ~ 'Mesenchymal_cells',
#     meta_HQ_sub[["hpca"]] == "Chondrocytes"  ~ 'Mesenchymal_cells',
#     meta_HQ_sub[["hpca"]] == "Osteoblasts"  ~ 'Mesenchymal_cells',
#     meta_HQ_sub[["hpca"]] == "T_cells"  ~ 'Hematopoeitic_cells',
#     meta_HQ_sub[["hpca"]] == "B_cell"  ~ 'Hematopoeitic_cells',
#     meta_HQ_sub[["hpca"]] == "Macrophage"  ~ 'Hematopoeitic_cells',
#     meta_HQ_sub[["hpca"]] == "Monocyte"  ~ 'Hematopoeitic_cells',
#     meta_HQ_sub[["hpca"]] == "NK_cell"  ~ 'Hematopoeitic_cells',
#     meta_HQ_sub[["hpca"]] == "Neutrophils"  ~ 'Hematopoeitic_cells',
#     meta_HQ_sub[["hpca"]] == "Platelets"  ~ 'Hematopoeitic_cells',
# # overwrite HQ cells according to JF annotations | Only HQ cells have "Celltype" annotaions
#     meta_HQ_sub[["Celltype"]] == "Cancer cells"  ~ 'Epithelial_cells',
#     meta_HQ_sub[["Celltype"]] == "Normal epithelial cells"  ~ 'Epithelial_cells',
#     meta_HQ_sub[["Celltype"]] == "Adipocytes"  ~ 'Mesenchymal_cells',
#     meta_HQ_sub[["Celltype"]] == "Fibroblasts"  ~ 'Mesenchymal_cells',
#     meta_HQ_sub[["Celltype"]] == "Endothelial cells"  ~ 'Mesenchymal_cells',
#     meta_HQ_sub[["Celltype"]] == "Pericytes"  ~ 'Mesenchymal_cells',
#     meta_HQ_sub[["Celltype"]] == "Macrophages"  ~ 'Hematopoeitic_cells',
#     meta_HQ_sub[["Celltype"]] == "T cells"  ~ 'Hematopoeitic_cells',
#     meta_HQ_sub[["Celltype"]] == "B cells"  ~ 'Hematopoeitic_cells'
# ))




## Subsample count table for matching cell ids ##
## Inefficient to run in a notebook. Run in base R or Radian
## creates and index of subsamples cells in metadata for each raw count table then subsample each table only for matching cells
system.time({
all.counts <- list.files(path = "/Users/ibishara/Desktop/FELINE_C1/raw/FELINE_cellranger_premRNA/", pattern = "*counts.txt", recursive = TRUE) # create a list of raw count tables from all batches
count.filter.list <- mclapply(paste('raw/FELINE_cellranger_premRNA/', all.counts, sep=''), FUN = function(x) {
count_sub <- fread(x, sep='auto', select= c("Gene Symbol",  meta_HQ_sub$Cell  )) # cell ids subsampled from the metadata
return(count_sub)  # returns warnings for every cell id NOT in each count table
}, mc.cores= numCores)

HQ_count_batches_merged <- as.data.frame(do.call(cbind, count.filter.list))  #join list of tables 
HQ_count_batches_merged <- HQ_count_batches_merged[,!duplicated(colnames(HQ_count_batches_merged))] # remove duplicate cell id
HQ_count_batches_merged <- HQ_count_batches_merged[!duplicated(HQ_count_batches_merged$"Gene Symbol"),] # remove duplicate genes
HQ_count_batches_merged <- HQ_count_batches_merged[HQ_count_batches_merged$'Gene Symbol' != '' ,] # remove null genes

# ordering counts and meta cells alphabetically 
temp <- HQ_count_batches_merged[,-1]
sorted_counts <- temp[, order(colnames(temp))] 
sorted_meta <- meta_HQ_sub[order(meta_HQ_sub$Cell) ,]

sorted_meta$nCount_RNA = colSums(sorted_counts)  # corrected nCount_RNA
sorted_counts[sorted_counts > 0] <- 1
sorted_meta$nFeature_RNA = colSums(sorted_counts)  # corrected nFeature_RNA
meta_HQ_sub <- sorted_meta

})

fwrite(HQ_count_batches_merged, 'raw_counts_subsample_HQ.txt', sep='\t') # Export subsampled raw counts
fwrite(meta_HQ_sub, 'metadata_subsample_HQ.txt', sep='\t') # Export subsampled metadata

# construct HQ seurat object
seu_HQ <- CreateSeuratObject(counts= HQ_count_batches_merged, min.features= 0, min.cells = 0, names.delim= "_", meta.data= meta_HQ_sub) 
qsave(seu_HQ, file="seu_HQ.qs")





## LQ metadata subsampling ##

meta_LQ <- meta.data[ !meta.data$Cell %in% meta_HQ$Cell, ] # remove HQ cells 
meta_LQ_sub <- meta_LQ[sample(nrow(meta_LQ), sample.ncell), ] #
meta_LQ_sub$V1 <- NULL

# Read corresponding cell counts to correct nCounts_RNA and nFeature_RNA 
all.counts <- list.files(path = "/Users/ibishara/Desktop/FELINE_C1/raw/FELINE_cellranger_premRNA/", pattern = "*counts.txt", recursive = TRUE) # create a list of raw count tables from all batches
count.filter.list <- mclapply(paste('raw/FELINE_cellranger_premRNA/', all.counts, sep=''), FUN = function(x) {
count_sub <- fread(x, sep='auto', select= c("Gene Symbol",  meta_LQ_sub$Cell  )) # cell ids subsampled from the metadata
return(count_sub)  # returns warnings for every cell id NOT in each count table
}, mc.cores= numCores)

LQ_count_batches_merged <- as.data.frame(do.call(cbind, count.filter.list))  #join list of tables 
LQ_count_batches_merged <- LQ_count_batches_merged[,!duplicated(colnames(LQ_count_batches_merged))] # remove duplicate cell id
LQ_count_batches_merged <- LQ_count_batches_merged[!duplicated(LQ_count_batches_merged$"Gene Symbol"),] # remove duplicate genes
LQ_count_batches_merged <- LQ_count_batches_merged[LQ_count_batches_merged$'Gene Symbol' != '' ,] # remove null genes

# ordering counts and meta cells alphabetically 
temp <- LQ_count_batches_merged[,-1]
sorted_counts <- temp[, order(colnames(temp))] 
sorted_meta <- meta_LQ_sub[order(meta_LQ_sub$Cell) ,]

sorted_meta$nCount_RNA = colSums(sorted_counts)  # corrected nCount_RNA
sorted_counts[sorted_counts > 0] <- 1
sorted_meta$nFeature_RNA = colSums(sorted_counts)  # corrected nFeature_RNA
meta_LQ_sub <- sorted_meta


# # Lineage annotations based off hpca annotations 
# meta_LQ_sub[["Lineage"]] <- meta_LQ_sub$hpca
# meta_LQ_sub <- meta_LQ_sub %>% mutate( Lineage = case_when(
#         meta_LQ_sub[["hpca"]] == "Epithelial_cells"  ~ 'Epithelial_cells',
#         meta_LQ_sub[["hpca"]] == "Fibroblasts"  ~ 'Mesenchymal_cells',
#         meta_LQ_sub[["hpca"]] == "Smooth_muscle_cells"  ~ 'Mesenchymal_cells',
#         meta_LQ_sub[["hpca"]] == "Endothelial_cells"  ~ 'Mesenchymal_cells',
#         meta_LQ_sub[["hpca"]] == "Chondrocytes"  ~ 'Mesenchymal_cells',
#         meta_LQ_sub[["hpca"]] == "Osteoblasts"  ~ 'Mesenchymal_cells',
#         meta_LQ_sub[["hpca"]] == "T_cells"  ~ 'Hematopoeitic_cells',
#         meta_LQ_sub[["hpca"]] == "B_cell"  ~ 'Hematopoeitic_cells',
#         meta_LQ_sub[["hpca"]] == "Macrophage"  ~ 'Hematopoeitic_cells',
#         meta_LQ_sub[["hpca"]] == "Monocyte"  ~ 'Hematopoeitic_cells',
#         meta_LQ_sub[["hpca"]] == "NK_cell"  ~ 'Hematopoeitic_cells',
#         meta_LQ_sub[["hpca"]] == "Neutrophils"  ~ 'Hematopoeitic_cells',
#         meta_LQ_sub[["hpca"]] == "Platelets"  ~ 'Hematopoeitic_cells'
# ))

names(meta_LQ_sub)[names(meta_LQ_sub) == 'Percent Mitochondria'] <- 'Percent.Mitochondria' # match names in HQ metadata 

fwrite(meta_LQ_sub, 'metadata_subsample_LQ.txt', sep='\t') # Export subsampled metadata


