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
ref1 <- ref1[ref1$Cohort == 'C1', c('Cell', 'nCount_RNA', 'nFeature_RNA', 'Percent.Mitochondria', 'Annotation')] # filter for FELINE C1 
rownames(ref1) <- ref1$Cell

ref2 <- as.data.frame(fread('post-filter/FEL001046_scRNA.metadata_JF.txt', sep='auto')) # post-filtered metadata as a reference for "high-quality" cells | 11 celltype classes
colnames(ref2)[which(names(ref2) == "Cell.ID")] <- "Cell"
rownames(ref2) <- ref2$Cell
common.cells <- intersect(ref1$Cell, ref2$Cell)
ref2 <- ref2[common.cells, ]
ref1 <- ref1[common.cells, ]

# FELINE C1 high quality metadata with cell type annotations 
meta_HQ <- left_join(ref1, ref2[, c('Cell','Celltype')], by = 'Cell' ) # combine both cell type annotations

sample.ncell <- 35000 # number of cells to sample for HQ and LQ


###############################
### HQ metadata subsampling ###
###############################
set.seed(100)
meta_HQ_sub <- meta_HQ[sample(nrow(meta_HQ), sample.ncell), ] 

# Generate "Lineage" annotations
meta_HQ_sub[["Lineage"]] <- meta_HQ_sub$Celltype
meta_HQ_sub <- meta_HQ_sub %>% mutate( Lineage = case_when(
# overwrite HQ cells according to JF annotations | Only HQ cells have "Celltype" annotaions
    meta_HQ_sub[["Celltype"]] == "Cancer cells"  ~ 'Epithelial_cells',
    meta_HQ_sub[["Celltype"]] == "Normal epithelial cells"  ~ 'Epithelial_cells',
    meta_HQ_sub[["Celltype"]] == "Adipocytes"  ~ 'Stromal_cells',
    meta_HQ_sub[["Celltype"]] == "Fibroblasts"  ~ 'Stromal_cells',
    meta_HQ_sub[["Celltype"]] == "Endothelial cells"  ~ 'Stromal_cells',
    meta_HQ_sub[["Celltype"]] == "Pericytes"  ~ 'Stromal_cells',
    meta_HQ_sub[["Celltype"]] == "Macrophages"  ~ 'Immune_cells',
    meta_HQ_sub[["Celltype"]] == "T cells"  ~ 'Immune_cells',
    meta_HQ_sub[["Celltype"]] == "B cells"  ~ 'Immune_cells'
))

## Subsample count table for matching cell ids ##
## Inefficient to run in a notebook. Run in base R or Radian
## creates and index of subsamples cells in metadata for each raw count table then subsample each table only for matching cells

all.counts <- list.files(path = "/Users/ibishara/Desktop/FELINE_C1/raw/FELINE_cellranger_premRNA/", pattern = "*counts.txt", recursive = TRUE) # create a list of raw count tables from all batches
count.filter.list <- mclapply(paste('raw/FELINE_cellranger_premRNA/', all.counts, sep=''), FUN = function(x) {
    count_sub <- fread(x, sep='auto', select= c("Gene Symbol",  meta_HQ_sub$Cell )) # cell ids subsampled from the metadata
    return(count_sub)  
}, mc.cores= numCores)

HQ_count_batches_merged <- as.data.frame(do.call(cbind, count.filter.list))  #join list of tables 
HQ_count_batches_merged <- HQ_count_batches_merged[,!duplicated(colnames(HQ_count_batches_merged))] # remove duplicate cell id
HQ_count_batches_merged <- HQ_count_batches_merged[!duplicated(HQ_count_batches_merged$"Gene Symbol"),] # remove duplicate genes
HQ_count_batches_merged <- HQ_count_batches_merged[HQ_count_batches_merged$'Gene Symbol' != '' ,] # remove missing genes
rownames(HQ_count_batches_merged) <- HQ_count_batches_merged$'Gene Symbol'
HQ_count_batches_merged$'Gene Symbol' <- NULL

fwrite(HQ_count_batches_merged, 'raw_counts_subsample_HQ.txt', sep='\t') # Export subsampled raw counts. includes cell ids
fwrite(meta_HQ_sub, 'metadata_subsample_HQ.txt', sep='\t') # Export subsampled metadata. includes cell ids

# construct HQ seurat object
rownames(meta_HQ_sub) <- meta_HQ_sub$Cell 
seu_HQ <- CreateSeuratObject(counts= HQ_count_batches_merged, min.features= 0, min.cells = 0, names.delim= "_", meta.data= meta_HQ_sub) 

# order counts and meta cells alphabetically to correct nCount_RNA &  nFeature_RNA
seu.HQ.counts <- GetAssayData(seu_HQ, assay = "RNA")
seu.HQ.counts <- seu.HQ.counts[, order(colnames(seu.HQ.counts))] # Edit
sorted_counts <- seu.HQ.counts
# sorted_counts <- seu.HQ.counts[, order(colnames(seu.HQ.counts))] # Edit

seu_HQ@meta.data <- seu_HQ@meta.data[order(rownames(seu_HQ@meta.data)) ,]

all(colnames(sorted_counts) ==  rownames(seu_HQ@meta.data)) # test

seu_HQ@meta.data$nCount_RNA = colSums(sorted_counts)  # corrected nCount_RNA
seu_HQ@meta.data$nFeature_RNA = apply(sorted_counts, 2, function (x) sum(x > 0))  # corrected nFeature_RNA

##  so far so good

seu_HQ@meta.data <- seu_HQ@meta.data[, c('Cell', 'Percent.Mitochondria', 'Celltype', 'nCount_RNA', 'nFeature_RNA', 'Lineage')]
qsave(seu_HQ, file="seu_HQ_no_id1.qs") # export seurat obj

all(colnames(seu.HQ.counts) ==  rownames(seu_HQ@meta.data)) # test

# Hide patient ids 
seu_HQ_rename <- RenameCells(seu_HQ, new.names =  paste('Cell', seq(1, ncol(sorted_counts), 1), sep='_') ) 
# ############################

# meta <- seu_HQ@meta.data
# seu.HQ.counts <- GetAssayData(seu_HQ, assay = "RNA")

# meta_rename <- seu_HQ_rename@meta.data
# seu.HQ.counts_rename <- GetAssayData(seu_HQ_rename, assay = "RNA")


# all(meta == meta_rename)
# all(seu.HQ.counts == seu.HQ.counts_rename)


# ############################
seu_HQ_rename@meta.data$Cell <- rownames(seu_HQ_rename@meta.data)

qsave(seu_HQ_rename, file="seu_HQ_no_id2.qs") # export seurat obj | out of subscript



meta <- seu_HQ@meta.data
seu.HQ.counts <- as.data.frame(GetAssayData(seu_HQ, assay = "RNA"))

fwrite(seu.HQ.counts, 'raw_counts_subsample_HQ_no_id.txt', sep='\t') # Export subsampled raw counts. no cell ids
fwrite(meta, 'metadata_subsample_HQ_no_id.txt', sep='\t') # Export subsampled metadata. no cell ids




#################################
#### LQ metadata subsampling ####
#################################
set.seed(200)
meta_LQ <- meta.data[ !meta.data$Cell %in% meta_HQ$Cell, ] # remove HQ cells 
meta_LQ_sub <- meta_LQ[sample(nrow(meta_LQ), sample.ncell), ] # random sampling
meta_LQ_sub$V1 <- NULL



# Read corresponding cell counts to correct nCounts_RNA and nFeature_RNA 
all.counts <- list.files(path = "/Users/ibishara/Desktop/FELINE_C1/raw/FELINE_cellranger_premRNA/", pattern = "*counts.txt", recursive = TRUE) # create a list of raw count tables from all batches
count.filter.list <- mclapply(paste('raw/FELINE_cellranger_premRNA/', all.counts, sep=''), FUN = function(x) {
count_sub <- fread(x, sep='auto', select= c("Gene Symbol",  meta_LQ_sub$Cell  )) # cell ids subsampled from the metadata
return(count_sub) 
}, mc.cores= numCores)

LQ_count_batches_merged <- as.data.frame(do.call(cbind, count.filter.list))  #join list of tables 
LQ_count_batches_merged <- LQ_count_batches_merged[,!duplicated(colnames(LQ_count_batches_merged))] # remove duplicate cell id
LQ_count_batches_merged <- LQ_count_batches_merged[!duplicated(LQ_count_batches_merged$"Gene Symbol"),] # remove duplicate genes
LQ_count_batches_merged <- LQ_count_batches_merged[LQ_count_batches_merged$'Gene Symbol' != '' ,] # remove null genes
rownames(LQ_count_batches_merged) <- LQ_count_batches_merged$'Gene Symbol'
LQ_count_batches_merged$'Gene Symbol' <- NULL

# ordering counts and meta cells alphabetically 
sorted_counts <- LQ_count_batches_merged[, order(colnames(LQ_count_batches_merged))]
sorted_meta <- meta_LQ_sub[order(meta_LQ_sub$Cell) ,]

all(colnames(sorted_counts) ==  sorted_meta$Cell) # test

sorted_meta$nCount_RNA = colSums(sorted_counts)  # corrected nCount_RNA
sorted_meta$nFeature_RNA = apply(sorted_counts, 2, function (x) sum(x > 0))  # corrected nFeature_RNA

meta_LQ_sub <- sorted_meta

# Lineage annotations based off hpca annotations 
meta_LQ_sub[["Lineage"]] <- meta_LQ_sub$hpca
meta_LQ_sub <- meta_LQ_sub %>% mutate( Lineage = case_when(
        meta_LQ_sub[["hpca"]] == "Epithelial_cells"  ~ 'Epithelial_cells',
        # meta_LQ_sub[["hpca"]] == "Neurons"  ~ 'Epithelial_cells',
        # meta_LQ_sub[["hpca"]] == "Chondrocytes"  ~ 'Epithelial_cells',
        meta_LQ_sub[["hpca"]] == "Fibroblasts"  ~ 'Stromal_cells',
        meta_LQ_sub[["hpca"]] == "Smooth_muscle_cells"  ~ 'Stromal_cells',
        meta_LQ_sub[["hpca"]] == "Endothelial_cells"  ~ 'Stromal_cells',
        meta_LQ_sub[["hpca"]] == "Chondrocytes"  ~ 'Stromal_cells',
        meta_LQ_sub[["hpca"]] == "Osteoblasts"  ~ 'Stromal_cells',
        meta_LQ_sub[["hpca"]] == "T_cells"  ~ 'Immune_cells',
        meta_LQ_sub[["hpca"]] == "B_cell"  ~ 'Immune_cells',
        meta_LQ_sub[["hpca"]] == "Macrophage"  ~ 'Immune_cells',
        meta_LQ_sub[["hpca"]] == "Monocyte"  ~ 'Immune_cells',
        meta_LQ_sub[["hpca"]] == "NK_cell"  ~ 'Immune_cells',
        meta_LQ_sub[["hpca"]] == "Neutrophils"  ~ 'Immune_cells'
))

names(meta_LQ_sub)[names(meta_LQ_sub) == 'Percent Mitochondria'] <- 'Percent.Mitochondria' # match names in HQ metadata 

meta_LQ_sub <- meta_LQ_sub[, c('Cell', 'Percent.Mitochondria', 'hpca', 'nCount_RNA', 'nFeature_RNA', 'Lineage')]

fwrite(meta_LQ_sub, 'metadata_subsample_LQ.txt', sep='\t') # Export subsampled metadata

## Randomize patient id
meta_LQ_sub$Cell <- paste('Cell', seq(1, nrow(meta_LQ_sub), 1), sep='_') # hide patient ids

fwrite(meta_LQ_sub, 'metadata_subsample_LQ_no_id.txt', sep='\t') # Export subsampled metadata


