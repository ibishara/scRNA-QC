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

## assign hpca labels to lineage annotations. JF's annotations "Celltype" are used for lineage annotation of high quality cells 
meta_sub[["Lineage"]] <- meta_sub$Celltype
meta_sub <- meta_sub %>% mutate( Lineage = case_when(
# relabel HQ cells according to JF annotations | Only HQ cells have "Celltype" annotaions
        meta_sub[["Celltype"]] == "Cancer cells"  ~ 'Epithelial_cells',
        meta_sub[["Celltype"]] == "Normal epithelial cells"  ~ 'Epithelial_cells',
        meta_sub[["Celltype"]] == "Adipocytes"  ~ 'Mesenchymal_cells',
        meta_sub[["Celltype"]] == "Fibroblasts"  ~ 'Mesenchymal_cells',
        meta_sub[["Celltype"]] == "Endothelial cells"  ~ 'Mesenchymal_cells',
        meta_sub[["Celltype"]] == "Pericytes"  ~ 'Mesenchymal_cells',
        meta_sub[["Celltype"]] == "Macrophages"  ~ 'Hematopoeitic_cells',
        meta_sub[["Celltype"]] == "T cells"  ~ 'Hematopoeitic_cells',
        meta_sub[["Celltype"]] == "B cells"  ~ 'Hematopoeitic_cells'
))


# export metadata
fwrite(meta_sub,'metadata_subsample_anno.txt', sep='\t', nThread = 16, row.names = TRUE) # metadata with singleR annotations (subsampled) | contains hpca and blueprint annotations

# construct seurat object
seu_HQ <- CreateSeuratObject(counts= counts, min.features= 0, min.cells = 0, names.delim= "_", meta.data= meta_sub) 
qsave(seu_HQ, file="seu_HQ.qs")




# seu_HQ <- qread(file = "seu_HQ.qs", nthreads = 16)

# # # seu_HQ <- SCTransform(seu_HQ, method = "glmGamPoi", verbose = TRUE) # "S", "G2M", removed.


# #dim reduction 
# seu_HQ <- RunPCA(seu_HQ, features = ) # "Features = NULL" to run on variable genes only


# #Cell clustering
# seu_HQ <- FindNeighbors(seu_HQ, dims = 1:20) #Matrix package version issue 
# seu_HQ <- FindClusters(seu_HQ, resolution = 0.4) # mod from resolution= 0.5 

# # Idents(seu_HQ) <- "Lineage"

# #UMAP
# seu_HQ <- RunUMAP(seu_HQ, dims = 1:20)
# DimPlot(seu_HQ, reduction = "umap", group.by = 'Lineage') # infercnv_annotation  FinalAnnotation


# ggplot(seu_HQ@meta.data, aes(seu_HQ@meta.data$nFeature_RNA, log10( seu_HQ@meta.data$nCount_RNA))) +
# geom_point()