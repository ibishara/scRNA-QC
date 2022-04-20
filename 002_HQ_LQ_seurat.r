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

# meta_singler <- as.data.frame(fread('raw/FELINE_C1_raw_singler_metadata.txt', sep='\t')) # full metadata with singleR annotations
# rownames(meta_singler) <- meta_singler$Cell
# meta_singler$V1 <- NULL

metadata_post_filter <- as.data.frame(fread('/Users/ibishara/Desktop/FELINE_C1/post-filter/FEL001046_scRNA.metadata_JF.txt', sep='auto')) # Jinfeng's high quality metadata 
rownames(metadata_post_filter) <- metadata_post_filter$Cell.ID

meta_sub <- metadata_post_filter[colnames(counts), ]

# # data clean-up
# meta_sub <- meta_singler[ colnames(counts), ] # filter metadata with singleR 
# meta_sub <- meta_sub[!duplicated(meta_sub$Cell),] # remove duplicates 

# celltype_anno_post_filter <- metadata_post_filter[,c('Cell.ID','Celltype')] # extract HQ cell type JF's annotations
# colnames(celltype_anno_post_filter)[1] <- 'Cell'
# rownames(celltype_anno_post_filter) <- celltype_anno_post_filter$Cell

# celltype_anno_post_filter <- celltype_anno_post_filter[rownames(meta_sub) , ]
# meta_sub <- left_join(meta_sub, celltype_anno_post_filter, by = 'Cell') # add JF annotations to subsample metadata
# rownames(meta_sub) <- meta_sub$Cell

## assign hpca labels to lineage annotations. JF's annotations "Celltype" are used for lineage annotation of high quality cells 
meta_sub[["Lineage"]] <- meta_sub$Celltype
meta_sub <- meta_sub %>% mutate( Lineage = case_when(
        # meta_sub[["hpca"]] == "Epithelial_cells"  ~ 'Epithelial_cells',
        # meta_sub[["hpca"]] == "Fibroblasts"  ~ 'Mesenchymal_cells',
        # meta_sub[["hpca"]] == "Smooth_muscle_cells"  ~ 'Mesenchymal_cells',
        # meta_sub[["hpca"]] == "Endothelial_cells"  ~ 'Mesenchymal_cells',
        # meta_sub[["hpca"]] == "Chondrocytes"  ~ 'Mesenchymal_cells',
        # meta_sub[["hpca"]] == "Osteoblasts"  ~ 'Mesenchymal_cells',
        # meta_sub[["hpca"]] == "T_cells"  ~ 'Hematopoeitic_cells',
        # meta_sub[["hpca"]] == "B_cell"  ~ 'Hematopoeitic_cells',
        # meta_sub[["hpca"]] == "Macrophage"  ~ 'Hematopoeitic_cells',
        # meta_sub[["hpca"]] == "Monocyte"  ~ 'Hematopoeitic_cells',
        # meta_sub[["hpca"]] == "NK_cell"  ~ 'Hematopoeitic_cells',
        # meta_sub[["hpca"]] == "Neutrophils"  ~ 'Hematopoeitic_cells',
        # meta_sub[["hpca"]] == "Platelets"  ~ 'Hematopoeitic_cells',
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




# # Add quality annotations 
# meta_sub[["Quality"]] <- meta_sub$Celltype
# meta_sub$Quality <- ifelse(is.na(meta_sub$Quality), 'Low', 'High')
    
# export metadata
fwrite(meta_sub,'metadata_subsample_anno.txt', sep='\t', nThread = 16, row.names = TRUE) # metadata with singleR annotations (subsampled) | contains hpca and blueprint annotations


# construct seurat object
seu_HQ <- CreateSeuratObject(counts= counts, min.features= 0, min.cells = 0, names.delim= "_", meta.data= meta_sub) 
qsave(seu_HQ, file="seu_HQ.qs")

# # Split into 2 seurat objects based on quality 
# Idents(seu) <- 'Quality'
# seu_HQ <- subset(seu, Quality == 'High')
# seu_LQ <- subset(seu, Quality == 'Low')
# qsave(seu_LQ, file="seu_LQ.qs")
# qsave(seu_HQ, file="seu_HQ.qs")

seu_HQ <- SCTransform(seu_HQ, method = "glmGamPoi", verbose = TRUE) # "S", "G2M", removed. Doesn't affect Archetype

#dim reduction 
seu_HQ <- RunPCA(seu_HQ, features = NULL) # "Features = NULL" to run on variable genes only


#Cell clustering
seu_HQ <- FindNeighbors(seu_HQ, dims = 1:20) #Matrix package version issue 
seu_HQ <- FindClusters(seu_HQ, resolution = 0.4) # mod from resolution= 0.5 

# Idents(seu_HQ) <- "Lineage"

#UMAP
seu_HQ <- RunUMAP(seu_HQ, dims = 1:20)
DimPlot(seu_HQ, reduction = "umap", group.by = 'Celltype') # infercnv_annotation  FinalAnnotation
