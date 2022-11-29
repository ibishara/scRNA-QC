# Author: "Isaac Bishara"
# identify marker genes for lineage and cell type levels 


setwd('/Users/ibishara/Desktop/FELINE_C1/')
#packages
library(Seurat)
library(qs)

# dataset <- c( 'feline', 'combes', 'pbmc' )
dataset <- 'feline' # select dataset to be analyzed 
dir.create(dataset, recursive = TRUE)

if (dataset == 'feline'){
    #########High quality FELINE C1 dataset###############
    seu <- qread(file = "seu_HQ_no_id.qs", nthreads = numCores)
    seu <- subset(seu, subset = Celltype != "Normal epithelial cells")   ## Removes normal epithelial cells. Genes unique to normal epi cells are removed from analysis downstream
    seu <- subset(seu, subset = nCount_RNA < 15000 ) # filter out cells with > 15k UMI 
   # Identify Lineage markers 
    Idents(seu) <- "Lineage"
    Annotation.lineage.markers <- FindAllMarkers(seu, test.use="negbinom", 
                                                    min.pct=0.5, max.cells.per.ident=2000, logfc.threshold=0.5)
    write.table(Annotation.lineage.markers, file=paste(dataset, "Annotation.lineage.markers.txt", sep='/'), sep = '\t')

}else if(dataset == 'combes'){
    #########Combes dataset###############
    seu <- qread('/Users/ibishara/Desktop/FELINE_C1/Additional_datasets/Combes_whole_blood/sobj_SCG1-2-3.qs', nthreads = numCores )
    seu <- subset(seu, subset = Celltype !='RBC' )

}else if(dataset == 'pbmc'){
    #########pbmc3k dataset###############
    # After running pbmc3k_tutorial.rmd #
    seu <- readRDS('/Users/ibishara/Desktop/FELINE_C1/Additional_datasets/PBMC3k/pbmc3k_final.rds')

}


# Identify Celltype markers 
Idents(seu) <- "Celltype"
Annotation.celltype.markers <- FindAllMarkers(seu, test.use="negbinom", 
                                                   min.pct=0.5, max.cells.per.ident=2000, logfc.threshold=0.5)
write.table(Annotation.celltype.markers, file=paste(dataset,"Annotation.celltype.markers.txt", sep='/'), sep = '\t')
