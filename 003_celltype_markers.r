


setwd('/Users/ibishara/Desktop/FELINE_C1/')
#packages
library(Seurat)

# data
seu_HQ <- qread(file = "seu_HQ.qs", nthreads = 16)

# Identify Lineage markers in HQ cells 
Idents(seu_HQ) <- "Lineage"

Annotation.lineage.latentCC.markers <- FindAllMarkers(seu_HQ, test.use="negbinom", 
                                                   min.pct=0.5, max.cells.per.ident=2000, logfc.threshold=0.5)

write.table(Annotation.lineage.latentCC.markers, file="Annotation.lineage.markers.txt", sep = '\t')

# Identify Celltype markers in HQ cells 
Idents(seu_HQ) <- "Celltype"
Annotation.celltype.latentCC.markers <- FindAllMarkers(seu_HQ, test.use="negbinom", 
                                                   min.pct=0.5, max.cells.per.ident=2000, logfc.threshold=0.5)

write.table(Annotation.celltype.latentCC.markers, file="Annotation.celltype.markers.txt", sep = '\t')
