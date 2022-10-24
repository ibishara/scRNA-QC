## A machine learning framework for reads threshold optimization and accurate classification of cell types using scRNA-seq

This repository contains all script and data used to produce figures in the 2022 Frontiers Genetics Publication.


### Contents:
- 001_marker_genes.r --    
- 002_SingleRclassifier.v1.0.r -- 
- 003_SCNclassifier.v1.0.r -- Train the SingleCellNet classifier using the HQ dataset and outputs model performance metrics.
- 004_figures.ipynb -- Produce paper figures using the output of the above.
- seu_HQ_no_id.qs -- Seurat object containing high quality cells.
- metadata_subsample_LQ_no_id.txt -- Metadata of low quality cells used in Figure 2. 
- Annotation.celltype.markers.txt -- Cell type marker genes. 
- Annotation.lineage.markers.txt -- Lineage marker genes. 

The script has been throughly annotated to describe its functions.  
Working dirctory paths can be modified as needed. All required sub-directories would be automatically created.

### Analysis pipeline
1. 001_marker_genes.r  -->  Identify marker genes for lineage and cell type levels. 
   - Output (included in the repository): 
     - Annotation.celltype.markers.txt
     - Annotation.lineage.markers.txt
2. SingleR or SingleCellNet classifiers:
   - 002_SingleRclassifier.v1.0.r  -->  Train the SingleR classifier using the HQ dataset and outputs model performance metrics.
