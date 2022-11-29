## A machine learning framework for reads threshold optimization and accurate classification of cell types using scRNA-seq

This repository contains all script and data used to produce figures in the 2022 Frontiers Genetics Publication.


### Contents:
- 001_marker_genes.r -- Identify marker genes
- 002_SingleRclassifier.v1.0.r -- SingleR classifier run
- 003_SCNclassifier.v1.0.r -- Single Cell Net classifier run
- 004_figures.ipynb -- Produce paper figures using the output of the above.
- seu_HQ_no_id.qs -- Seurat object containing high quality cells.
- metadata_subsample_LQ_no_id.txt -- Metadata of low quality cells used in Figure 2. 
- Annotation.celltype.markers.txt -- Cell type marker genes. 
- Annotation.lineage.markers.txt -- Lineage marker genes. 

The script has been throughly annotated to describe its functions.  
Working dirctory paths can be modified as needed. All required sub-directories would be automatically created.

### Analysis pipeline
1. 001_marker_genes.r  -->  Identify marker genes for lineage and cell type levels for SingleCellNet classifier. 
   - Output (included in the repository): 
     - Annotation.celltype.markers.txt
     - Annotation.lineage.markers.txt
2. Run SingleR and/or SingleCellNet classifiers:
   - 002_SingleRclassifier.v1.0.r  -->  Train the SingleR classifier using the HQ dataset and outputs model performance metrics.
   - 003_SCNclassifier.v1.0.r  -->  Train the SingleCellNet classifier using the HQ dataset and outputs model performance metrics.
     - Output: SR/SCN directory containng performance, UMI, gene distribution metrics. Used as input to generate plots. 
3. 004_figures.ipynb  -->  Generate plots. Jupyter notebook cells are annotated by figure. 
   - Output: Plots in PDF format. Can be modified to show plots directly in the notebook.


### Contact us: 
ibishara@coh.org

