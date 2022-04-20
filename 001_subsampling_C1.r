## This script randomly subsamples 10% cells from each batch pre- and post-filtering from metadata and raw counts (FELINE cohort 1)


# Author: "Isaac Bishara"

setwd('/Users/ibishara/Desktop/FELINE_C1/')

# library
library(data.table)
library(stringr)


# data
ref <- as.data.frame(fread('/Users/ibishara/Desktop/FELINE_C1/post-filter/FEL001046_scRNA.metadata_JF.txt', sep='auto')) # post-filtered metadata as a reference for "high-quality" cells
filtered_HQ_cellID <- ref[ref$Platform == '10x','Cell.ID'] # cell ids filtered for 10X technology 
dim(ref)

set.seed(100)
## metadata subsampling
all.meta <- list.files(path = "/Users/ibishara/Desktop/FELINE_C1/raw/FELINE_cellranger_premRNA/", pattern = "*metadata.txt", recursive = TRUE) # create a list of raw count tables from all batches
## This function subsamples 5% of all high quality and low quality cells, using post-filtered metadata (Jin Feng) as a reference for HQ, adds quality annotation then merge meta tables 
meta.filter.list <- lapply(paste('raw/FELINE_cellranger_premRNA/', all.meta, sep=''), FUN = function(x) {
meta.data <- as.data.frame(fread(x, sep='auto'))
# subset cells by quality, subsample and annotate metadata
meta.data_HQ <- meta.data[ meta.data$Cell %in% filtered_HQ_cellID,]
meta_HQ_sub <- meta.data_HQ[sample(nrow(meta.data_HQ), round(0.1 * nrow(meta.data_HQ))), ] # subsample 10% cells 
# meta_HQ_sub['Quality'] <- rep('high', nrow(meta_HQ_sub) )

# meta.data_LQ <- meta.data[! meta.data$Cell %in% filtered_HQ_cellID,]
# meta_LQ_sub <- meta.data_LQ[sample(nrow(meta.data_LQ), round(0.1 * nrow(meta.data_LQ))), ] # subsample 10% cells
# meta_LQ_sub['Quality'] <- rep('low', nrow(meta_LQ_sub) )

# meta_subsample <- rbind(meta_HQ_sub, meta_LQ_sub)
return(meta_HQ_sub)
})

meta_batches_merged <- do.call(rbind, meta.filter.list)
meta_batches_merged <- meta_batches_merged[!duplicated(meta_batches_merged$Cell),]
rownames(meta_batches_merged) <- meta_batches_merged$Cell

write.table(meta_batches_merged, 'metadata_subsample_HQ.txt', sep='\t') # Export subsampled metadata



## Subsample count table for matching cell ids 
## Takes forever to run in the notebook. Run in standard R or Radian
## creates and index of subsamples cells in metadata for each raw count table then subsample each table only for matching cells
all.counts <- list.files(path = "/Users/ibishara/Desktop/FELINE_C1/raw/FELINE_cellranger_premRNA/", pattern = "*counts.txt", recursive = TRUE) # create a list of raw count tables from all batches
count.filter.list <- lapply(paste('raw/FELINE_cellranger_premRNA/', all.counts, sep=''), FUN = function(x) {
count_sub <- fread(x, sep='auto', select= c("Gene Symbol",  meta_batches_merged$Cell  )) # cell ids subsampled from the metadata
return(count_sub)  # returns warnings for every cell id NOT in each count table
})

count_batches_merged <- do.call(cbind, count.filter.list)  #join list of tables 
count_batches_merged <- as.data.frame(count_batches_merged)
count_batches_merged <- count_batches_merged[,!duplicated(colnames(count_batches_merged))] # remove duplicate cell id
count_batches_merged <- count_batches_merged[!duplicated(count_batches_merged$"Gene Symbol"),] # remove duplicate genes
count_batches_merged <- count_batches_merged[count_batches_merged$'Gene Symbol' != '' ,] # remove null genes

fwrite(count_batches_merged, 'raw_counts_subsample_HQ.txt', sep='\t') # Export subsampled raw counts