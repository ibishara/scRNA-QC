

# packages 
library(Seurat)
library(qs)
library(singleCellNet)
library(ggpubr)
library(data.table)
library(stringr)
setwd('/Users/ibishara/Desktop/FELINE_C1/')

# data
markers <- read.table('Annotation.lineage.latentCC.markers.txt', sep = '\t' )

# High quality FELINE C1 data
seu_HQ <- qread(file="seu_HQ.qs", nthreads=16)
# seu_HQ <- subset(x = seu_HQ, subset = Celltype != "Normal epithelial cells")   ## Removes normal epithelial cells 
meta <- seu_HQ@meta.data
seu.HQ.counts <- GetAssayData(seu_HQ, assay = "RNA")

# Genes
counts <- as.data.frame(fread( 'reads_downsample/reads_down_1.5.txt', sep='\t')) # downsamples genes 
rownames(counts) <- counts$V1
counts$V1 <- NULL


total.mito.reads <- colSums(counts[grepl('MT-', rownames(counts)), ])  # mitochondrial reads per cell
total.non.mito.reads <- colSums(counts[!grepl('MT-', rownames(counts)),])
total.reads <- colSums(counts)  # mitochondrial reads per cell
perc.mito <-  total.mito.reads / total.reads
hist( perc.mito)



mito.pos <- which(grepl('MT-', rownames(counts)))  # mitochondrial genes position
non.mito.pos <- which(!grepl('MT-', rownames(counts)))  # non-mitochondrial genes position

# Function to reduce the number of counts by a fraction 
foo = function(x, y){
    # x <- counts[,2]
    # y <- 0.4
    total.mito.reads <- sum(x[mito.pos])  # mitochondrial reads per cell
    total.non.mito.reads <- sum(x[non.mito.pos]) 
    total.reads = sum(x)
    perc.mito <-  total.mito.reads / total.reads # calculate cell's mitochondrial percentage
# perc.mito
    if(perc.mito < y & perc.mito > 0 ){
        ## approach # 1 
        # x[non.mito.pos] <- x[non.mito.pos] * (1 - y)/ (1 - perc.mito) 
        
        ## approach # 2
        new.total.reads <- (perc.mito * total.reads)/y
        new.total.non.mito.reads <- new.total.reads - total.mito.reads
        x[non.mito.pos] <- x[non.mito.pos]   * new.total.non.mito.reads/ total.reads    
    }
    return(x)
}

mito_down_0.4  <- as.data.frame(apply(counts, MARGIN=2, FUN=foo, y = 0.4))


total.mito.reads <- colSums(mito_down_0.4[mito.pos, ])  # mitochondrial reads per cell
total.non.mito.reads <- colSums(mito_down_0.4[non.mito.pos,])
total.reads <- colSums(mito_down_0.4)  # mitochondrial reads per cell
perc.mito <-  total.mito.reads / total.reads
hist(perc.mito)
hist(total.reads)



# model validation
stTestList = splitCommon(sampTab=meta, ncells= ncells, dLevel= var2) #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = mito_down_0.4 [ ,rownames(stTest)]

classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, 
                                dLevelSID = "Cell", classTrain = var2, classQuery = var2, nRand = 50)
plot_PRs(tm_heldoutassessment)










# classifier 
var2 <- 'Lineage'
nGenes <- 25
ncells <- 400
   
#Split for training and assessment, and transform training data
set.seed(100) 
stList = splitCommon(sampTab = meta, ncells = ncells, dLevel = var2)
stTrain = stList[[1]]
expTrain = seu.HQ.counts[, rownames(stTrain)]

# model training 
system.time(class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, 
                        nTopGenes = nGenes, nRand = 100, nTrees = 1000, nTopGenePairs = 2*nGenes, 
                        dLevel = var2, colName_samp = "Cell"))

# model validation
stTestList = splitCommon(sampTab=meta, ncells= ncells, dLevel= var2) #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = mito_down_0.4 [ ,rownames(stTest)]

classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, 
                                dLevelSID = "Cell", classTrain = var2, classQuery = var2, nRand = 50)
plot_PRs(tm_heldoutassessment)




plot_attr(classRes = classRes_val_all, sampTab=stTest, nrand=50, dLevel=var2, sid="Cell")
plot_metrics(tm_heldoutassessment)






x <- seu.HQ.counts[-grep('LINC', rownames(seu.HQ.counts)),]
x[x > 0] <- 1 
y <- as.data.frame(cbind(meta,colSums(x)))
names(y)[names(y) == "colSums(x)"] <- "LincRNA_genes" 
names(y)[names(y) == "LincRNA_genes"] <- "non_LincRNA_genes" 



hist(colSums(x), breaks = 50) # number of lncRNA expressed/captured 

plot(colSums(x), log10(y[rownames(y) %in% colnames(x), ]$nCount_RNA) , xlab= 'LincRNA genes / cell', ylab = 'nCount_RNA')
# hist(colSums(x)/y[rownames(y) %in% colnames(x), ]$nCount_RNA, breaks = 50)


# cor(colSums(x), log10(y[rownames(y) %in% colnames(x), ]$nCount_RNA ))



ggplot(y, aes(non_LincRNA_genes, log10(nCount_RNA), fill = )) + 
geom_point (aes( color = orig.ident))


colSums(seu.HQ.counts[!grepl('LINC', rownames(seu.HQ.counts)),]) # non.mito counts/cell

colSums(seu.HQ.counts) # total counts / cell 

hist(colSums(seu.HQ.counts[grepl('MT-', rownames(seu.HQ.counts)),]) / colSums(seu.HQ.counts) )



# filter for genes identified from FindAllMarkers
common.genes <- intersect(rownames(seu.HQ.counts), markers$gene)
seu.HQ.counts <- seu.HQ.counts[common.genes, ]

#Split for training and assessment, and transform training data
set.seed(100) 
stList = splitCommon(sampTab=meta, ncells=400, dLevel="Lineage")
stTrain = stList[[1]]
expTrain = seu.HQ.counts[, rownames(stTrain)]




# model training 

system.time(class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, 
                            nTopGenes = 25, nRand = 100, nTrees = 1000, nTopGenePairs = 50, 
                            dLevel = "Lineage", colName_samp = "Cell"))



# model validation
stTestList = splitCommon(sampTab=meta, ncells=400, dLevel="Lineage") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = counts[ ,rownames(stTest)]




#predict
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, 
                                   dLevelSID = "Cell", classTrain = "Lineage", classQuery = "Lineage", nRand = 50)


plot_PRs(tm_heldoutassessment)

plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=50, dLevel="Lineage", sid="Cell")


plot_metrics(tm_heldoutassessment)