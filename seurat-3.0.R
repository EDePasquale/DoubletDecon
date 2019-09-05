library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(base)
library(methods)
library(utils)
library(stats)
library(gdata)
library(graphics)
library(grDevices)

args<- commandArgs(trailingOnly=TRUE)
path <- args[1]
name <- args[2]
# Load the dataset
setwd(path)
M.data <- Read10X(data.dir = path)
#Examine the memory savings between regular and sparse matrices
dense.size <- object.size(as.matrix(M.data))
dense.size
sparse.size <- object.size(M.data,sparseMatrixClass='Matrix')
sparse.size
dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data)
# Note that this is slightly different than the older Seurat workflow, where log-normalized values were passed in directly.
# You can continue to pass in log-normalized values, just set do.logNormalize=F in the next step.
#M <- new ("seurat", raw.data = M.data)

# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
# Initialize the Seurat object with the raw (non-normalized data).
M <- CreateSeuratObject(counts = M.data, project = name, min.cells = 3, min.features = 200)
M

#nGene and nUMI are automatically calculated for every object by Seurat. For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
M[["percent.mt"]] <- PercentageFeatureSet(object = M, pattern = "^MT-")
pdf("QCstats.pdf", width = 12, height = 15)
VlnPlot(object = M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(object = M, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = M, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("Features-plot.pdf", width = 15, height = 15)
CombinePlots(plots = list(plot1, plot2))
dev.off()
#We filter out cells that have unique gene counts over 2,500
#Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nGene but on anything in the object (as in GenePlot above)
# M <- SubsetData(M, subset.name = "nGene", accept.high = 4500)
# M <- SubsetData(M, subset.name = "percent.mito", accept.high = 0.6)

M <- subset(x = M, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt <40)
###Normalizing the data
###After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
M <- NormalizeData(object = M, normalization.method = "LogNormalize", scale.factor = 10000)

# ## Detection of variable genes across the single cells
# #M <- MeanVarPlot(M,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.5, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
# M <- FindVariableGenes(object = M, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
# length(M@var.genes)


# Identification of highly variable features (feature selection)
# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

M <- FindVariableFeatures(object = M, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = M), 10)
write.table(top10,file="top10-high-variable-genes.txt",sep="\t",col.names= NA)

# plot variable features with and without labels
pdf("Variable-Feature-Plot.pdf", width = 13, height = 15)
plot1 <- VariableFeaturePlot(object = M)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


##Scaling the data and removing unwanted sources of variation
all.genes <- rownames(x = M)
M <- ScaleData(object = M, features = all.genes)

## Perform linear dimensional reduction
M <- RunPCA(object = M, features = VariableFeatures(object = M))
M

# Examine  and visualize PCA results a few different ways
pdf("PCA.pdf", width = 8, height = 15)
VizDimLoadings(object = M, dims = 1:2, reduction = "pca")
dev.off()

pdf("PCA-plot.pdf", width = 8, height = 15)
DimPlot(object = M, reduction = "pca")
dev.off()

pdf("PC-Heatmap.pdf", width = 8, height = 15)
DimHeatmap(object = M, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = M, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


####Determine statistically significant principal components
M <- JackStraw(object = M, num.replicate = 100)
M <- ScoreJackStraw(object = M, dims = 1:20)
png("JackStraw.png", width = 8, height = 15, units = 'in', res = 600)
pdf("JackStraw.pdf", width = 8, height = 15)
JackStrawPlot(object = M, dims = 1:20)
dev.off()

png("PCE-bowl-plot.png", width = 8, height = 15, units = 'in', res = 600)
pdf("PCE-bowl-plot.pdf", width = 8, height = 15)
ElbowPlot(object = M)
dev.off()

## Cluster the cells
#############GIVE THE NUMBER OF PCS TO USE FROM JACKSTRAW AND ELBOWPLOT############### IN THE CLUSTERS AND UMAP OPTION##################
#FindNeighbors function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).
M <- FindNeighbors(object = M, dims = 1:20)
M <- FindClusters(object = M, resolution = 0.5)
head(x = Idents(object = M), 5)

####Run non-linear dimensional reduction (UMAP/tSNE)####
pdf("UMAP-plot.pdf", width = 8, height = 15)
M <- RunUMAP(object = M, dims = 1:20)
DimPlot(object = M, reduction = "umap",label= TRUE)
dev.off()

saveRDS(M, file = "Seurat.rds")


## Finding differentially expressed features (cluster biomarkers)

M.markers <- FindAllMarkers(object = M, only.pos = TRUE, min.pct=0.25, logfc.threshold = 0.25)
#M.markers %>% group_by(cluster) %>% top_n(n=2,avg_logFC)

cluster1.markers<- FindMarkers(object = M, ident.1 =0, thresh.use = 0.25, test.use = "roc", only.pos=TRUE)

top50 <- M.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
###color manipulration of scale DoHeatmap(object = pbmc_small) + scale_fill_gradientn(colors = c("blue", "white", "red"))#######
pdf("Heatmap-clusters.pdf", width = 8, height = 15)
DoHeatmap(object = M, features = top50$gene,label=TRUE,size=3, hjust=3) + NoLegend() + FontSize(y.text=3)
dev.off()
## Heatmap with Legend##
pdf("Heatmap-clusters-Legend.pdf", width = 8, height = 15)
DoHeatmap(object = M, features = top50$gene,label=TRUE,size=3, hjust=3) + FontSize(y.text=3)
dev.off()

write.table(top50,file="Top50Genes.txt",sep="\t",col.names= NA)
#write.table(M@ident,"Cluster.txt",sep="\t",col.names= NA)
write.table(x = Idents(object = M),"Cluster.txt",sep="\t",col.names= NA)
dev.off()

