#Create a folder in your Documents folder named: "scRNA-workshop"
#Download the raw data, unzip it and put it in that folder
#https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
#This tutorial is based on https://satijalab.org/seurat/articles/pbmc3k_tutorial
#The order of certain analysis was changed for didactic purposes
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

###### Loading and processing a Seurat object #####
#set work directory
setwd("C:/Users/dzsid/Documents/scRNA-workshop/")

#load dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# create seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3", min.cells = 3,
                           min.features = 200)

#prepare for PCA by normalization, selecting variable features and scaling
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#First 5 principal components
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#Plot PCA
DimPlot(pbmc, reduction = "pca") + NoLegend()

ElbowPlot(pbmc)

# Cluster samples
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Run and plot UMAP (non-linear dimensionality reduction method)
pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap")
################################################################################

##### Quality control #####
#calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#plot number of genes, number of reads and mitochondrial percentage
VlnPlot(pbmc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

FeaturePlot(pbmc, "percent.mt", max.cutoff = 7,)
FeaturePlot(pbmc, "nFeature_RNA")

# Create a filtered dataset
pbmc_filt <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3", 
                                min.cells = 3, min.features = 200)
pbmc_filt[["percent.mt"]] <- PercentageFeatureSet(pbmc_filt, pattern = "^MT-")

pbmc_filt <- subset(pbmc_filt, subset = nFeature_RNA > 200 & 
                      nFeature_RNA < 2500 & percent.mt < 5)

# Repeat the same preprocessing with the filtered object
pbmc_filt <- NormalizeData(pbmc_filt, normalization.method = "LogNormalize",
                           scale.factor = 10000)
pbmc_filt <- FindVariableFeatures(pbmc_filt, selection.method = "vst",
                                  nfeatures = 2000)
all.genes <- rownames(pbmc_filt)
pbmc_filt <- ScaleData(pbmc_filt, features = all.genes)
pbmc_filt <- RunPCA(pbmc_filt, features = VariableFeatures(object = pbmc_filt))
pbmc_filt <- FindNeighbors(pbmc_filt, dims = 1:10)
pbmc_filt <- FindClusters(pbmc_filt, resolution = 0.5)
pbmc_filt <- RunUMAP(pbmc_filt, dims = 1:10)

VlnPlot(pbmc_filt, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

#Plot UMAP
DimPlot(pbmc_filt, reduction = "umap")

#Check mitochondrial percentage
FeaturePlot(pbmc_filt, "percent.mt")
FeaturePlot(pbmc_filt, "nFeature_RNA")
FeaturePlot(pbmc_filt, "nCount_RNA")
################################################################################

##### Differential gene expression #####
#To find marker genes for cluster2
cluster2.markers <- FindMarkers(pbmc_filt, ident.1 = 2)
head(cluster2.markers, n = 5)

#Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc_filt, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#Take the top10 makers in each cluster
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc_filt, features = top10$gene) + NoLegend()

##### Cell typing #####
#Show selected genes on a UMAP
FeaturePlot(pbmc_filt, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
                               "FCGR3A", "LYZ", "PPBP", "CD8A"))
DimPlot(pbmc_filt, reduction = "umap", label = T)

#Set cell-types based on https://satijalab.org/seurat/articles/pbmc3k_tutorial
#Normally, you would run a cell-typing algorithm such as SingleR (2019-2022) or nowadays Azimuth
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc_filt)
pbmc_filt <- RenameIdents(pbmc_filt, new.cluster.ids)

#Plot the cell-types on a UMAP
DimPlot(pbmc_filt, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(pbmc_filt, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + NoAxes()
#Plot them on a PCA
DimPlot(pbmc_filt, reduction = "pca", label = TRUE) + NoLegend()

##### Data structure #####
#see how the count matrix is stored
View(pbmc_filt@assays$RNA$data)

#View metadata
View(pbmc_filt@meta.data)

#Create a "cell_type" column in the metadata
pbmc_filt@meta.data$cell_type <- Idents(pbmc_filt)

#Count each cell type
table(pbmc_filt@meta.data$cell_type)

################################################################################

##### Data visualisation #####
#Visualize cell type distribution
p2 <- ggplot() + 
  geom_bar(data = pbmc_filt@meta.data, aes(x = orig.ident, fill = cell_type),
           position = "fill") +  
  theme_classic() + 
  xlab("sample") + 
  ylab("percentage")
p2

#Visualize single genes
VlnPlot(pbmc_filt, features = c("PDCD1"))
FeaturePlot(pbmc_filt, "PDCD1")

# Find differentially expressed genes between Naive and Memory CD4+ cells
T4markers <- FindMarkers(pbmc_filt, ident.1 = "Naive CD4 T", 
                         ident.2 = "Memory CD4 T")

# Prepare a dataframe for a volcano plot
T4markers$diffexpressed <- "NO"
T4markers$diffexpressed[T4markers$avg_log2FC>2 & T4markers$p_val<0.001]<-"Naive"
T4markers$diffexpressed[T4markers$avg_log2FC<(-2)&T4markers$p_val<0.001]<-"Memory"
T4markers$direction <- NA
T4markers$direction[T4markers$diffexpressed != "NO"] <- rownames(T4markers[
  T4markers$diffexpressed != "NO",])

#Plot a volcano plot
ggplot(data=T4markers, aes(x=avg_log2FC, y=-log10(p_val), 
                           col=diffexpressed, label=direction)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  ggtitle("Differentially expressed genes between CD4+ T-cells") +
  scale_color_manual(values=c("blue", "red", "black")) +
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept=-log10(0.001), col="red")

################################################################################
