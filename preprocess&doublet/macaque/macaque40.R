library(Seurat)
library(magrittr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(metap)
library(biomaRt)
library(hrbrthemes)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(slingshot)
library(scater)
library(ggbreak)
library(stringr)
library(sva)
library(ggrepel)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(ggpubr)
library(reshape2)
library(knitr)
library(RColorBrewer)
suppressMessages(library(gplots))
library(amap) 
library(BiocParallel)
library(tximport)
library(readr)
library(BuenColors)
rm(list = ls())

setwd('~/Documents/Doctor/lab/develoment_of_early_ge/split analysis/macaque')

mac40.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mac40_filtered_feature_bc_matrix")
mac40 <- CreateSeuratObject(counts = mac40.data, 
                           project = 'mac40', 
                           min.cells = 13, 
                           min.features = 800)
gc()
mac40
##
mac40[["percent.mt"]] <- PercentageFeatureSet(mac40, features = mt.genes)
VlnPlot(mac40, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mac40, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mac40, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
mac40<- subset(mac40, 
               subset = nFeature_RNA > 800 & 
                 nFeature_RNA < 6000 & 
                 percent.mt < 20)

mac40 <- NormalizeData(mac40,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
mac40<- FindVariableFeatures(mac40, selection.method = "vst", nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mac40), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mac40)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(mac40)
mac40 <- ScaleData(mac40, features = all.genes)

#Perform linear dimensional reduction
mac40 <- RunPCA(mac40, features = VariableFeatures(object = mac40))

#Examine and visualize PCA results a few different ways
print(mac40[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(mac40, dims = 1:2, reduction = "pca")
DimPlot(mac40, reduction = "pca")
DimHeatmap(mac40, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
mac40 <- JackStraw(mac40, num.replicate = 100)
mac40 <- ScoreJackStraw(mac40, dims = 1:50)
JackStrawPlot(mac40, dims = 1:40)
ElbowPlot(object = mac40, ndims = 40, reduction = "pca") 

#FindNeighbors
mac40 <- FindNeighbors(mac40, 
                          reduction = "pca", 
                          dims = 1:11)
mac40 <- FindClusters(mac40, 
                         resolution = 1.7)

#Look at cluster IDs of the first 5 cells
head(Idents(mac40), 5)

#deDoulblets 
sweep.mac40 <- paramSweep_v3(mac40, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.mac40 , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(mac40@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(mac40@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
mac40 <- doubletFinder_v3(mac40, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###mac40 <- doubletFinder_v3(mac40, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
mac40@meta.data[,"DF_hi.lo"] <- mac40@meta.data$DF.classifications_0.25_0.09_643
mac40@meta.data$DF_hi.lo[which(mac40@meta.data$DF_hi.lo == "Doublet" & mac40@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(mac40, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
mac40_1<-subset(mac40,subset=DF_hi.lo=="Singlet")

#UMAP
mac40 <- RunUMAP(mac40, 
                 dims = 1:30, 
                 label = T)
head(mac40@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(mac40, reduction = "umap")
LabelClusters(DimPlot(mac40, reduction = "umap"),id = 'ident')
DimPlot(mac40, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
mac40 <- RunTSNE(mac40, dims = 1:30)
head(mac40@reductions$tsne@cell.embeddings)
p2 <- DimPlot(mac40, reduction = "tsne")
LabelClusters(DimPlot(mac40, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(mac40, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(mac40, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
mac40.markers <- FindAllMarkers(mac40, only.pos = TRUE, min.pct = 0.25)
mac40.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(mac40.markers,file="allmarkers_mac40.csv")

#DoHeatmap
top10 <- mac40.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(mac40, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(mac40, c("ASCL1","HES5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(mac40, c(" NKX2-1","LHX6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(mac40, c(" MEIS2", "ISL1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(mac40, c("CALB2"," NR2F2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(mac40, c("PDGFRA")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(mac40, c("CX3CR1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(mac40, c("COL3A1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(mac40,file="mac40_.Rdata")
