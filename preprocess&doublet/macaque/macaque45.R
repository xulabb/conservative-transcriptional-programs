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

mac45.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mac45_filtered_feature_bc_matrix")
mac45 <- CreateSeuratObject(counts = mac45.data, 
                            project = 'mac45', 
                            min.cells = 13, 
                            min.features = 800)
gc()
mac45
##
mac45[["percent.mt"]] <- PercentageFeatureSet(mac45, features = mt.genes)
VlnPlot(mac45, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mac45, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mac45, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
mac45<- subset(mac45, 
               subset = nFeature_RNA > 800 & 
                 nFeature_RNA < 6000 & 
                 percent.mt < 20)

mac45 <- NormalizeData(mac45,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
mac45<- FindVariableFeatures(mac45, selection.method = "vst", nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mac45), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mac45)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(mac45)
mac45 <- ScaleData(mac45, features = all.genes)

#Perform linear dimensional reduction
mac45 <- RunPCA(mac45, features = VariableFeatures(object = mac45))

#Examine and visualize PCA results a few different ways
print(mac45[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(mac45, dims = 1:2, reduction = "pca")
DimPlot(mac45, reduction = "pca")
DimHeatmap(mac45, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
mac45 <- JackStraw(mac45, num.replicate = 100)
mac45 <- ScoreJackStraw(mac45, dims = 1:50)
JackStrawPlot(mac45, dims = 1:40)
ElbowPlot(object = mac45, ndims = 40, reduction = "pca") 

#FindNeighbors
mac45 <- FindNeighbors(mac45, 
                       reduction = "pca", 
                       dims = 1:11)
mac45 <- FindClusters(mac45, 
                      resolution = 1.7)

#Look at cluster IDs of the first 5 cells
head(Idents(mac45), 5)

#deDoulblets 
sweep.mac45 <- paramSweep_v3(mac45, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.mac45 , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(mac45@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(mac45@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
mac45 <- doubletFinder_v3(mac45, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###mac45 <- doubletFinder_v3(mac45, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
mac45@meta.data[,"DF_hi.lo"] <- mac45@meta.data$DF.classifications_0.25_0.09_643
mac45@meta.data$DF_hi.lo[which(mac45@meta.data$DF_hi.lo == "Doublet" & mac45@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(mac45, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
mac45_1<-subset(mac45,subset=DF_hi.lo=="Singlet")

#UMAP
mac45 <- RunUMAP(mac45, 
                 dims = 1:30, 
                 label = T)
head(mac45@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(mac45, reduction = "umap")
LabelClusters(DimPlot(mac45, reduction = "umap"),id = 'ident')
DimPlot(mac45, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
mac45 <- RunTSNE(mac45, dims = 1:30)
head(mac45@reductions$tsne@cell.embeddings)
p2 <- DimPlot(mac45, reduction = "tsne")
LabelClusters(DimPlot(mac45, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(mac45, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(mac45, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
mac45.markers <- FindAllMarkers(mac45, only.pos = TRUE, min.pct = 0.25)
mac45.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(mac45.markers,file="allmarkers_mac45.csv")

#DoHeatmap
top10 <- mac45.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(mac45, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(mac45, c("ASCL1","HES5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(mac45, c(" NKX2-1","LHX6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(mac45, c(" MEIS2", "ISL1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(mac45, c("CALB2"," NR2F2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(mac45, c("PDGFRA")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(mac45, c("CX3CR1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(mac45, c("COL3A1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(mac45,file="mac45_.Rdata")
