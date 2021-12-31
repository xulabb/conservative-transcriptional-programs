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

setwd('~/Documents/Doctor/lab/develoment_of_early_ge/split analysis/human')

hsa13.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/hsa13_filtered_feature_bc_matrix")
hsa13 <- CreateSeuratObject(counts = hsa13.data, 
                           project = 'hsa13', 
                           min.cells = 13, 
                           min.features = 800)
gc()
hsa13
##
hsa13[["percent.mt"]] <- PercentageFeatureSet(hsa13, features = "^MT")
VlnPlot(hsa13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(hsa13, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(hsa13, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
hsa13<- subset(hsa13, 
              subset = nFeature_RNA > 800 & 
                nFeature_RNA < 6000 & 
                percent.mt < 20)

hsa13 <- NormalizeData(hsa13,
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
hsa13<- FindVariableFeatures(hsa13, selection.method = "vst", nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hsa13), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hsa13)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(hsa13)
hsa13 <- ScaleData(hsa13, features = all.genes)

#Perform linear dimensional reduction
hsa13 <- RunPCA(hsa13, features = VariableFeatures(object = hsa13))

#Examine and visualize PCA results a few different ways
print(hsa13[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(hsa13, dims = 1:2, reduction = "pca")
DimPlot(hsa13, reduction = "pca")
DimHeatmap(hsa13, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
hsa13 <- JackStraw(hsa13, num.replicate = 100)
hsa13 <- ScoreJackStraw(hsa13, dims = 1:50)
JackStrawPlot(hsa13, dims = 1:40)
ElbowPlot(object = hsa13, ndims = 40, reduction = "pca") 

#FindNeighbors
hsa13 <- FindNeighbors(hsa13, 
                      reduction = "pca", 
                      dims = 1:11)
hsa13 <- FindClusters(hsa13, 
                     resolution = 1.7)

#Look at cluster IDs of the first 5 cells
head(Idents(hsa13), 5)

#deDoulblets 
sweep.hsa13 <- paramSweep_v3(hsa13, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.hsa13 , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(hsa13@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(hsa13@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
hsa13 <- doubletFinder_v3(hsa13, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###hsa13 <- doubletFinder_v3(hsa13, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
hsa13@meta.data[,"DF_hi.lo"] <- hsa13@meta.data$DF.classifications_0.25_0.09_643
hsa13@meta.data$DF_hi.lo[which(hsa13@meta.data$DF_hi.lo == "Doublet" & hsa13@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(hsa13, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
hsa13_1<-subset(hsa13,subset=DF_hi.lo=="Singlet")

#UMAP
hsa13 <- RunUMAP(hsa13, 
                dims = 1:30, 
                label = T)
head(hsa13@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(hsa13, reduction = "umap")
LabelClusters(DimPlot(hsa13, reduction = "umap"),id = 'ident')
DimPlot(hsa13, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
hsa13 <- RunTSNE(hsa13, dims = 1:30)
head(hsa13@reductions$tsne@cell.embeddings)
p2 <- DimPlot(hsa13, reduction = "tsne")
LabelClusters(DimPlot(hsa13, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(hsa13, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(hsa13, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
hsa13.markers <- FindAllMarkers(hsa13, only.pos = TRUE, min.pct = 0.25)
hsa13.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(hsa13.markers,file="allmarkers_hsa13.csv")

#DoHeatmap
top10 <- hsa13.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(hsa13, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(hsa13, c("ASCL1","HES5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(hsa13, c(" NKX2-1","LHX6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(hsa13, c(" MEIS2", "ISL1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(hsa13, c("CALB2"," NR2F2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(hsa13, c("PDGFRA")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(hsa13, c("CX3CR1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(hsa13, c("COL3A1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(hsa13,file="hsa13_.Rdata")
