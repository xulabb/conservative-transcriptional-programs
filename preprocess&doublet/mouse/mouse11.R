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

setwd('~/Documents/Doctor/lab/develoment_of_early_ge/split analysis/mouse')

mou11.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mou11_filtered_feature_bc_matrix")
mou11 <- CreateSeuratObject(counts = mou11.data, 
                                  project = 'mou11', 
                                  min.cells = 13, 
                                  min.features = 800)
gc()
mou11
##
mou11[["percent.mt"]] <- PercentageFeatureSet(mou11, features = "^mt")
VlnPlot(mou11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mou11, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mou11, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
mou11<- subset(mou11, 
                     subset = nFeature_RNA > 800 & 
                       nFeature_RNA < 6000 & 
                       percent.mt < 20)

mou11 <- NormalizeData(mou11,
                             normalization.method = "LogNormalize", 
                             scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
mou11<- FindVariableFeatures(mou11, selection.method = "vst", nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mou11), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mou11)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(mou11)
mou11 <- ScaleData(mou11, features = all.genes)

#Perform linear dimensional reduction
mou11 <- RunPCA(mou11, features = VariableFeatures(object = mou11))

#Examine and visualize PCA results a few different ways
print(mou11[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(mou11, dims = 1:2, reduction = "pca")
DimPlot(mou11, reduction = "pca")
DimHeatmap(mou11, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
mou11 <- JackStraw(mou11, num.replicate = 100)
mou11 <- ScoreJackStraw(mou11, dims = 1:50)
JackStrawPlot(mou11, dims = 1:40)
ElbowPlot(object = mou11, ndims = 40, reduction = "pca") 

#FindNeighbors
mou11 <- FindNeighbors(mou11, 
                             reduction = "pca", 
                             dims = 1:11)
mou11 <- FindClusters(mou11, 
                            resolution = 1.7)

#Look at cluster IDs of the first 5 cells
head(Idents(mou11), 5)

#deDoulblets 
sweep.mou11 <- paramSweep_v3(mou11, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.mou11 , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(mou11@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(mou11@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
mou11 <- doubletFinder_v3(mou11, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###mou11 <- doubletFinder_v3(mou11, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
mou11@meta.data[,"DF_hi.lo"] <- mou11@meta.data$DF.classifications_0.25_0.09_643
mou11@meta.data$DF_hi.lo[which(mou11@meta.data$DF_hi.lo == "Doublet" & mou11@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(mou11, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
mou11_1<-subset(mou11,subset=DF_hi.lo=="Singlet")

#UMAP
mou11 <- RunUMAP(mou11, 
                       dims = 1:30, 
                       label = T)
head(mou11@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(mou11, reduction = "umap")
LabelClusters(DimPlot(mou11, reduction = "umap"),id = 'ident')
DimPlot(mou11, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
mou11 <- RunTSNE(mou11, dims = 1:30)
head(mou11@reductions$tsne@cell.embeddings)
p2 <- DimPlot(mou11, reduction = "tsne")
LabelClusters(DimPlot(mou11, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(mou11, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(mou11, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
mou11.markers <- FindAllMarkers(mou11, only.pos = TRUE, min.pct = 0.25)
mou11.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(mou11.markers,file="allmarkers_mou11.csv")

#DoHeatmap
top10 <- mou11.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(mou11, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(mou11, c("Ascl1","Hes5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(mou11, c(" Nkx2-1","Lhx6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(mou11, c(" Mesi2", "Isl1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(mou11, c("Calb2"," Nr2f2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(mou11, c("Pdgfra")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(mou11, c("Cx3cr1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(mou11, c("Col3a1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(mou11,file="mou11_.Rdata")
