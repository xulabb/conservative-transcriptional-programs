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

mou12.5_MGE.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mou12.5_MGE_filtered_feature_bc_matrix")
mou12.5_MGE <- CreateSeuratObject(counts = mou12.5_MGE.data, 
                                  project = 'mou12.5_MGE', 
                                  min.cells = 13, 
                                  min.features = 800)
gc()
mou12.5_MGE
##
mou12.5_MGE[["percent.mt"]] <- PercentageFeatureSet(mou12.5_MGE, features = "^mt")
VlnPlot(mou12.5_MGE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mou12.5_MGE, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mou12.5_MGE, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
mou12.5_MGE<- subset(mou12.5_MGE, 
                     subset = nFeature_RNA > 800 & 
                       nFeature_RNA < 6000 & 
                       percent.mt < 20)

mou12.5_MGE <- NormalizeData(mou12.5_MGE,
                             normalization.method = "LogNormalize", 
                             scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
mou12.5_MGE<- FindVariableFeatures(mou12.5_MGE, selection.method = "vst", nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mou12.5_MGE), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mou12.5_MGE)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(mou12.5_MGE)
mou12.5_MGE <- ScaleData(mou12.5_MGE, features = all.genes)

#Perform linear dimensional reduction
mou12.5_MGE <- RunPCA(mou12.5_MGE, features = VariableFeatures(object = mou12.5_MGE))

#Examine and visualize PCA results a few different ways
print(mou12.5_MGE[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(mou12.5_MGE, dims = 1:2, reduction = "pca")
DimPlot(mou12.5_MGE, reduction = "pca")
DimHeatmap(mou12.5_MGE, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
mou12.5_MGE <- JackStraw(mou12.5_MGE, num.replicate = 100)
mou12.5_MGE <- ScoreJackStraw(mou12.5_MGE, dims = 1:50)
JackStrawPlot(mou12.5_MGE, dims = 1:40)
ElbowPlot(object = mou12.5_MGE, ndims = 40, reduction = "pca") 

#FindNeighbors
mou12.5_MGE <- FindNeighbors(mou12.5_MGE, 
                             reduction = "pca", 
                             dims = 1:11)
mou12.5_MGE <- FindClusters(mou12.5_MGE, 
                            resolution = 1.7)

#Look at cluster IDs of the first 5 cells
head(Idents(mou12.5_MGE), 5)

#deDoulblets 
sweep.mou12.5_MGE <- paramSweep_v3(mou12.5_MGE, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.mou12.5_MGE , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(mou12.5_MGE@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(mou12.5_MGE@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
mou12.5_MGE <- doubletFinder_v3(mou12.5_MGE, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###mou12.5_MGE <- doubletFinder_v3(mou12.5_MGE, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
mou12.5_MGE@meta.data[,"DF_hi.lo"] <- mou12.5_MGE@meta.data$DF.classifications_0.25_0.09_643
mou12.5_MGE@meta.data$DF_hi.lo[which(mou12.5_MGE@meta.data$DF_hi.lo == "Doublet" & mou12.5_MGE@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(mou12.5_MGE, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
mou12.5_MGE_1<-subset(mou12.5_MGE,subset=DF_hi.lo=="Singlet")

#UMAP
mou12.5_MGE <- RunUMAP(mou12.5_MGE, 
                       dims = 1:30, 
                       label = T)
head(mou12.5_MGE@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(mou12.5_MGE, reduction = "umap")
LabelClusters(DimPlot(mou12.5_MGE, reduction = "umap"),id = 'ident')
DimPlot(mou12.5_MGE, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
mou12.5_MGE <- RunTSNE(mou12.5_MGE, dims = 1:30)
head(mou12.5_MGE@reductions$tsne@cell.embeddings)
p2 <- DimPlot(mou12.5_MGE, reduction = "tsne")
LabelClusters(DimPlot(mou12.5_MGE, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(mou12.5_MGE, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(mou12.5_MGE, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
mou12.5_MGE.markers <- FindAllMarkers(mou12.5_MGE, only.pos = TRUE, min.pct = 0.25)
mou12.5_MGE.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(mou12.5_MGE.markers,file="allmarkers_mou12.5_MGE.csv")

#DoHeatmap
top10 <- mou12.5_MGE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(mou12.5_MGE, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(mou12.5_LGE, c("Ascl1","Hes5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(mou12.5_LGE, c(" Nkx2-1","Lhx6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(mou12.5_LGE, c(" Mesi2", "Isl1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(mou12.5_LGE, c("Calb2"," Nr2f2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(mou12.5_LGE, c("Pdgfra")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(mou12.5_LGE, c("Cx3cr1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(mou12.5_LGE, c("Col3a1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(mou12.5_MGE,file="mou12.5_MGE_.Rdata")
