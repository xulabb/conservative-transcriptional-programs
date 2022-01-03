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

setwd('~/Documents/Doctor/lab/develoment_of_early_ge/split analysis/mouse1')

mou12.5_LGE.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mou12.5_LGE_filtered_feature_bc_matrix")
mou12.5_LGE <- CreateSeuratObject(counts = mou12.5_LGE.data, 
                                  project = 'mou12.5_LGE', 
                                  min.cells = 13, 
                                  min.features = 800)
mou12.5_LGE <- NormalizeData(mou12.5_LGE,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
mou12.5_LGE<- FindVariableFeatures(mou12.5_LGE, 
                             selection.method = "vst", 
                             nfeatures = 3000)

mou12.5_MGE.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mou12.5_MGE_filtered_feature_bc_matrix")
mou12.5_MGE <- CreateSeuratObject(counts = mou12.5_MGE.data, 
                                  project = 'mou12.5_MGE', 
                                  min.cells = 13, 
                                  min.features = 800)
mou12.5_MGE <- NormalizeData(mou12.5_MGE,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
mou12.5_MGE<- FindVariableFeatures(mou12.5_MGE, 
                             selection.method = "vst", 
                             nfeatures = 3000)

gc()
mou12.5_MGE
mou12.5_LGE

##Integrate
mouse1.anc <- FindIntegrationAnchors(object.list = list(mou12.5_MGE,mou12.5_LGE), 
                                           anchor.features = 2000, 
                                           scale = T, 
                                           reduction = 'cca', 
                                           dims = 1:30, 
                                           k.filter = 200)
mouse1 <- IntegrateData(anchorset = mouse1.anc, dims = 1:30)
mouse1
gc()

##
mouse1[["percent.mt"]] <- PercentageFeatureSet(mouse1, features = "^mt")
VlnPlot(mouse1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mouse1, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse1, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
mouse1<- subset(mouse1, 
                     subset = nFeature_RNA > 800 & 
                       nFeature_RNA < 6000 & 
                       percent.mt < 20)

mouse1 <- NormalizeData(mouse1,
                             normalization.method = "LogNormalize", 
                             scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
mouse1<- FindVariableFeatures(mouse1, selection.method = "vst", nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mouse1), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mouse1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(mouse1)
mouse1 <- ScaleData(mouse1, features = all.genes)

#Perform linear dimensional reduction
mouse1 <- RunPCA(mouse1, 
                features = VariableFeatures(object = mouse1),
                npcs=50)

#Examine and visualize PCA results a few different ways
print(mouse1[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(mouse1, dims = 1:2, reduction = "pca")
DimPlot(mouse1, reduction = "pca")
DimHeatmap(mouse1, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
mouse1 <- JackStraw(mouse1, num.replicate = 100)
mouse1 <- ScoreJackStraw(mouse1, dims = 1:50)
JackStrawPlot(mouse1, dims = 1:40)
ElbowPlot(object = mouse1, ndims = 40, reduction = "pca") 

#FindNeighbors
mouse1 <- FindNeighbors(mouse1, 
                             reduction = "pca", 
                             dims = 1:11)
mouse1 <- FindClusters(mouse1, 
                            resolution = 0.5)

#Look at cluster IDs of the first 5 cells
head(Idents(mouse1), 5)

#deDoulblets 
sweep.mouse1 <- paramSweep_v3(mouse1, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.mouse1 , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(mouse1@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(mouse1@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
mouse1 <- doubletFinder_v3(mouse1, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###mouse1 <- doubletFinder_v3(mouse1, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
mouse1@meta.data[,"DF_hi.lo"] <- mouse1@meta.data$DF.classifications_0.25_0.09_643
mouse1@meta.data$DF_hi.lo[which(mouse1@meta.data$DF_hi.lo == "Doublet" & mouse1@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(mouse1, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
mouse1_1<-subset(mouse1,subset=DF_hi.lo=="Singlet")

#UMAP
mouse1 <- RunUMAP(mouse1, 
                       dims = 1:30, 
                       label = T)
head(mouse1@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(mouse1, reduction = "umap")
LabelClusters(DimPlot(mouse1, reduction = "umap"),id = 'ident')
DimPlot(mouse1, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
mouse1 <- RunTSNE(mouse1, dims = 1:30)
head(mouse1@reductions$tsne@cell.embeddings)
p2 <- DimPlot(mouse1, reduction = "tsne")
LabelClusters(DimPlot(mouse1, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(mouse1, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(mouse1, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
mouse1.markers <- FindAllMarkers(mouse1, only.pos = TRUE, min.pct = 0.25)
mouse1.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(mouse1.markers,file="allmarkers_mouse1.csv")

#DoHeatmap
top10 <- mouse1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(mouse1, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(mouse1, c("Ascl1","Hes5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(mouse1, c(" Nkx2-1","Lhx6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(mouse1, c(" Mesi2", "Isl1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(mouse1, c("Calb2"," Nr2f2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(mouse1, c("Pdgfra")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(mouse1, c("Cx3cr1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(mouse1, c("Col3a1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(mouse1,file="mouse1_.Rdata")
