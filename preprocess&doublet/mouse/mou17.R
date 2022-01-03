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

mou17.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mou17_filtered_feature_bc_matrix")
mou17 <- CreateSeuratObject(counts = mou17.data, 
                            project = 'mou17', 
                            min.cells = 13, 
                            min.features = 800)
gc()
mou17
##
mou17[["percent.mt"]] <- PercentageFeatureSet(mou17, features = "^mt")
VlnPlot(mou17, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mou17, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mou17, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
mou17<- subset(mou17, 
               subset = nFeature_RNA > 800 & 
                 nFeature_RNA < 6000 & 
                 percent.mt < 20)

mou17 <- NormalizeData(mou17,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
mou17<- FindVariableFeatures(mou17, selection.method = "vst", nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mou17), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mou17)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(mou17)
mou17 <- ScaleData(mou17, features = all.genes)

#Perform linear dimensional reduction
mou17 <- RunPCA(mou17, features = VariableFeatures(object = mou17))

#Examine and visualize PCA results a few different ways
print(mou17[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(mou17, dims = 1:2, reduction = "pca")
DimPlot(mou17, reduction = "pca")
DimHeatmap(mou17, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
mou17 <- JackStraw(mou17, num.replicate = 100)
mou17 <- ScoreJackStraw(mou17, dims = 1:50)
JackStrawPlot(mou17, dims = 1:40)
ElbowPlot(object = mou17, ndims = 40, reduction = "pca") 

#FindNeighbors
mou17 <- FindNeighbors(mou17, 
                       reduction = "pca", 
                       dims = 1:11)
mou17 <- FindClusters(mou17, 
                      resolution = 1.7)

#Look at cluster IDs of the first 5 cells
head(Idents(mou17), 5)

#deDoulblets 
sweep.mou17 <- paramSweep_v3(mou17, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.mou17 , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(mou17@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(mou17@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
mou17 <- doubletFinder_v3(mou17, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###mou17 <- doubletFinder_v3(mou17, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
mou17@meta.data[,"DF_hi.lo"] <- mou17@meta.data$DF.classifications_0.25_0.09_643
mou17@meta.data$DF_hi.lo[which(mou17@meta.data$DF_hi.lo == "Doublet" & mou17@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(mou17, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
mou17_1<-subset(mou17,subset=DF_hi.lo=="Singlet")

#UMAP
mou17 <- RunUMAP(mou17, 
                 dims = 1:30, 
                 label = T)
head(mou17@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(mou17, reduction = "umap")
LabelClusters(DimPlot(mou17, reduction = "umap"),id = 'ident')
DimPlot(mou17, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
mou17 <- RunTSNE(mou17, dims = 1:30)
head(mou17@reductions$tsne@cell.embeddings)
p2 <- DimPlot(mou17, reduction = "tsne")
LabelClusters(DimPlot(mou17, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(mou17, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(mou17, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
mou17.markers <- FindAllMarkers(mou17, only.pos = TRUE, min.pct = 0.25)
mou17.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(mou17.markers,file="allmarkers_mou17.csv")

#DoHeatmap
top10 <- mou17.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(mou17, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(mou17, c("Ascl1","Hes5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(mou17, c(" Nkx2-1","Lhx6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(mou17, c(" Mesi2", "Isl1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(mou17, c("Calb2"," Nr2f2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(mou17, c("Pdgfra")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(mou17, c("Cx3cr1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(mou17, c("Col3a1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(mou17,file="mou17_.Rdata")
