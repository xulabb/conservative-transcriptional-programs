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

mou15.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mou15_filtered_feature_bc_matrix")
mou15 <- CreateSeuratObject(counts = mou15.data, 
                            project = 'mou15', 
                            min.cells = 13, 
                            min.features = 800)
gc()
mou15
##
mou15[["percent.mt"]] <- PercentageFeatureSet(mou15, features = "^mt")
VlnPlot(mou15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mou15, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mou15, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
mou15<- subset(mou15, 
               subset = nFeature_RNA > 800 & 
                 nFeature_RNA < 6000 & 
                 percent.mt < 20)

mou15 <- NormalizeData(mou15,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
mou15<- FindVariableFeatures(mou15, selection.method = "vst", nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mou15), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mou15)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(mou15)
mou15 <- ScaleData(mou15, features = all.genes)

#Perform linear dimensional reduction
mou15 <- RunPCA(mou15, features = VariableFeatures(object = mou15))

#Examine and visualize PCA results a few different ways
print(mou15[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(mou15, dims = 1:2, reduction = "pca")
DimPlot(mou15, reduction = "pca")
DimHeatmap(mou15, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
mou15 <- JackStraw(mou15, num.replicate = 100)
mou15 <- ScoreJackStraw(mou15, dims = 1:50)
JackStrawPlot(mou15, dims = 1:40)
ElbowPlot(object = mou15, ndims = 40, reduction = "pca") 

#FindNeighbors
mou15 <- FindNeighbors(mou15, 
                       reduction = "pca", 
                       dims = 1:11)
mou15 <- FindClusters(mou15, 
                      resolution = 1.7)

#Look at cluster IDs of the first 5 cells
head(Idents(mou15), 5)

#deDoulblets 
sweep.mou15 <- paramSweep_v3(mou15, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.mou15 , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(mou15@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(mou15@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
mou15 <- doubletFinder_v3(mou15, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###mou15 <- doubletFinder_v3(mou15, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
mou15@meta.data[,"DF_hi.lo"] <- mou15@meta.data$DF.classifications_0.25_0.09_643
mou15@meta.data$DF_hi.lo[which(mou15@meta.data$DF_hi.lo == "Doublet" & mou15@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(mou15, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
mou15_1<-subset(mou15,subset=DF_hi.lo=="Singlet")

#UMAP
mou15 <- RunUMAP(mou15, 
                 dims = 1:30, 
                 label = T)
head(mou15@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(mou15, reduction = "umap")
LabelClusters(DimPlot(mou15, reduction = "umap"),id = 'ident')
DimPlot(mou15, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
mou15 <- RunTSNE(mou15, dims = 1:30)
head(mou15@reductions$tsne@cell.embeddings)
p2 <- DimPlot(mou15, reduction = "tsne")
LabelClusters(DimPlot(mou15, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(mou15, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(mou15, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
mou15.markers <- FindAllMarkers(mou15, only.pos = TRUE, min.pct = 0.25)
mou15.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(mou15.markers,file="allmarkers_mou15.csv")

#DoHeatmap
top10 <- mou15.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(mou15, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(mou15, c("Ascl1","Hes5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(mou15, c(" Nkx2-1","Lhx6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(mou15, c(" Mesi2", "Isl1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(mou15, c("Calb2"," Nr2f2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(mou15, c("Pdgfra")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(mou15, c("Cx3cr1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(mou15, c("Col3a1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(mou15,file="mou15_.Rdata")
