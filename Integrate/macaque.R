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

setwd('~/Documents/Doctor/lab/develoment_of_early_ge/split analysis/macintegrated')

mac40.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mac40_filtered_feature_bc_matrix")
mac40 <- CreateSeuratObject(counts = mac40.data, 
                                    project = 'mac40', 
                                    min.cells = 13, 
                                    min.features = 800)
mac40 <- NormalizeData(mac40,
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)
mac40<- FindVariableFeatures(mac40, 
                            selection.method = "vst", 
                            nfeatures = 3000)
mac45.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mac45_filtered_feature_bc_matrix")
mac45 <- CreateSeuratObject(counts = mac45.data, 
                            project = 'mac45', 
                            min.cells = 13, 
                            min.features = 800)
mac45 <- NormalizeData(mac45,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
mac45<- FindVariableFeatures(mac45, 
                             selection.method = "vst", 
                             nfeatures = 3000)
mac53.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mac53_filtered_feature_bc_matrix")
mac53 <- CreateSeuratObject(counts = mac53.data, 
                            project = 'mac53', 
                            min.cells = 13, 
                            min.features = 800)
mac53 <- NormalizeData(mac53,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
mac53<- FindVariableFeatures(mac53, 
                             selection.method = "vst", 
                             nfeatures = 3000)
gc()
mac53
mac45
mac40

##Integrate
macintegrated.anc <- FindIntegrationAnchors(object.list = list(mac45,mac40,mac53), 
                                            anchor.features = 2000, 
                                            scale = T, 
                                            reduction = 'cca', 
                                            dims = 1:30, 
                                            k.filter = 200)
macintegrated <- IntegrateData(anchorset = macintegrated.anc, dims = 1:30)
macintegrated
gc()

##
macintegrated[["percent.mt"]] <- PercentageFeatureSet(macintegrated, features = mt.genes)
VlnPlot(macintegrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(macintegrated, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(macintegrated, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
macintegrated<- subset(macintegrated, 
                       subset = nFeature_RNA > 800 & 
                         nFeature_RNA < 6000 & 
                         percent.mt < 20)

macintegrated <- NormalizeData(macintegrated,
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
macintegrated<- FindVariableFeatures(macintegrated, 
                                     selection.method = "vst", 
                                     nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(macintegrated), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(macintegrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(macintegrated)
macintegrated <- ScaleData(macintegrated, features = all.genes)

#Perform linear dimensional reduction
macintegrated <- RunPCA(macintegrated, 
                        features = VariableFeatures(object = macintegrated),
                        npcs=50)

#Examine and visualize PCA results a few different ways
print(macintegrated[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(macintegrated, dims = 1:2, reduction = "pca")
DimPlot(macintegrated, reduction = "pca")
DimHeatmap(macintegrated, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
macintegrated <- JackStraw(macintegrated, num.replicate = 100)
macintegrated <- ScoreJackStraw(macintegrated, dims = 1:50)
JackStrawPlot(macintegrated, dims = 1:40)
ElbowPlot(object = macintegrated, ndims = 40, reduction = "pca") 

#FindNeighbors
macintegrated <- FindNeighbors(macintegrated, 
                               reduction = "pca", 
                               dims = 1:11)
macintegrated <- FindClusters(macintegrated, 
                              resolution = 1.7)

#Look at cluster IDs of the first 5 cells
head(Idents(macintegrated), 5)

#deDoulblets 
sweep.macintegrated <- paramSweep_v3(macintegrated, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.macintegrated , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(macintegrated@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(macintegrated@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
macintegrated <- doubletFinder_v3(macintegrated, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###macintegrated <- doubletFinder_v3(macintegrated, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
macintegrated@meta.data[,"DF_hi.lo"] <- macintegrated@meta.data$DF.classifications_0.25_0.09_643
macintegrated@meta.data$DF_hi.lo[which(macintegrated@meta.data$DF_hi.lo == "Doublet" & macintegrated@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(macintegrated, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
macintegrated_1<-subset(macintegrated,subset=DF_hi.lo=="Singlet")

#UMAP
macintegrated <- RunUMAP(macintegrated, 
                         dims = 1:30, 
                         label = T)
head(macintegrated@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(macintegrated, reduction = "umap")
LabelClusters(DimPlot(macintegrated, reduction = "umap"),id = 'ident')
DimPlot(macintegrated, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
macintegrated <- RunTSNE(macintegrated, dims = 1:30)
head(macintegrated@reductions$tsne@cell.embeddings)
p2 <- DimPlot(macintegrated, reduction = "tsne")
LabelClusters(DimPlot(macintegrated, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(macintegrated, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(macintegrated, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
macintegrated.markers <- FindAllMarkers(macintegrated, only.pos = TRUE, min.pct = 0.25)
macintegrated.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(macintegrated.markers,file="allmarkers_macintegrated.csv")

#DoHeatmap
top10 <- macintegrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(macintegrated, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(macintegrated, c("ASCL1","HES5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(macintegrated, c(" NKX2-1","LHX6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(macintegrated, c(" MEIS2", "ISL1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(macintegrated, c("CALB2"," NR2F2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(macintegrated, c("PDGFRA")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(macintegrated, c("CX3CR1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(macintegrated, c("COL3A1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(macintegrated,file="macintegrated_.Rdata")
