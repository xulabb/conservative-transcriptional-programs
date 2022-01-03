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

setwd('~/Documents/Doctor/lab/develoment_of_early_ge/split analysis/hsaintegrated')

hsa9.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/hsa9_filtered_feature_bc_matrix")
hsa9 <- CreateSeuratObject(counts = hsa9.data, 
                            project = 'hsa9', 
                            min.cells = 13, 
                            min.features = 800)
hsa9 <- NormalizeData(hsa9,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
hsa9<- FindVariableFeatures(hsa9, 
                             selection.method = "vst", 
                             nfeatures = 3000)
hsa13.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/hsa13_filtered_feature_bc_matrix")
hsa13 <- CreateSeuratObject(counts = hsa13.data, 
                            project = 'hsa13', 
                            min.cells = 13, 
                            min.features = 800)
hsa13 <- NormalizeData(hsa13,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
hsa13<- FindVariableFeatures(hsa13, 
                             selection.method = "vst", 
                             nfeatures = 3000)

gc()
hsa13
hsa9

##Integrate
hsaintegrated.anc <- FindIntegrationAnchors(object.list = list(hsa13,hsa9), 
                                            anchor.features = 2000, 
                                            scale = T, 
                                            reduction = 'cca', 
                                            dims = 1:30, 
                                            k.filter = 200)
hsaintegrated <- IntegrateData(anchorset = hsaintegrated.anc, dims = 1:30)
hsaintegrated
gc()

##
hsaintegrated[["percent.mt"]] <- PercentageFeatureSet(hsaintegrated, features = "^MT")
VlnPlot(hsaintegrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(hsaintegrated, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(hsaintegrated, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
hsaintegrated<- subset(hsaintegrated, 
                       subset = nFeature_RNA > 800 & 
                         nFeature_RNA < 6000 & 
                         percent.mt < 20)

hsaintegrated <- NormalizeData(hsaintegrated,
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
hsaintegrated<- FindVariableFeatures(hsaintegrated, 
                                     selection.method = "vst", 
                                     nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hsaintegrated), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hsaintegrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(hsaintegrated)
hsaintegrated <- ScaleData(hsaintegrated, features = all.genes)

#Perform linear dimensional reduction
hsaintegrated <- RunPCA(hsaintegrated, 
                        features = VariableFeatures(object = hsaintegrated),
                        npcs=50)

#Examine and visualize PCA results a few different ways
print(hsaintegrated[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(hsaintegrated, dims = 1:2, reduction = "pca")
DimPlot(hsaintegrated, reduction = "pca")
DimHeatmap(hsaintegrated, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
hsaintegrated <- JackStraw(hsaintegrated, num.replicate = 100)
hsaintegrated <- ScoreJackStraw(hsaintegrated, dims = 1:50)
JackStrawPlot(hsaintegrated, dims = 1:40)
ElbowPlot(object = hsaintegrated, ndims = 40, reduction = "pca") 

#FindNeighbors
hsaintegrated <- FindNeighbors(hsaintegrated, 
                               reduction = "pca", 
                               dims = 1:18)
hsaintegrated <- FindClusters(hsaintegrated, 
                              resolution = 0.5)

#Look at cluster IDs of the first 5 cells
head(Idents(hsaintegrated), 5)

#deDoulblets 
sweep.hsaintegrated <- paramSweep_v3(hsaintegrated, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.hsaintegrated , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(hsaintegrated@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(hsaintegrated@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
hsaintegrated <- doubletFinder_v3(hsaintegrated, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###hsaintegrated <- doubletFinder_v3(hsaintegrated, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
hsaintegrated@meta.data[,"DF_hi.lo"] <- hsaintegrated@meta.data$DF.classifications_0.25_0.09_643
hsaintegrated@meta.data$DF_hi.lo[which(hsaintegrated@meta.data$DF_hi.lo == "Doublet" & hsaintegrated@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(hsaintegrated, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
hsaintegrated_1<-subset(hsaintegrated,subset=DF_hi.lo=="Singlet")

#UMAP
hsaintegrated <- RunUMAP(hsaintegrated, 
                         dims = 1:30, 
                         label = T)
head(hsaintegrated@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(hsaintegrated, reduction = "umap")
LabelClusters(DimPlot(hsaintegrated, reduction = "umap"),id = 'ident')
DimPlot(hsaintegrated, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
hsaintegrated <- RunTSNE(hsaintegrated, dims = 1:30)
head(hsaintegrated@reductions$tsne@cell.embeddings)
p2 <- DimPlot(hsaintegrated, reduction = "tsne")
LabelClusters(DimPlot(hsaintegrated, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(hsaintegrated, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(hsaintegrated, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
hsaintegrated.markers <- FindAllMarkers(hsaintegrated, only.pos = TRUE, min.pct = 0.25)
hsaintegrated.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(hsaintegrated.markers,file="allmarkers_hsaintegrated.csv")

#DoHeatmap
top10 <- hsaintegrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(hsaintegrated, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(hsaintegrated, c("ASCL1","HES5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(hsaintegrated, c(" NKX2-1","LHX6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(hsaintegrated, c(" MEIS2", "ISL1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(hsaintegrated, c("CALB2"," NR2F2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(hsaintegrated, c("PDGFRA")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(hsaintegrated, c("CX3CR1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(hsaintegrated, c("COL3A1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(hsaintegrated,file="hsaintegrated_.Rdata")
