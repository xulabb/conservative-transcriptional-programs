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

setwd('~/Documents/Doctor/lab/develoment_of_early_ge/split analysis/mouse2')

mou11.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mou11_filtered_feature_bc_matrix")
mou11 <- CreateSeuratObject(counts = mou11.data, 
                                  project = 'mou11', 
                                  min.cells = 13, 
                                  min.features = 800)
mou11 <- NormalizeData(mou11,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
mou11<- FindVariableFeatures(mou11, 
                             selection.method = "vst", 
                             nfeatures = 3000)

mou15.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mou15_filtered_feature_bc_matrix")
mou15 <- CreateSeuratObject(counts = mou15.data, 
                                  project = 'mou15', 
                                  min.cells = 13, 
                                  min.features = 800)
mou15 <- NormalizeData(mou15,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
mou15<- FindVariableFeatures(mou15, 
                             selection.method = "vst", 
                             nfeatures = 3000)

mou17.data <- Read10X(data.dir = "~/Documents/Doctor/lab/develoment_of_early_ge/matrix/mou17_filtered_feature_bc_matrix")
mou17 <- CreateSeuratObject(counts = mou17.data, 
                                  project = 'mou17', 
                                  min.cells = 13, 
                                  min.features = 800)
mou17 <- NormalizeData(mou17,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)
mou17<- FindVariableFeatures(mou17, 
                             selection.method = "vst", 
                             nfeatures = 3000)

gc()
mou17
mou15
mou11

##Integrate
mouse2.anc <- FindIntegrationAnchors(object.list = list(mou15,mou11,mou17), 
                                     anchor.features = 2000, 
                                     scale = T, 
                                     reduction = 'cca', 
                                     dims = 1:30, 
                                     k.filter = 200)
mouse2 <- IntegrateData(anchorset = mouse2.anc, dims = 1:30)
mouse2
gc()

##
mouse2[["percent.mt"]] <- PercentageFeatureSet(mouse2, features = "^mt")
VlnPlot(mouse2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(mouse2, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse2, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
mouse2<- subset(mouse2, 
                subset = nFeature_RNA > 800 & 
                  nFeature_RNA < 6000 & 
                  percent.mt < 20)

mouse2 <- NormalizeData(mouse2,
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
mouse2<- FindVariableFeatures(mouse2, selection.method = "vst", nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mouse2), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mouse2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(mouse2)
mouse2 <- ScaleData(mouse2, features = all.genes)

#Perform linear dimensional reduction
mouse2 <- RunPCA(mouse2, 
                 features = VariableFeatures(object = mouse2),
                 npcs=50)

#Examine and visualize PCA results a few different ways
print(mouse2[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(mouse2, dims = 1:2, reduction = "pca")
DimPlot(mouse2, reduction = "pca")
DimHeatmap(mouse2, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
mouse2 <- JackStraw(mouse2, num.replicate = 100)
mouse2 <- ScoreJackStraw(mouse2, dims = 1:50)
JackStrawPlot(mouse2, dims = 1:40)
ElbowPlot(object = mouse2, ndims = 40, reduction = "pca") 

#FindNeighbors
mouse2 <- FindNeighbors(mouse2, 
                        reduction = "pca", 
                        dims = 1:16)
mouse2 <- FindClusters(mouse2, 
                       resolution = 0.5)

#Look at cluster IDs of the first 5 cells
head(Idents(mouse2), 5)

#deDoulblets 
sweep.mouse2 <- paramSweep_v3(mouse2, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.mouse2 , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(mouse2@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(mouse2@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
mouse2 <- doubletFinder_v3(mouse2, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###mouse2 <- doubletFinder_v3(mouse2, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
mouse2@meta.data[,"DF_hi.lo"] <- mouse2@meta.data$DF.classifications_0.25_0.09_643
mouse2@meta.data$DF_hi.lo[which(mouse2@meta.data$DF_hi.lo == "Doublet" & mouse2@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(mouse2, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
mouse2_1<-subset(mouse2,subset=DF_hi.lo=="Singlet")

#UMAP
mouse2 <- RunUMAP(mouse2, 
                  dims = 1:16, 
                  label = T)
head(mouse2@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(mouse2, reduction = "umap")
LabelClusters(DimPlot(mouse2, reduction = "umap"),id = 'ident')
DimPlot(mouse2, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
mouse2 <- RunTSNE(mouse2, dims = 1:30)
head(mouse2@reductions$tsne@cell.embeddings)
p2 <- DimPlot(mouse2, reduction = "tsne")
LabelClusters(DimPlot(mouse2, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(mouse2, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(mouse2, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
mouse2.markers <- FindAllMarkers(mouse2, only.pos = TRUE, min.pct = 0.25)
mouse2.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(mouse2.markers,file="allmarkers_mouse2.csv")

#DoHeatmap
top10 <- mouse2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(mouse2, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(mouse2, c("Ascl1","Hes5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(mouse2, c(" Nkx2-1","Lhx6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(mouse2, c(" Mesi2", "Isl1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(mouse2, c("Calb2"," Nr2f2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(mouse2, c("Pdgfra")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(mouse2, c("Cx3cr1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(mouse2, c("Col3a1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(mouse2,file="mouse2_.Rdata")
