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

setwd('~/Documents/Doctor/lab/develoment_of_early_ge/split analysis/primate')

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
hsa13
hsa9
mac53
mac45
mac40
mou12.5_MGE
mou12.5_LGE

##Integrate
primate.anc <- FindIntegrationAnchors(object.list = list(mac45,mac40,mac53,hsa13,hsa9,mou12.5_MGE,mou12.5_LGE), 
                                      anchor.features = 2000, 
                                      scale = T, 
                                      reduction = 'cca', 
                                      dims = 1:30, 
                                      k.filter = 200)
primate <- IntegrateData(anchorset = primate.anc, dims = 1:30)
primate
gc()

##
primate <- NormalizeData(primate,
                         normalization.method = "LogNormalize", 
                         scale.factor = 10000)
#Identification of genes with high expression variation between cells
# The goal of this step is to identify genes that vary greatly from cell to cell for subsequent identification of cell types,
# We use the default parameter "VST" to select 2000 highly variable genes.
primate<- FindVariableFeatures(primate, 
                               selection.method = "vst", 
                               nfeatures = 3000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(primate), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(primate)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(primate)
primate <- ScaleData(primate, features = all.genes)

#Perform linear dimensional reduction
primate <- RunPCA(primate, 
                  features = VariableFeatures(object = primate),
                  npcs=50)

#Examine and visualize PCA results a few different ways
print(primate[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(primate, dims = 1:2, reduction = "pca")
DimPlot(primate, reduction = "pca")
DimHeatmap(primate, dims = 1, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
primate <- JackStraw(primate, num.replicate = 100)
primate <- ScoreJackStraw(primate, dims = 1:50)
JackStrawPlot(primate, dims = 1:40)
ElbowPlot(object = primate, ndims = 40, reduction = "pca") 

#FindNeighbors
primate <- FindNeighbors(primate, 
                         reduction = "pca", 
                         dims = 1:11)
primate <- FindClusters(primate, 
                        resolution = 0.9)

#Look at cluster IDs of the first 5 cells
head(Idents(primate), 5)

#deDoulblets 
sweep.primate <- paramSweep_v3(primate, PCs = 1:30, sct = FALSE)
sweep.stats_ge77 <- summarizeSweep(sweep.primate , GT = FALSE)
bcmvn_ge77 <- find.pK(sweep.stats_ge77)
## 
homotypic.prop <- modelHomotypic(primate@meta.data$ClusteringResults)         
nExp_poi <- round(0.075*nrow(primate@meta.data)) 
## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 
primate <- doubletFinder_v3(primate, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
###primate <- doubletFinder_v3(primate, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

## Plot results --------------------------------------------------------------------------------------------------------------
primate@meta.data[,"DF_hi.lo"] <- primate@meta.data$DF.classifications_0.25_0.09_643
primate@meta.data$DF_hi.lo[which(primate@meta.data$DF_hi.lo == "Doublet" & primate@meta.data$DF.classifications_0.25_0.09_643 == "Singlet")] <- "Doublet_lo"
TSNEPlot(primate, group.by="DF_hi.lo", plot.order=c("Doublet_hi","Doublet_lo","Singlet"), label.color=c("black","gold","red"))
###########
primate_1<-subset(primate,subset=DF_hi.lo=="Singlet")

#UMAP
primate <- RunUMAP(primate, 
                   dims = 1:11, 
                   label = T)
head(primate@reductions$umap@cell.embeddings) 
# Extract UMAP coordinates.
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(primate, reduction = "umap")
LabelClusters(DimPlot(primate, reduction = "umap"),id = 'ident')
DimPlot(primate, cols = umapcol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#T-SNE
primate <- RunTSNE(primate, dims = 1:30)
head(primate@reductions$tsne@cell.embeddings)
p2 <- DimPlot(primate, reduction = "tsne")
LabelClusters(DimPlot(primate, reduction = "tsne"),id = 'ident')
p1 + p2
DimPlot(primate, cols = tsnecol) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) & labs(title = '') & 
  NoLegend()

#find all markers of cluster 1
cluster1.markers <- FindMarkers(primate, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 30)

#find all marker genes
primate.markers <- FindAllMarkers(primate, only.pos = TRUE, min.pct = 0.25)
primate.markers %>% group_by(cluster) %>% top_n(n= 2, wt = "avg_logFC")
write.csv(primate.markers,file="allmarkers_primate.csv")

#DoHeatmap
top10 <- primate.markers %>% group_by(cluster) %>% top_n(n = 10, wt = "avg_logFC")
DoHeatmap(primate, features = top10$gene) + NoLegend()

#marker gene
#RGC&IPC
FeaturePlot(primate, c("ASCL1","HES5")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#MGE
FeaturePlot(primate, c(" NKX2-1","LHX6")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#LGE
FeaturePlot(primate, c(" MEIS2", "ISL1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#CGE
FeaturePlot(primate, c("CALB2"," NR2F2 ")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#OPC
FeaturePlot(primate, c("PDGFRA")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#microglia 
FeaturePlot(primate, c("CX3CR1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

#pericytes 
FeaturePlot(primate, c("COL3A1")) & 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) &
  NoLegend()

saveRDS(primate,file="primate_.Rdata")
