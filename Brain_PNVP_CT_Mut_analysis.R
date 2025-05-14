library(Seurat)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(viridis)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(GGally)

hg_to_mg <- function(y) {
  genes <- tolower(y)
  for(i in 1:length(genes)){
  	if(grepl("^mt-", genes[i])){
  		genes[i] <- paste0("mt-", toupper(substring(genes[i], 4, 4)), substring(genes[i], 5))
  	}else{
  		genes[i] <- paste0(toupper(substring(genes[i], 1, 1)), substring(genes[i], 2))
  	}
  }
  return(genes)
}

setwd("~/projects/def-dubraca/jhowa105/singlecell/AVM/integrated/Brain_PNVP_CT_Mut")

brain <- readRDS("../Brain_CT_Mut/brain.rds")
brain$type <- "brain"
brain_lst <- SplitObject(brain, split.by = "condition")
pnvp <- readRDS("../PNVP_CT_Mut/brain.rds")
pnvp$type <- "pnvp"
pnvp_lst <- SplitObject(pnvp, split.by = "condition")

obj_lst <- c(brain_lst, pnvp_lst)
obj_lst <- lapply(obj_lst, function(x){
    DefaultAssay(x) <- "RNA"
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    return(x)
})

features <- SelectIntegrationFeatures(object.list = obj_lst)
anchors <- FindIntegrationAnchors(object.list = obj_lst, anchor.features = features)
seur <- IntegrateData(anchorset = anchors)

DefaultAssay(seur) <- "integrated"

seur <- ScaleData(seur, verbose = FALSE)
seur <- RunPCA(seur, npcs = 50, verbose = FALSE)
png("elbowplot.png")
ElbowPlot(seur, ndims = 50)
dev.off()
seur <- RunUMAP(seur, reduction = "pca", dims = 1:30)
seur <- FindNeighbors(seur, reduction = "pca", dims = 1:30)
seur <- FindClusters(seur, resolution = 0.5)

png("umap.png", units = "in", height = 5, width = 5, res = 150)
DimPlot(seur, group.by = "seurat_clusters", label = T) + NoLegend()
dev.off()

png("umap_condition.png", units = "in", height = 5, width = 5, res = 150)
DimPlot(seur, group.by = "condition")
dev.off()

png("umap_type.png", units = "in", height = 5, width = 5, res = 150)
DimPlot(seur, group.by = "type")
dev.off()

png("umap_annotation1.png", units = "in", height = 5, width = 10, res = 150)
DimPlot(seur, group.by = "annotation1")
dev.off()

png("umap_annotation2.png", units = "in", height = 5, width = 10, res = 150)
DimPlot(seur, group.by = "annotation2")
dev.off()

saveRDS(seur, "seur.rds")