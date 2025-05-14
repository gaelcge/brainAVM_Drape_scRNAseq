rm(list = ls())
library(plotly)
library(ggrepel)
library(gtools)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(useful)
library(Matrix)
library(Seurat)
#library(DoubletFinder)
#library(tidyverse)
#library(scico)
#library(DropletUtils)
theme_set(theme_bw())
library(future)
#library(SoupX)
library(SingleR)
library(pals)
plan("multicore", workers = (availableCores()-2))
options(future.globals.maxSize = 30000000000 * 1024^3)



###Set directory




setwd("/home/gaelcge/projects/ctb-lavvin00/gaelcge/Pigu/Dubrac_AVM_Elise/AVM_Elise/SeuratObject/")

Seurat_object <- readRDS(paste("brain_pnvp_ct_mut_2.rds", sep="."))


setwd("/home/gaelcge/projects/ctb-dubraca/shared/projects/Dubrac_AVM_Elise/AVM_Elise/SeuratObject/New_annotation_2")



# Seurat_object <- merge(x=Seurat_object_1, 
#   y=c(Seurat_object_2,Seurat_object_3,Seurat_object_4,Seurat_object_5,
#     Seurat_object_6,Seurat_object_7,Seurat_object_8,Seurat_object_9,
#     Seurat_object_10,Seurat_object_11,Seurat_object_12))

# Seurat_object <- JoinLayers(Seurat_object)

# Seurat_object[["RNA"]] <- split(Seurat_object[["RNA"]], f = Seurat_object$Rep)

# Seurat_object <- SCTransform(Seurat_object, verbose = TRUE, vars.to.regress = c("nFeature_RNA", "percent.mt"))

# Seurat_object <- RunPCA(Seurat_object, features = VariableFeatures(object = Seurat_object))


# Seurat_object <- IntegrateLayers(object = Seurat_object, method = CCAIntegration, normalization.method = "SCT", verbose = F)


#Seurat_object$timepoint <- factor(Seurat_object$timepoint, levels = c("D1","D3",  "D7" , "D21", "D56"))

#Seurat_object$celltype_timepoint <- factor(paste(Seurat_object$celltype,Seurat_object$timepoint, sep="_"))

plot1 <- FeatureScatter(Seurat_object, feature1 = "percent.mt", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="nCount_RNA.nGenePlot.png", width=2300, height=600, bg = "white", res = 150)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#We filter out cells that have unique gene counts over 6000
#Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nGene but on anything in the object (as in GenePlot above)


#Seurat_object <- subset(Seurat_object, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & nCount_RNA < 30000 & percent.mt < 10)



png(filename="VlnPlotQC.orig.ident.png", width=1500, height=1000, bg = "white", res = 150)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=-1, group.by="orig.ident")
dev.off()



png(filename="VlnPlotQC.condition.png", width=3500, height=1000, bg = "white", res = 150)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=-1, group.by="condition")
dev.off()


png(filename="VlnPlotQC.is_doublet.png", width=1500, height=1000, bg = "white", res = 150)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=-1, group.by="is_doublet")
dev.off()


png(filename="VlnPlotQC.type.png", width=1500, height=1000, bg = "white", res = 150)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=-1, group.by="type")
dev.off()



##Normalize data

Seurat_object <- SCTransform(Seurat_object, verbose = TRUE, vars.to.regress = c("nFeature_RNA", "percent.mt"))


#Cell Cycle scoring

#cc.genes <- readLines(con="/home/gaelcge/projects/def-jsjoyal/gaelcge/Seurat_ressource/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

#s.genes <- cc.genes[1:43]
#g2m.genes <- cc.genes[44:98]


#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes

#Seurat_object <- CellCycleScoring(Seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)




plot1 <- FeatureScatter(Seurat_object, feature1 = "percent.mt", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="nFeature_RNA.nGenePlot.SCT.png", width=900, height=600, bg = "white", res = 50)
CombinePlots(plots = list(plot1, plot2))
dev.off()


# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(Seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(filename="VariableFeaturePlot.png", width=2500, height=600, bg = "white", res = 150)
CombinePlots(plots = list(plot1, plot2))
dev.off()


Seurat_object <- RunPCA(Seurat_object, features = VariableFeatures(object = Seurat_object))

# Examine and visualize PCA results a few different ways
print(Seurat_object[["pca"]], dims = 1:5, nfeatures = 5)


png(filename="VizDimLoadings.sct.png", width=1500, height=900, bg = "white", res = 150)
VizDimLoadings(Seurat_object, dims = 1:2, reduction = "pca")
dev.off()

png(filename="DimHeatmap.sct.png", width=1500, height=3000, bg = "white", res = 150)
DimHeatmap(Seurat_object, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

png(filename="ElbowPlot.sct.png", width=1500, height=1000, bg = "white", res = 150)
ElbowPlot(Seurat_object)
dev.off()


png(filename="DimPlot.cellcycle.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "pca", group.by ="Phase")
dev.off()


DefaultAssay(Seurat_object) <- "SCT"


##Run non-linear dimensional reduction (UMAP/tSNE)


Seurat_object <- RunUMAP(Seurat_object, dims = 1:15,
                            n.components = 2L)

##Cluster the cells

Seurat_object <- FindNeighbors(Seurat_object, 
                dims = 1:2, 
                reduction = "umap"
                )

Seurat_object <- FindClusters(Seurat_object, 
                resolution = 1.2, 
                reduction = "umap"
                )


png(filename="umap.res1.2.sct.png", width=1100, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE)
dev.off()

png(filename="umap.Phase.sct.png", width=1000, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, group.by="Phase")
dev.off()

##
####################
#doublet detection

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
#sweep.res.list_Seurat_object <- paramSweep_v3(Seurat_object, PCs = 1:10, sct = FALSE)
#sweep.stats_Seurat_object <- summarizeSweep(sweep.res.list_Seurat_object, GT = FALSE)
#bcmvn_Seurat_object <- find.pK(sweep.stats_Seurat_object)

#


## Initial run, nExp set to 0.075 Poisson loading estimate (e.g., 1000 total doublet predictions)
#nExp_poi <- round(0.075*length(rownames(Seurat_object@meta.data)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset

## Run DoubletFinder with varying classification stringencies ---------------------

#Seurat_object <- doubletFinder_v3(Seurat_object, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------

#homotypic.prop <- modelHomotypic(Seurat_object@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#Seurat_object <- doubletFinder(Seurat_object, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct=TRUE)


summary(Seurat_object@meta.data)

#nbcol <- length(colnames(Seurat_object@meta.data))

#colnames(Seurat_object@meta.data)[(nbcol-1)] <- "pANN"
#colnames(Seurat_object@meta.data)[(nbcol)] <- "DF.classifications"

png(filename="FeaturePlot_QC.sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

png(filename="FeaturePlot_ImmuneCells.sct.png", width=1500, height=1500, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("Cx3cr1", "Aif1", "Tmem119", "P2ry12", "Olfml3", "Sall1", "Ptprc", "Itgam"))
dev.off()


png(filename="umap.type.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "type")
dev.off()


png(filename="umap.annotation1_type.png", width=3600, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, group.by = "annotation1_type", repel=TRUE, split.by="type_condition")
dev.off()


png(filename="umap.annotation2.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, group.by = "annotation2", repel=TRUE)
dev.off()

png(filename="umap.annotation1.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, group.by = "annotation1", repel=TRUE)
dev.off()

png(filename="umap.cell_type.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, group.by = "cell_type", repel=TRUE)
dev.off()


tab = prop.table(table(Seurat_object$annotation1_type, Seurat_object$seurat_clusters), margin = 2)*100
library(pals)
pdf(paste("prop.table","annotation1_type_seurat_clusters","plot.pdf", sep="."), width=12, height=8)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=as.vector(alphabet2(18)))
dev.off()

tab = prop.table(table(Seurat_object$annotation2, Seurat_object$seurat_clusters), margin = 2)*100
library(pals)
pdf(paste("prop.table","annotation2_seurat_clusters","plot.pdf", sep="."), width=12, height=8)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=as.vector(alphabet2(18)))
dev.off()

tab = prop.table(table(Seurat_object$cell_type, Seurat_object$seurat_clusters), margin = 2)*100
library(pals)
pdf(paste("prop.table","cell_type_seurat_clusters","plot.pdf", sep="."), width=12, height=8)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=as.vector(alphabet2(18)))
dev.off()


Seurat_object$type_condition <- paste(Seurat_object@meta.data$type, Seurat_object@meta.data$condition, sep="_")



png(filename="DotPlot_markers_Alex.seuratcluster.type_condition.sct.png", res = 150, width=800, height=3500)
DotPlot(
  Seurat_object,
  assay = "SCT",
  c("Kit", "Pgf", "Slc16a3", "Stc1"),
  cols = c("blue", "red", "orange", "purple"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by=NULL,
  split.by = "type_condition",
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()

DefaultAssay(Seurat_object) <- "SCT"

png(filename="FeaturePlot_markers_Alex.seuratcluster.type_condition.sct.png", res = 150, width=1500, height=1500)
FeaturePlot(
  Seurat_object,
  c("Kit", "Pgf", "Slc16a3", "Stc1"),
  cols = c("lightblue", "red"),
  split.by = "type_condition",
  reduction = "umap", order=T
) 
dev.off()


tab = prop.table(table(Seurat_object$type_condition, Seurat_object$seurat_clusters), margin = 2)*100
library(pals)
pdf(paste("prop.table","type_condition_seurat_clusters","plot.pdf", sep="."), width=12, height=8)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=as.vector(alphabet2(18)))
dev.off()

# Seurat_object <- SubsetData(
#   Seurat_object,
#   assay = NULL,
#   cells = NULL,
#   subset.name = "DF.classifications",
#   ident.use = NULL,
#   ident.remove = NULL,
#   low.threshold = -Inf,
#   high.threshold = Inf,
#   accept.value = "Singlet",
#   max.cells.per.ident = Inf,
#   random.seed = 1)

# png(filename="umap.doublet_filtered.sct.png", width=1500, height=900, bg = "white", res = 150)
# DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "DF.classifications")
# dev.off()



# find markers for every cluster compared to all remaining cells, report only the positive ones

DefaultAssay(Seurat_object) <- "RNA"

Seurat_object <- JoinLayers(Seurat_object)

Seurat_object <- NormalizeData(Seurat_object, verbose = FALSE)

Seurat_object <- ScaleData(Seurat_object, verbose = FALSE, vars.to.regress = c("nFeature_RNA", "percent.mt"))


markers.dotplot <- rev(c("Cldn5", "Pecam1", "Kcnj8", "Rgs5", "Gdf2", "C1qa"))


png(filename="DotPlot_markers.rna.seuratcluster.png", res = 150, width=1200, height=1500)
DotPlot(
  Seurat_object,
  assay = NULL,
  markers.dotplot,
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by=NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()


Seurat_object.markers <- FindAllMarkers(Seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Seurat_object.markers, "Seurat_object.cluster.markers.txt")

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

png(filename="VlnPlot_top2Markers.sct.png", width=3500, height=4500, bg = "white", res = 150)
VlnPlot(Seurat_object, features = unique(top2$gene))
dev.off()

png(filename="FeaturePlot_top2Markers.sct.png", width=3000, height=8000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c(unique(top2$gene)))
dev.off()

png(filename="DoHeatmap_top10Markers.sct.png", width=2000, height=3000, bg = "white", res = 150)
DoHeatmap(Seurat_object, features = top10$gene) + NoLegend()
dev.off()

png(filename="DoHeatmap_top2Markers.sct.png", width=900, height=700, bg = "white", res = 150)
DoHeatmap(Seurat_object, features = top2$gene) + NoLegend()
dev.off()

png(filename="DotPlot_markers.top10.sct.png", res = 150, width=4500, height=1000)
DotPlot(
  Seurat_object,
  assay = NULL,
  unique(top10$gene),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()

png(filename="DotPlot_markers.top2.sct.png", res = 150, width=3500, height=1000)
DotPlot(
  Seurat_object,
  assay = NULL,
  unique(top2$gene),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()


####


##Using scMAGIC (Not working)

#list.Ref <- readRDS('/home/gaelcge/projects/ctb-lavvin00/gaelcge/Pigu/SingleCellRefdb/brain/mouse/Tasic.Rdata')
#ref.mtx <- list.Ref$mat_exp
#ref.labels <-list.Ref$label[, 1]
# load target dataset

#exp_sc_mat <- Seurat_object@assays$RNA$counts

#Cell type classification by scMAGIC
#output.scMAGIC <- scMAGIC(exp_sc_mat, ref.mtx, ref.labels, atlas = 'MCA', num_threads = 4)

#Error because of Array5
#pred.scMAGIC <- output.scMAGIC$scMAGIC.tag



### Using singleR to determine annotation

# transform into sce

Seurat_object.se <- as.SingleCellExperiment(Seurat_object, assay = "SCT")

#ref_Seurat <- readRDS("/home/gaelcge/projects/def-jsjoyal/gaelcge/Dipolo/SeuratV4/Rod_depleted/Clustering/Seurat_object.CTL_OHT_integrated.rds")


#Predict with Tasic dataset
list.Ref <- readRDS('/home/gaelcge/projects/ctb-lavvin00/gaelcge/Pigu/SingleCellRefdb/brain/mouse/Tasic.Rdata')
ref.mtx <- list.Ref$mat_exp
ref.labels <-list.Ref$label[, 1]

ref_Seurat <- CreateSeuratObject(counts = ref.mtx, project = "Tasic_brain_ref", 
  min.cells = 3, min.features = 100)
ref_Seurat$labels <- ref.labels

ref_Seurat[["percent.mt"]] <- PercentageFeatureSet(ref_Seurat, pattern = "^mt-")
ref_Seurat <- SCTransform(ref_Seurat, verbose = TRUE, vars.to.regress = c("nFeature_RNA", "percent.mt"))

ref.se <- as.SingleCellExperiment(ref_Seurat, assay = "SCT")


pred.Seurat_object <- SingleR(test=Seurat_object.se, ref=ref.se, labels=ref.se$labels, de.method="wilcox")


png(filename="plotScoreHeatmap_SingleR.Tasic.png", width=1500, height=900, bg = "white", res = 150)
plotScoreHeatmap(pred.Seurat_object)
dev.off()

Seurat_object$SingleR_annot_Tasic <- pred.Seurat_object$pruned.labels

png(filename="umap.sct_filtered.SingleR_annot_Tasic.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label = F, pt.size = 0.5, group.by="SingleR_annot_Tasic") #+ NoLegend()
dev.off()


tab = prop.table(table(Seurat_object$SingleR_annot_Tasic, Seurat_object$seurat_clusters), margin = 2)*100
library(pals)
pdf(paste("prop.table","SingleR_annot_Tasic_seurat_clusters","plot.pdf", sep="."), width=8, height=4)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=as.vector(alphabet2(13)))
dev.off()



#Predict with Celldex dataset

library(celldex)

MouseRNAseqData <- readRDS("/home/gaelcge/projects/ctb-lavvin00/gaelcge/Pigu/SingleCellRefdb/MouseRNAseqData.rds")

pred.Seurat_object <- SingleR(test=Seurat_object.se, ref=MouseRNAseqData, labels=MouseRNAseqData$label.main, de.method="wilcox")


Seurat_object$SingleR_annot_MouseRNAseqData <- pred.Seurat_object$pruned.labels

png(filename="umap.sct_filtered.SingleR_annot_MouseRNAseqData.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label = F, pt.size = 0.5, group.by="SingleR_annot_MouseRNAseqData") #+ NoLegend()
dev.off()


tab = prop.table(table(Seurat_object$type_condition, Seurat_object$seurat_clusters), margin = 2)*100
library(pals)
pdf(paste("prop.table","type_condition_seurat_clusters","plot.pdf", sep="."), width=8, height=4)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=as.vector(alphabet2(13)))
dev.off()




Idents(Seurat_object) <- "Cell_Type_ECs"

library("pals")
tab = prop.table(table(Seurat_object$Cell_Type_ECs, Seurat_object$seurat_clusters), margin = 2)*100
pdf(paste("prop.table","Cell_Type_ECs.seurat_clusters","plot.pdf", sep="."), width=8, height=4)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=as.vector(tol.rainbow(length(unique(Seurat_object@active.ident)))))
dev.off()

### Assigning annotation
Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")

DefaultAssay(Seurat_object) <- "SCT"
Seurat_object <- BuildClusterTree(Seurat_object)
pdf("PlotClusterTree.pdf", width=10, height=5)
PlotClusterTree(Seurat_object)
dev.off()

DefaultAssay(Seurat_object) <- "RNA"

Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")
new.cluster.ids <- c("Venous_ECs_like",      # cluster 0
                      "Prolif_ECs",     # cluster 1
                      "CapVenous_ECs",     # cluster 2
                      "CapArterial_ECs",     # cluster 3
                      "Cap_ECs",     # cluster 4
                      "Venous_ECs_like",     # cluster 5
                      "Tip_ECs_like",     # cluster 6
                      "Tip_ECs_like",     # cluster 7
                      "CapArterial_ECs",     # cluster 8
                      "Tip_ECs_like",     # cluster 9
                      "Tip_ECs",     # cluster 10
                      "Cap_ECs",     # cluster 11
                      "Arteriole_ECs",     # cluster 12
                      "Venular_ECs_like", # cluster 13
                      "Tip_ECs_like", # cluster 14
                      "Cap_ECs",     # cluster 15
                      "Tip_ECs",     # cluster 16
                      "Venous_ECs",  # cluster 17
                      "CapVenous_ECs",     # cluster 18
                      "CapVenous_ECs",     # cluster 19
                      "Arterial_ECs", # cluster 20
                      "CapArterial_ECs", # cluster 21
                      "Choroidal_ECs",      # cluster 22
                      "Slc16a3_ECs", #cluster 23
                      "Prolif_ECs",    # cluster 24
                      "CapArterial_ECs",     # cluster 25
                       "Tip_ECs",     # cluster 26
                       "Cap_ECs", # cluster 27
                       "Venular_ECs",      # cluster 28
                       "Prolif_ECs",     # cluster 29
                       "Tip_ECs_like",     # cluster 30
                       "Venous_ECs",  # cluster 31
                       "Arterial_ECs_like",     # cluster 32
                       "Tip_ECs"     # cluster 33
                      # "Neuron", # cluster 34
                      # "Neuron", # cluster 35
                      # "Neuron", # cluster 36
                      # "Mixed_1", # cluster 37
                      # "Neuron", # cluster 38
                      # "Neuron", # cluster 39
                      # "Neuron", # cluster 40
                      # "Neuron", # cluster 41
                      # "Mixed_2", # cluster 42
                      # "Neuron", # cluster 43
                      # "Immune_2", # cluster 44
                      # "Neuron", # cluster 45
                      # "Neuron", # cluster 46
                      # "Neuron" # cluster 47
                      )     



names(new.cluster.ids) <- levels(Seurat_object)
Seurat_object <- RenameIdents(Seurat_object, new.cluster.ids)

Seurat_object[["Cell_Type_ECs"]] <- Idents(object = Seurat_object)

Seurat_object$Cell_Type_ECs_Condition <- paste(Seurat_object$Cell_Type_ECs, Seurat_object$condition, sep="_")
Seurat_object$Cell_Type_ECs_type_condition <- paste(Seurat_object$Cell_Type_ECs, Seurat_object$type_condition, sep="_")


Seurat_object <- subset(Seurat_object, 
  subset=Cell_Type_ECs_type_condition %in% c("Cap_ECs_pnvp_MUT", "CapArterial_ECs_pnvp_MUT", "CapVenous_ECs_pnvp_MUT",
    "Slc16a3_ECs_brain_CT", "Slc16a3_ECs_brain_MUT", "Slc16a3_ECs_pnvp_CT", "Venular_ECs_brain_MUT"), 
  invert=T)

pdf("umap.sct_filtered.annotated.pdf", width=6, height=6)
DimPlot(Seurat_object, reduction = "umap", label = TRUE, pt.size = 1,
  cols=as.vector(tol.rainbow(length(unique(Seurat_object$Cell_Type_ECs))))) + NoLegend()
dev.off()

pdf("umap.sct_filtered.split.by.type_condition.pdf", width=24, height=8)
DimPlot(Seurat_object, reduction = "umap", label = T, pt.size = 1, shuffle = T, split.by="type_condition",
  cols=as.vector(tol.rainbow(length(unique(Seurat_object$Cell_Type_ECs)))))
dev.off()


pdf("umap.sct_filtered.split.by.condition.pdf", width=14, height=8)
DimPlot(Seurat_object, reduction = "umap", label = T, pt.size = 1, shuffle = T, split.by="condition",
  cols=as.vector(tol.rainbow(length(unique(Seurat_object$Cell_Type_ECs)))))
dev.off()




Idents(Seurat_object) <- "Cell_Type_ECs"

library("pals")
tab = prop.table(table(Seurat_object$Cell_Type_ECs, Seurat_object$type_condition), margin = 2)*100
pdf(paste("prop.table","Cell_Type_ECs.type_condition","plot.pdf", sep="."), width=4, height=6)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=as.vector(tol.rainbow(length(unique(Seurat_object@active.ident)))))
dev.off()


write.table(tab, "Prop_Cell_Type_ECs_type_condition.txt", sep="\t", quote=F, col.name=NA, row.name=T)




##Find markers of cell types

df.Conserved.markers <- data.frame()
listcells <- c("Cap_ECs", "CapArterial_ECs", "Choroidal_ECs", 
  "Prolif_ECs", "Tip_ECs", "Arteriole_ECs", "CapVenous_ECs", "Venous_ECs", "Arterial_ECs")

for(i in listcells) {
  Seurat_object.Conserved.markers <- FindConservedMarkers(Seurat_object, ident.1=i, grouping.var = "condition",
                    only.pos = TRUE,
                    min.diff.pct = 0.25,
                    min.pct = 0.25,
                    logfc.threshold = 0.25)
  Seurat_object.Conserved.markers$celltype <- i
  Seurat_object.Conserved.markers$genes <- rownames(Seurat_object.Conserved.markers)
  df.Conserved.markers <- rbind(df.Conserved.markers, Seurat_object.Conserved.markers)
}

df.Conserved.markers$ALL_avg_log2FC <- (df.Conserved.markers$MUT_avg_log2FC+df.Conserved.markers$CT_avg_log2FC)/2

write.table(df.Conserved.markers, "Seurat_object.celltypes.conserved.markers.txt", sep="\t", quote=F, col.name=NA, row.name=T)

top2 <- df.Conserved.markers %>% group_by(celltype) %>% top_n(n = 2, wt = ALL_avg_log2FC)


pdf("DotPlot_Conserved.markers.Cell_Type_ECs_type_condition.pdf", width=10, height=14)
DotPlot(
  Seurat_object,
  assay = NULL,
  unique(top2$genes),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_Type_ECs_type_condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()



top10 <- df.Conserved.markers %>% group_by(celltype) %>% top_n(n = 10, wt = ALL_avg_log2FC)


pdf("DotPlot_Conserved.markers.top10.Cell_Type_ECs_type_condition.pdf", width=20, height=14)
DotPlot(
  Seurat_object,
  assay = NULL,
  unique(top10$genes),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_Type_ECs_type_condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()



Seurat_object.markers <- FindAllMarkers(Seurat_object, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

write.table(Seurat_object.markers, "Seurat_object.celltypes.ECs.markers.txt", sep="\t", quote=F, col.name=NA, row.name=T)


Seurat_object.markers <- FindAllMarkers(Seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

png(filename="VlnPlot_top2Markers.celltypes.png", width=3500, height=4500, bg = "white", res = 150)
VlnPlot(Seurat_object, features = unique(top2$gene))
dev.off()

png(filename="FeaturePlot_top2Markers.celltypes.png", width=3000, height=3000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c(unique(top2$gene)))
dev.off()

png(filename="DoHeatmap_top10Markers.celltypes.png", width=2000, height=3000, bg = "white", res = 150)
DoHeatmap(Seurat_object, features = top10$gene) + NoLegend()
dev.off()

png(filename="DoHeatmap_top2Markers.celltypes.png", width=900, height=700, bg = "white", res = 150)
DoHeatmap(Seurat_object, features = top2$gene) + NoLegend()
dev.off()

pdf("DotPlot_markers.top10.celltypes.pdf", width=35, height=10)
DotPlot(
  Seurat_object,
  assay = NULL,
  unique(top10$gene),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_Type_ECs_type_condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf("DotPlot_markers.celltypes.pdf", width=10, height=6)
DotPlot(
  Seurat_object,
  assay = NULL,
  unique(top2$gene),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()



## Calculate distance caused by Mutant

## Regroup

Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")

new.cluster.ids <- c("Venous_ECs",      # cluster 0
                      "Prolif_ECs",     # cluster 1
                      "CapVenous_ECs",     # cluster 2
                      "CapArterial_ECs",     # cluster 3
                      "Cap_ECs",     # cluster 4
                      "Venous_ECs",     # cluster 5
                      "Tip_ECs",     # cluster 6
                      "Tip_ECs",     # cluster 7
                      "CapArterial_ECs",     # cluster 8
                      "Tip_ECs",     # cluster 9
                      "Tip_ECs",     # cluster 10
                      "Cap_ECs",     # cluster 11
                      "Arteriole_ECs",     # cluster 12
                      "Venular_ECs", # cluster 13
                      "Tip_ECs", # cluster 14
                      "Cap_ECs",     # cluster 15
                      "Tip_ECs",     # cluster 16
                      "Venous_ECs",  # cluster 17
                      "CapVenous_ECs",     # cluster 18
                      "CapVenous_ECs",     # cluster 19
                      "Arterial_ECs", # cluster 20
                      "CapArterial_ECs", # cluster 21
                      "Choroid_ECs",      # cluster 22
                      "Slc16a3_ECs", #cluster 23
                      "Prolif_ECs",    # cluster 24
                      "CapArterial_ECs",     # cluster 25
                       "Tip_ECs",     # cluster 26
                       "Cap_ECs", # cluster 27
                       "Venular_ECs",      # cluster 28
                       "Prolif_ECs",     # cluster 29
                       "Tip_ECs",     # cluster 30
                       "Venous_ECs",  # cluster 31
                       "Arterial_ECs",     # cluster 32
                       "Tip_ECs"     # cluster 33
                      # "Neuron", # cluster 34
                      # "Neuron", # cluster 35
                      # "Neuron", # cluster 36
                      # "Mixed_1", # cluster 37
                      # "Neuron", # cluster 38
                      # "Neuron", # cluster 39
                      # "Neuron", # cluster 40
                      # "Neuron", # cluster 41
                      # "Mixed_2", # cluster 42
                      # "Neuron", # cluster 43
                      # "Immune_2", # cluster 44
                      # "Neuron", # cluster 45
                      # "Neuron", # cluster 46
                      # "Neuron" # cluster 47
                      )     




names(new.cluster.ids) <- levels(Seurat_object)
Seurat_object <- RenameIdents(Seurat_object, new.cluster.ids)

Seurat_object[["Cell_Type_ECs_grouped"]] <- Idents(object = Seurat_object)



pdf("umap.sct_filtered.grouped.split.by.type_condition.pdf", width=24, height=8)
DimPlot(Seurat_object, reduction = "umap", label = T, pt.size = 1, shuffle = T, split.by="type_condition",
  cols=as.vector(tol.rainbow(length(unique(Seurat_object$Cell_Type_ECs_grouped)))))
dev.off()

saveRDS(Seurat_object, paste("brain_pnvp_ct_mut_no_integration", "Seurat_object.rds", sep="."))

##

setwd("/home/gaelcge/projects/ctb-dubraca/shared/projects/Dubrac_AVM_Elise/AVM_Elise/SeuratObject/New_annotation_2")

Seurat_object <- readRDS(paste("brain_pnvp_ct_mut_no_integration", "Seurat_object.rds", sep="."))



Seurat_object_scDist <- subset(Seurat_object, ident=c("Slc16a3_ECs"), invert=T)

table(Idents(object = Seurat_object_scDist))


library(scDist)


out <- scDist(as.matrix(Seurat_object_scDist@assays$SCT$data),Seurat_object_scDist@meta.data,fixed.effects = "condition",
              random.effects="type",
              clusters="Cell_Type_ECs_grouped")

pdf("scDist_Cell_Type_ECs_grouped_condition_random_type.pdf", width = 5, height = 5)
DistPlot(out)
dev.off()


out$results

pdf("distGenes_Venular_ECs_condition_random_type.pdf", width = 5, height = 5)
distGenes(out, cluster = "Venular_ECs")
dev.off()

pdf("distGenes_Arterial_ECs_condition_random_type.pdf", width = 5, height = 5)
distGenes(out, cluster = "Arterial_ECs")
dev.off()



q("no")
