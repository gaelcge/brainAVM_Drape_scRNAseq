library(Seurat)
library(ggplot2)
library(DropletUtils)
library(loomR)
library(SeuratDisk)
library(reticulate)
library(gtools)

hg_to_mg <- function(y) {
  genes <- tolower(y)
  for(i in 1:length(genes)){
  	if(grepl("^mt-", genes[i])){
  		genes[i] <- paste0("mt-", toupper(substring(genes[i], 4, 4)), substrin(genes[i, 5]))
  	}else{
  		genes[i] <- paste0(toupper(substring(genes[i], 1, 1)), substring(genes[i], 2))
  	}
  }
  return(genes)
}

setwd("/home/jhowa105/projects/def-dubraca/jhowa105/singlecell/AVM/BrainCT/analysis")

# # emptyDrops - Remove empty droplets from data #
# data_dir <- "/home/jhowa105/projects/def-dubraca/jhowa105/singlecell/AVM/BrainCT/cellranger_outs/raw_feature_bc_matrix"
# sce <- read10xCounts(data_dir)
# sce

# br.out <- barcodeRanks(counts(sce))
# # Making barcode rank plot
# png("barcode_rank.png")
# plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
# o <- order(br.out$rank)
# lines(br.out$rank[o], br.out$fitted[o], col="red")
# abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
# abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
# legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#     legend=c("knee", "inflection"))
# dev.off()

# e.out <- emptyDrops(counts(sce))
# is.cell <- e.out$FDR <= 0.01
# sum(is.cell, na.rm=TRUE)

# table(Limited=e.out$Limited, Significant=is.cell)

# png("emptydrops.png")
# plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),xlab="Total UMI count", ylab="-Log Probability")
# dev.off()

# sce <- sce[,which(is.cell)]
# #-----------------------------------------------#
# #Convert sce to andata for doublet detection #
# sce_mat <- as.matrix(counts(sce))
# colnames(sce_mat) <- sce@colData@listData$Barcode
# rownames(sce_mat) <- sce@rowRanges@elementMetadata@listData$Symbol
# seur <- CreateSeuratObject(counts = sce_mat)
# seur <- FindVariableFeatures(seur) # Needed for as.loom to work
# loom <- as.loom(seur, filename = "for_solo.loom")
# loom$close()
# closeAllConnections()
# #---------------------------------------------#
# #saveRDS(sce, "sce.rds")

sce <- readRDS("sce.rds")

# Remove Doublets from data #
np <- import("numpy")
is_doublet <- np$load("./solo_28.06.22_14_51_22/is_doublet.npy")
preds <- np$load("./solo_28.06.22_14_51_22/preds.npy")
sce_mat <- as.matrix(counts(sce))
colnames(sce_mat) <- sce@colData@listData$Barcode
rownames(sce_mat) <- sce@rowRanges@elementMetadata@listData$Symbol

tmp_mat <- sce_mat
tmp_mat[which(tmp_mat>0)] <- 1
cell_has_gene <- rowSums(tmp_mat)
gene_in_cell <- colSums(tmp_mat)
min_cells <- 10
min_features <- 100
keep_gene <- which(cell_has_gene > min_cells)
bcs <- colnames(tmp_mat)[which(gene_in_cell > min_features)]
sce_mat <- sce_mat[keep_gene,]

#seur <- CreateSeuratObject(counts = sce_mat, min.cells = 10, min.features = 100)
seur <- CreateSeuratObject(counts = sce_mat)
seur[["percent.mt"]] <- PercentageFeatureSet(seur, "^mt-")
seur$is_doublet <- is_doublet

seur <- subset(seur, cells = bcs)

png("qc_scatter.png")
ggplot(seur@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point(aes(colour = percent.mt))
dev.off()
png("qc_scatter2.png")
ggplot(seur@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point(aes(colour = is_doublet))
dev.off()

png("qc_vln.png", units = "in", height = 5, width = 15, res = 150)
VlnPlot(seur, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
dev.off()

seur <- subset(seur, subset = nCount_RNA < 45000 & nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 10)

png("qc_scatter_post_filt.png")
ggplot(seur@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point(aes(colour = percent.mt))
dev.off()

# Seurat Work FLow
seur <- NormalizeData(seur)
seur <- CellCycleScoring(seur, g2m.features = hg_to_mg(cc.genes.updated.2019$g2m.genes), s.features = hg_to_mg(cc.genes.updated.2019$s.genes))
#seur <- FindVariableFeatures(seur)
#seur <- ScaleData(seur, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "G2M.Score", "S.Score"))
seur <- SCTransform(seur)
seur <- RunPCA(seur)

# seur <- JackStraw(seur, num.replicate = 100)
# seur <- ScoreJackStraw(seur, dims = 1:20)

# png("jackstraw.png", units = "in", height = 10, width = 10, res = 150)
# JackStrawPlot(seur, dims = 1:20)
# dev.off()

png("elbowplot.png")
ElbowPlot(seur, ndim = 50)
dev.off()

seur <- RunUMAP(seur, dims = 1:15, n.neighbors = 20)
seur <- FindNeighbors(seur, dims = 1:15)
seur <- FindClusters(seur)

png("umap.png", units = "in", height = 5, width = 5, res = 150)
DimPlot(seur, label = T)
dev.off()

png("qc_feature.png", units = "in", height = 10, width = 10, res = 150)
FeaturePlot(seur, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "is_doublet"), ncol = 2)
dev.off()

ec_genes <- c("Bmx", "Fbln5", "Gja4", "Unc5b", "Dll4", "Notch1", "Esm1", "Apod", "Apln", "Mcam", "Chst1", "Kcne3", "Nrp2", "Nr2f2", "Aplnr", "Prcp", "Mfsd2a", "Spock2", "Mki67", "Pcna", "Birc5")
DefaultAssay(seur) <- "RNA"
png("ec_markers_featureplot.png", units = "in", height = 25, width = 25, res = 150)
FeaturePlot(seur, features = ec_genes, ncol = 5)
dev.off()

DefaultAssay(seur) <- "RNA"
markers <- FindAllMarkers(seur)
markers$cluster <- as.character(markers$cluster)
markers <- markers[which(markers$p_val_adj < 0.05),]
markers <- markers[order(markers$avg_log2FC, decreasing = TRUE),]
grp_lst <- list()
for(i in unique(markers$cluster)){
    grp_lst[[i]] <- NULL
}
for(i in 1:nrow(markers)){
    if(markers$gene[i] %in% unlist(grp_lst)) next
    if(length(unlist(grp_lst)) == length(unique(markers$cluster)) * 2) break
    if(length(grp_lst[[markers$cluster[i]]]) < 2){
        grp_lst[[markers$cluster[i]]] <- c(grp_lst[[markers$cluster[i]]], markers$gene[i])
    }
}
grp_lst <- unlist(grp_lst)[mixedorder(names(unlist(grp_lst)))]
png(paste0("vln_plot_top_2.png"), units = "in", height = 10, width = 10, res = 150)
print({VlnPlot(seur, features = grp_lst, pt.size = 0, assay = "RNA", stack = TRUE, fill.by = "ident", flip = TRUE)})
dev.off()

grp_lst <- list()
for(i in unique(markers$cluster)){
    grp_lst[[i]] <- NULL
}
for(i in 1:nrow(markers)){
    if(markers$gene[i] %in% unlist(grp_lst)) next
    if(length(unlist(grp_lst)) == length(unique(markers$cluster)) * 5) break
    if(length(grp_lst[[markers$cluster[i]]]) < 5){
        grp_lst[[markers$cluster[i]]] <- c(grp_lst[[markers$cluster[i]]], markers$gene[i])
    }
}
grp_lst <- unlist(grp_lst)[mixedorder(names(unlist(grp_lst)))]
png(paste0("dot_plot_top_5.png"), units = "in", height = 10, width = 30, res = 150)
DotPlot(seur, features = grp_lst, assay = "RNA") + RotatedAxis()
dev.off()

seur$Phase <- factor(seur$Phase, levels = rev(c("G1", "S", "G2M")))
seur$cluster_condition <- paste(seur$seurat_clusters, seur$condition)
seur$cluster_condition <- factor(seur$cluster_condition, levels = unlist(lapply(sort(unique(as.numeric(as.character(seur$seurat_clusters)))), function(x) paste(x, c('Control', 'Alk5iEKO')))))
png("phase_barplot.png", units = "in", height = 5, width = 20, res = 150)
ggplot(seur@meta.data) + 
    geom_bar(color = 'black', position = 'fill') + 
    aes(x = cluster_condition, fill = Phase, by = cluster_condition) +
    geom_text(stat = "prop", position = position_fill(.5), size = 3) +
    theme_base() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_blank()) 
dev.off()

seur$condition <- "CT"
seur <- subset(seur, subset = seurat_clusters == 3, invert = TRUE)
seur$cell_type <- "Proliferative ECs"
seur$cell_type[which(seur$seurat_clusters == 0)] <- "Venous Capillary ECs"
seur$cell_type[which(seur$seurat_clusters == 1)] <- "Capillary ECs"
seur$cell_type[which(seur$seurat_clusters == 2)] <- "Arterial Capillary ECs"
seur$cell_type[which(seur$seurat_clusters == 4)] <- "Gm26917 Pos Capillary ECs"
seur$cell_type[which(seur$seurat_clusters == 5)] <- "Capillary ECs"
seur$cell_type[which(seur$seurat_clusters == 6)] <- "Tip Cells"
seur$cell_type[which(seur$seurat_clusters == 7)] <- "Arterial ECs"
seur$cell_type[which(seur$seurat_clusters == 9)] <- "Venous ECs"
seur$cell_type[which(seur$seurat_clusters == 11)] <- "Choroids"

png("umap_celltype.png", units = "in", height = 5, width = 7, res = 150)
DimPlot(seur, group.by = "cell_type")
dev.off()

saveRDS(seur, "seur.rds")

# for(d in c(10, 15, 20, 30)){
#     for(n in c(10, 20, 30)){
#         DefaultAssay(seur) <- "SCT"
#         seur <- RunUMAP(seur, dims = 1:d, n.neighbors = n)
#         seur <- FindNeighbors(seur, dims = 1:d)
#         seur <- FindClusters(seur)

#         png(paste0("./umaps/dims_", d, "_neighbors_", n, "_umap.png"), units = "in", height = 5, width = 5, res = 150)
#         print({DimPlot(seur, label = T) + NoLegend()})
#         dev.off()

#         DefaultAssay(seur) <- "RNA"
#         markers <- FindAllMarkers(seur)
#         markers$cluster <- as.character(markers$cluster)
#         markers <- markers[which(markers$p_val_adj < 0.05),]
#         markers <- markers[order(markers$avg_log2FC, decreasing = TRUE),]
#         grp_lst <- list()
#         for(i in unique(markers$cluster)){
#             grp_lst[[i]] <- NULL
#         }
#         for(i in 1:nrow(markers)){
#             if(markers$gene[i] %in% unlist(grp_lst)) next
#             if(length(unlist(grp_lst)) == length(unique(markers$cluster)) * 2) break
#             if(length(grp_lst[[markers$cluster[i]]]) < 2){
#                 grp_lst[[markers$cluster[i]]] <- c(grp_lst[[markers$cluster[i]]], markers$gene[i])
#             }
#         }
#         grp_lst <- unlist(grp_lst)[mixedorder(names(unlist(grp_lst)))]
#         png(paste0("./umaps/dims_", d, "_neighbors_", n, "_vln.png"), units = "in", height = 10, width = 10, res = 150)
#         print({VlnPlot(seur, features = grp_lst, pt.size = 0, assay = "RNA", stack = TRUE, fill.by = "ident", flip = TRUE)})
#         dev.off()
#     }
# }
