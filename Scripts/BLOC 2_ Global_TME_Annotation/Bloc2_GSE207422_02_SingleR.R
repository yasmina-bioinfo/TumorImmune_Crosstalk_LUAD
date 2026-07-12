#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc2 Script 02: SingleR automated annotation
# NOTE: matrix conversion via BPCells::convert_matrix_type() (used
# for GSE243013) is not needed here, since GSE207422 was built from
# a plain text matrix, not a BPCells-backed object. Standard dgCMatrix
# extraction is sufficient.
# Input:  Objects/Bloc1_GSE207422_04_seu_clustered.rds
# Output: Results/Tables/Bloc2_GSE207422_SingleR_main_per_cluster.csv
#         Results/Tables/Bloc2_GSE207422_SingleR_fine_per_cluster.csv
#         Results/Figures/Annotations_TME_GSE207422/Bloc2_GSE207422_SingleR_*.png
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc1_GSE207422_04_seu_clustered.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/Annotations_TME_GSE207422")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

message("Loading clustered Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

message("Loading HumanPrimaryCellAtlas reference...")
ref <- HumanPrimaryCellAtlasData()

seu <- JoinLayers(seu)
mat <- GetAssayData(seu, layer = "data")
mat <- as(mat, "dgCMatrix")

message("Running SingleR (label.main)...")
pred_main <- SingleR(
  test    = mat,
  ref     = ref,
  labels  = ref$label.main,
  BPPARAM = BiocParallel::SnowParam(1)
)
message("SingleR main annotation done.")

message("Running SingleR (label.fine)...")
pred_fine <- SingleR(
  test    = mat,
  ref     = ref,
  labels  = ref$label.fine,
  BPPARAM = BiocParallel::SnowParam(1)
)
message("SingleR fine annotation done.")

seu$SingleR_main <- pred_main$labels
seu$SingleR_fine <- pred_fine$labels

message("Generating summary tables...")

summary_main <- as.data.frame(table(
  Cluster = seu$seurat_clusters,
  SingleR = seu$SingleR_main
))
summary_main <- summary_main[summary_main$Freq > 0, ]
fwrite(summary_main, file.path(OUT_TAB, "Bloc2_GSE207422_SingleR_main_per_cluster.csv"))

summary_fine <- as.data.frame(table(
  Cluster = seu$seurat_clusters,
  SingleR = seu$SingleR_fine
))
summary_fine <- summary_fine[summary_fine$Freq > 0, ]
fwrite(summary_fine, file.path(OUT_TAB, "Bloc2_GSE207422_SingleR_fine_per_cluster.csv"))

message("Saved: SingleR summary tables")

message("Generating heatmaps...")

png(file.path(OUT_FIG, "Bloc2_GSE207422_SingleR_heatmap_main.png"),
    width = 12, height = 8, units = "in", res = 300)
plotScoreHeatmap(pred_main)
dev.off()
message("Saved: Bloc2_GSE207422_SingleR_heatmap_main.png")

png(file.path(OUT_FIG, "Bloc2_GSE207422_SingleR_heatmap_fine.png"),
    width = 16, height = 10, units = "in", res = 300)
plotScoreHeatmap(pred_fine)
dev.off()
message("Saved: Bloc2_GSE207422_SingleR_heatmap_fine.png")

message("Generating UMAP figures...")

png(file.path(OUT_FIG, "Bloc2_GSE207422_UMAP_SingleR_main.png"),
    width = 14, height = 8, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "SingleR_main",
              label     = FALSE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "SingleR annotation (main) — GSE207422 NSCLC"))
dev.off()
message("Saved: Bloc2_GSE207422_UMAP_SingleR_main.png")

png(file.path(OUT_FIG, "Bloc2_GSE207422_UMAP_SingleR_main_labeled.png"),
    width = 14, height = 8, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "SingleR_main",
              label     = FALSE,
              raster    = FALSE) +
        geom_text(data = data.frame(
          Embeddings(seu, "umap"),
          cluster = seu$seurat_clusters) %>%
            group_by(cluster) %>%
            summarise(umap_1 = median(umap_1), umap_2 = median(umap_2)),
          aes(x = umap_1, y = umap_2, label = cluster),
          size = 3, fontface = "bold", color = "black") +
        theme_bw() +
        labs(title = "SingleR annotation (main) with cluster numbers — GSE207422 NSCLC"))
dev.off()
message("Saved: Bloc2_GSE207422_UMAP_SingleR_main_labeled.png")

png(file.path(OUT_FIG, "Bloc2_GSE207422_UMAP_SingleR_fine.png"),
    width = 18, height = 10, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "SingleR_fine",
              label     = FALSE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "SingleR annotation (fine) — GSE207422 NSCLC") +
        theme(legend.text = element_text(size = 6)))
dev.off()
message("Saved: Bloc2_GSE207422_UMAP_SingleR_fine.png")

message("Saving updated Seurat object...")
saveRDS(seu, file.path(DATA_DIR, "Objects/Bloc2_GSE207422_02_seu_singler.rds"))
message("Saved: Objects/Bloc2_GSE207422_02_seu_singler.rds")
message("DONE Bloc2 GSE207422 Script 02")