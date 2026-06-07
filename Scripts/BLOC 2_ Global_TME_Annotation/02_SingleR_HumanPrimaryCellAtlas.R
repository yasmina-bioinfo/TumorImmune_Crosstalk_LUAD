#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc2 Script 02: SingleR automated annotation
# Input:  Objects/Bloc1_04_seu_clustered.rds
# Output: Results/Tables/Bloc2_SingleR_annotation.csv
#         Results/Figures/Bloc2_SingleR_heatmap_main.png
#         Results/Figures/Bloc2_SingleR_heatmap_fine.png
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc1_04_seu_clustered.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# -----------------------------
# 1) Load clustered object
# -----------------------------
message("Loading clustered Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

# -----------------------------
# 2) Load reference : HumanPrimaryCellAtlas
# -----------------------------
# HumanPrimaryCellAtlas covers both immune and non-immune cell types
# Suitable for TME datasets with epithelial + immune compartments
message("Loading HumanPrimaryCellAtlas reference...")
ref <- HumanPrimaryCellAtlasData()

# -----------------------------
# 3) Extract normalized expression matrix
# -----------------------------
# NOTE: JoinLayers() required for Seurat v5 multi-layer objects after merge()
seu <- JoinLayers(seu)
mat <- GetAssayData(seu, layer = "data")

# NOTE: BPCells matrix not directly compatible with SingleR
# Two options depending on available RAM:

# Option A — Server use (high RAM): converts to full sparse matrix
# mat <- as.matrix(mat)
# mat <- as(mat, "sparseMatrix")

# Option B — 16GB RAM (used here): lightweight float conversion
# convert_matrix_type(float) reduces memory footprint significantly
mat <- BPCells::convert_matrix_type(mat, type = "float")
mat <- as(mat, "dgCMatrix")

# -----------------------------
# 4) Run SingleR : label.main (broad cell types)
# -----------------------------
message("Running SingleR (label.main)...")
pred_main <- SingleR(
  test      = mat,
  ref       = ref,
  labels    = ref$label.main,
  BPPARAM   = BiocParallel::SnowParam(1)  # single core for stability and MulticoreParam not supported by Windows
)
message("SingleR main annotation done.")

# -----------------------------
# 5) Run SingleR : label.fine (fine-grained subtypes)
# -----------------------------
# NOTE: label.fine retained for reference but not primary analysis tool
# Too many categories for robust interpretation at this stage
# Will be revisited if specific sub-populations require finer resolution
message("Running SingleR (label.fine)...")
pred_fine <- SingleR(
  test    = mat,
  ref     = ref,
  labels  = ref$label.fine,
  BPPARAM = BiocParallel::SnowParam(1)
)
message("SingleR fine annotation done.")

# -----------------------------
# 6) Add annotations to Seurat object
# -----------------------------
seu$SingleR_main <- pred_main$labels
seu$SingleR_fine <- pred_fine$labels

# Primary analysis tool: compare SingleR main labels vs seurat_clusters
# More reliable than UMAP visualization for concordance assessment
message("Generating summary tables...")

summary_main <- as.data.frame(table(
  Cluster = seu$seurat_clusters,
  SingleR = seu$SingleR_main
))
summary_main <- summary_main[summary_main$Freq > 0, ]
fwrite(summary_main, file.path(OUT_TAB, "Bloc2_SingleR_main_per_cluster.csv"))

summary_fine <- as.data.frame(table(
  Cluster = seu$seurat_clusters,
  SingleR = seu$SingleR_fine
))
summary_fine <- summary_fine[summary_fine$Freq > 0, ]
fwrite(summary_fine, file.path(OUT_TAB, "Bloc2_SingleR_fine_per_cluster.csv"))

message("Saved: SingleR summary tables")

# -----------------------------
# 8) Heatmaps
# -----------------------------
message("Generating heatmaps...")

png(file.path(OUT_FIG, "Bloc2_SingleR_heatmap_main.png"),
    width = 12, height = 8, units = "in", res = 300)
plotScoreHeatmap(pred_main)
dev.off()
message("Saved: Bloc2_SingleR_heatmap_main.png")

png(file.path(OUT_FIG, "Bloc2_SingleR_heatmap_fine.png"),
    width = 16, height = 10, units = "in", res = 300)
plotScoreHeatmap(pred_fine)
dev.off()
message("Saved: Bloc2_SingleR_heatmap_fine.png")

# -----------------------------
# 9) UMAP figures
# -----------------------------
message("Generating UMAP figures...")

# NOTE: raster = FALSE required for large datasets (>100,000 points)
# label = FALSE for fine annotation ; too many categories to label legibly

# UMAP main — no labels (colors too similar for reliable reading)
png(file.path(OUT_FIG, "Bloc2_UMAP_SingleR_main.png"),
    width = 14, height = 8, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "SingleR_main",
              label     = FALSE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "SingleR annotation (main) — GSE243013 LUAD"))
dev.off()
message("Saved: Bloc2_UMAP_SingleR_main.png")

# UMAP main, with Seurat cluster numbers overlaid
# NOTE: cluster numbers overlaid on SingleR colors for cross-reference
# Primary visual tool for concordance check between manual and SingleR annotation
png(file.path(OUT_FIG, "Bloc2_UMAP_SingleR_main_labeled.png"),
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
        labs(title = "SingleR annotation (main) with cluster numbers — GSE243013 LUAD"))
dev.off()
message("Saved: Bloc2_UMAP_SingleR_main_labeled.png")

# UMAP fine — no labels (too many categories, legend only)
png(file.path(OUT_FIG, "Bloc2_UMAP_SingleR_fine.png"),
    width = 18, height = 10, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "SingleR_fine",
              label     = FALSE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "SingleR annotation (fine) — GSE243013 LUAD") +
        theme(legend.text = element_text(size = 6)))
dev.off()
message("Saved: Bloc2_UMAP_SingleR_fine.png")

# -----------------------------
# 10) Save updated Seurat object
# -----------------------------
message("Saving updated Seurat object...")
saveRDS(seu, file.path(DATA_DIR, "Objects/Bloc2_02_seu_singler.rds"))
message("Saved: Objects/Bloc2_02_seu_singler.rds")
message("DONE Bloc2 Script 02")