#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc2 Script 03: sctype automated annotation
# NOTE: Azimuth was initially planned but could not be installed
# on Windows due to heavy genomic dependencies:
# BSgenome.Hsapiens.UCSC.hg38, EnsDb.Hsapiens.v86
# sctype selected as lightweight alternative —
# marker-based automated annotation, no external reference required
# Input:  Objects/Bloc2_02_seu_singler.rds
# Output: Results/Tables/Bloc2_sctype_per_cluster.csv
#         Results/Figures/Annotations/Bloc2_UMAP_sctype.png
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(HGNChelper)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc2_02_seu_singler.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/Annotations")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# -----------------------------
# 1) Load sctype functions from GitHub
# -----------------------------
# sctype does not require installation, loaded directly from source
# Reference: Ianevski et al., Nature Communications 2022
message("Loading sctype functions...")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# -----------------------------
# 2) Load Seurat object
# -----------------------------
message("Loading Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

# -----------------------------
# 3) Load sctype marker database
# -----------------------------
# Tissue = "Lung" covers both immune and epithelial cell types
# Database downloaded directly from sctype GitHub repository
message("Loading sctype marker database (Lung tissue)...")
db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
db_path <- file.path(DATA_DIR, "Data/ScTypeDB_full.xlsx")

# Download only if not already present
if (!file.exists(db_path)) {
  download.file(db_url, db_path, mode = "wb")
  message("Database downloaded.")
} else {
  message("Database already present, skipping download.")
}

gs_list <- gene_sets_prepare(db_path, "Lung")

# -----------------------------
# 4) Compute sctype scores
# -----------------------------
# sctype scores each cluster based on expression of positive/negative markers
# Input: scaled data matrix (genes x cells)
message("Computing sctype scores...")

# NOTE: JoinLayers() required for Seurat v5 before extracting scaled data
seu <- JoinLayers(seu)
# NOTE: BPCells matrix not compatible with sctype, convert to standard matrix
# as.matrix() loads into RAM, feasible here as scale.data covers only 3000 HVGs
scaled_mat <- as.matrix(GetAssayData(seu, layer = "scale.data"))

sctype_scores <- sctype_score(
  scRNAseqData = scaled_mat,
  scaled       = TRUE,
  gs           = gs_list$gs_positive,
  gs2          = gs_list$gs_negative
)

# -----------------------------
# 5) Assign cell type per cluster
# -----------------------------
# For each cluster, the cell type with the highest cumulative score is assigned
message("Assigning cell types per cluster...")

cluster_results <- do.call("rbind", lapply(unique(seu$seurat_clusters), function(cl) {
  cl_cells <- which(seu$seurat_clusters == cl)
  cl_scores <- sort(rowSums(sctype_scores[, cl_cells, drop = FALSE]), decreasing = TRUE)
  data.frame(
    cluster   = cl,
    cell_type = names(cl_scores)[1],
    score     = cl_scores[1],
    ncells    = length(cl_cells)
  )
}))

cluster_results <- cluster_results[order(as.numeric(as.character(cluster_results$cluster))), ]
print(cluster_results)

# Save table
write.csv(cluster_results,
          file.path(OUT_TAB, "Bloc2_sctype_per_cluster.csv"),
          row.names = FALSE)
message("Saved: Bloc2_sctype_per_cluster.csv")

# -----------------------------
# 6) Add sctype annotation to Seurat object
# -----------------------------
# NOTE: names of sctype_map are cluster numbers, not barcodes
# Must reassign names to cell barcodes before adding to Seurat object
sctype_map <- setNames(as.character(cluster_results$cell_type), 
                       as.character(cluster_results$cluster))
test <- sctype_map[as.character(seu$seurat_clusters)]
names(test) <- colnames(seu)
seu$sctype <- test

# -----------------------------
# 7) UMAP colored by sctype annotation
# -----------------------------
message("Generating UMAP...")

png(file.path(OUT_FIG, "Bloc2_UMAP_sctype.png"),
    width = 14, height = 8, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "sctype",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "sctype annotation for GSE243013 LUAD"))
dev.off()
message("Saved: Bloc2_UMAP_sctype.png")

# UMAP with cluster numbers overlaid
png(file.path(OUT_FIG, "Bloc2_UMAP_sctype_labeled.png"),
    width = 14, height = 8, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "sctype",
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
        labs(title = "sctype annotation with cluster numbers — GSE243013 LUAD"))
dev.off()
message("Saved: Bloc2_UMAP_sctype_labeled.png")

# -----------------------------
# 8) Save updated Seurat object
# -----------------------------
message("Saving updated Seurat object...")
saveRDS(seu, file.path(DATA_DIR, "Objects/Bloc2_03_seu_sctype.rds"))
message("Saved: Objects/Bloc2_03_seu_sctype.rds")
message("DONE Bloc2 Script 03")