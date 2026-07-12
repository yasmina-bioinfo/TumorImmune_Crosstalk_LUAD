#!/usr/bin/env Rscript
# ============================================================
# GSE207422 — Bloc2 Script 03: sctype automated annotation
# Input:  Objects/Bloc2_GSE207422_02_seu_singler.rds
# Output: Results/Tables/Bloc2_GSE207422_sctype_per_cluster.csv
#         Results/Figures/Annotations_TME_GSE207422/Bloc2_GSE207422_UMAP_sctype*.png
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(HGNChelper)
  library(ggplot2)
  library(dplyr)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc2_GSE207422_02_seu_singler.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/Annotations_TME_GSE207422")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

message("Loading sctype functions...")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

message("Loading Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

message("Loading sctype marker database (Lung tissue)...")
db_path <- file.path(DATA_DIR, "Data/ScTypeDB_full.xlsx")

if (!file.exists(db_path)) {
  db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  download.file(db_url, db_path, mode = "wb")
  message("Database downloaded.")
} else {
  message("Database already present, skipping download.")
}

gs_list <- gene_sets_prepare(db_path, "Lung")

message("Computing sctype scores...")
seu <- JoinLayers(seu)
scaled_mat <- as.matrix(GetAssayData(seu, layer = "scale.data"))

sctype_scores <- sctype_score(
  scRNAseqData = scaled_mat,
  scaled       = TRUE,
  gs           = gs_list$gs_positive,
  gs2          = gs_list$gs_negative
)

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

write.csv(cluster_results,
          file.path(OUT_TAB, "Bloc2_GSE207422_sctype_per_cluster.csv"),
          row.names = FALSE)
message("Saved: Bloc2_GSE207422_sctype_per_cluster.csv")

sctype_map <- setNames(as.character(cluster_results$cell_type), 
                       as.character(cluster_results$cluster))
test <- sctype_map[as.character(seu$seurat_clusters)]
names(test) <- colnames(seu)
seu$sctype <- test

message("Generating UMAP...")

png(file.path(OUT_FIG, "Bloc2_GSE207422_UMAP_sctype.png"),
    width = 14, height = 8, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "sctype",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "sctype annotation for GSE207422 NSCLC"))
dev.off()
message("Saved: Bloc2_GSE207422_UMAP_sctype.png")

png(file.path(OUT_FIG, "Bloc2_GSE207422_UMAP_sctype_labeled.png"),
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
        labs(title = "sctype annotation with cluster numbers — GSE207422 NSCLC"))
dev.off()
message("Saved: Bloc2_GSE207422_UMAP_sctype_labeled.png")

message("Saving updated Seurat object...")
saveRDS(seu, file.path(DATA_DIR, "Objects/Bloc2_GSE207422_03_seu_sctype.rds"))
message("Saved: Objects/Bloc2_GSE207422_03_seu_sctype.rds")
message("DONE Bloc2 GSE207422 Script 03")