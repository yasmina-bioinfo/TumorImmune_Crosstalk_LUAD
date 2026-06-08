#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc3 Script 02: Harmony + Clustering + UMAP (T cells)
# Input:  Objects/Bloc3_01_seu_Tcells_pca.rds
# Output: Objects/Bloc3_02_seu_Tcells_clustered.rds
#         Results/Figures/Tcells/Bloc3_UMAP_Tcells_*.png
# NOTE: dims = 1:30 based on ElbowPlot inspection (Script 01)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_01_seu_Tcells_pca.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/Tcells")

# 1) Load T cells PCA object
message("Loading T cells PCA object...")
seu_T <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_T))

# 2) Harmony batch correction
# NOTE: dims.use = 1:30, based on ElbowPlot inspection (Script 01)
# T cells subset shows stabilization around PC30
message("Running Harmony...")
seu_T <- RunHarmony(seu_T,
                    group.by.vars = "sampleID",
                    dims.use      = 1:30,
                    verbose       = FALSE)

# 3) Find neighbors
message("Finding neighbors...")
seu_T <- FindNeighbors(seu_T,
                       reduction = "harmony",
                       dims      = 1:30,
                       verbose   = FALSE)

# 4) Clustering
# NOTE: resolution = 0.4, lower than TME global (0.5)
# T cells subset is more homogeneous, conservative resolution for first pass
# Will be refined after ProjecTILs annotation (Script 03)
message("Clustering (resolution = 0.4)...")
seu_T <- FindClusters(seu_T,
                      resolution = 0.4,
                      verbose    = FALSE)
message("Number of T cell clusters: ", length(levels(seu_T$seurat_clusters)))

# 5) UMAP
message("Running UMAP...")
seu_T <- RunUMAP(seu_T,
                 reduction = "harmony",
                 dims      = 1:30,
                 verbose   = FALSE)

# 6) UMAP figures
message("Generating UMAP figures...")

# By new cluster
png(file.path(OUT_FIG, "Bloc3_UMAP_Tcells_clusters.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_T,
              reduction = "umap",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "T cells reclustering — GSE243013 LUAD"))
dev.off()

# By pathological response
png(file.path(OUT_FIG, "Bloc3_UMAP_Tcells_response.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_T,
              reduction = "umap",
              group.by  = "pathological_response",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "T cells by pathological response — GSE243013 LUAD"))
dev.off()

# By original TME cluster for traceability
png(file.path(OUT_FIG, "Bloc3_UMAP_Tcells_original_clusters.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_T,
              reduction = "umap",
              group.by  = "seurat_clusters",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "T cells : original TME clusters - GSE243013 LUAD"))
dev.off()

# By patient : Harmony check
png(file.path(OUT_FIG, "Bloc3_UMAP_Tcells_patient.png"),
    width = 12, height = 7, units = "in", res = 300)
print(DimPlot(seu_T,
              reduction = "umap",
              group.by  = "sampleID",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "T cells by patient (Harmony check) — GSE243013 LUAD"))
dev.off()

message("Saved: UMAP figures")

# 7) Distribution per response group
message("T cell distribution per response group:")
print(table(seu_T$seurat_clusters, seu_T$pathological_response))

# 8) Save clustered object
message("Saving T cells clustered object...")
saveRDS(seu_T, file.path(OUT_OBJ, "Bloc3_02_seu_Tcells_clustered.rds"))
message("Saved: Objects/Bloc3_02_seu_Tcells_clustered.rds")
message("DONE Bloc3 Script 02")