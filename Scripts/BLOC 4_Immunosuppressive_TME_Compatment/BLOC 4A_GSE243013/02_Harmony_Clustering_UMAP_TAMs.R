#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc4A Script 02: Harmony + Clustering + UMAP (T cells)
# Input:  Objects/Bloc4A_01_seu_TAM_pca.rds
# Output: Objects/Bloc4A_02_seu_TAM_clustered.rds
#         Results/Figures/BLOC4A_TAMs/Bloc4A_UMAP_TAMs_*.png
# NOTE: dims = 1:20 based on ElbowPlot inspection (Script 01)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4A_01_seu_TAMs_pca.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs")

# 1) Load TAMs PCA object
message("Loading TAMs PCA object...")
seu_T <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_T))

# 2) Harmony batch correction
# NOTE: dims.use = 1:20, based on ElbowPlot inspection (Script 01)
# TAMs subset shows stabilization around PC25
message("Running Harmony...")
seu_TAM <- RunHarmony(seu_TAM,
                    group.by.vars = "sampleID",
                    dims.use      = 1:20,
                    verbose       = FALSE)

# 3) Find neighbors
message("Finding neighbors...")
seu_TAM <- FindNeighbors(seu_TAM,
                       reduction = "harmony",
                       dims      = 1:20,
                       verbose   = FALSE)

# 4) Clustering
# NOTE: resolution = 0.1, lower than TME global (0.5) and lower than Tcells (0.4)
# resolution = 0.4 gave 19 clusters, too much, so lowering it was a better choice
# TAMs subset is more homogeneous, conservative resolution for first pass
message("Clustering (resolution = 0.1)...")
seu_TAM <- FindClusters(seu_TAM,
                      resolution = 0.1,
                      verbose    = FALSE)
message("Number of TAMs clusters: ", length(levels(seu_TAM$seurat_clusters)))

# 5) UMAP
message("Running UMAP...")
seu_TAM <- RunUMAP(seu_TAM,
                 reduction = "harmony",
                 dims      = 1:20,
                 verbose   = FALSE)

# 6) UMAP figures
message("Generating UMAP figures...")

# By new cluster
png(file.path(OUT_FIG, "Bloc4A_UMAP_TAMs_clusters.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_TAM,
              reduction = "umap",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "TAMs reclustering — GSE243013 LUAD"))
dev.off()

# By pathological response
png(file.path(OUT_FIG, "Bloc4_UMAP_TAMs_response.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_TAM,
              reduction = "umap",
              group.by  = "pathological_response",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "TAMs by pathological response — GSE243013 LUAD"))
dev.off()

# By original TME cluster for traceability
png(file.path(OUT_FIG, "Bloc4A_UMAP_TAMs_original_clusters.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_TAM,
              reduction = "umap",
              group.by  = "seurat_clusters",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "TAMs : original TME clusters - GSE243013 LUAD"))
dev.off()

# By patient : Harmony check
png(file.path(OUT_FIG, "Bloc4A_UMAP_TAMs_patient.png"),
    width = 12, height = 7, units = "in", res = 300)
print(DimPlot(seu_TAM,
              reduction = "umap",
              group.by  = "sampleID",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "TAMs by patient (Harmony check) — GSE243013 LUAD"))
dev.off()

message("Saved: UMAP figures")

# 7) Distribution per response group
message("TAMs distribution per response group:")
print(table(seu_TAM$seurat_clusters, seu_TAM$pathological_response))

# 8) Save clustered object
message("Saving TAMs clustered object...")
saveRDS(seu_TAM, file.path(OUT_OBJ, "Bloc4A_02_seu_TAM_clustered.rds"))
message("Saved: Objects/Bloc4_02_seu_TAM_clustered.rds")
message("DONE Bloc4A Script 02")