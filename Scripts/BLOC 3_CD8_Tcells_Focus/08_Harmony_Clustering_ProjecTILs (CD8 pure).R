#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 08: Harmony + Clustering + ProjecTILs (CD8 pure)
# Corrected workflow: ProjecTILs on confirmed CD8 clusters only
# dims = 1:25 based on ElbowPlot inspection (Script 07)
# Input:  Objects/Bloc3_07_seu_CD8_pca.rds
# Output: Objects/Bloc3_08_seu_CD8_ProjecTILs.rds
#         Results/Figures/CD8/Bloc3_UMAP_CD8_*.png
#         Results/Tables/Bloc3_CD8_ProjecTILs_per_cluster.csv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ProjecTILs)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_07_seu_CD8_pca.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load CD8 PCA object
message("Loading CD8 PCA object...")
seu_CD8 <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_CD8))

# 2) Harmony batch correction
# NOTE: dims.use = 1:25 based on ElbowPlot inspection (Script 07)
# CD8 subset more homogeneous than full T cell subset, stabilization at PC25
message("Running Harmony...")
seu_CD8 <- RunHarmony(seu_CD8,
                      group.by.vars = "sampleID",
                      dims.use      = 1:25,
                      verbose       = FALSE)

# 3) Find neighbors
message("Finding neighbors...")
seu_CD8 <- FindNeighbors(seu_CD8,
                         reduction = "harmony",
                         dims      = 1:25,
                         verbose   = FALSE)

# 4) Clustering
# NOTE: resolution = 0.3, CD8 subset is the most homogeneous population
# Lower resolution to avoid over-fragmenting continuous exhaustion gradient
message("Clustering (resolution = 0.3)...")
seu_CD8 <- FindClusters(seu_CD8,
                        resolution = 0.3,
                        verbose    = FALSE)
message("Number of CD8 clusters: ", length(levels(seu_CD8$seurat_clusters)))

# 5) UMAP
message("Running UMAP...")
seu_CD8 <- RunUMAP(seu_CD8,
                   reduction = "harmony",
                   dims      = 1:25,
                   verbose   = FALSE)

# 6) UMAP figures before ProjecTILs
message("Generating UMAP figures...")

png(file.path(OUT_FIG, "Bloc3_UMAP_CD8_clusters.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "CD8 T cells clustering — GSE243013 LUAD"))
dev.off()

png(file.path(OUT_FIG, "Bloc3_UMAP_CD8_response.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              group.by  = "pathological_response",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "CD8 T cells by pathological response"))
dev.off()

png(file.path(OUT_FIG, "Bloc3_UMAP_CD8_patient.png"),
    width = 12, height = 7, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              group.by  = "sampleID",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "CD8 T cells by patient (Harmony check)"))
dev.off()

message("Saved: UMAP CD8 figures")

# 7) ProjecTILs, corrected workflow
# NOTE: ProjecTILs now applied to confirmed CD8 cells only
# Expected: cleaner annotations, fewer artifactual labels
# options() set to avoid STACAS alignment failure (see Script 03)
message("Loading ProjecTILs reference...")
options(future.globals.maxSize = 1000 * 1024^2)
ref <- get.reference.maps()$human$CD8

message("Projecting CD8 cells onto reference atlas...")
seu_CD8 <- Run.ProjecTILs(seu_CD8,
                          ref          = ref,
                          filter.cells = TRUE)

message("ProjecTILs annotation done.")
message("Functional state distribution:")
print(table(seu_CD8$functional.cluster))

# 8) UMAP colored by ProjecTILs
png(file.path(OUT_FIG, "Bloc3_UMAP_CD8_ProjecTILs.png"),
    width = 12, height = 7, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              group.by  = "functional.cluster",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "ProjecTILs annotation — CD8 T cells GSE243013 LUAD"))
dev.off()

png(file.path(OUT_FIG, "Bloc3_UMAP_CD8_ProjecTILs_split_response.png"),
    width = 16, height = 6, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              group.by  = "functional.cluster",
              split.by  = "pathological_response",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "ProjecTILs CD8 — split by pathological response"))
dev.off()

message("Saved: UMAP ProjecTILs CD8 figures")

# 9) Summary table per cluster
proj_CD8 <- data.frame(
  cluster    = seu_CD8$seurat_clusters,
  ProjecTILs = seu_CD8$functional.cluster
) %>%
  group_by(cluster) %>%
  count(ProjecTILs) %>%
  slice_max(n, n = 1) %>%
  summarise(ProjecTILs_majority = first(ProjecTILs)) %>%
  ungroup()

proj_CD8$cluster <- as.integer(as.character(proj_CD8$cluster))
proj_CD8 <- proj_CD8[order(proj_CD8$cluster), ]
print(proj_CD8)

fwrite(proj_CD8, file.path(OUT_TAB, "Bloc3_CD8_ProjecTILs_per_cluster.csv"))
message("Saved: Bloc3_CD8_ProjecTILs_per_cluster.csv")

# 10) Save final CD8 object
message("Saving CD8 ProjecTILs object...")
saveRDS(seu_CD8, file.path(OUT_OBJ, "Bloc3_08_seu_CD8_ProjecTILs.rds"))
message("Saved: Objects/Bloc3_08_seu_CD8_ProjecTILs.rds")
message("DONE Bloc3 Script 08")