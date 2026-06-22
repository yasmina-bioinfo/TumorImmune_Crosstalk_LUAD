#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc3_GSE207422 Script 01: ProjecTILs (CD8 annotated)
# Input:  C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy/objects/08_CD8_MPR_NMPR.rds
# Output: Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds
#         Results/Figures/CD8/Bloc3_GSE207422_UMAP_CD8_*.png
#         Results/Tables/Bloc3_GSE207422_CD8_ProjecTILs_per_cluster.csv
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
IN_OBJ   <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy/objects/08_CD8_MPR_NMPR.rds"
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load CD8 PCA object
message("Loading CD8 PCA object...")
seu_CD8 <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_CD8))

# 2) UMAP figures before ProjecTILs
message("Generating UMAP figures...")

png(file.path(OUT_FIG, "Bloc3_GSE207422_UMAP_CD8_clusters.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "CD8 T cells clustering — GSE207422 NSCLC"))
dev.off()

png(file.path(OUT_FIG, "Bloc3_GSE207422_UMAP_CD8_response.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              group.by  = "PathResponse",
              raster    = FALSE) +
        theme_bw() +
        scale_color_manual(values = c("MPR" = "#4393C3", "NMPR" = "#D73027")) +
        labs(title = "CD8 T cells by pathological response"))
dev.off()

png(file.path(OUT_FIG, "Bloc3_GSE207422_UMAP_CD8_patient.png"),
    width = 12, height = 7, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              group.by  = "Sample",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "CD8 T cells by patient (Harmony check)"))
dev.off()

message("Saved: UMAP CD8 figures")

# 3) ProjecTILs, corrected workflow
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

# 4) UMAP colored by ProjecTILs
png(file.path(OUT_FIG, "Bloc3_GSE207422_UMAP_CD8_ProjecTILs.png"),
    width = 12, height = 7, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              group.by  = "functional.cluster",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "ProjecTILs annotation — CD8 T cells GSE207422 NSCLC"))
dev.off()

png(file.path(OUT_FIG, "Bloc3_GSE207422_UMAP_CD8_ProjecTILs_split_response.png"),
    width = 16, height = 6, units = "in", res = 300)
print(DimPlot(seu_CD8,
              reduction = "umap",
              group.by  = "functional.cluster",
              split.by  = "PathResponse",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "ProjecTILs CD8 — split by pathological response"))
dev.off()

message("Saved: UMAP ProjecTILs CD8 figures")

# 5) Save final CD8 object
message("Saving CD8 ProjecTILs object...")
saveRDS(seu_CD8, file.path(OUT_OBJ, "Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds"))
message("Saved: Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds")
message("DONE Bloc3_GSE207422 script 01")