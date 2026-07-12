#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc1 Script 04: Harmony + Clustering + UMAP
# NOTE: dims.use = 30, based on quantitative PC selection
# (cumulative variance >90% at PC30 = 92.22%, individual PC
# contribution <0.6% at PC30 = 0.56%), not visual elbow inspection.
# Harmony batch correction is applied here (unlike an earlier,
# now-superseded exploration of this dataset that skipped it),
# following confirmation of patient-driven batch effects in
# preliminary cluster composition checks.
# Input:  Objects/Bloc1_GSE207422_03_seu_pca.rds
# Output: Objects/Bloc1_GSE207422_04_seu_clustered.rds
#         Results/Figures/QC and Clustering/Bloc1_GSE207422_UMAP_*.png
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc1_GSE207422_03_seu_pca.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/QC and Clustering")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

message("Loading PCA Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

message("Running Harmony batch correction...")
seu <- RunHarmony(seu,
                  group.by.vars = "Sample",
                  dims.use      = 1:30,
                  verbose       = FALSE)

message("Finding neighbors...")
seu <- FindNeighbors(seu,
                     reduction = "harmony",
                     dims      = 1:30,
                     verbose   = FALSE)

message("Clustering (resolution = 0.5)...")
seu <- FindClusters(seu,
                    resolution = 0.5,
                    verbose    = FALSE)
message("Number of clusters: ", length(levels(seu$seurat_clusters)))

message("Running UMAP...")
seu <- RunUMAP(seu,
               reduction = "harmony",
               dims      = 1:30,
               verbose   = FALSE)

message("Generating UMAP plots...")

p_clust <- DimPlot(seu,
                   reduction = "umap",
                   label     = TRUE,
                   repel     = TRUE) +
  theme_bw() +
  labs(title = "Global TME clustering — GSE207422 NSCLC")

p_response <- DimPlot(seu,
                      reduction = "umap",
                      group.by  = "PathResponse") +
  theme_bw() +
  labs(title = "By pathological response")

p_split <- DimPlot(seu,
                   reduction = "umap",
                   split.by  = "PathResponse",
                   ncol      = 2) +
  theme_bw()

p_patient <- DimPlot(seu,
                     reduction = "umap",
                     group.by  = "Sample") +
  theme_bw() +
  labs(title = "By patient (Harmony check)")

ggsave(file.path(OUT_FIG, "Bloc1_GSE207422_UMAP_clusters.png"),
       p_clust, width = 7, height = 5, dpi = 300, bg = "white")

ggsave(file.path(OUT_FIG, "Bloc1_GSE207422_UMAP_response.png"),
       p_response, width = 7, height = 5, dpi = 300, bg = "white")

ggsave(file.path(OUT_FIG, "Bloc1_GSE207422_UMAP_split_response.png"),
       p_split, width = 10, height = 5, dpi = 300, bg = "white")

ggsave(file.path(OUT_FIG, "Bloc1_GSE207422_UMAP_patient.png"),
       p_patient, width = 8, height = 5, dpi = 300, bg = "white")

message("Saved: UMAP figures")

message("Saving clustered object...")
saveRDS(seu, file.path(OUT_OBJ, "Bloc1_GSE207422_04_seu_clustered.rds"))
message("Saved: Objects/Bloc1_GSE207422_04_seu_clustered.rds")
message("DONE Bloc1 GSE207422 Script 04")