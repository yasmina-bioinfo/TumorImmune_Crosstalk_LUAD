#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc1 Script 04: Harmony + Clustering + UMAP
# Input:  Objects/Bloc1_03_seu_pca.rds
# Output: Objects/Bloc1_04_seu_clustered.rds
#         Results/Figures/Bloc1_UMAP_*.png
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc1_03_seu_pca.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures")

# -----------------------------
# 1) Load PCA object
# -----------------------------
message("Loading PCA Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

# -----------------------------
# 2) Harmony batch correction
# -----------------------------
# Harmony corrects batch effects across patients in PCA space
# group.by.vars = "sampleID": correct for inter-patient technical variation
# dims_use = 1:40: based on ElbowPlot inspection (gradual decline until ~PC40)
# Dataset complexity (immune + epithelial TME, 298567 cells) justifies 40 PCs
message("Running Harmony batch correction...")
seu <- RunHarmony(seu,
                  group.by.vars = "sampleID",
                  dims.use      = 1:40,
                  verbose       = FALSE)

# -----------------------------
# 3) Find neighbors (on Harmony embedding)
# -----------------------------
# NOTE: reduction = "harmony" — not "pca"
# All downstream steps use Harmony-corrected embedding
message("Finding neighbors...")
seu <- FindNeighbors(seu,
                     reduction = "harmony",
                     dims      = 1:40,
                     verbose   = FALSE)

# -----------------------------
# 4) Clustering
# -----------------------------
# resolution = 0.5 : conservative starting point for global TME annotation
# Lower resolution = fewer, broader clusters (good for first annotation pass)
# We will refine resolution after inspecting UMAP and marker genes
message("Clustering (resolution = 0.5)...")
seu <- FindClusters(seu,
                    resolution = 0.5,
                    verbose    = FALSE)
message("Number of clusters: ", length(levels(seu$seurat_clusters)))

# -----------------------------
# 5) UMAP (on Harmony embedding)
# -----------------------------
message("Running UMAP...")
seu <- RunUMAP(seu,
               reduction = "harmony",
               dims      = 1:40,
               verbose   = FALSE)

# -----------------------------
# 6) UMAP plots
# -----------------------------
message("Generating UMAP plots...")

# By cluster
p_clust <- DimPlot(seu,
                   reduction = "umap",
                   label     = TRUE,
                   repel     = TRUE) +
  theme_bw() +
  labs(title = "Global TME clustering — GSE243013 LUAD")

# By pathological response
p_response <- DimPlot(seu,
                      reduction = "umap",
                      group.by  = "pathological_response") +
  theme_bw() +
  labs(title = "By pathological response")

# Split by pathological response
p_split <- DimPlot(seu,
                   reduction = "umap",
                   split.by  = "pathological_response",
                   ncol      = 3) +
  theme_bw()

# By patient (to check Harmony integration)
p_patient <- DimPlot(seu,
                     reduction = "umap",
                     group.by  = "sampleID") +
  theme_bw() +
  labs(title = "By patient (Harmony check)")

# Save all plots
ggsave(file.path(OUT_FIG, "Bloc1_UMAP_clusters.png"),
       p_clust, width = 7, height = 5, dpi = 300, bg = "white")

ggsave(file.path(OUT_FIG, "Bloc1_UMAP_response.png"),
       p_response, width = 7, height = 5, dpi = 300, bg = "white")

ggsave(file.path(OUT_FIG, "Bloc1_UMAP_split_response.png"),
       p_split, width = 14, height = 5, dpi = 300, bg = "white")

ggsave(file.path(OUT_FIG, "Bloc1_UMAP_patient.png"),
       p_patient, width = 8, height = 5, dpi = 300, bg = "white")

message("Saved: UMAP figures")

# -----------------------------
# 7) Save clustered object
# -----------------------------
message("Saving clustered object...")
saveRDS(seu, file.path(OUT_OBJ, "Bloc1_04_seu_clustered.rds"))
message("Saved: Objects/Bloc1_04_seu_clustered.rds")
message("DONE Bloc1 Script 04")