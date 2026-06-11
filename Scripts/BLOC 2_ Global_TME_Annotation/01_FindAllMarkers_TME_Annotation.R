#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc2 Script 01: FindAllMarkers (TME annotation)
# Input:  Objects/Bloc1_04_seu_clustered.rds
# Output: Results/Markers/Bloc2_markers_all_clusters.tsv
#         Results/Markers/Bloc2_top50_per_cluster.tsv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc1_04_seu_clustered.rds")
OUT_DIR  <- file.path(DATA_DIR, "Results/Markers")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load clustered object
# -----------------------------
message("Loading clustered Seurat object...")
seu <- readRDS(IN_OBJ)

if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("seurat_clusters absent. Run FindNeighbors + FindClusters first.")
}

Idents(seu) <- "seurat_clusters"
DefaultAssay(seu) <- "RNA"
# NOTE: after merge(), Seurat v5 creates multiple layers (one per sample)
# JoinLayers() merges them into a single "data" layer for FindAllMarkers()
seu <- JoinLayers(seu)
message("Layers joined.")
message("Cells: ", ncol(seu), 
        " | Genes: ", nrow(seu), 
        " | Clusters: ", length(levels(Idents(seu))))

# -----------------------------
# 2) FindAllMarkers
# -----------------------------
# only.pos = TRUE: keep only upregulated markers per cluster
# min.pct = 0.25: gene must be detected in at least 25% of cells in the cluster
# logfc.threshold = 0.25: minimum log2 fold-change
# test.use = "wilcox": Wilcoxon rank-sum test (standard, RAM-friendly)
# # NOTE: JoinLayers() required before FindAllMarkers() in Seurat v5
# Runtime: ~15-20 min on 298567 cells with 16GB RAM after joining layers
message("Running FindAllMarkers (this may take a while)...")
markers <- FindAllMarkers(
  object          = seu,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

# Harmonize gene column name
if (!"gene" %in% colnames(markers)) {
  markers <- markers %>% tibble::rownames_to_column("gene")
}

message("Total markers found: ", nrow(markers))

# Save all markers
write_tsv(markers, file.path(OUT_DIR, "Bloc2_markers_all_clusters.tsv"))
message("Saved: Bloc2_markers_all_clusters.tsv")

# -----------------------------
# 3) Top 50 markers per cluster
# -----------------------------
# Ranked by avg_log2FC (highest fold-change = most specific marker)
# 50 markers chosen empirically — sufficient signal, avoids noise
top50 <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  ungroup()

write_tsv(top50, file.path(OUT_DIR, "Bloc2_top50_per_cluster.tsv"))
message("Saved: Bloc2_top50_per_cluster.tsv")

message("DONE Bloc2 Script 01")