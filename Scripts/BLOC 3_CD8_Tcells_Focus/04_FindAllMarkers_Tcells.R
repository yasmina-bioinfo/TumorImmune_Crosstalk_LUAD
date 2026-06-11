#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 04: FindAllMarkers T cells
# Canonical marker-based annotation of T cell subclusters
# Objective: validate CD8 state annotations against canonical markers
# and cross-validate with ProjecTILs results (Script 05)
# This reinforces the biological credibility of portfolio findings (GSE207422)
# Input:  Objects/Bloc3_02_seu_Tcells_clustered.rds
# Output: Results/Markers/Bloc3_markers_Tcells_all.tsv
#         Results/Markers/Bloc3_top50_Tcells_per_cluster.tsv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_02_seu_Tcells_clustered.rds")
OUT_MRK  <- file.path(DATA_DIR, "Results/Markers")

# 1) Load T cells clustered object
message("Loading T cells clustered object...")
seu_T <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_T), " | Clusters: ", length(levels(seu_T$seurat_clusters)))

Idents(seu_T) <- "seurat_clusters"
DefaultAssay(seu_T) <- "RNA"

# 2) JoinLayers, required for Seurat v5 before FindAllMarkers
# NOTE: same requirement as Bloc2 Script 01
seu_T <- JoinLayers(seu_T)
message("Layers joined.")

# 3) FindAllMarkers
# NOTE: same parameters as Bloc2 Script 01 for methodological consistency
# only.pos = TRUE: upregulated markers only
# min.pct = 0.25: gene detected in at least 25% of cluster cells
# logfc.threshold = 0.25: minimum log2 fold-change
# test.use = "wilcox": Wilcoxon rank-sum test (standard, RAM-friendly)
message("Running FindAllMarkers on T cells...")
markers_T <- FindAllMarkers(
  object          = seu_T,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

if (!"gene" %in% colnames(markers_T)) {
  markers_T <- markers_T %>% tibble::rownames_to_column("gene")
}

message("Total markers found: ", nrow(markers_T))

# Save all markers
write_tsv(markers_T, file.path(OUT_MRK, "Bloc3_markers_Tcells_all.tsv"))
message("Saved: Bloc3_markers_Tcells_all.tsv")

# 4) Top 50 markers per cluster
# Ranked by avg_log2FC, most specific markers first
top50_T <- markers_T %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  ungroup()

write_tsv(top50_T, file.path(OUT_MRK, "Bloc3_top50_Tcells_per_cluster.tsv"))
message("Saved: Bloc3_top50_Tcells_per_cluster.tsv")

message("DONE Bloc3 Script 04 — inspect top50 markers before running Script 05")