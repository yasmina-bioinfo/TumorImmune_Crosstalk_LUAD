#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc4A Script 03: FindAllMarkers TAMs
# Canonical marker-based annotation of TAMs subclusters
# Objective: extract data-driven markers to profile TAMs heterogeneity
# and identify specific immunosuppressive mechanisms before CD8 crosstalk
# Input:  Objects/Bloc4A_02_seu_TAM_clustered.rds
# Output: Results/Markers/Bloc4A_markers_TAMs_all.tsv
#         Results/Markers/Bloc4A_top20_TAMs_per_cluster.tsv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4A_02_seu_TAM_clustered.rds")
OUT_MRK  <- file.path(DATA_DIR, "Results/Markers")

# 1) Load TAMs clustered object
message("Loading TAMs clustered object...")
seu_TAM <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_TAM), " | Clusters: ", length(levels(seu_TAM$seurat_clusters)))

Idents(seu_TAM) <- "seurat_clusters"
DefaultAssay(seu_TAM) <- "RNA"

# 2) JoinLayers, required for Seurat v5 before FindAllMarkers
# Ensures counts across layers are integrated for differential expression
seu_TAM <- JoinLayers(seu_TAM)
message("Layers joined.")

# 3) FindAllMarkers
# Non-biased approach to capture the distinct signatures of the 8 clusters
# only.pos = TRUE: upregulated markers only
# min.pct = 0.25: gene detected in at least 25% of cluster cells
# logfc.threshold = 0.25: minimum log2 fold-change
# test.use = "wilcox": Wilcoxon rank-sum test
message("Running FindAllMarkers on TAMs...")
markers_TAM <- FindAllMarkers(
  object          = seu_TAM,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

# Fix row names to column if needed (depending on Seurat output format)
if (!"gene" %in% colnames(markers_TAM)) {
  markers_TAM <- markers_TAM %>% tibble::rownames_to_column("gene")
}

message("Total markers found: ", nrow(markers_TAM))

# Save all markers
write_tsv(markers_TAM, file.path(OUT_MRK, "Bloc4A_markers_TAMs_all.tsv"))
message("Saved: Bloc4A_markers_TAMs_all.tsv")

# 4) Top 20 markers per cluster
# Ranked by avg_log2FC, most specific markers first for data-driven annotation
top20_TAM <- markers_TAM %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  ungroup()

write_tsv(top20_TAM, file.path(OUT_MRK, "Bloc4A_top20_TAMs_per_cluster.tsv"))
message("Saved: Bloc4A_top20_TAMs_per_cluster.tsv")

message("DONE Bloc4A Script 03 — ready to profile your 8 TAMs clusters!")