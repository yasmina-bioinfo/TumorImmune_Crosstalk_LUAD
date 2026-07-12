#!/usr/bin/env Rscript
# ============================================================
# GSE207422 — Bloc2 Script 04: Consensus annotation
# Crosses manual (top15/50 markers) and SingleR annotations.
# sctype excluded from consensus for immune cell types: its Lung
# reference lacked sufficient resolution, broadly mislabeling
# immune clusters with epithelial identities. Retained only as a
# reference column for structural/epithelial clusters, where all
# three methods agreed (clusters 11, 12, 16).
# Input:  Objects/Bloc2_GSE207422_03_seu_sctype.rds
# Output: Results/Tables/Bloc2_GSE207422_consensus_annotation.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc2_GSE207422_03_seu_sctype.rds")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

message("Loading Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

get_majority_label <- function(seu, column) {
  df <- data.frame(
    cluster = seu$seurat_clusters,
    label   = seu@meta.data[[column]]
  )
  df %>%
    group_by(cluster) %>%
    count(label) %>%
    slice_max(n, n = 1) %>%
    summarise(label = first(label)) %>%
    ungroup()
}

# Manual annotation, from top15 markers analysis, cross-checked against
# SingleR fine for clusters 15 and 18 (see Progress log)
manual <- data.frame(
  cluster = 0:18,
  manual_annotation = c(
    "CD8 T cells (effector/cytotoxic)",              # 0
    "Neutrophils",                                    # 1
    "Monocytes (inflammatory)",                       # 2
    "B cells",                                        # 3
    "CD8 T cells (exhausted/cytotoxic)",              # 4
    "CD4 T cells (naive/memory)",                     # 5
    "Tregs (tumor-infiltrating)",                     # 6
    "Tumor epithelial cells (squamous)",              # 7
    "NK cells (cytotoxic)",                           # 8
    "Proliferating cells (cycling)",                  # 9
    "Plasma cells",                                   # 10
    "Epithelial cells (AT2, normal)",                 # 11
    "Alveolar macrophages (tissue-resident)",         # 12
    "Dendritic cells (cDC2)",                         # 13
    "Mast cells",                                     # 14
    "Stromal cells (mesenchymal/endothelial, unresolved)", # 15
    "Ciliated epithelial cells",                      # 16
    "Plasmacytoid dendritic cells (pDCs)",            # 17
    "T cells (mixed CD4/CD8, low purity)"             # 18
  ),
  manual_confidence = c(
    "high", "high", "high", "high", "high",
    "high", "high", "high", "high", "high",
    "high", "high", "high", "high", "high",
    "medium", "high", "high", "medium"
  )
)

singler <- get_majority_label(seu, "SingleR_main")
colnames(singler) <- c("cluster", "SingleR_main")

sctype_ann <- get_majority_label(seu, "sctype")
colnames(sctype_ann) <- c("cluster", "sctype")

singler$cluster    <- as.integer(as.character(singler$cluster))
sctype_ann$cluster <- as.integer(as.character(sctype_ann$cluster))

message("Building consensus table...")

consensus <- data.frame(cluster = 0:18) %>%
  left_join(manual,      by = "cluster") %>%
  left_join(singler,     by = "cluster") %>%
  left_join(sctype_ann,  by = "cluster")

consensus$final_annotation <- consensus$manual_annotation
consensus$final_confidence <- consensus$manual_confidence
consensus$notes <- NA
consensus$notes[consensus$cluster == 18] <- "Revised from initial epithelial call after SingleR fine cross-check showed majority CD4/CD8 lymphocyte composition"
consensus$notes[consensus$cluster == 15] <- "sctype and SingleR main disagree (Fibroblasts vs Endothelial); SingleR fine shows mixed mesenchymal/endothelial composition, no single dominant subtype"

print(consensus)

fwrite(consensus, file.path(OUT_TAB, "Bloc2_GSE207422_consensus_annotation.csv"))
message("Saved: Bloc2_GSE207422_consensus_annotation.csv")
message("DONE Bloc2 GSE207422 Script 04")