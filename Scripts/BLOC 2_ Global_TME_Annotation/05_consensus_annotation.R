#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc2 Script 05: Consensus annotation
# Crosses manual (top50 markers), SingleR and sctype annotations
# Input:  Objects/Bloc2_03_seu_sctype.rds
# Output: Results/Tables/Bloc2_consensus_annotation.csv
# NOTE: Final annotation column to be filled manually after reviewing
#       the three methods. Azimuth (Script 04) pending server execution.
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc2_03_seu_sctype.rds")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load Seurat object
message("Loading Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

# 2) Define helper function: extract majority label per cluster
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

# 3) Manual annotation, hardcoded from top50 markers analysis
# NOTE: cluster_annotation.csv was created in Excel but had encoding issues
# (invalid multibyte characters, Microsoft Excel 2007+ format saved as .csv)
# Hardcoded here for full reproducibility
manual <- data.frame(
  cluster = 0:19,
  manual_annotation = c(
    "B cells (memory/activated)",           # 0
    "CD4 T cells (naive/memory activated)", # 1
    "CD8 T cells (exhausted/cytotoxic TRM)",# 2
    "CD8 T cells (effector memory)",        # 3
    "Tregs (tumor-infiltrating)",           # 4
    "CD4 Tfh / exhausted CD4",             # 5
    "Monocytes (inflammatory)",             # 6
    "TAMs (M2-like/immunosuppressive)",     # 7
    "NK cells (cytotoxic)",                # 8
    "NKT / gamma-delta T cells",           # 9
    "CD4 T cells (naive/resting)",         # 10
    "Proliferating cells (cycling)",       # 11
    "B cells (naive/transitional)",        # 12
    "Mast cells",                          # 13
    "Dendritic cells (mregDC)",            # 14
    "Plasma cells",                        # 15
    "B cells (memory)",                    # 16
    "IFN-stimulated cells (ISG-high)",     # 17
    "Plasmacytoid Dendritic cells (pDCs)", # 18
    "Monocytes (classical CD14+)"          # 19
  ),
  manual_confidence = c(
    "high", "high", "high", "high", "high",
    "medium", "high", "high", "high", "medium",
    "high", "high", "high", "high", "high",
    "high", "high", "high", "high", "high"
  )
)

# 4) SingleR main: majority label per cluster
singler <- get_majority_label(seu, "SingleR_main")
colnames(singler) <- c("cluster", "SingleR_main")

# 5) sctype: majority label per cluster
sctype_ann <- get_majority_label(seu, "sctype")
colnames(sctype_ann) <- c("cluster", "sctype")

# 6) Convert cluster column from factor to integer for join compatibility
# NOTE: seurat_clusters is stored as factor. Must convert to integer
# to match manual annotation cluster column (integer 0:19)
singler$cluster    <- as.integer(as.character(singler$cluster))
sctype_ann$cluster <- as.integer(as.character(sctype_ann$cluster))

# 7) Build consensus table
message("Building consensus table...")

consensus <- data.frame(cluster = 0:19) %>%
  left_join(manual,      by = "cluster") %>%
  left_join(singler,     by = "cluster") %>%
  left_join(sctype_ann,  by = "cluster")

# Add empty columns for final annotation (to be filled manually in Excel)
consensus$final_annotation <- NA
consensus$final_confidence  <- NA
consensus$notes             <- NA

print(consensus)

# 8) Save
fwrite(consensus, file.path(OUT_TAB, "Bloc2_consensus_annotation.csv"))
message("Saved: Bloc2_consensus_annotation.csv")
message("NOTE: Fill final_annotation column manually in Excel before running next step.")
message("DONE Bloc2 Script 05")