#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 05: T cells consensus annotation
# Crosses manual (top50 markers) and ProjecTILs annotations
# Objective: validate CD8 state annotations and reinforce
# biological credibility of portfolio findings (GSE207422)
# Input:  Objects/Bloc3_03_seu_Tcells_ProjecTILs.rds
# Output: Results/Tables/Bloc3_Tcells_canonical_annotation.csv
#         Results/Tables/Bloc3_Tcells_consensus_annotation.csv
# NOTE: Clusters 14 (B cells contamination) and 15 (artifact)
#       flagged for exclusion from downstream analyses
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_03_seu_Tcells_ProjecTILs.rds")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load ProjecTILs annotated object
message("Loading ProjecTILs annotated object...")
seu_T <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_T), " | Clusters: ", length(levels(seu_T$seurat_clusters)))

# 2) Manual annotation, hardcoded from top50 markers analysis
# NOTE: hardcoded to avoid Excel encoding issues (as in Bloc2 Script 05)
manual_T <- data.frame(
  cluster = 0:15,
  manual_annotation = c(
    "CD4 T cells (naive/central memory)",          # 0
    "CD8 T cells (effector memory/GZMK+)",         # 1
    "CD8 T cells (exhausted/cytotoxic TRM)",       # 2
    "Tregs (tumor-infiltrating)",                   # 3
    "CD4 Tfh / exhausted CD4",                     # 4
    "CD8 T cells (TRM/transitional exhausted)",    # 5
    "CD8 T cells (TRM/quiescent)",                 # 6
    "NKT / gamma-delta T cells",                   # 7
    "CD4 T cells (naive/resting)",                 # 8
    "CD4 T cells (activated/metabolically active)",# 9
    "Proliferating T cells (cycling)",             # 10
    "MAIT cells",                                   # 11
    "IFN-stimulated T cells (ISG-high)",           # 12
    "gamma-delta T cells / innate T cells",        # 13
    "B cells (contamination - to exclude)",        # 14
    "Artifact (2 cells - to exclude)"              # 15
  ),
  manual_confidence = c(
    "high", "high", "high", "high", "medium",
    "high", "high", "high", "high", "medium",
    "high", "high", "high", "medium", "high",
    "high"
  )
)

# Save canonical annotation table
fwrite(manual_T, file.path(OUT_TAB, "Bloc3_Tcells_canonical_annotation.csv"))
message("Saved: Bloc3_Tcells_canonical_annotation.csv")

# 3) Extract ProjecTILs majority label per cluster
message("Extracting ProjecTILs majority label per cluster...")

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

proj <- get_majority_label(seu_T, "functional.cluster")
colnames(proj) <- c("cluster", "ProjecTILs")

# NOTE: seurat_clusters is factor, convert to integer for join
proj$cluster <- as.integer(as.character(proj$cluster))

# 4) Build consensus table
message("Building consensus table...")

consensus_T <- data.frame(cluster = 0:15) %>%
  left_join(manual_T, by = "cluster") %>%
  left_join(proj, by = "cluster")

# Add columns for final annotation
consensus_T$final_annotation <- NA
consensus_T$final_confidence  <- NA
consensus_T$notes             <- NA

# NOTE: clusters 14 and 15 flagged for exclusion
consensus_T$notes[consensus_T$cluster == 14] <- "B cells contamination — exclude from CD8 analyses"
consensus_T$notes[consensus_T$cluster == 15] <- "Artifact — 2 cells only — exclude"

print(consensus_T)

# Save consensus table
fwrite(consensus_T, file.path(OUT_TAB, "Bloc3_Tcells_consensus_annotation.csv"))
message("Saved: Bloc3_Tcells_consensus_annotation.csv")
message("NOTE: Fill final_annotation column after reviewing manual vs ProjecTILs concordance")
message("DONE Bloc3 Script 05")