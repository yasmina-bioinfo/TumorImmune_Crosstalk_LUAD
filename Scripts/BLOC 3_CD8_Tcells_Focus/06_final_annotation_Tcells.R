#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 06: T cells final annotation
# Fills final_annotation based on manual vs ProjecTILs concordance
# Input:  Results/Tables/Bloc3_Tcells_consensus_annotation.csv
# Output: Results/Tables/Bloc3_Tcells_final_annotation.csv
# NOTE: ProjecTILs labels validated only for confirmed CD8 clusters
#       Non-CD8 clusters rely on canonical marker annotation
#       See README for methodological note on workflow order
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load consensus table
message("Loading consensus annotation table...")
consensus_T <- fread(file.path(OUT_TAB, "Bloc3_Tcells_consensus_annotation.csv"))

# 2) Fill final annotation
# NOTE: ProjecTILs run before canonical marker annotation (methodological limitation)
# ProjecTILs scGate incompletely filtered non-CD8 cells, artifactual CD8 labels
# for non-CD8 clusters (CD4, Tregs, NKT). Final annotation relies on canonical
# markers for these clusters. ProjecTILs validated only for confirmed CD8 clusters.
# Recommended future workflow: subset confirmed CD8 clusters first, then run ProjecTILs

consensus_T$final_annotation <- c(
  "CD4 T cells (naive/central memory)",           # 0  — manual (ProjecTILs CD8.CM artifactual)
  "CD8 T cells (effector memory/GZMK+)",          # 1  — manual + ProjecTILs CD8.EM confirmed
  "CD8 T cells (exhausted/cytotoxic TRM)",        # 2  — manual + ProjecTILs CD8.TEX confirmed
  "Tregs (tumor-infiltrating)",                    # 3  — manual (ProjecTILs CD8.CM artifactual)
  "CD4 Tfh / exhausted CD4",                      # 4  — manual (ProjecTILs CD8.CM artifactual)
  "CD8 T cells (TRM/transitional exhausted)",     # 5  — manual + ProjecTILs CD8.EM confirmed
  "CD8 T cells (TRM/quiescent)",                  # 6  — manual (ProjecTILs CD8.CM partial)
  "NKT / gamma-delta T cells",                    # 7  — manual (ProjecTILs CD8.CM artifactual)
  "CD4 T cells (naive/resting)",                  # 8  — manual (ProjecTILs CD8.CM artifactual)
  "CD4 T cells (activated/metabolically active)", # 9  — manual (ProjecTILs CD8.CM artifactual)
  "Proliferating T cells (cycling)",              # 10 — manual (ProjecTILs CD8.TEX mixed)
  "MAIT cells",                                    # 11 — manual + ProjecTILs CD8.MAIT confirmed
  "IFN-stimulated T cells (ISG-high)",            # 12 — manual (ProjecTILs CD8.EM partial)
  "gamma-delta T cells / innate T cells",         # 13 — manual (ProjecTILs CD8.CM artifactual)
  "EXCLUDE — B cells contamination",              # 14 — exclude
  "EXCLUDE — Artifact (2 cells)"                  # 15 — exclude
)

consensus_T$final_confidence <- c(
  "high", "high", "high", "high", "medium",
  "high", "medium", "high", "high", "medium",
  "high", "high", "high", "medium", "high",
  "high"
)

# 3) Add concordance notes
consensus_T$notes[consensus_T$cluster %in% c(1, 2, 5, 11)] <-
  "Manual + ProjecTILs concordant"
consensus_T$notes[consensus_T$cluster %in% c(0, 3, 4, 7, 8, 9, 13)] <-
  "Manual annotation — ProjecTILs label artifactual (non-CD8 passed scGate filter)"
consensus_T$notes[consensus_T$cluster %in% c(6, 10, 12)] <-
  "Manual annotation — ProjecTILs partial concordance"
consensus_T$notes[consensus_T$cluster == 14] <-
  "B cells contamination — exclude from CD8 analyses"
consensus_T$notes[consensus_T$cluster == 15] <-
  "Artifact — 2 cells only — exclude"

print(consensus_T)

# 4) Save final annotation
fwrite(consensus_T, file.path(OUT_TAB, "Bloc3_Tcells_final_annotation.csv"))
message("Saved: Bloc3_Tcells_final_annotation.csv")
message("DONE Bloc3 Script 06")