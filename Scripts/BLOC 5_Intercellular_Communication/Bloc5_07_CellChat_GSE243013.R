#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc5 Script 07: CellChat v2 intercellular communication
# Differential communication analysis MPR vs non-MPR
# No epithelial compartment available for this dataset
# Input:  Objects/Bloc5_03_seu_TME_GSE243013.rds
# Output: Results/Figures/BLOC5_Communication/GSE243013/
#         Results/Tables/BLOC5/GSE243013/
# Reference: Jin et al., Nature Communications 2021 (CellChat v2)
# Parameters: CellChatDB.human (full), nboot=100, min.cells=5
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(patchwork)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE243013")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5/GSE243013")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) Load TME object
# ============================================================
message("Loading TME object...")
seu_TME <- readRDS(file.path(DATA_DIR, "Objects/Bloc5_03_seu_TME_GSE243013.rds"))
message("Total cells: ", ncol(seu_TME))
message("Conditions: ", paste(unique(seu_TME$pathological_response), collapse=", "))

# Apply short labels to TAM types
short_tam_labels <- c(
  "Proliferating TAMs (cycling/MKI67+)"                                = "Proliferating",
  "Inflammatory monocyte-derived TAMs (FCN1+/S100A8+)"                 = "Inflammatory Mono-derived",
  "Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like)"  = "Resident M2",
  "TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs)"        = "Lipid-associated",
  "IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+)"      = "IFN-stimulated",
  "Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high)"    = "Stress-response",
  "Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+)"           = "Classical Mono-derived"
)

seu_TME$cell_type <- ifelse(seu_TME$cell_type %in% names(short_tam_labels),
                            short_tam_labels[seu_TME$cell_type],
                            seu_TME$cell_type)
message("Cell types: ", paste(unique(seu_TME$cell_type), collapse=", "))

# ============================================================
# 2) Load CellChatDB
# ============================================================
CellChatDB <- CellChatDB.human
message("CellChatDB interactions: ", nrow(CellChatDB$interaction))

# ============================================================
# 3) Function to create CellChat object per condition
# ============================================================
create_cellchat <- function(seu_obj, condition,
                            cell_type_col = "cell_type",
                            condition_col = "pathological_response") {
  message("Creating CellChat object for: ", condition)
  
  seu_sub <- subset(seu_obj, subset = pathological_response == condition)
  
  cellchat <- createCellChat(object   = seu_sub,
                             group.by = cell_type_col,
                             assay    = "RNA")
  
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = 100)
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  message("Done: ", condition)
  return(cellchat)
}

# ============================================================
# 4) Create CellChat objects — MPR, non-MPR, pCR
# ============================================================
cellchat_MPR    <- create_cellchat(seu_TME, "MPR")
cellchat_nonMPR <- create_cellchat(seu_TME, "non-MPR")
cellchat_pCR    <- create_cellchat(seu_TME, "pCR")

# Save individual objects
saveRDS(cellchat_MPR,
        file.path(DATA_DIR, "Objects/Bloc5_07_cellchat_MPR_GSE243013.rds"))
saveRDS(cellchat_nonMPR,
        file.path(DATA_DIR, "Objects/Bloc5_07_cellchat_nonMPR_GSE243013.rds"))
saveRDS(cellchat_pCR,
        file.path(DATA_DIR, "Objects/Bloc5_07_cellchat_pCR_GSE243013.rds"))
message("Saved: individual CellChat objects")

# ============================================================
# 5) Merge MPR vs non-MPR (main narrative)
# ============================================================
message("Merging MPR vs non-MPR...")
object.list_main <- list(MPR = cellchat_MPR, nonMPR = cellchat_nonMPR)
cellchat_main    <- mergeCellChat(object.list_main, 
                                  add.names = names(object.list_main))

saveRDS(cellchat_main,
        file.path(DATA_DIR, "Objects/Bloc5_07_cellchat_merged_MPR_nonMPR_GSE243013.rds"))

# Merge MPR vs pCR (discussion)
message("Merging MPR vs pCR...")
object.list_pcr <- list(MPR = cellchat_MPR, pCR = cellchat_pCR)
cellchat_pcr    <- mergeCellChat(object.list_pcr,
                                 add.names = names(object.list_pcr))

saveRDS(cellchat_pcr,
        file.path(DATA_DIR, "Objects/Bloc5_07_cellchat_merged_MPR_pCR_GSE243013.rds"))
message("Saved: merged CellChat objects")

# ============================================================
# 6) Save interaction tables
# ============================================================
df_MPR    <- subsetCommunication(cellchat_MPR)
df_nonMPR <- subsetCommunication(cellchat_nonMPR)
df_pCR    <- subsetCommunication(cellchat_pCR)

fwrite(df_MPR,
       file.path(OUT_TAB, "Bloc5_07_CellChat_GSE243013_MPR_interactions.csv"))
fwrite(df_nonMPR,
       file.path(OUT_TAB, "Bloc5_07_CellChat_GSE243013_nonMPR_interactions.csv"))
fwrite(df_pCR,
       file.path(OUT_TAB, "Bloc5_07_CellChat_GSE243013_pCR_interactions.csv"))
message("Saved: interaction tables")

message("DONE CellChat GSE243013")