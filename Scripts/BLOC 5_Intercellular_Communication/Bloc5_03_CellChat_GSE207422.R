#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 03: CellChat v2 intercellular communication
# Differential communication analysis MPR vs NMPR
# Input:  Objects/Bloc5_01_seu_TME_GSE207422.rds
# Output: Results/Figures/BLOC5_Communication/GSE207422/
#         Results/Tables/BLOC5/
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
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE207422")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5")

# ============================================================
# 1) Load TME object
# ============================================================
message("Loading TME object...")
seu_TME <- readRDS(file.path(DATA_DIR, "Objects/Bloc5_01_seu_TME_GSE207422.rds"))
message("Total cells: ", ncol(seu_TME))

seu_TME <- subset(seu_TME, subset = cell_type %in% c(
  "CD8.TEX", "CD8.TPEX", "CD8.EM", "CD8.CM",
               "CD8.TEMRA", "CD8.NaiveLike",
  "MRC1+ M2-like", "SPP1+ immunosuppressive", "IFN-stimulated",
  "M2-SIGLEC8+", "Resident M2", "Monocyte-derived", "Lipid-associated",
  "Stress-response", "Regulatory",
  "tumor_MPR", "normal_MPR", "tumor_NMPR", "normal_NMPR", "Ciliated"
))

# ============================================================
# 2) Load CellChatDB
# ============================================================
CellChatDB <- CellChatDB.human
message("CellChatDB interactions: ", nrow(CellChatDB$interaction))

# ============================================================
# 3) Function to create CellChat object per condition
# ============================================================
create_cellchat <- function(seu_obj, condition, cell_type_col = "cell_type") {
  message("Creating CellChat object for: ", condition)
  
  seu_sub <- subset(seu_obj, subset = PathResponse == condition)
  message("Cells: ", ncol(seu_sub))
  
  # Create CellChat object
  cellchat <- createCellChat(object   = seu_sub,
                             group.by = cell_type_col,
                             assay    = "RNA")
  
  # Set database
  cellchat@DB <- CellChatDB
  
  # Preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probabilities
  cellchat <- computeCommunProb(cellchat,
                                type    = "triMean",
                                nboot   = 100)
  
  # Filter min cells
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  # Compute pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Aggregate network
  cellchat <- aggregateNet(cellchat)
  
  message("Done: ", condition)
  return(cellchat)
}

# ============================================================
# 4) Create CellChat objects per condition
# ============================================================
cellchat_MPR  <- create_cellchat(seu_TME, "MPR")
cellchat_NMPR <- create_cellchat(seu_TME, "NMPR")

# Save individual objects
saveRDS(cellchat_MPR,
        file.path(DATA_DIR, "Objects/Bloc5_03_cellchat_MPR.rds"))
saveRDS(cellchat_NMPR,
        file.path(DATA_DIR, "Objects/Bloc5_03_cellchat_NMPR.rds"))
message("Saved: individual CellChat objects")

# ============================================================
# 5) Merge and compare MPR vs NMPR
# ============================================================
message("Merging CellChat objects...")
object.list <- list(MPR = cellchat_MPR, NMPR = cellchat_NMPR)
cellchat    <- mergeCellChat(object.list, add.names = names(object.list))

# Save merged object
saveRDS(cellchat,
        file.path(DATA_DIR, "Objects/Bloc5_03_cellchat_merged_GSE207422.rds"))
message("Saved: merged CellChat object")

# ============================================================
# 6) Save interaction tables
# ============================================================
df_MPR  <- subsetCommunication(cellchat_MPR)
df_NMPR <- subsetCommunication(cellchat_NMPR)

fwrite(df_MPR,  file.path(OUT_TAB, "Bloc5_03_CellChat_GSE207422_MPR_interactions.csv"))
fwrite(df_NMPR, file.path(OUT_TAB, "Bloc5_03_CellChat_GSE207422_NMPR_interactions.csv"))
message("Saved: interaction tables")

# ============================================================
# 7) Define cell type groups
# ============================================================
cd8_types <- c("CD8.TEX", "CD8.TPEX", "CD8.EM", "CD8.CM",
               "CD8.TEMRA", "CD8.NaiveLike")
tam_types <- c("MRC1+ M2-like", "SPP1+ immunosuppressive", "Resident M2",
               "IFN-stimulated", "Monocyte-derived", "Lipid-associated",
               "Stress-response", "Regulatory", "M2-SIGLEC8+")
epi_types  <- c("tumor_MPR", "normal_MPR", "tumor_NMPR", "normal_NMPR", "Ciliated")

# ============================================================
# 8) Bubble plots — 3 main biological axes
# ============================================================
message("Generating bubble plots...")

# Get common cell types between MPR and NMPR
ct_MPR  <- levels(cellchat@idents$MPR)
ct_NMPR <- levels(cellchat@idents$NMPR)
ct_common <- intersect(ct_MPR, ct_NMPR)

# Filter groups to common cell types only
cd8_types_common <- intersect(cd8_types, ct_common)
tam_types_common <- intersect(tam_types, ct_common)
epi_types_common <- intersect(epi_types, ct_common)

message("Common CD8: ", paste(cd8_types_common, collapse=", "))
message("Common TAMs: ", paste(tam_types_common, collapse=", "))
message("Common Epi: ", paste(epi_types_common, collapse=", "))

# Figure 1 — CD8 -> TAMs
png(file.path(OUT_FIG, "Bloc5_03_CellChat_bubbleplot_CD8_TAMs.png"),
    width = 16, height = 10, units = "in", res = 300)
netVisual_bubble(cellchat,
                 sources.use = cd8_types_common,
                 targets.use = tam_types_common,
                 comparison  = c(1, 2),
                 angle.x     = 45,
                 title.name  = "CD8 → TAMs — MPR vs NMPR (GSE207422)")
dev.off()
message("Saved: CD8 -> TAMs")

# Figure 2 — CD8 -> Epithelial
png(file.path(OUT_FIG, "Bloc5_03_CellChat_bubbleplot_CD8_Epithelial.png"),
    width = 16, height = 10, units = "in", res = 300)
netVisual_bubble(cellchat,
                 sources.use = cd8_types_common,
                 targets.use = epi_types_common,
                 comparison  = c(1, 2),
                 angle.x     = 45,
                 title.name  = "CD8 → Epithelial — MPR vs NMPR (GSE207422)")
dev.off()
message("Saved: CD8 -> Epithelial")

# Figure 3 — TAMs -> Epithelial
png(file.path(OUT_FIG, "Bloc5_03_CellChat_bubbleplot_TAMs_Epithelial.png"),
    width = 16, height = 10, units = "in", res = 300)
netVisual_bubble(cellchat,
                 sources.use = tam_types_common,
                 targets.use = epi_types_common,
                 comparison  = c(1, 2),
                 angle.x     = 45,
                 title.name  = "TAMs → Epithelial — MPR vs NMPR (GSE207422)")
dev.off()
message("Saved: TAMs -> Epithelial")

# ============================================================
# 9) Supplementary figures
# ============================================================
message("Generating supplementary figures...")

# Barplot — number of interactions
png(file.path(OUT_FIG, "Bloc5_03_CellChat_barplot_nInteractions.png"),
    width = 10, height = 6, units = "in", res = 300)
compareInteractions(cellchat,
                    show.legend = FALSE,
                    group       = c(1, 2),
                    measure     = "count")
dev.off()
message("Saved: barplot nInteractions")

# Barplot — interaction weights
png(file.path(OUT_FIG, "Bloc5_03_CellChat_barplot_weight.png"),
    width = 10, height = 6, units = "in", res = 300)
compareInteractions(cellchat,
                    show.legend = FALSE,
                    group       = c(1, 2),
                    measure     = "weight")
dev.off()
message("Saved: barplot weight")

# Heatmap MPR
png(file.path(OUT_FIG, "Bloc5_03_CellChat_heatmap_MPR.png"),
    width = 14, height = 10, units = "in", res = 300)
netVisual_heatmap(cellchat_MPR,
                  color.heatmap = "Blues",
                  title.name    = "Number of interactions — MPR")
dev.off()
message("Saved: heatmap MPR")

# Heatmap NMPR
png(file.path(OUT_FIG, "Bloc5_03_CellChat_heatmap_NMPR.png"),
    width = 14, height = 10, units = "in", res = 300)
netVisual_heatmap(cellchat_NMPR,
                  color.heatmap = "Blues",
                  title.name    = "Number of interactions — NMPR")
dev.off()
message("Saved: heatmap NMPR")

message("DONE CellChat GSE207422")