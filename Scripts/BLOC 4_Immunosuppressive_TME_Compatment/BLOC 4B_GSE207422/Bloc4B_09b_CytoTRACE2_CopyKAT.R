#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 09b: CytoTRACE2 on epithelial cells
# Labels from CopyKAT predictions (cross-validation CNV tool)
# 6 groups: aneuploid_MPR, diploid_MPR, aneuploid_NMPR, diploid_NMPR,
#           EMT_NMPR (independent group), Ciliated (reference)
# Input:  Objects/Bloc4B_09a_seu_Epithelial_CytoTRACE2.rds
#         Results/Tables/Bloc4B_08_CopyKAT_predictions_MPR.csv
#         Results/Tables/Bloc4B_08h_CopyKAT_predictions_NMPR_combined.csv
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/CytoTRACE2_CopyKAT/
#         Results/Tables/Bloc4B_09b_CytoTRACE2_CopyKAT_scores.csv
# Reference: Kang et al., Nature Methods 2025
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CytoTRACE2)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/CytoTRACE2_CopyKAT")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# 1) Load epithelial object with CytoTRACE2 scores already computed
message("Loading epithelial object...")
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_09a_seu_Epithelial_CytoTRACE2.rds"))
message("Total epithelial cells: ", ncol(seu_Epi))

# 2) Load CopyKAT predictions
message("Loading CopyKAT predictions...")
pred_MPR  <- fread(file.path(OUT_TAB, "Bloc4B_08_CopyKAT_predictions_MPR.csv")) %>%
  filter(!cell.names %in% colnames(seu_Epi)[seu_Epi$TME_cell_type == "Ciliated_epithelial"])
pred_NMPR <- fread(file.path(OUT_TAB, "Bloc4B_08h_CopyKAT_predictions_NMPR_combined.csv"))

# 3) Create copykat_group column
message("Creating CopyKAT group labels...")
seu_Epi$copykat_group <- NA

# MPR predictions
idx_mpr <- match(pred_MPR$cell.names, colnames(seu_Epi))
idx_mpr <- idx_mpr[!is.na(idx_mpr)]
for (i in seq_along(idx_mpr)) {
  cell_idx  <- idx_mpr[i]
  pred_val  <- pred_MPR$copykat.pred[i]
  cell_type <- seu_Epi$TME_cell_type[cell_idx]
  if (!is.na(pred_val) && pred_val != "not.defined") {
    seu_Epi$copykat_group[cell_idx] <- paste0(pred_val, "_MPR")
  }
}

# NMPR predictions
idx_nmpr <- match(pred_NMPR$cell.names, colnames(seu_Epi))
idx_nmpr <- idx_nmpr[!is.na(idx_nmpr)]
for (i in seq_along(idx_nmpr)) {
  cell_idx  <- idx_nmpr[i]
  cell_type <- seu_Epi$TME_cell_type[cell_idx]
  pred_val  <- pred_NMPR$copykat.pred[i]
  if (cell_type == "Tumor_epithelial_EMT") {
    seu_Epi$copykat_group[cell_idx] <- "EMT_NMPR"
  } else if (!is.na(pred_val) && pred_val != "not.defined") {
    seu_Epi$copykat_group[cell_idx] <- paste0(pred_val, "_NMPR")
  }
}

# Ciliated cells
seu_Epi$copykat_group[seu_Epi$TME_cell_type == "Ciliated_epithelial"] <- "Ciliated"

message("Group distribution:")
print(table(seu_Epi$copykat_group, useNA = "always"))

# 4) CytoTRACE2 scores already in object — no need to rerun
# Use scores from Script 09a
message("CytoTRACE2 scores already computed in Script 09a — reusing.")
message("Score NA count: ", sum(is.na(seu_Epi$CytoTRACE2_Score)))

# 5) Save scores with CopyKAT labels
scores_df <- data.frame(
  cell               = colnames(seu_Epi),
  CytoTRACE2_Score   = seu_Epi$CytoTRACE2_Score,
  CytoTRACE2_Potency = seu_Epi$CytoTRACE2_Potency,
  copykat_group      = seu_Epi$copykat_group,
  response           = seu_Epi$PathResponse,
  cell_type          = seu_Epi$TME_cell_type
)
fwrite(scores_df, file.path(OUT_TAB, "Bloc4B_09b_CytoTRACE2_CopyKAT_scores.csv"))
message("Saved: Bloc4B_09b_CytoTRACE2_CopyKAT_scores.csv")

# 6) Visualizations
group_colors <- c(
  "aneuploid_MPR"  = "#E63946",
  "diploid_MPR"    = "#457B9D",
  "aneuploid_NMPR" = "#F4A261",
  "diploid_NMPR"   = "#2A9D8F",
  "EMT_NMPR"       = "#9B2226",
  "Ciliated"       = "#ADB5BD"
)

# UMAP
tryCatch({
  seu_plot <- subset(seu_Epi, cells = colnames(seu_Epi)[!is.na(seu_Epi$CytoTRACE2_Score)])
  p1 <- FeaturePlot(seu_plot, features = "CytoTRACE2_Score",
                    cols = c("darkblue", "white", "red")) +
    ggtitle("CytoTRACE2 Score — Epithelial cells GSE207422 (CopyKAT labels)")
  ggsave(file.path(OUT_FIG, "Bloc4B_09b_UMAP_CytoTRACE2_Score.png"),
         p1, width = 8, height = 6, dpi = 150)
  message("Saved: UMAP")
}, error = function(e) message("UMAP failed: ", e$message))

# Violin
tryCatch({
  seu_vln <- subset(seu_Epi, cells = colnames(seu_Epi)[!is.na(seu_Epi$copykat_group) &
                                                        !is.na(seu_Epi$CytoTRACE2_Score)])
  p2 <- VlnPlot(seu_vln, features = "CytoTRACE2_Score",
                group.by = "copykat_group",
                cols = group_colors, pt.size = 0) +
    ggtitle("CytoTRACE2 Score by CopyKAT group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_FIG, "Bloc4B_09b_Violin_CytoTRACE2_CopyKAT.png"),
         p2, width = 10, height = 6, dpi = 150)
  message("Saved: Violin")
}, error = function(e) message("Violin failed: ", e$message))

# Save updated Seurat object
saveRDS(seu_Epi, file.path(DATA_DIR, "Objects/Bloc4B_09b_seu_Epithelial_CytoTRACE2_CopyKAT.rds"))
message("Saved: Bloc4B_09b_seu_Epithelial_CytoTRACE2_CopyKAT.rds")

message("DONE CytoTRACE2 CopyKAT")