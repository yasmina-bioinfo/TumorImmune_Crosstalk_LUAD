#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08f: CopyKAT on Tumor_epithelial_basal — NMPR
# Methodological revision: rerun per subtype x condition
# Normal reference: Ciliated_epithelial
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Tables/Bloc4B_08f_CopyKAT_Basal_NMPR_predictions.csv
# Reference: Gao et al., Nature Genetics 2021
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
  library(data.table)
})

DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4B_07_seu_Epithelial.rds")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/CopyKAT_NMPR")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

seu_Epi <- readRDS(IN_OBJ)
DefaultAssay(seu_Epi) <- "RNA"

cells_keep <- colnames(seu_Epi)[
  (seu_Epi$TME_cell_type == "Tumor_epithelial_basal" & seu_Epi$PathResponse == "NMPR") |
   seu_Epi$TME_cell_type == "Ciliated_epithelial"]
seu_sub <- subset(seu_Epi, cells = cells_keep)
message("Basal NMPR: ", sum(seu_sub$TME_cell_type == "Tumor_epithelial_basal"))
message("Ciliated reference: ", sum(seu_sub$TME_cell_type == "Ciliated_epithelial"))
message("Total: ", ncol(seu_sub))

raw_counts <- as.matrix(seu_sub[["RNA"]]$counts)
normal_cells <- colnames(seu_sub)[seu_sub$TME_cell_type == "Ciliated_epithelial"]

setwd(OUT_FIG)
copykat_result <- copykat(
  rawmat          = raw_counts,
  id.type         = "S",
  ngene.chr       = 5,
  win.size        = 15,
  KS.cut          = 0.1,
  sam.name        = "CopyKAT_Basal_NMPR",
  distance        = "euclidean",
  norm.cell.names = normal_cells,
  output.seg      = FALSE,
  plot.genes      = FALSE,
  LOW.DR          = 0.05,
  UP.DR           = 0.1,
  n.cores         = 1
)

predictions <- as.data.frame(copykat_result$prediction)
predictions$cell_type <- seu_sub$TME_cell_type[match(predictions$cell.names, colnames(seu_sub))]
predictions$response  <- seu_sub$PathResponse[match(predictions$cell.names, colnames(seu_sub))]
print(table(predictions$copykat.pred))
fwrite(predictions, file.path(OUT_TAB, "Bloc4B_08f_CopyKAT_Basal_NMPR_predictions.csv"))
message("DONE CopyKAT Basal NMPR")