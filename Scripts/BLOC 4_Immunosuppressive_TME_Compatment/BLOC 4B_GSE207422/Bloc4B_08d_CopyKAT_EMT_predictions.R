#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08d: CopyKAT validation on EMT epithelial cells
# Tests whether Tumor_epithelial_EMT cells have aneuploid CNV profile
# Addresses SCEVAN classification of EMT as "normal"
# Normal reference: Ciliated_epithelial
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Tables/Bloc4B_08d_CopyKAT_EMT_predictions.csv
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
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs")

# 1) Load and subset EMT + Ciliated
message("Loading epithelial object...")
seu_Epi <- readRDS(IN_OBJ)
DefaultAssay(seu_Epi) <- "RNA"

cells_keep <- colnames(seu_Epi)[seu_Epi$TME_cell_type %in% 
                                  c("Tumor_epithelial_EMT", "Ciliated_epithelial")]
seu_sub <- subset(seu_Epi, cells = cells_keep)
message("EMT cells: ", sum(seu_sub$TME_cell_type == "Tumor_epithelial_EMT"))
message("Ciliated reference: ", sum(seu_sub$TME_cell_type == "Ciliated_epithelial"))
message("Total: ", ncol(seu_sub))

# 2) Prepare count matrix
message("Preparing count matrix...")
raw_counts <- as.matrix(GetAssayData(seu_sub, layer = "counts"))

# 3) Normal reference
normal_cells <- colnames(seu_sub)[seu_sub$TME_cell_type == "Ciliated_epithelial"]
message("Normal reference cells: ", length(normal_cells))

# 4) Run CopyKAT
message("Running CopyKAT on EMT cells...")
setwd(OUT_FIG)
copykat_result <- copykat(
  rawmat          = raw_counts,
  id.type         = "S",
  ngene.chr       = 5,
  win.size        = 15,
  KS.cut          = 0.1,
  sam.name        = "Bloc4B_CopyKAT_EMT",
  distance        = "euclidean",
  norm.cell.names = normal_cells,
  output.seg      = FALSE,
  plot.genes      = FALSE,
  LOW.DR          = 0.05,
  UP.DR           = 0.1,
  n.cores         = 1
)
message("CopyKAT done.")

# 5) Save predictions
predictions <- as.data.frame(copykat_result$prediction)
message("CopyKAT predictions:")
print(table(predictions$copykat.pred))

fwrite(predictions, file.path(OUT_TAB, "Bloc4B_08d_CopyKAT_EMT_predictions.csv"))
message("Saved: Bloc4B_08d_CopyKAT_EMT_predictions.csv")
message("DONE CopyKAT EMT validation")