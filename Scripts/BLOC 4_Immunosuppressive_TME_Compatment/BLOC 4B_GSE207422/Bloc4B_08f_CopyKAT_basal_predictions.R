#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08f: CopyKAT on Tumor_epithelial_basal
# Both MPR and NMPR analyzed together
# Normal reference: Ciliated_epithelial
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Tables/Bloc4B_08f_CopyKAT_basal_predictions.csv
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

# 1) Load and subset basal + Ciliated
message("Loading epithelial object...")
seu_Epi <- readRDS(IN_OBJ)
DefaultAssay(seu_Epi) <- "RNA"

cells_keep <- colnames(seu_Epi)[seu_Epi$TME_cell_type %in% 
                                  c("Tumor_epithelial_basal", "Ciliated_epithelial")]
seu_sub <- subset(seu_Epi, cells = cells_keep)
message("Basal MPR: ", sum(seu_sub$TME_cell_type == "Tumor_epithelial_basal" & seu_sub$PathResponse == "MPR"))
message("Basal NMPR: ", sum(seu_sub$TME_cell_type == "Tumor_epithelial_basal" & seu_sub$PathResponse == "NMPR"))
message("Ciliated reference: ", sum(seu_sub$TME_cell_type == "Ciliated_epithelial"))
message("Total: ", ncol(seu_sub))

# 2) Prepare count matrix
message("Preparing count matrix...")
raw_counts <- as.matrix(GetAssayData(seu_sub, layer = "counts"))

# 3) Normal reference
normal_cells <- colnames(seu_sub)[seu_sub$TME_cell_type == "Ciliated_epithelial"]
message("Normal reference cells: ", length(normal_cells))

# 4) Run CopyKAT
message("Running CopyKAT on basal cells...")
setwd(OUT_FIG)
copykat_result <- copykat(
  rawmat          = raw_counts,
  id.type         = "S",
  ngene.chr       = 5,
  win.size        = 15,
  KS.cut          = 0.1,
  sam.name        = "Bloc4B_CopyKAT_basal",
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
predictions$cell_type <- seu_sub$TME_cell_type[match(predictions$cell.names, colnames(seu_sub))]
predictions$response  <- seu_sub$PathResponse[match(predictions$cell.names, colnames(seu_sub))]

message("CopyKAT predictions:")
print(table(predictions$copykat.pred, predictions$response))

fwrite(predictions, file.path(OUT_TAB, "Bloc4B_08f_CopyKAT_basal_predictions.csv"))
message("Saved: Bloc4B_08f_CopyKAT_basal_predictions.csv")
message("DONE CopyKAT Basal")