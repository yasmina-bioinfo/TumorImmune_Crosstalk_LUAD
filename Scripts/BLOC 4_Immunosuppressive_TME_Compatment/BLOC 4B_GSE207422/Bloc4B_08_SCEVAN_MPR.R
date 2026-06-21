#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08a: SCEVAN, CNV inference MPR epithelial cells
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Tables/Bloc4B_08_SCEVAN_predictions_MPR.csv
# Reference: De Falco et al., Nature Communications 2023
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SCEVAN)
  library(dplyr)
  library(data.table)
})

DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4B_07_seu_Epithelial.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load and subset MPR + Ciliated
message("Loading epithelial object...")
seu_Epi <- readRDS(IN_OBJ)
DefaultAssay(seu_Epi) <- "RNA"

cells_MPR <- colnames(seu_Epi)[seu_Epi$PathResponse == "MPR" | 
                                seu_Epi$TME_cell_type == "Ciliated_epithelial"]
seu_sub <- subset(seu_Epi, cells = cells_MPR)
message("MPR cells: ", sum(seu_sub$PathResponse == "MPR"))
message("Ciliated reference: ", sum(seu_sub$TME_cell_type == "Ciliated_epithelial"))

# 2) Count matrix
count_matrix <- GetAssayData(seu_sub, layer = "counts")

# 3) Normal reference
normal_cells <- colnames(seu_sub)[seu_sub$TME_cell_type == "Ciliated_epithelial"]
message("Normal reference: ", length(normal_cells))

# 4) Run SCEVAN
message("Running SCEVAN MPR...")
setwd(file.path(OUT_FIG))
results <- pipelineCNA(
  count_matrix,
  norm_cell  = normal_cells,
  par_cores  = 1,
  SUBCLONES  = FALSE,
  plotTree   = FALSE,
  ClonalCN   = FALSE,
  output_dir = file.path(OUT_FIG, "SCEVAN_MPR")
)
message("SCEVAN MPR done.")

# 5) Save predictions
predictions <- data.frame(
  cell         = rownames(results),
  scevan_class = results$class,
  response     = "MPR"
)
print(table(predictions$scevan_class))
fwrite(predictions, file.path(OUT_TAB, "Bloc4B_08_SCEVAN_predictions_MPR.csv"))
message("Saved: Bloc4B_08_SCEVAN_predictions_MPR.csv")
message("DONE SCEVAN MPR")