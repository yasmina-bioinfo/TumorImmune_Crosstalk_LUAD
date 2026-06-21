#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08: CopyKAT, CNV inference on epithelial cells
# Identifies malignant vs normal epithelial cells based on copy number variations
# Ciliated_epithelial used as normal reference (diploid, non-malignant)
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_08_CopyKAT_*.png
#         Results/Tables/Bloc4B_08_CopyKAT_predictions.csv
#         Objects/Bloc4B_08_seu_Epithelial_CopyKAT.rds
# NOTE: CopyKAT killed at step 4/7 due to RAM constraints (WSL 13GB, Colab 12GB)
# Replaced by SCEVAN, see Bloc4B_08_SCEVAN.R
# Reference: Gao et al., Nature Genetics 2021, doi:10.1038/s41588-021-00913-z
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# Paths
# Before: "C:/Users/yasmi/OneDrive/Desktop/..." in RStudio on Windows
# After:  "/mnt/c/Users/yasmi/OneDrive/Desktop/..." in WSL2 Ubuntu
DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4B_07_seu_Epithelial.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load epithelial object
message("Loading epithelial object...")
seu_Epi <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_Epi))

# 2) Prepare count matrix for CopyKAT
# NOTE: CopyKAT requires raw counts (not normalized)
message("Preparing count matrix...")
DefaultAssay(seu_Epi) <- "RNA"
raw_counts <- GetAssayData(seu_Epi, layer = "counts")
# FIX: Removed as.matrix() to preserve memory by keeping the native dgCMatrix format

# 3) Define normal reference cells
# Ciliated_epithelial = normal bronchial cells, expected diploid
normal_cells <- colnames(seu_Epi)[seu_Epi$TME_cell_type == "Ciliated_epithelial"]
message("Normal reference cells (Ciliated_epithelial): ", length(normal_cells))

# 4) Run CopyKAT
# NOTE: CopyKAT is computationally intensive — may take 1-2 hours
# ngene.chr = 5 : minimum genes per chromosome arm
# LOW.DR = 0.05, UP.DR = 0.1 : dropout rate thresholds
# win.size = 15 : smoothing window
# distance = "euclidean" : clustering distance
message("Running CopyKAT...")
message("This may take 1-2 hours...")

# FIX: CopyKAT automatically writes output plots to the current working directory.
# Temporarily switch to OUT_FIG so CopyKAT saves its heatmaps directly there.
old_wd <- getwd()
setwd(OUT_FIG)

copykat_result <- copykat(
  rawmat          = raw_counts,
  id.type         = "S",
  ngene.chr       = 5,
  win.size        = 15,
  KS.cut          = 0.1,
  sam.name        = "Bloc4B_GSE207422_Epithelial",
  distance        = "euclidean",
  norm.cell.names = normal_cells,
  output.seg      = FALSE,
  plot.genes      = TRUE,
  LOW.DR          = 0.05,
  UP.DR           = 0.1,
  n.cores         = 1
  # FIX: Must be 1 to avoid the mclapply error
)

# Restore the original working directory after execution
setwd(old_wd)
message("CopyKAT done.")

# 5) Extract predictions
# NOTE: copykat prediction: "aneuploid" = malignant, "diploid" = normal
predictions <- as.data.frame(copykat_result$prediction)
message("CopyKAT predictions:")
print(table(predictions$copykat.pred))

# 6) Add predictions to Seurat object
seu_Epi$copykat_pred <- predictions$copykat.pred[match(
  colnames(seu_Epi), predictions$cell.names)]
message("Predictions added to Seurat object")
print(table(seu_Epi$copykat_pred, seu_Epi$TME_cell_type))

# 7) Save predictions table
fwrite(predictions, file.path(OUT_TAB, "Bloc4B_08_CopyKAT_predictions.csv"))
message("Saved: Bloc4B_08_CopyKAT_predictions.csv")

# 8) UMAP colored by CopyKAT prediction
png(file.path(OUT_FIG, "Bloc4B_08_UMAP_CopyKAT_pred.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_Epi,
              reduction = "umap",
              group.by  = "copykat_pred",
              raster    = FALSE) +
        scale_color_manual(values = c("aneuploid" = "#D73027",
                                      "diploid"   = "#4393C3",
                                      "not.defined" = "grey80")) +
        theme_bw() +
        theme(legend.text = element_text(size = 12)) +
        labs(title = "CopyKAT predictions — Epithelial cells GSE207422"))
dev.off()
message("Saved: Bloc4B_08_UMAP_CopyKAT_pred.png")

# 9) Barplot malignant vs normal by response
df_pred <- seu_Epi@meta.data %>%
  filter(!is.na(copykat_pred)) %>%
  group_by(PathResponse, copykat_pred) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PathResponse) %>%
  mutate(prop = n / sum(n))

png(file.path(OUT_FIG, "Bloc4B_08_Barplot_CopyKAT_response.png"),
    width = 8, height = 6, units = "in", res = 300)
print(ggplot(df_pred, aes(x = PathResponse, y = prop, fill = copykat_pred)) +
        geom_col(width = 0.8, color = "white") +
        scale_fill_manual(values = c("aneuploid"   = "#D73027",
                                     "diploid"     = "#4393C3",
                                     "not.defined" = "grey80")) +
        scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
        theme_classic() +
        theme(axis.title.x  = element_blank(),
              axis.text.x   = element_text(size = 13, face = "bold"),
              legend.title  = element_blank()) +
        ylab("Proportion") +
        labs(title = "Malignant vs normal epithelial cells by response — GSE207422"))
dev.off()
message("Saved: Bloc4B_08_Barplot_CopyKAT_response.png")

# 10) Save updated object
message("Saving updated epithelial object...")
saveRDS(seu_Epi, file.path(OUT_OBJ, "Bloc4B_08_seu_Epithelial_CopyKAT.rds"))
message("Saved: Objects/Bloc4B_08_seu_Epithelial_CopyKAT.rds")
message("DONE Bloc4B Script 08")