#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08: SCEVAN, CNV inference on epithelial cells
# Identifies malignant vs normal epithelial cells based on copy number variations
# Ciliated_epithelial used as normal reference (diploid, non-malignant)
# Replaces CopyKAT (killed at step 4 due to RAM constraints)
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_08_SCEVAN_*.png
#         Results/Tables/Bloc4B_08_SCEVAN_predictions.csv
#         Objects/Bloc4B_08_seu_Epithelial_SCEVAN.rds
# Reference: De Falco et al., Nature Communications 2023
#            doi:10.1038/s41467-023-36790-9
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SCEVAN)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# Paths
#DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4B_07_seu_Epithelial.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load epithelial object
message("Loading epithelial object...")
seu_Epi <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_Epi))
DefaultAssay(seu_Epi) <- "RNA"

# 2) Prepare count matrix
# SCEVAN works with sparse matrices, no dense conversion required
message("Preparing count matrix...")
count_matrix <- GetAssayData(seu_Epi, layer = "counts")
message("Matrix dimensions: ", nrow(count_matrix), " genes x ", ncol(count_matrix), " cells")

# 3) Define normal reference cells
# Ciliated_epithelial = normal bronchial cells, expected diploid
normal_cells <- colnames(seu_Epi)[seu_Epi$TME_cell_type == "Ciliated_epithelial"]
message("Normal reference cells (Ciliated_epithelial): ", length(normal_cells))

# 4) Run SCEVAN
# norm_cell: vector of normal cell barcodes
# assay: count matrix (genes x cells)
# SCEVAN does not require dense matrix conversion
message("Running SCEVAN...")
results <- pipelineCNA(
  count_matrix,
  norm_cell    = normal_cells,
  par_cores    = 1,
  SUBCLONES    = FALSE,
  plotTree     = FALSE
)
message("SCEVAN done.")

# 5) Extract predictions
# SCEVAN output: "tumor" = malignant, "normal" = diploid
message("SCEVAN predictions:")
print(table(results$class))

# 6) Add predictions to Seurat object
seu_Epi$scevan_class <- results$class[match(colnames(seu_Epi), rownames(results))]
message("Predictions by cell type:")
print(table(seu_Epi$scevan_class, seu_Epi$TME_cell_type))
print(table(seu_Epi$scevan_class, seu_Epi$PathResponse))

# 7) Save predictions table
pred_df <- data.frame(
  cell = colnames(seu_Epi),
  scevan_class = seu_Epi$scevan_class,
  cell_type = seu_Epi$TME_cell_type,
  response = seu_Epi$PathResponse
)
fwrite(pred_df, file.path(OUT_TAB, "Bloc4B_08_SCEVAN_predictions.csv"))
message("Saved: Bloc4B_08_SCEVAN_predictions.csv")

# 8) UMAP colored by SCEVAN prediction
png(file.path(OUT_FIG, "Bloc4B_08_UMAP_SCEVAN_class.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_Epi,
              reduction = "umap",
              group.by  = "scevan_class",
              raster    = FALSE) +
        scale_color_manual(values = c("tumor"  = "#D73027",
                                      "normal" = "#4393C3")) +
        theme_bw() +
        theme(legend.text = element_text(size = 12)) +
        labs(title = "SCEVAN predictions — Epithelial cells GSE207422"))
dev.off()
message("Saved: Bloc4B_08_UMAP_SCEVAN_class.png")

# 9) Barplot malignant vs normal by response
df_pred <- seu_Epi@meta.data %>%
  filter(!is.na(scevan_class)) %>%
  group_by(PathResponse, scevan_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PathResponse) %>%
  mutate(prop = n / sum(n))

png(file.path(OUT_FIG, "Bloc4B_08_Barplot_SCEVAN_response.png"),
    width = 8, height = 6, units = "in", res = 300)
print(ggplot(df_pred, aes(x = PathResponse, y = prop, fill = scevan_class)) +
        geom_col(width = 0.8, color = "white") +
        scale_fill_manual(values = c("tumor"  = "#D73027",
                                     "normal" = "#4393C3")) +
        scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
        theme_classic() +
        theme(axis.title.x  = element_blank(),
              axis.text.x   = element_text(size = 13, face = "bold"),
              legend.title  = element_blank()) +
        ylab("Proportion") +
        labs(title = "Malignant vs normal epithelial cells by response — GSE207422"))
dev.off()
message("Saved: Bloc4B_08_Barplot_SCEVAN_response.png")

# 10) Save updated object
message("Saving updated epithelial object...")
saveRDS(seu_Epi, file.path(OUT_OBJ, "Bloc4B_08_seu_Epithelial_SCEVAN.rds"))
message("Saved: Objects/Bloc4B_08_seu_Epithelial_SCEVAN.rds")
message("DONE Bloc4B Script 08 SCEVAN")