#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08h: CopyKAT , Combine predictions and visualize
# Combines all CopyKAT predictions by cell type, filters Ciliated
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
#         Results/Tables/Bloc4B_08d_CopyKAT_EMT_predictions.csv
#         Results/Tables/Bloc4B_08e_CopyKAT_AT2_predictions.csv
#         Results/Tables/Bloc4B_08f_CopyKAT_basal_predictions.csv
#         Results/Tables/Bloc4B_08g_CopyKAT_Tumor_epithelial_predictions.csv
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_08h_UMAP_CopyKAT.png
#         Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_08h_Barplot_CopyKAT.png
#         Results/Tables/Bloc4B_08h_CopyKAT_predictions_combined.csv
#         Objects/Bloc4B_08h_seu_Epithelial_CopyKAT.rds
# Reference: Gao et al., Nature Genetics 2021
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4B_07_seu_Epithelial.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")

# 1) Load epithelial object
message("Loading epithelial object...")
seu_Epi <- readRDS(IN_OBJ)
message("Total cells: ", ncol(seu_Epi))

# 2) Load all predictions
pred_EMT   <- fread(file.path(OUT_TAB, "Bloc4B_08d_CopyKAT_EMT_predictions.csv"))
pred_AT2   <- fread(file.path(OUT_TAB, "Bloc4B_08e_CopyKAT_AT2_predictions.csv"))
pred_basal <- fread(file.path(OUT_TAB, "Bloc4B_08f_CopyKAT_basal_predictions.csv"))
pred_tumor <- fread(file.path(OUT_TAB, "Bloc4B_08g_CopyKAT_epithelial_predictions.csv"))

# Add missing columns to EMT and Tumor_epithelial
pred_EMT$cell_type <- "Tumor_epithelial_EMT"
pred_EMT$response  <- "NMPR"
pred_tumor$cell_type <- "Tumor_epithelial"
pred_tumor$response  <- "NMPR"

# 3) Combine and filter Ciliated
pred_all <- rbind(pred_EMT, pred_AT2, pred_basal, pred_tumor) %>%
  filter(cell_type != "Ciliated_epithelial")

# 4) Save combined predictions
fwrite(pred_all, file.path(OUT_TAB, "Bloc4B_08h_CopyKAT_predictions_combined.csv"))
message("Saved: Bloc4B_08h_CopyKAT_predictions_combined.csv")

# 5) Add to Seurat object
seu_Epi$copykat_pred <- pred_all$copykat.pred[match(
  colnames(seu_Epi), pred_all$cell.names)]
message("Predictions added to Seurat object")
message("By cell type:")
print(table(seu_Epi$copykat_pred, seu_Epi$TME_cell_type))

# 6) UMAP by CopyKAT prediction
png(file.path(OUT_FIG, "Bloc4B_08h_UMAP_CopyKAT.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_Epi,
              reduction = "umap",
              group.by  = "copykat_pred",
              raster    = FALSE) +
        scale_color_manual(values = c("aneuploid"    = "#D73027",
                                      "diploid"      = "#4393C3",
                                      "not.defined"  = "grey80")) +
        theme_bw() +
        theme(legend.text = element_text(size = 12)) +
        labs(title = "CopyKAT predictions — Epithelial cells GSE207422"))
dev.off()
message("Saved: Bloc4B_08h_UMAP_CopyKAT.png")

# 7) Barplot by response
df_bar <- pred_all %>%
  filter(copykat.pred %in% c("aneuploid", "diploid")) %>%
  group_by(response, copykat.pred) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(response) %>%
  mutate(prop = n / sum(n))

png(file.path(OUT_FIG, "Bloc4B_08h_Barplot_CopyKAT_response.png"),
    width = 8, height = 6, units = "in", res = 300)
print(ggplot(df_bar, aes(x = response, y = prop, fill = copykat.pred)) +
        geom_col(width = 0.8, color = "white") +
        scale_fill_manual(values = c("aneuploid" = "#D73027",
                                     "diploid"   = "#4393C3")) +
        scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
        theme_classic() +
        theme(axis.title.x = element_blank(),
              axis.text.x  = element_text(size = 13, face = "bold"),
              legend.title = element_blank()) +
        ylab("Proportion") +
        labs(title = "Malignant vs normal epithelial cells by response — GSE207422 (CopyKAT)"))
dev.off()
message("Saved: Bloc4B_08h_Barplot_CopyKAT_response.png")

# 8) Save updated object
saveRDS(seu_Epi, file.path(OUT_OBJ, "Bloc4B_08h_seu_Epithelial_CopyKAT.rds"))
message("Saved: Objects/Bloc4B_08h_seu_Epithelial_CopyKAT.rds")

# 9) Session info
writeLines(capture.output(sessionInfo()),
           file.path(DATA_DIR, "session_info_linux.txt"))
message("DONE Bloc4B Script 08h")