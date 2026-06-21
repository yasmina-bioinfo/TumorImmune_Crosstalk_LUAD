#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08c: SCEVAN , Combine predictions and visualize
# Combines MPR and NMPR SCEVAN predictions, generates UMAP and barplot
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
#         Results/Tables/Bloc4B_08_SCEVAN_predictions_MPR.csv
#         Results/Tables/Bloc4B_08_SCEVAN_predictions_NMPR.csv
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_08_UMAP_SCEVAN.png
#         Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_08_Barplot_SCEVAN.png
#         Results/Tables/Bloc4B_08_SCEVAN_predictions_combined.csv
#         Objects/Bloc4B_08_seu_Epithelial_SCEVAN.rds
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

# 2) Load and combine predictions
message("Loading SCEVAN predictions...")
pred_MPR  <- fread(file.path(OUT_TAB, "Bloc4B_08_SCEVAN_predictions_MPR.csv"))
pred_NMPR <- fread(file.path(OUT_TAB, "Bloc4B_08_SCEVAN_predictions_NMPR.csv"))
pred_all  <- rbind(pred_MPR, pred_NMPR)

message("Combined predictions:")
print(table(pred_all$scevan_class, pred_all$response))

# 3) Save combined predictions
fwrite(pred_all, file.path(OUT_TAB, "Bloc4B_08_SCEVAN_predictions_combined.csv"))
message("Saved: Bloc4B_08_SCEVAN_predictions_combined.csv")

# 4) Add to Seurat object
seu_Epi$scevan_class <- pred_all$scevan_class[match(
  colnames(seu_Epi), pred_all$cell)]
message("Predictions added to Seurat object")
message("By cell type:")
print(table(seu_Epi$scevan_class, seu_Epi$TME_cell_type))

# 5) UMAP by SCEVAN prediction
png(file.path(OUT_FIG, "Bloc4B_08_UMAP_SCEVAN_class.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_Epi,
              reduction = "umap",
              group.by  = "scevan_class",
              raster    = FALSE) +
        scale_color_manual(values = c("tumor"    = "#D73027",
                                      "normal"   = "#4393C3",
                                      "filtered" = "grey80")) +
        theme_bw() +
        theme(legend.text = element_text(size = 12)) +
        labs(title = "SCEVAN predictions — Epithelial cells GSE207422"))
dev.off()
message("Saved: Bloc4B_08_UMAP_SCEVAN_class.png")

# 6) Barplot by response
df_bar <- pred_all %>%
  filter(scevan_class %in% c("tumor", "normal")) %>%
  group_by(response, scevan_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(response) %>%
  mutate(prop = n / sum(n))

png(file.path(OUT_FIG, "Bloc4B_08_Barplot_SCEVAN_response.png"),
    width = 8, height = 6, units = "in", res = 300)
print(ggplot(df_bar, aes(x = response, y = prop, fill = scevan_class)) +
        geom_col(width = 0.8, color = "white") +
        scale_fill_manual(values = c("tumor"  = "#D73027",
                                     "normal" = "#4393C3")) +
        scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
        theme_classic() +
        theme(axis.title.x = element_blank(),
              axis.text.x  = element_text(size = 13, face = "bold"),
              legend.title = element_blank()) +
        ylab("Proportion") +
        labs(title = "Malignant vs normal epithelial cells by response — GSE207422"))
dev.off()
message("Saved: Bloc4B_08_Barplot_SCEVAN_response.png")

# 7) Save updated object
saveRDS(seu_Epi, file.path(OUT_OBJ, "Bloc4B_08_seu_Epithelial_SCEVAN.rds"))
message("Saved: Objects/Bloc4B_08_seu_Epithelial_SCEVAN.rds")

# 8) Session info
writeLines(capture.output(sessionInfo()),
           file.path(DATA_DIR, "session_info_linux.txt"))
message("DONE Bloc4B Script 08c")