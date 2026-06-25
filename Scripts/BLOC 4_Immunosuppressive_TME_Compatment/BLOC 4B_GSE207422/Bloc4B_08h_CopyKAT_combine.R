#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08h: CopyKAT combine
# Combines predictions from all epithelial subtypes MPR + NMPR
# Produces: Barplot NMPR, Barplot MPR+NMPR, UMAP MPR+NMPR
# Input:  Results/Tables/Bloc4B_08*_CopyKAT_*_predictions.csv
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/
# Reference: Gao et al., Nature Genetics 2021
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs")

# 0) Ciliated barcodes for filtering
message("Loading Seurat object to extract ciliated barcodes...")
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_07_seu_Epithelial.rds"))
ciliated_barcodes <- colnames(seu_Epi)[seu_Epi$TME_cell_type == "Ciliated_epithelial"]
message("Ciliated barcodes to exclude: ", length(ciliated_barcodes))

# 1) Load and combine all NMPR predictions
message("Loading NMPR predictions...")
pred_AT2   <- fread(file.path(OUT_TAB, "Bloc4B_08e_CopyKAT_AT2_NMPR_predictions.csv"))
pred_Basal <- fread(file.path(OUT_TAB, "Bloc4B_08f_CopyKAT_Basal_NMPR_predictions.csv"))
pred_EMT   <- fread(file.path(OUT_TAB, "Bloc4B_08d_CopyKAT_EMT_NMPR_predictions.csv"))
pred_Tumor <- fread(file.path(OUT_TAB, "Bloc4B_08g_CopyKAT_TumorEpi_NMPR_predictions.csv"))

pred_NMPR <- bind_rows(pred_AT2, pred_Basal, pred_EMT, pred_Tumor) %>%
  filter(response == "NMPR") %>%
  filter(!cell.names %in% ciliated_barcodes)

message("NMPR cells after filtering: ", nrow(pred_NMPR))
print(table(pred_NMPR$copykat.pred))
fwrite(pred_NMPR, file.path(OUT_TAB, "Bloc4B_08h_CopyKAT_predictions_NMPR_combined.csv"))

# 2) Load MPR predictions
pred_MPR <- fread(file.path(OUT_TAB, "Bloc4B_08_CopyKAT_predictions_MPR.csv")) %>%
  filter(!cell.names %in% ciliated_barcodes)
message("MPR cells after filtering: ", nrow(pred_MPR))

# 3) Barplot NMPR only
bar_nmpr <- pred_NMPR %>%
  filter(!is.na(copykat.pred)) %>%
  group_by(copykat.pred) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop = n / sum(n), response = "NMPR")

p_nmpr <- ggplot(bar_nmpr, aes(x = response, y = prop, fill = copykat.pred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("aneuploid" = "#E63946", "diploid" = "#457B9D", "not.defined" = "grey80")) +
  labs(title = "Malignant vs normal epithelial cells — NMPR (CopyKAT)",
       x = "", y = "Proportion") +
  theme_classic()
ggsave(file.path(OUT_FIG, "CopyKAT_NMPR", "CopyKAT_NMPR_Barplot.png"),
       p_nmpr, width = 5, height = 5, dpi = 300)
message("Saved: CopyKAT_NMPR_Barplot.png")

# 4) Barplot MPR + NMPR combiné
pred_all <- bind_rows(pred_MPR, pred_NMPR)

bar_all <- pred_all %>%
  filter(!is.na(copykat.pred)) %>%
  group_by(response, copykat.pred) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(response) %>%
  mutate(prop = n / sum(n))

p_all <- ggplot(bar_all, aes(x = response, y = prop, fill = copykat.pred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("aneuploid" = "#E63946", "diploid" = "#457B9D", "not.defined" = "grey80")) +
  labs(title = "Malignant vs normal epithelial cells by response — GSE207422 (CopyKAT)",
       x = "", y = "Proportion") +
  theme_classic()
ggsave(file.path(OUT_FIG, "Bloc4B_08h_Barplot_CopyKAT_response.png"),
       p_all, width = 7, height = 5, dpi = 300)
message("Saved: Bloc4B_08h_Barplot_CopyKAT_response.png")

# 5) UMAP MPR + NMPR — exclude cells without CopyKAT prediction
pred_umap <- pred_all %>%
  dplyr::select(cell.names, copykat.pred) %>%
  dplyr::rename(CopyKAT_pred = copykat.pred)

seu_Epi$CopyKAT_pred <- pred_umap$CopyKAT_pred[match(colnames(seu_Epi), pred_umap$cell.names)]

# Keep only cells with a CopyKAT prediction
seu_sub_umap <- subset(seu_Epi, cells = pred_umap$cell.names)

p_umap <- DimPlot(seu_sub_umap, group.by = "CopyKAT_pred",
                  cols = c("aneuploid" = "#E63946", "diploid" = "#457B9D", "not.defined" = "grey80")) +
  ggtitle("CopyKAT predictions — Epithelial cells GSE207422")
ggsave(file.path(OUT_FIG, "Bloc4B_08h_UMAP_CopyKAT.png"),
       p_umap, width = 8, height = 6, dpi = 300)
message("Saved: Bloc4B_08h_UMAP_CopyKAT.png")