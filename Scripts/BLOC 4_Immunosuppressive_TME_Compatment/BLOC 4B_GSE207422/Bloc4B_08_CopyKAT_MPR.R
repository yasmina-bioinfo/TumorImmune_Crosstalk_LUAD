#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08: CopyKAT, CNV inference MPR epithelial cells
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Tables/Bloc4B_08_CopyKAT_predictions_MPR.csv
# Reference: Gao et al., Nature Genetics 2021
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4B_07_seu_Epithelial.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/CopyKAT_MPR")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# 1) Load and subset MPR + Ciliated
message("Loading epithelial object...")
seu_Epi <- readRDS(IN_OBJ)
DefaultAssay(seu_Epi) <- "RNA"

cells_MPR <- colnames(seu_Epi)[seu_Epi$PathResponse == "MPR" | 
                                seu_Epi$TME_cell_type == "Ciliated_epithelial"]
seu_sub <- subset(seu_Epi, cells = cells_MPR)
message("MPR cells: ", sum(seu_sub$PathResponse == "MPR"))
message("Ciliated reference: ", sum(seu_sub$TME_cell_type == "Ciliated_epithelial"))

# 2) Count matrix (raw counts, genes x cells)
count_matrix <- as.matrix(seu_sub[["RNA"]]$counts)

# 3) Normal reference barcodes
normal_cells <- colnames(seu_sub)[seu_sub$TME_cell_type == "Ciliated_epithelial"]
message("Normal reference: ", length(normal_cells))

# 4) Run CopyKAT
message("Running CopyKAT MPR...")
setwd(OUT_FIG)
copykat_result <- copykat(
  rawmat        = count_matrix,
  id.type       = "S",           # single cell
  cell.line     = "no",
  norm.cell.names = normal_cells,
  sam.name      = "CopyKAT_MPR",
  n.cores       = 1
)
message("CopyKAT MPR done.")

# 5) Save predictions
predictions <- copykat_result$prediction
predictions$response <- "MPR"
print(table(predictions$copykat.pred))
fwrite(predictions, file.path(OUT_TAB, "Bloc4B_08_CopyKAT_predictions_MPR.csv"))
message("Saved: Bloc4B_08_CopyKAT_predictions_MPR.csv")

# 6) Add predictions to Seurat object and produce UMAP + Barplot
pred_df <- predictions %>%
  dplyr::select(cell.names, copykat.pred) %>%
  dplyr::rename(CopyKAT_pred = copykat.pred)

seu_sub$CopyKAT_pred <- pred_df$CopyKAT_pred[match(colnames(seu_sub), pred_df$cell.names)]

# UMAP
p1 <- DimPlot(seu_sub, group.by = "CopyKAT_pred",
              cols = c("aneuploid" = "#E63946", "diploid" = "#457B9D", "not.defined" = "grey80")) +
  ggtitle("CopyKAT predictions — Epithelial MPR GSE207422")
ggsave(file.path(OUT_FIG, "CopyKAT_MPR_UMAP.png"), p1, width = 8, height = 6, dpi = 150)

# Barplot
bar_df <- seu_sub@meta.data %>%
  filter(!is.na(CopyKAT_pred), PathResponse == "MPR") %>%
  group_by(PathResponse, CopyKAT_pred) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PathResponse) %>%
  mutate(prop = n / sum(n))

p2 <- ggplot(bar_df, aes(x = PathResponse, y = prop, fill = CopyKAT_pred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("aneuploid" = "#E63946", "diploid" = "#457B9D", "not.defined" = "grey80")) +
  labs(title = "Malignant vs normal epithelial cells — MPR (CopyKAT)",
       x = "", y = "Proportion") +
  theme_classic()
ggsave(file.path(OUT_FIG, "CopyKAT_MPR_Barplot.png"), p2, width = 6, height = 5, dpi = 150)

message("DONE CopyKAT MPR")