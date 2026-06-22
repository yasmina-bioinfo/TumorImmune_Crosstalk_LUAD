#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08c: CopyKAT, CNV inference NMPR epithelial cells
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Tables/Bloc4B_08_CopyKAT_predictions_NMPR.csv
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
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/CopyKAT_NMPR")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# 1) Load and subset NMPR + Ciliated
message("Loading epithelial object...")
seu_Epi <- readRDS(IN_OBJ)
DefaultAssay(seu_Epi) <- "RNA"

cells_NMPR <- colnames(seu_Epi)[seu_Epi$PathResponse == "NMPR" | 
                                 seu_Epi$TME_cell_type == "Ciliated_epithelial"]
seu_sub <- subset(seu_Epi, cells = cells_NMPR)
message("NMPR cells: ", sum(seu_sub$PathResponse == "NMPR"))
message("Ciliated reference: ", sum(seu_sub$TME_cell_type == "Ciliated_epithelial"))

# 2) Count matrix
count_matrix <- as.matrix(seu_sub[["RNA"]]$counts)

# 3) Normal reference
normal_cells <- colnames(seu_sub)[seu_sub$TME_cell_type == "Ciliated_epithelial"]
message("Normal reference: ", length(normal_cells))

# 4) Run CopyKAT
message("Running CopyKAT NMPR...")
setwd(OUT_FIG)
copykat_result <- copykat(
  rawmat          = count_matrix,
  id.type         = "S",
  cell.line       = "no",
  norm.cell.names = normal_cells,
  sam.name        = "CopyKAT_NMPR",
  n.cores         = 1
)
message("CopyKAT NMPR done.")

# 5) Save predictions
predictions <- copykat_result$prediction
predictions$response <- "NMPR"
print(table(predictions$copykat.pred))
fwrite(predictions, file.path(OUT_TAB, "Bloc4B_08_CopyKAT_predictions_NMPR.csv"))
message("Saved: Bloc4B_08_CopyKAT_predictions_NMPR.csv")

# 6) Add predictions to Seurat and produce UMAP + Barplot
pred_df <- predictions %>%
  dplyr::select(cell.names, copykat.pred) %>%
  dplyr::rename(CopyKAT_pred = copykat.pred)

seu_sub$CopyKAT_pred <- pred_df$CopyKAT_pred[match(colnames(seu_sub), pred_df$cell.names)]

# UMAP
p1 <- DimPlot(seu_sub, group.by = "CopyKAT_pred",
              cols = c("aneuploid" = "#E63946", "diploid" = "#457B9D", "not.defined" = "grey80")) +
  ggtitle("CopyKAT predictions — Epithelial NMPR GSE207422")
ggsave(file.path(OUT_FIG, "CopyKAT_NMPR_UMAP.png"), p1, width = 8, height = 6, dpi = 150)

# Barplot NMPR uniquement
bar_df <- seu_sub@meta.data %>%
  filter(!is.na(CopyKAT_pred), PathResponse == "NMPR") %>%
  group_by(CopyKAT_pred) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop = n / sum(n),
         response = "NMPR")

p2 <- ggplot(bar_df, aes(x = response, y = prop, fill = CopyKAT_pred)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("aneuploid" = "#E63946", "diploid" = "#457B9D", "not.defined" = "grey80")) +
  labs(title = "Malignant vs normal epithelial cells — NMPR (CopyKAT)",
       x = "", y = "Proportion") +
  theme_classic()
ggsave(file.path(OUT_FIG, "CopyKAT_NMPR_Barplot.png"), p2, width = 5, height = 5, dpi = 150)

message("DONE CopyKAT NMPR")