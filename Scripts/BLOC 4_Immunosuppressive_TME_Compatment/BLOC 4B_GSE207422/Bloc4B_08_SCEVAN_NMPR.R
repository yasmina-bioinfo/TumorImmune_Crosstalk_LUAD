#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08b: SCEVAN , CNV inference NMPR epithelial cells
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Tables/Bloc4B_08_SCEVAN_predictions_NMPR.csv
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
count_matrix <- GetAssayData(seu_sub, layer = "counts")

# 3) Normal reference
normal_cells <- colnames(seu_sub)[seu_sub$TME_cell_type == "Ciliated_epithelial"]
message("Normal reference: ", length(normal_cells))

# 4) Run SCEVAN
message("Running SCEVAN NMPR...")
setwd(file.path(OUT_FIG))
results <- pipelineCNA(
  count_matrix,
  norm_cell  = normal_cells,
  par_cores  = 1,
  SUBCLONES  = FALSE,
  plotTree   = FALSE,
  ClonalCN   = FALSE,
  output_dir = file.path(OUT_FIG, "SCEVAN_NMPR")
)
message("SCEVAN NMPR done.")

# 5) Save predictions
predictions <- data.frame(
  cell         = rownames(results),
  scevan_class = results$class,
  response     = "NMPR"
)

# 6) Barplot SCEVAN NMPR
library(ggplot2)

bar_df <- predictions %>%
  filter(scevan_class != "filtered") %>%
  group_by(scevan_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop = n / sum(n),
         response = "NMPR")

p <- ggplot(bar_df, aes(x = response, y = prop, fill = scevan_class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("tumor" = "#E63946", "normal" = "#457B9D")) +
  labs(title = "Malignant vs normal epithelial cells — NMPR (SCEVAN)",
       x = "", y = "Proportion") +
  theme_classic()

ggsave(file.path(OUT_FIG, "SCEVAN_NMPR", "Bloc4B_08_Barplot_SCEVAN_NMPR.png"), 
       p, width = 5, height = 5, dpi = 300)
message("Barplot SCEVAN NMPR saved.")

print(table(predictions$scevan_class))
fwrite(predictions, file.path(OUT_TAB, "Bloc4B_08_SCEVAN_predictions_NMPR.csv"))
message("Saved: Bloc4B_08_SCEVAN_predictions_NMPR.csv")
message("DONE SCEVAN NMPR")