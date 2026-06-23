#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 09a: CytoTRACE2 on epithelial cells
# Labels from SCEVAN predictions (primary CNV tool)
# 6 groups: tumor_MPR, normal_MPR, tumor_NMPR, normal_NMPR, 
#           EMT_NMPR (independent group), Ciliated (reference)
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
#         Results/Tables/Bloc4B_08_SCEVAN_predictions_MPR.csv
#         Results/Tables/Bloc4B_08_SCEVAN_predictions_NMPR.csv
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/CytoTRACE2_SCEVAN/
#         Results/Tables/Bloc4B_09a_CytoTRACE2_SCEVAN_scores.csv
# Reference: Kang et al., Nature Methods 2025
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CytoTRACE2)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/CytoTRACE2_SCEVAN")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# 1) Load epithelial object
message("Loading epithelial object...")
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_07_seu_Epithelial.rds"))
DefaultAssay(seu_Epi) <- "RNA"
message("Total epithelial cells: ", ncol(seu_Epi))

# 2) Load SCEVAN predictions
message("Loading SCEVAN predictions...")
pred_MPR  <- fread(file.path(OUT_TAB, "Bloc4B_08_SCEVAN_predictions_MPR.csv"))
pred_NMPR <- fread(file.path(OUT_TAB, "Bloc4B_08_SCEVAN_predictions_NMPR.csv"))

pred_all <- bind_rows(pred_MPR, pred_NMPR)

# 3) Create scevan_group column
message("Creating SCEVAN group labels...")

seu_Epi$scevan_group <- NA

# Vectorized assignment
pred_lookup <- pred_all %>%
  mutate(scevan_group = case_when(
    seu_Epi$TME_cell_type[cell] == "Tumor_epithelial_EMT" ~ "EMT_NMPR",
    scevan_class == "tumor" ~ paste0("tumor_", response),
    scevan_class == "normal" ~ paste0("normal_", response),
    TRUE ~ NA_character_
  ))

seu_Epi$scevan_group[pred_lookup$cell] <- pred_lookup$scevan_group

# Ciliated cells
seu_Epi$scevan_group[seu_Epi$TME_cell_type == "Ciliated_epithelial"] <- "Ciliated"

message("Group distribution:")
print(table(seu_Epi$scevan_group, useNA = "always"))

# 4) Prepare count matrix for CytoTRACE2
message("Preparing count matrix...")
count_matrix <- as.data.frame(as.matrix(seu_Epi[["RNA"]]$counts))

# 5) Run CytoTRACE2 with batching to reduce RAM
message("Running CytoTRACE2...")
cytotrace2_result <- cytotrace2(
  input             = count_matrix,
  species           = "human",
  ncores            = 1,
  seed              = 42,
  batch_size        = 3000,
  smooth_batch_size = 1000
)
message("CytoTRACE2 done.")

# Save raw result immediately
saveRDS(cytotrace2_result, file.path(DATA_DIR, "Objects/Bloc4B_09a_cytotrace2_raw_result.rds"))
message("Raw result saved.")

# 6) Add scores to Seurat object
message("Matching scores to Seurat object...")
scores_match <- cytotrace2_result[colnames(seu_Epi), ]
seu_Epi$CytoTRACE2_Score   <- scores_match$CytoTRACE2_Score
seu_Epi$CytoTRACE2_Potency <- scores_match$CytoTRACE2_Potency
message("Score NA count: ", sum(is.na(seu_Epi$CytoTRACE2_Score)))

# 7) Save scores
scores_df <- data.frame(
  cell               = colnames(seu_Epi),
  CytoTRACE2_Score   = seu_Epi$CytoTRACE2_Score,
  CytoTRACE2_Potency = seu_Epi$CytoTRACE2_Potency,
  scevan_group       = seu_Epi$scevan_group,
  response           = seu_Epi$PathResponse,
  cell_type          = seu_Epi$TME_cell_type
)
fwrite(scores_df, file.path(OUT_TAB, "Bloc4B_09a_CytoTRACE2_SCEVAN_scores.csv"))
message("Saved: Bloc4B_09a_CytoTRACE2_SCEVAN_scores.csv")

# 8) Visualizations
group_colors <- c(
  "tumor_MPR"   = "#E63946",
  "normal_MPR"  = "#457B9D",
  "tumor_NMPR"  = "#F4A261",
  "normal_NMPR" = "#2A9D8F",
  "EMT_NMPR"    = "#9B2226",
  "Ciliated"    = "#ADB5BD"
)

# Save updated Seurat object
saveRDS(seu_Epi, file.path(DATA_DIR, "Objects/Bloc4B_09a_seu_Epithelial_CytoTRACE2.rds"))
message("Saved: Bloc4B_09a_seu_Epithelial_CytoTRACE2.rds")

message("DONE CytoTRACE2 SCEVAN")