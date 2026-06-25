#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 09c: Final epithelial group classification
# EMT_NMPR reclassified as normal_NMPR based on:
# 1. SCEVAN highest specificity (0.75) in published benchmarks
#    (Lanucara et al. 2024, Biomedicines; De Falco et al. 2023)
# 2. CytoTRACE2 low stemness scores for EMT cells (both SCEVAN and CopyKAT labels)
# 3. CopyKAT tendency to overestimate tumor fractions in benchmarks
# Final groups: tumor_MPR, normal_MPR, tumor_NMPR, normal_NMPR, Ciliated
# Input:  Objects/Bloc4B_09a_seu_Epithelial_CytoTRACE2.rds
# Output: Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds
#         Results/Tables/Bloc4B_09c_Final_groups.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/CytoTRACE2_SCEVAN")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load epithelial object with CytoTRACE2 scores and SCEVAN groups
message("Loading epithelial object...")
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_09a_seu_Epithelial_CytoTRACE2.rds"))
message("Total epithelial cells: ", ncol(seu_Epi))

# 2) Reclassify EMT_NMPR as normal_NMPR
message("Reclassifying EMT_NMPR as normal_NMPR...")
seu_Epi$final_group <- seu_Epi$scevan_group
seu_Epi$final_group[seu_Epi$final_group == "EMT_NMPR"] <- "normal_NMPR"

message("Final group distribution:")
print(table(seu_Epi$final_group, useNA = "always"))

# 3) Violin plot — final groups
group_colors <- c(
  "tumor_MPR"   = "#E63946",
  "normal_MPR"  = "#457B9D",
  "tumor_NMPR"  = "#F4A261",
  "normal_NMPR" = "#2A9D8F",
  "Ciliated"    = "#ADB5BD"
)

seu_vln <- subset(seu_Epi, cells = colnames(seu_Epi)[!is.na(seu_Epi$final_group) &
                                                      !is.na(seu_Epi$CytoTRACE2_Score)])
p1 <- VlnPlot(seu_vln, features = "CytoTRACE2_Score",
              group.by = "final_group",
              cols = group_colors, pt.size = 0) +
  ggtitle("CytoTRACE2 Score — Final epithelial groups (EMT reclassified as normal_NMPR)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_FIG, "Bloc4B_09c_Violin_CytoTRACE2_FinalGroups.png"),
       p1, width = 10, height = 6, dpi = 300)
message("Saved: Violin final groups")

# 4) Save final groups CSV
scores_df <- data.frame(
  cell          = colnames(seu_Epi),
  final_group   = seu_Epi$final_group,
  scevan_group  = seu_Epi$scevan_group,
  CytoTRACE2_Score   = seu_Epi$CytoTRACE2_Score,
  CytoTRACE2_Potency = seu_Epi$CytoTRACE2_Potency,
  response      = seu_Epi$PathResponse,
  cell_type     = seu_Epi$TME_cell_type
)
fwrite(scores_df, file.path(OUT_TAB, "Bloc4B_09c_Final_groups.csv"))
message("Saved: Bloc4B_09c_Final_groups.csv")

# 5) Save final Seurat object
saveRDS(seu_Epi, file.path(DATA_DIR, "Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds"))
message("Saved: Bloc4B_09c_seu_Epithelial_FinalGroups.rds")

message("DONE Final Groups")