#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 01: LIANA+ intercellular communication
# Infers ligand-receptor interactions across CD8, TAMs, and Epithelial compartments
# Compares interactions across MPR and NMPR conditions
# Methods: sca, natmi, connectome (aggregated consensus)
# Input:  Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds
#         Objects/Bloc4B_04_seu_TAMs_combined.rds
#         Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds
# Output: Results/Figures/BLOC5_Communication/GSE207422/
#         Results/Tables/BLOC5/
# Reference: Dimitrov et al., Nature Communications 2022 (LIANA+)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(liana)
  library(SingleCellExperiment)
  library(tidyverse)
  library(ggplot2)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE207422")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) Load objects
# ============================================================
message("Loading objects...")
seu_CD8 <- readRDS(file.path(DATA_DIR, "Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds"))
seu_TAM <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_04_seu_TAMs_combined.rds"))
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds"))

# ============================================================
# 2) Unify cell_type column
# ============================================================
seu_CD8$cell_type <- seu_CD8$functional.cluster
seu_TAM$cell_type <- seu_TAM$combined_annotation
seu_Epi$cell_type <- seu_Epi$final_group

# ============================================================
# 3) Merge objects
# ============================================================
message("Merging objects...")
seu_TME <- merge(seu_CD8, y = list(seu_TAM, seu_Epi),
                 add.cell.ids = c("CD8", "TAM", "Epi"),
                 merge.data = TRUE)

# ============================================================
# 4) Remove NA cell types
# ============================================================
seu_TME <- subset(seu_TME, cells = colnames(seu_TME)[!is.na(seu_TME$cell_type)])
message("Total cells after QC: ", ncol(seu_TME))
message("Cell types: ", length(unique(seu_TME$cell_type)))
message("Conditions: ", paste(unique(seu_TME$PathResponse), collapse=", "))

# ============================================================
# 5) Normalize and prepare object
# ============================================================
DefaultAssay(seu_TME) <- "RNA"
seu_TME <- JoinLayers(seu_TME)
seu_TME <- NormalizeData(seu_TME)
Idents(seu_TME) <- "cell_type"

# ============================================================
# 6) Convert to SCE and run LIANA+ - 3 methods aggregated
# ============================================================
# Convert to SingleCellExperiment for LIANA+ compatibility with Seurat v5
sce_TME <- as.SingleCellExperiment(seu_TME, assay = "RNA")

message("Running LIANA+ (sca + natmi + connectome)...")
liana_results <- liana_wrap(sce_TME,
                            method = c("sca", "natmi", "connectome"),
                            resource = "Consensus",
                            idents_col = "cell_type")

# Aggregate results across methods
liana_agg <- liana_aggregate(liana_results)
message("LIANA aggregated — ", nrow(liana_agg), " interactions")

# Save results
fwrite(as.data.frame(liana_agg),
       file.path(OUT_TAB, "Bloc5_01_LIANA_GSE207422_aggregated.csv"))
message("Saved: CSV aggregated interactions")

# ============================================================
# 7) Run LIANA per condition - MPR and NMPR separately
# ============================================================
message("Running LIANA per condition...")

for (cond in c("MPR", "NMPR")) {
  sce_sub <- sce_TME[, sce_TME$PathResponse == cond]
  message("Running LIANA for ", cond, " — ", ncol(sce_sub), " cells")
  
  liana_cond <- liana_wrap(sce_sub,
                           method = c("sca", "natmi", "connectome"),
                           resource = "Consensus",
                           idents_col = "cell_type")
  
  liana_cond_agg <- liana_aggregate(liana_cond)
  
  fwrite(as.data.frame(liana_cond_agg),
         file.path(OUT_TAB, paste0("Bloc5_01_LIANA_GSE207422_", cond, "_aggregated.csv")))
  message("Saved: LIANA ", cond)
}

# ============================================================
# 8) Visualizations - dot plot top interactions
# ============================================================
message("Generating figures...")

p1 <- liana_agg %>%
  filter(aggregate_rank <= 0.05) %>%
  liana_dotplot(source_groups = c("CD8.TEX", "CD8.TPEX"),
                target_groups = c("TAM_like_IFN (PD-L1+/IDO1+/CXCL9+)",
                                  "TAM_like_stress (HSP-high/M1-like)",
                                  "tumor_NMPR", "tumor_MPR",
                                  "normal_MPR", "normal_NMPR"),
                ntop = 20) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 17, face = "bold", hjust = 0.5)) +
  ggtitle("Top CD8 → TAM/Epithelial interactions — GSE207422")

ggsave(file.path(OUT_FIG, "Bloc5_01_LIANA_dotplot_CD8_to_TME.png"),
        p1, width = 18, height = 11, dpi = 300, bg = "white")

message("Saved: Dotplot CD8 to TME")

# ============================================================
# 9) Save TME object
# ============================================================
saveRDS(seu_TME, file.path(DATA_DIR, "Objects/Bloc5_01_seu_TME_GSE207422.rds"))
message("Saved: Bloc5_01_seu_TME_GSE207422.rds")

message("DONE LIANA+ GSE207422")