#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 03: MultiNicheNet differential communication
# Identifies ligand-receptor interactions differentially active
# between MPR and NMPR at the sample level (pseudobulk framework)
# Input:  Objects/Bloc5_01_seu_TME_GSE207422.rds
# Output: Results/Figures/BLOC5_Communication/GSE207422/
#         Results/Tables/BLOC5/
# Reference: Browaeys et al., bioRxiv 2023 (MultiNicheNet)
# ============================================================
# ============================================================
# NOTE — MultiNicheNet not used in final analysis
# MultiNicheNet was evaluated but not retained for the following reasons:
# 1. GSE207422 MPR group contains only n=3 patients — insufficient for 
#    robust pseudobulk differential analysis (recommended minimum: n=4-5)
# 2. Key cell types excluded at min_cells=10: tumor_MPR, tumor_NMPR, 
#    CD8.TPEX — central to the biological narrative
# 3. Methodological consistency across cohorts: applying MultiNicheNet 
#    to GSE207422 (n=3 MPR) but not GSE243013 would create asymmetry
# Decision: CellChat v2 retained as primary differential communication 
# tool — robust permutation-based statistics, no minimum patient requirement
# Script retained for methodological transparency and reproducibility
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(multinichenetr)
  library(SingleCellExperiment)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE207422")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5")

# ============================================================
# 1) Load TME object
# ============================================================
message("Loading TME object...")
seu_TME <- readRDS(file.path(DATA_DIR, "Objects/Bloc5_01_seu_TME_GSE207422.rds"))
message("Total cells: ", ncol(seu_TME))

# ============================================================
# 2) Convert to SingleCellExperiment
# ============================================================
message("Converting to SCE...")
seu_TME <- JoinLayers(seu_TME)
sce_TME <- as.SingleCellExperiment(seu_TME, assay = "RNA")

# Clean cell type names for MultiNicheNet compatibility
sce_TME$cell_type_clean <- make.names(sce_TME$cell_type)
message("Cell types cleaned: ", length(unique(sce_TME$cell_type_clean)), " types")

# ============================================================
# 3) Define MultiNicheNet parameters
# ============================================================
sample_id   <- "Sample"
group_id    <- "PathResponse"
celltype_id <- "cell_type_clean"   # use the cleaned version
min_cells   <- 5

message("Sample distribution:")
print(table(seu_TME$Sample, seu_TME$PathResponse))

# ============================================================
# 4) Run MultiNicheNet
# ============================================================
message("Running MultiNicheNet...")
multinichenet_output <- multi_nichenet_analysis(
  sce                    = sce_TME,
  celltype_id            = celltype_id,
  sample_id              = sample_id,
  group_id               = group_id,
  batches                = NA,
  covariates             = NA,
  min_cells              = 10,
  contrasts_oi           = c("MPR-NMPR"),
  contrast_tbl           = tibble(
    contrast = c("MPR-NMPR", "NMPR-MPR"),
    group    = c("MPR", "NMPR")
  ),
  ligand_target_matrix   = ligand_target_matrix,
  lr_network             = lr_network
)

message("MultiNicheNet done.")

# ============================================================
# 5) Save results
# ============================================================
message("Saving results...")
saveRDS(multinichenet_output,
        file.path(DATA_DIR, "Objects/Bloc5_03_multinichenet_GSE207422.rds"))

# Save prioritized interactions table
fwrite(as.data.frame(multinichenet_output$prioritization_tables$group_prioritization_tbl),
       file.path(OUT_TAB, "Bloc5_03_MultiNicheNet_GSE207422_prioritized.csv"))
message("Saved: prioritized interactions CSV")

# ============================================================
# 6) Figures — top prioritized interactions
# ============================================================
message("Generating figures...")

# Bubble plot — top interactions MPR enriched
p1 <- make_sample_lr_prod_activity_plots(
  multinichenet_output$prioritization_tables,
  prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
    filter(group == "MPR") %>%
    arrange(desc(prioritization_score)) %>%
    head(25)
)

ggsave(file.path(OUT_FIG, "Bloc5_03_MultiNicheNet_bubbleplot_MPR.png"),
       p1, width = 16, height = 10, dpi = 300, bg = "white")
message("Saved: bubbleplot MPR")

# Bubble plot — top interactions NMPR enriched
p2 <- make_sample_lr_prod_activity_plots(
  multinichenet_output$prioritization_tables,
  prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
    filter(group == "NMPR") %>%
    arrange(desc(prioritization_score)) %>%
    head(25)
)

ggsave(file.path(OUT_FIG, "Bloc5_03_MultiNicheNet_bubbleplot_NMPR.png"),
       p2, width = 16, height = 10, dpi = 300, bg = "white")
message("Saved: bubbleplot NMPR")

message("DONE MultiNicheNet GSE207422")