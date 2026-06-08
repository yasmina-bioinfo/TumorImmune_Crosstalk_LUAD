#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 01: T cells subsetting + ElbowPlot
# Input:  Objects/Bloc2_03_seu_sctype.rds
# Output: Objects/Bloc3_01_seu_Tcells_pca.rds
#         Results/Figures/Tcells/Bloc3_ElbowPlot_Tcells.png
# NOTE: Inspect ElbowPlot before running Script 02 (Harmony + clustering)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(BPCells)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc2_03_seu_sctype.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/Tcells")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# 1) Load TME object
message("Loading TME Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

# 2) Subset T cells
# Clusters identified as T cell populations in Bloc2 consensus annotation:
# 1  = CD4 T cells (naive/memory activated)
# 2  = CD8 T cells (exhausted/cytotoxic TRM)
# 3  = CD8 T cells (effector memory)
# 4  = Tregs (tumor-infiltrating)
# 5  = CD4 Tfh / exhausted CD4
# 9  = NKT / gamma-delta T cells
# 10 = CD4 T cells (naive/resting)
# 11 = Proliferating cells (cycling, majority T cells)
# 17 = IFN-stimulated cells (ISG-high, transversal, majority T cells)
# NOTE: cluster 8 (NK cells) excluded — distinct lineage from T cells
t_clusters <- c(1, 2, 3, 4, 5, 9, 10, 11, 17)

seu_T <- subset(seu, subset = seurat_clusters %in% t_clusters)
message("T cells subset: ", ncol(seu_T), " cells")
message("Original cluster composition:")
print(table(seu_T$seurat_clusters))

# Free memory
rm(seu); gc()

# 3) Reset reductions and graphs
# NOTE: old reductions from TME object are no longer relevant for T cells subset
seu_T@reductions <- list()
seu_T@graphs     <- list()
DefaultAssay(seu_T) <- "RNA"

# 4) JoinLayers : required for Seurat v5 after subsetting
seu_T <- JoinLayers(seu_T)

# 5) Normalization
message("Normalizing...")
seu_T <- NormalizeData(seu_T,
                       normalization.method = "LogNormalize",
                       scale.factor = 1e4,
                       verbose = FALSE)

# 6) HVG selection
# NOTE: 2000 HVGs for T cells subset,  more focused population than full TME
# 2000 sufficient here vs 3000 for full TME (epithelial + immune compartments)
message("Selecting HVGs (nfeatures = 2000)...")
seu_T <- FindVariableFeatures(seu_T,
                              selection.method = "vst",
                              nfeatures = 2000,
                              verbose = FALSE)

# 7) Scaling
message("Scaling HVGs...")
seu_T <- ScaleData(seu_T,
                   features = VariableFeatures(seu_T),
                   verbose = FALSE)

# 8) PCA
# NOTE: npcs = 40, same as TME global, will adjust after ElbowPlot inspection
message("Running PCA...")
seu_T <- RunPCA(seu_T,
                features = VariableFeatures(seu_T),
                npcs = 40,
                verbose = FALSE)

# 9) ElbowPlot : inspect before running Script 02
message("Generating ElbowPlot...")
p_elbow <- ElbowPlot(seu_T, ndims = 40) +
  theme_bw() +
  labs(title = "PCA Elbow Plot: T cells subset GSE243013 LUAD",
       x = "Principal Component",
       y = "Standard Deviation") +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(file.path(OUT_FIG, "Bloc3_ElbowPlot_Tcells.png"),
       p_elbow, width = 6, height = 4, dpi = 300, bg = "white")
message("Saved: Bloc3_ElbowPlot_Tcells.png")

# 10) Save object with PCA
message("Saving T cells PCA object...")
saveRDS(seu_T, file.path(OUT_OBJ, "Bloc3_01_seu_Tcells_pca.rds"))
message("Saved: Objects/Bloc3_01_seu_Tcells_pca.rds")
message("DONE Bloc3 Script 01 — inspect ElbowPlot before running Script 02")