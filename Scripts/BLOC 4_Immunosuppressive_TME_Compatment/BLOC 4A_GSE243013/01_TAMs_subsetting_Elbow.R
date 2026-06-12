#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc4A Script 01: TAMs subsetting and ElbowPlot
# Input:  Objects/Bloc2_06_seu_final_annotated.rds
# Output: Objects/Bloc4A_01_seu_TAMs_pca.rds
#         Results/Figures/BLOC4A_TAMs/Bloc4A_ElbowPlot_TAMs.png
# NOTE: TAMs = cluster 7 (M2-like/immunosuppressive): 14,901 cells
#       Inspect ElbowPlot before running Script 02
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(BPCells)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc2_06_seu_final_annotated.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# 1) Load TME annotated object
message("Loading TME annotated object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu))

# 2) Subset TAMs : cluster 7
# NOTE: TAMs (M2-like/immunosuppressive) = cluster 7 in Bloc2 consensus annotation
# 14,901 cells, single cluster, high confidence annotation
message("Subsetting TAMs (cluster 7)...")
seu_TAM <- subset(seu, subset = seurat_clusters == 7)
message("TAMs subset: ", ncol(seu_TAM), " cells")
print(table(seu_TAM$pathological_response))

# Free memory
rm(seu); gc()

# 3) Reset reductions and graphs
# NOTE: reductions from TME object not relevant for TAMs subset
seu_TAM@reductions <- list()
seu_TAM@graphs     <- list()
DefaultAssay(seu_TAM) <- "RNA"

# 4) JoinLayers — required for Seurat v5 after subsetting
seu_TAM <- JoinLayers(seu_TAM)

# 5) Normalization
message("Normalizing...")
seu_TAM <- NormalizeData(seu_TAM,
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4,
                         verbose = FALSE)

# 6) HVG selection
# NOTE: 2000 HVGs, TAMs subset is a single cell type, more homogeneous than TME
message("Selecting HVGs (nfeatures = 2000)...")
seu_TAM <- FindVariableFeatures(seu_TAM,
                                selection.method = "vst",
                                nfeatures = 2000,
                                verbose = FALSE)

# 7) Scaling
message("Scaling HVGs...")
seu_TAM <- ScaleData(seu_TAM,
                     features = VariableFeatures(seu_TAM),
                     verbose = FALSE)

# 8) PCA
message("Running PCA...")
seu_TAM <- RunPCA(seu_TAM,
                  features = VariableFeatures(seu_TAM),
                  npcs = 30,
                  verbose = FALSE)

# 9) ElbowPlot
message("Generating ElbowPlot...")
p_elbow <- ElbowPlot(seu_TAM, ndims = 30) +
  theme_bw() +
  labs(title = "PCA Elbow Plot — TAMs subset GSE243013 LUAD",
       x = "Principal Component",
       y = "Standard Deviation") +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(file.path(OUT_FIG, "Bloc4A_ElbowPlot_TAMs.png"),
       p_elbow, width = 6, height = 4, dpi = 300, bg = "white")
message("Saved: Bloc4A_ElbowPlot_TAMs.png")

# 10) Save TAMs PCA object
message("Saving TAMs PCA object...")
saveRDS(seu_TAM, file.path(OUT_OBJ, "Bloc4A_01_seu_TAMs_pca.rds"))
message("Saved: Objects/Bloc4A_01_seu_TAMs_pca.rds")
message("DONE Bloc4A Script 01 — inspect ElbowPlot before running Script 02")