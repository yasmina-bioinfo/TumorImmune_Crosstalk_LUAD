#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc1 Script 03: Normalization + HVG + Scaling + PCA
# Input:  Objects/Bloc1_02_seu_qc.rds
# Output: Objects/Bloc1_03_seu_pca.rds
#         Results/Figures/Bloc1_ElbowPlot.png
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc1_02_seu_qc.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures")

# -----------------------------
# 1) Load QC object
# -----------------------------
message("Loading QC Seurat object...")
seu <- readRDS(IN_OBJ)
DefaultAssay(seu) <- "RNA"
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

# -----------------------------
# 2) Normalization
# -----------------------------
# LogNormalize: divides each cell by total UMIs x scale.factor, then log1p
# scale.factor = 10000 is standard (counts per 10k)
message("Normalizing...")
seu <- NormalizeData(seu,
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4,
                     verbose = FALSE)

# -----------------------------
# 3) HVG selection
# -----------------------------
# 3000 HVGs selected (instead of default 2000) to better capture
# diversity across immune AND epithelial compartments in the TME
# Method: vst (variance stabilizing transformation)
message("Selecting highly variable genes (nfeatures = 3000)...")
seu <- FindVariableFeatures(seu,
                            selection.method = "vst",
                            nfeatures = 3000,
                            verbose = FALSE)
message("HVGs selected: ", length(VariableFeatures(seu)))

# -----------------------------
# 4) Scaling (HVGs only — RAM efficient)
# -----------------------------
# ScaleData centers and scales each gene (mean=0, variance=1)
# Applied only to HVGs to avoid RAM overload on all 26686 genes
# This ensures all genes contribute equally to PCA
message("Scaling HVGs...")
seu <- ScaleData(seu,
                 features = VariableFeatures(seu),
                 verbose = FALSE)

# -----------------------------
# 5) PCA
# -----------------------------
# npcs = 50 to have enough dimensions to evaluate the ElbowPlot
# We will select the optimal number of PCs after visual inspection
message("Running PCA...")
seu <- RunPCA(seu,
              features = VariableFeatures(seu),
              npcs = 50,
              verbose = FALSE)

# -----------------------------
# 6) ElbowPlot — to determine optimal number of PCs
# -----------------------------
# Read this plot before running Script 04 (Harmony + clustering)
# Choose the number of PCs where the curve flattens (the "elbow")
message("Generating ElbowPlot...")
p_elbow <- ElbowPlot(seu, ndims = 50) +
  theme_bw() +
  labs(title = "PCA Elbow Plot — GSE243013 LUAD",
       x = "Principal Component",
       y = "Standard Deviation") +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(
  filename = file.path(OUT_FIG, "Bloc1_ElbowPlot.png"),
  plot     = p_elbow,
  width = 6, height = 4, dpi = 300, bg = "white"
)
message("Saved: Bloc1_ElbowPlot.png")

# -----------------------------
# 7) Save object
# -----------------------------
message("Saving PCA object...")
saveRDS(seu, file.path(OUT_OBJ, "Bloc1_03_seu_pca.rds"))
message("Saved: Objects/Bloc1_03_seu_pca.rds")
message("DONE Bloc1 Script 03")