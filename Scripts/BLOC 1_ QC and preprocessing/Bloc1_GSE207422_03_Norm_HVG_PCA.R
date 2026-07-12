#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc1 Script 03: Normalization + HVG + Scaling + PCA
# Input:  Objects/Bloc1_GSE207422_02_seu_qc.rds
# Output: Objects/Bloc1_GSE207422_03_seu_pca.rds
#         Results/Figures/QC and Clustering/Bloc1_GSE207422_ElbowPlot.png
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc1_GSE207422_02_seu_qc.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/QC and Clustering")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

message("Loading QC Seurat object...")
seu <- readRDS(IN_OBJ)
DefaultAssay(seu) <- "RNA"
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

message("Normalizing...")
seu <- NormalizeData(seu,
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4,
                     verbose = FALSE)

message("Selecting highly variable genes (nfeatures = 3000)...")
seu <- FindVariableFeatures(seu,
                            selection.method = "vst",
                            nfeatures = 3000,
                            verbose = FALSE)
message("HVGs selected: ", length(VariableFeatures(seu)))

message("Scaling HVGs...")
seu <- ScaleData(seu,
                 features = VariableFeatures(seu),
                 verbose = FALSE)

message("Running PCA...")
seu <- RunPCA(seu,
              features = VariableFeatures(seu),
              npcs = 50,
              verbose = FALSE)

message("Generating ElbowPlot...")
p_elbow <- ElbowPlot(seu, ndims = 50) +
  theme_bw() +
  labs(title = "PCA Elbow Plot — GSE207422 NSCLC",
       x = "Principal Component",
       y = "Standard Deviation") +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(
  filename = file.path(OUT_FIG, "Bloc1_GSE207422_ElbowPlot.png"),
  plot     = p_elbow,
  width = 6, height = 4, dpi = 300, bg = "white"
)
message("Saved: Bloc1_GSE207422_ElbowPlot.png")

message("Saving PCA object...")
saveRDS(seu, file.path(OUT_OBJ, "Bloc1_GSE207422_03_seu_pca.rds"))
message("Saved: Objects/Bloc1_GSE207422_03_seu_pca.rds")
message("DONE Bloc1 GSE207422 Script 03")