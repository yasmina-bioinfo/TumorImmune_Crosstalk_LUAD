#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 07: CD8 T cells subsetting
# Corrected workflow: canonical annotation first, then ProjecTILs
# Extracts confirmed CD8 clusters only before ProjecTILs projection
# This corrects the methodological order error from Script 03
# Input:  Objects/Bloc3_02_seu_Tcells_clustered.rds
# Output: Objects/Bloc3_07_seu_CD8_pca.rds
#         Results/Figures/CD8/Bloc3_ElbowPlot_CD8.png
# CD8 confirmed clusters:
#   1  = CD8 T cells (effector memory/GZMK+)
#   2  = CD8 T cells (exhausted/cytotoxic TRM)
#   5  = CD8 T cells (TRM/transitional exhausted)
#   6  = CD8 T cells (TRM/quiescent)
#   10 = Proliferating T cells (cycling — majority CD8)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_02_seu_Tcells_clustered.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# 1) Load T cells clustered object
message("Loading T cells clustered object...")
seu_T <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_T), " | Clusters: ", length(levels(seu_T$seurat_clusters)))

# 2) Subset confirmed CD8 clusters
# NOTE: only clusters confirmed by canonical marker annotation
# Cluster 10 (Proliferating) included, majority CD8, to be validated by ProjecTILs
cd8_clusters <- c(1, 2, 5, 6, 10)

seu_CD8 <- subset(seu_T, subset = seurat_clusters %in% cd8_clusters)
message("CD8 subset: ", ncol(seu_CD8), " cells")
message("Cluster composition:")
print(table(seu_CD8$seurat_clusters))

# Free memory
rm(seu_T); gc()

# 3) Reset reductions and graphs
# NOTE: reductions from T cells object not relevant for CD8 subset
seu_CD8@reductions <- list()
seu_CD8@graphs     <- list()
DefaultAssay(seu_CD8) <- "RNA"

# 4) JoinLayers, required for Seurat v5 after subsetting
seu_CD8 <- JoinLayers(seu_CD8)

# 5) Normalization
message("Normalizing...")
seu_CD8 <- NormalizeData(seu_CD8,
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4,
                         verbose = FALSE)

# 6) HVG selection
# NOTE: 2000 HVGs, CD8 subset is more homogeneous than full T cell subset
message("Selecting HVGs (nfeatures = 2000)...")
seu_CD8 <- FindVariableFeatures(seu_CD8,
                                selection.method = "vst",
                                nfeatures = 2000,
                                verbose = FALSE)

# 7) Scaling
message("Scaling HVGs...")
seu_CD8 <- ScaleData(seu_CD8,
                     features = VariableFeatures(seu_CD8),
                     verbose = FALSE)

# 8) PCA
message("Running PCA...")
seu_CD8 <- RunPCA(seu_CD8,
                  features = VariableFeatures(seu_CD8),
                  npcs = 30,
                  verbose = FALSE)

# 9) ElbowPlot
message("Generating ElbowPlot...")
p_elbow <- ElbowPlot(seu_CD8, ndims = 30) +
  theme_bw() +
  labs(title = "PCA Elbow Plot — CD8 T cells subset GSE243013 LUAD",
       x = "Principal Component",
       y = "Standard Deviation") +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(file.path(OUT_FIG, "Bloc3_ElbowPlot_CD8.png"),
       p_elbow, width = 6, height = 4, dpi = 300, bg = "white")
message("Saved: Bloc3_ElbowPlot_CD8.png")

# 10) Save CD8 PCA object
message("Saving CD8 PCA object...")
saveRDS(seu_CD8, file.path(OUT_OBJ, "Bloc3_07_seu_CD8_pca.rds"))
message("Saved: Objects/Bloc3_07_seu_CD8_pca.rds")
message("DONE Bloc3 Script 07 — inspect ElbowPlot before running Script 08")