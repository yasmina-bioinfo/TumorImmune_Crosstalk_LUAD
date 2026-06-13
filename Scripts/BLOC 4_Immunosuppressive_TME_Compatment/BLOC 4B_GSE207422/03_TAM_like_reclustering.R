#!/usr/bin/env Rscript
# ============================================================
# GSE207422  Bloc4B Script 03: TAM_like reclustering and annotation
# TAM_like (n=7,969) reclustered to identify hidden subpopulations
# NOTE: Initial annotation was less granular, reclustering reveals
#       biologically distinct subtypes within TAM_like population
# Input:  Objects/Bloc4B_01_seu_TAMs_annotated.rds
# Output: Objects/Bloc4B_03_seu_TAMlike_annotated.rds
#         Results/Tables/BLOC4B/Bloc4B_TAMlike_markers_top20.csv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

# Paths
DATA_DIR_OUTPUT <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ  <- file.path(DATA_DIR_OUTPUT, "Objects/Bloc4B_01_seu_TAMs_annotated.rds")
OUT_OBJ <- file.path(DATA_DIR_OUTPUT, "Objects")
OUT_TAB <- file.path(DATA_DIR_OUTPUT, "Results/Tables/BLOC4B")

# 1) Load TAMs object and subset TAM_like
message("Loading TAMs object...")
seu_TAM <- readRDS(IN_OBJ)
seu_TAMlike <- subset(seu_TAM, subset = final_annotation == "TAM_like")
message("TAM_like cells: ", ncol(seu_TAMlike))
rm(seu_TAM); gc()

# 2) Reset reductions
seu_TAMlike@reductions <- list()
seu_TAMlike@graphs     <- list()
DefaultAssay(seu_TAMlike) <- "RNA"
seu_TAMlike <- JoinLayers(seu_TAMlike)

# 3) Preprocessing
message("Preprocessing...")
seu_TAMlike <- NormalizeData(seu_TAMlike, verbose = FALSE)
seu_TAMlike <- FindVariableFeatures(seu_TAMlike, nfeatures = 2000, verbose = FALSE)
seu_TAMlike <- ScaleData(seu_TAMlike, features = VariableFeatures(seu_TAMlike), verbose = FALSE)
seu_TAMlike <- RunPCA(seu_TAMlike, npcs = 20, verbose = FALSE)

# 4) PCA stdev , ElbowPlot equivalent
message("PCA stdev:")
print(seu_TAMlike[["pca"]]@stdev)
# NOTE: 15 PCs retained based on stdev inspection — plateau after PC 15

# 5) Clustering
# NOTE: resolution = 0.3 gives 8 biologically meaningful clusters
message("Clustering...")
seu_TAMlike <- FindNeighbors(seu_TAMlike, dims = 1:15, verbose = FALSE)
seu_TAMlike <- FindClusters(seu_TAMlike, resolution = 0.3, verbose = FALSE)
seu_TAMlike <- RunUMAP(seu_TAMlike, dims = 1:15, verbose = FALSE)
message("Clusters: ", length(unique(seu_TAMlike$seurat_clusters)))

# 6) FindAllMarkers, top 20 per cluster
message("Finding markers...")
seu_TAMlike <- JoinLayers(seu_TAMlike)
markers_TAMlike <- FindAllMarkers(seu_TAMlike,
                                  only.pos       = TRUE,
                                  min.pct        = 0.25,
                                  logfc.threshold = 0.25,
                                  test.use       = "wilcox",
                                  verbose        = FALSE)

top20 <- markers_TAMlike %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20)

fwrite(as.data.frame(top20),
       file.path(OUT_TAB, "Bloc4B_TAMlike_markers_top20.csv"))
message("Saved: Bloc4B_TAMlike_markers_top20.csv")

# 7) Manual annotation , hardcoded from top20 markers analysis
# Cluster 7 excluded , T cell contamination (CD3G, CD3D, TRAC, GZMA)
tamlike_annotation <- data.frame(
  cluster = 0:7,
  final_annotation = c(
    "TAM_like_resident_M2 (iron metabolism/anti-inflammatory)",
    "TAM_like_IFN (PD-L1+/IDO1+/CXCL9+)",
    "TAM_like_monocyte (classical inflammatory)",
    "TAM_like_lipid (CCL18+/AKR+)",
    "TAM_like_stress (HSP-high/M1-like)",
    "TAM_like_regulatory (glucocorticoid-responsive)",
    "TAM_like_M2 (SIGLEC8+/CCL18+)",
    "EXCLUDE — T cell contamination"
  )
)

annot_map <- setNames(tamlike_annotation$final_annotation, tamlike_annotation$cluster)
test <- annot_map[as.character(seu_TAMlike$seurat_clusters)]
names(test) <- colnames(seu_TAMlike)
seu_TAMlike$final_annotation_detailed <- test

message("Annotation distribution:")
print(table(seu_TAMlike$final_annotation_detailed))

# 8) Save
message("Saving...")
saveRDS(seu_TAMlike,
        file.path(OUT_OBJ, "Bloc4B_03_seu_TAMlike_annotated.rds"))
message("Saved Objects/Bloc4B_03_seu_TAMlike_annotated.rds")
message("DONE Bloc4B Script 03")