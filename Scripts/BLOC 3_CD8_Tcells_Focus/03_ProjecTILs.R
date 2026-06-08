#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 03: ProjecTILs annotation
# Projects T cells onto reference TIL atlas for precise state annotation
# Reference: Andreatta et al., Nature Communications 2021
# Input:  Objects/Bloc3_02_seu_Tcells_clustered.rds
# Output: Objects/Bloc3_03_seu_Tcells_ProjecTILs.rds
#         Results/Figures/Tcells/Bloc3_UMAP_ProjecTILs.png
#         Results/Tables/Bloc3_ProjecTILs_per_cluster.csv

# NOTE: get.reference.maps() downloads all ProjecTILs references (~2GB total)
# Human CD8 reference: ~500MB, download time ~2h on standard connection
# Run once, references cached locally in ProjecTILs_references/ folder

# NOTE: STACAS alignment failed due to memory limit (649MB > 500MB future default)
# ProjecTILs automatically switched to direct projection mode
# Direct projection is less precise than STACAS alignment but still valid
# To fix: options(future.globals.maxSize = 1000 * 1024^2) before Run.ProjecTILs()
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ProjecTILs)
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_02_seu_Tcells_clustered.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/Tcells")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load T cells clustered object
message("Loading T cells clustered object...")
seu_T <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_T), " | Clusters: ", length(levels(seu_T$seurat_clusters)))

# 2) Load ProjecTILs reference
# NOTE: ProjecTILs uses a pre-built TIL atlas reference
# Default reference: human CD8 TIL atlas (Andreatta et al. 2021)
# Downloaded automatically on first use
# NOTE: load.reference.map() only supports default mouse reference
# Use get.reference.maps() to access human references
# First run downloads ~2h — references cached locally after
message("Loading ProjecTILs reference...")
ref <- get.reference.maps()$human$CD8


# 3) Project T cells onto reference
# NOTE: ProjecTILs projects cells geometrically into the reference space
# Each cell receives the label of the nearest T cell state in the reference
# filter.cells = TRUE: removes cells too distant from reference (low quality projection)

message("Projecting T cells onto reference atlas...")
# NOTE: increase memory limit for future globals to avoid STACAS alignment failure
# Default 500MB limit causes alignment to fail — direct projection used as fallback
options(future.globals.maxSize = 1000 * 1024^2)

# NOTE: Run.ProjecTILs() with filter.cells=TRUE removed 108,456/172,110 cells (63%)
# Non-CD8 pure cells (CD4, Tregs, NKT, ISG-high) excluded by scGate
# 63,654 cells projected onto human CD8 reference
# STACAS alignment failed (649MB > 500MB limit) — direct projection used

seu_T <- Run.ProjecTILs(seu_T,
                        ref = ref,
                        filter.cells = TRUE)

message("ProjecTILs annotation done.")
message("Functional state distribution:")
print(table(seu_T$functional.cluster))

# 4) Summary table per cluster
message("Generating summary table...")

# Majority ProjecTILs label per new T cell cluster
proj_summary <- data.frame(
  cluster = seu_T$seurat_clusters,
  ProjecTILs = seu_T$functional.cluster
) %>%
  group_by(cluster) %>%
  count(ProjecTILs) %>%
  slice_max(n, n = 1) %>%
  summarise(ProjecTILs_majority = first(ProjecTILs)) %>%
  ungroup()

proj_summary$cluster <- as.integer(as.character(proj_summary$cluster))
proj_summary <- proj_summary[order(proj_summary$cluster), ]
print(proj_summary)

fwrite(proj_summary, file.path(OUT_TAB, "Bloc3_ProjecTILs_per_cluster.csv"))
message("Saved: Bloc3_ProjecTILs_per_cluster.csv")

# 5) UMAP colored by ProjecTILs annotation
message("Generating UMAP figures...")

png(file.path(OUT_FIG, "Bloc3_UMAP_ProjecTILs.png"),
    width = 12, height = 7, units = "in", res = 300)
print(DimPlot(seu_T,
              reduction = "umap",
              group.by  = "functional.cluster",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        labs(title = "ProjecTILs annotation : T cells GSE243013 LUAD"))
dev.off()

# Split by pathological response
png(file.path(OUT_FIG, "Bloc3_UMAP_ProjecTILs_split_response.png"),
    width = 16, height = 6, units = "in", res = 300)
print(DimPlot(seu_T,
              reduction = "umap",
              group.by  = "functional.cluster",
              split.by  = "pathological_response",
              raster    = FALSE) +
        theme_bw() +
        labs(title = "ProjecTILs : split by pathological response"))
dev.off()

message("Saved: UMAP ProjecTILs figures")

# 6) Save updated object
message("Saving ProjecTILs annotated object...")
saveRDS(seu_T, file.path(OUT_OBJ, "Bloc3_03_seu_Tcells_ProjecTILs.rds"))
message("Saved: Objects/Bloc3_03_seu_Tcells_ProjecTILs.rds")
message("DONE Bloc3 Script 03")