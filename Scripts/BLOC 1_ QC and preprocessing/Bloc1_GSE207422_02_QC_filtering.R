#!/usr/bin/env Rscript
# ============================================================
# GSE207422 — Bloc1 Script 02: QC filtering
# NOTE: nFeature_RNA max threshold differs by group, based on
# observed distributions for this dataset specifically:
# MPR → 6000 (broader distribution, q95=5346)
# NMPR → 5000 (narrower distribution, q95=4691)
# This is the opposite pattern from GSE243013, where non-MPR had
# the broader distribution — thresholds were set independently
# for each dataset based on its own observed QC plots, not copied
# across datasets.
# Input:  Objects/Bloc1_GSE207422_01_seu_raw.rds
# Output: Objects/Bloc1_GSE207422_02_seu_qc.rds
#         Results/Figures/QC and Clustering/Bloc1_GSE207422_QC_violin_post_filter.png
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc1_GSE207422_01_seu_raw.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/QC and Clustering")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

message("Loading raw Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells before filtering: ", ncol(seu))

MIN_FEATURES        <- 200
MAX_FEATURES_MPR    <- 6000
MAX_FEATURES_NONMPR <- 5000
MAX_COUNTS <- 12000
MAX_MT     <- 20
MIN_COUNTS <- 500

message("Applying QC filters...")
mpr_cells   <- colnames(seu)[seu$PathResponse == "MPR"]
nmpr_cells  <- colnames(seu)[seu$PathResponse == "NMPR"]

seu_mpr <- subset(seu, cells = mpr_cells)
seu_mpr <- subset(seu_mpr,
                  subset = nFeature_RNA > MIN_FEATURES &
                    nFeature_RNA < MAX_FEATURES_MPR &
                    nCount_RNA   > MIN_COUNTS &
                    nCount_RNA   < MAX_COUNTS &
                    percent.mt   < MAX_MT)

seu_nmpr <- subset(seu, cells = nmpr_cells)
seu_nmpr <- subset(seu_nmpr,
                   subset = nFeature_RNA > MIN_FEATURES &
                     nFeature_RNA < MAX_FEATURES_NONMPR &
                     nCount_RNA   > MIN_COUNTS &
                     nCount_RNA   < MAX_COUNTS &
                     percent.mt   < MAX_MT)

seu_qc <- merge(seu_mpr, seu_nmpr)
message("Cells after filtering: ", ncol(seu_qc))
message("Cells removed: ", ncol(seu) - ncol(seu_qc))

message("PathResponse distribution after filtering:")
print(table(seu_qc$PathResponse, useNA = "ifany"))

message("Generating post-filter violin plots...")
p_vln <- VlnPlot(
  seu_qc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "PathResponse",
  ncol     = 3,
  pt.size  = 0
) & theme_bw()
ggsave(
  filename = file.path(OUT_FIG, "Bloc1_GSE207422_QC_violin_post_filter.png"),
  plot     = p_vln,
  width = 12, height = 4, dpi = 300, bg = "white"
)
message("Saved: Bloc1_GSE207422_QC_violin_post_filter.png")

message("Saving filtered Seurat object...")
saveRDS(seu_qc, file.path(OUT_OBJ, "Bloc1_GSE207422_02_seu_qc.rds"))
message("Saved: Objects/Bloc1_GSE207422_02_seu_qc.rds")
message("DONE Bloc1 GSE207422 Script 02")