#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc1 Script 02: QC filtering
# Input:  Objects/Bloc1_01_seu_raw.rds
# Output: Objects/Bloc1_02_seu_qc.rds
#         Results/Figures/Bloc1_QC_violin_post_filter.png
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc1_01_seu_raw.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures")

# -----------------------------
# 1) Load raw Seurat object
# -----------------------------
message("Loading raw Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells before filtering: ", ncol(seu))

# -----------------------------
# 2) Define QC thresholds
# -----------------------------
# NOTE: nFeature_RNA max threshold differs by group:
# non-MPR → 4000 (broader natural distribution observed in QC plots)
# MPR and pCR → 3000
# This avoids artificially standardizing distributions that are biologically distinct
# All other thresholds are uniform across groups

MIN_FEATURES <- 200
MAX_FEATURES_DEFAULT <- 3000   # MPR and pCR
MAX_FEATURES_NONMPR  <- 4000   # non-MPR
MAX_COUNTS   <- 12000
MAX_MT       <- 7
MIN_COUNTS   <- 500

# -----------------------------
# 3) Apply filtering
# -----------------------------
message("Applying QC filters...")

# Identify non-MPR cells
nonmpr_cells <- colnames(seu)[seu$pathological_response == "non-MPR"]
other_cells  <- colnames(seu)[seu$pathological_response != "non-MPR"]

# Filter non-MPR cells
seu_nonmpr <- subset(seu, cells = nonmpr_cells)
seu_nonmpr <- subset(seu_nonmpr,
                     subset = nFeature_RNA > MIN_FEATURES &
                       nFeature_RNA < MAX_FEATURES_NONMPR &
                       nCount_RNA   > MIN_COUNTS &
                       nCount_RNA   < MAX_COUNTS &
                       percent.mt   < MAX_MT)

# Filter MPR and pCR cells
seu_other <- subset(seu, cells = other_cells)
seu_other <- subset(seu_other,
                    subset = nFeature_RNA > MIN_FEATURES &
                      nFeature_RNA < MAX_FEATURES_DEFAULT &
                      nCount_RNA   > MIN_COUNTS &
                      nCount_RNA   < MAX_COUNTS &
                      percent.mt   < MAX_MT)

# Merge filtered objects
seu_qc <- merge(seu_nonmpr, seu_other)

message("Cells after filtering: ", ncol(seu_qc))
message("Cells removed: ", ncol(seu) - ncol(seu_qc))

# -----------------------------
# 4) Distribution after filtering
# -----------------------------
message("Pathological response distribution after filtering:")
print(table(seu_qc$pathological_response, useNA = "ifany"))

# -----------------------------
# 5) Violin plots post-filtering (for comparison)
# -----------------------------
message("Generating post-filter violin plots...")
p_vln <- VlnPlot(
  seu_qc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "pathological_response",
  ncol     = 3,
  pt.size  = 0
) & theme_bw()

ggsave(
  filename = file.path(OUT_FIG, "Bloc1_QC_violin_post_filter.png"),
  plot     = p_vln,
  width = 12, height = 4, dpi = 300, bg = "white"
)
message("Saved: Bloc1_QC_violin_post_filter.png")

# -----------------------------
# 6) Save filtered Seurat object
# -----------------------------
message("Saving filtered Seurat object...")
saveRDS(seu_qc, file.path(OUT_OBJ, "Bloc1_02_seu_qc.rds"))
message("Saved: Objects/Bloc1_02_seu_qc.rds")
message("DONE Bloc1 Script 02")