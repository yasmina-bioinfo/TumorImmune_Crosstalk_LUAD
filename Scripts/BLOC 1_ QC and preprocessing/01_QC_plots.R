#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc1 Script 01: Create Seurat object + QC plots
# Input:  Objects/01_LUAD_matrix_BPCells/
#         Results/Tables/02_clinical_metadata_LUAD.csv
# Output: Objects/Bloc1_01_seu_raw.rds
#         Results/Figures/Bloc1_QC_violin.png
#         Results/Figures/Bloc1_QC_scatter.png
# ============================================================
suppressPackageStartupMessages({
  library(BPCells)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR    <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
MAT_DIR     <- file.path(DATA_DIR, "Objects/01_LUAD_matrix_BPCells")
META_FILE   <- file.path(DATA_DIR, "Results/Tables/02_clinical_metadata_LUAD.csv")
OUT_OBJ     <- file.path(DATA_DIR, "Objects")
OUT_FIG     <- file.path(DATA_DIR, "Results/Figures")

dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Open BPCells matrix from disk
# -----------------------------
message("Opening BPCells matrix...")
mat <- open_matrix_dir(MAT_DIR)
message("Matrix dim (genes x cells): ", paste(dim(mat), collapse = " x "))

# -----------------------------
# 2) Create Seurat object
# -----------------------------
message("Creating Seurat object...")
seu <- CreateSeuratObject(
  counts      = mat,
  project     = "GSE243013_LUAD",
  min.cells   = 3,
  min.features = 200
)
message("Seurat object created: ", ncol(seu), " cells, ", nrow(seu), " genes")

# -----------------------------
# 3) Add clinical metadata
# -----------------------------
message("Adding clinical metadata...")

# Load full metadata to get cell-level sampleID
meta_full <- fread(file.path(DATA_DIR, "Data/GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"))

# Filter to LUAD cells present in the Seurat object
meta_luad <- meta_full[cellID %in% colnames(seu), .(cellID, sampleID)]

# Load clinical metadata
meta_clin <- fread(META_FILE)

# Join: cellID → sampleID → pathological_response
meta_luad <- merge(meta_luad, meta_clin, by = "sampleID", all.x = TRUE)

# Map to Seurat object (ordered by cell barcode)
meta_luad <- meta_luad[match(colnames(seu), meta_luad$cellID), ]

seu$sampleID <- meta_luad$sampleID
seu$pathological_response <- meta_luad$pathological_response

message("Pathological response distribution:")
print(table(seu$pathological_response, useNA = "ifany"))

# -----------------------------
# 4) QC metrics
# -----------------------------
message("Computing QC metrics...")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# -----------------------------
# 5) Violin plots (by pathological_response)
# -----------------------------
message("Generating violin plots...")
p_vln <- VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "pathological_response",
  ncol     = 3,
  pt.size  = 0
) & theme_bw()

ggsave(
  filename = file.path(OUT_FIG, "Bloc1_QC_violin.png"),
  plot     = p_vln,
  width = 12, height = 4, dpi = 300, bg = "white"
)
message("Saved: Bloc1_QC_violin.png")

# -----------------------------
# 6) Scatter plots
# -----------------------------
message("Generating scatter plots...")
p_sc1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme_bw()
p_sc2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme_bw()
p_sc  <- p_sc1 + p_sc2

ggsave(
  filename = file.path(OUT_FIG, "Bloc1_QC_scatter.png"),
  plot     = p_sc,
  width    = 10, height = 4, dpi = 300, bg = "white"
)
message("Saved: Bloc1_QC_scatter.png")

# -----------------------------
# 7) Save raw Seurat object
# -----------------------------
message("Saving raw Seurat object...")
saveRDS(seu, file.path(OUT_OBJ, "Bloc1_01_seu_raw.rds"))
message("Saved: Objects/Bloc1_01_seu_raw.rds")
message("DONE Bloc1 Script 01")