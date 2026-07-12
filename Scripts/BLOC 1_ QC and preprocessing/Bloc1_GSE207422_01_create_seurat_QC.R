#!/usr/bin/env Rscript
# ============================================================
# GSE207422 Bloc1 Script 01: Create Seurat object + QC plots
# Cohort correction applied at the earliest possible step:
# BD_immune01 (NE, pre-treatment), BD_immune05 and BD_immune08
# (pre-treatment biopsies mislabeled with a post-treatment
# PathResponse outcome in some downstream uses), and BD_immune06
# (pCR, outside the binary MPR/NMPR comparison framework used in
# this project) are excluded here, before any clustering or
# annotation, per official GEO sample metadata
# (GSE207422_NSCLC_scRNAseq_metadata.xlsx).
# Input:  Data/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz
#         Data/GSE207422_NSCLC_scRNAseq_metadata.xlsx
# Output: Objects/Bloc1_GSE207422_01_seu_raw.rds
#         Results/Figures/Bloc1_GSE207422_QC_violin.png
#         Results/Figures/Bloc1_GSE207422_QC_scatter.png
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(readxl)
  library(dplyr)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR    <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
MAT_FILE    <- file.path(DATA_DIR, "Data/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz")
META_FILE   <- file.path(DATA_DIR, "Data/GSE207422_NSCLC_scRNAseq_metadata.xlsx")
OUT_OBJ     <- file.path(DATA_DIR, "Objects")
OUT_FIG     <- file.path(DATA_DIR, "Results/Figures")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load raw UMI matrix (text format, not BPCells for this dataset)
# -----------------------------
message("Loading UMI matrix (this may take a few minutes given file size)...")
mat_dt <- fread(MAT_FILE)
gene_names <- mat_dt[[1]]
mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- gene_names
message("Matrix dim (genes x cells): ", paste(dim(mat), collapse = " x "))

# -----------------------------
# 2) Load official GEO sample metadata and apply cohort correction
# -----------------------------
message("Loading GEO sample metadata...")
meta_samples <- read_excel(META_FILE)
message("Full cohort: ", nrow(meta_samples), " samples")
print(table(meta_samples$Resource, meta_samples$`Pathologic Response`))

samples_to_exclude <- c("BD_immune01", "BD_immune05", "BD_immune06", "BD_immune08")
message("\nExcluding: ", paste(samples_to_exclude, collapse = ", "),
        " (pre-treatment biopsies mislabeled with post-treatment outcome, or non-binary NE/pCR categories)")

meta_samples_corrected <- meta_samples %>% filter(!Sample %in% samples_to_exclude)
message("Corrected cohort: ", nrow(meta_samples_corrected), " samples")
print(table(meta_samples_corrected$`Pathologic Response`))

# -----------------------------
# 3) Identify cell barcodes belonging to excluded samples and remove
#    them from the matrix BEFORE creating the Seurat object
#    (assumes cell barcodes are prefixed with the Sample ID, e.g.
#    "BD_immune02_AAACCTG..."; adjust the separator/pattern below
#    if your barcode format differs)
# -----------------------------
cell_sample <- sapply(strsplit(colnames(mat), "_"), function(x) paste(x[1:2], collapse = "_"))
cells_to_keep <- colnames(mat)[!cell_sample %in% samples_to_exclude]
message("\nCells before exclusion: ", ncol(mat))
mat <- mat[, cells_to_keep]
message("Cells after excluding samples: ", ncol(mat))

# -----------------------------
# 4) Create Seurat object
# -----------------------------
message("Creating Seurat object...")
seu <- CreateSeuratObject(
  counts       = mat,
  project      = "GSE207422_NSCLC",
  min.cells    = 3,
  min.features = 200
)
message("Seurat object created: ", ncol(seu), " cells, ", nrow(seu), " genes")

# -----------------------------
# 5) Add clinical metadata (cell-level Sample assignment + response)
# -----------------------------
message("Adding clinical metadata...")
cell_sample_final <- sapply(strsplit(colnames(seu), "_"), function(x) paste(x[1:2], collapse = "_"))
seu$Sample <- cell_sample_final

meta_map <- meta_samples_corrected %>% select(Sample, `Pathologic Response`, `Clinical Stage`, Pathology, Sex, Age)
cell_meta <- data.frame(Sample = seu$Sample) %>% left_join(meta_map, by = "Sample")
seu$PathResponse <- cell_meta$`Pathologic Response`

message("PathResponse distribution:")
print(table(seu$PathResponse, useNA = "ifany"))
message("\nPer-patient cell counts:")
print(table(seu$Sample, seu$PathResponse))

# -----------------------------
# 6) QC metrics
# -----------------------------
message("Computing QC metrics...")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# -----------------------------
# 7) Violin plots (by PathResponse)
# -----------------------------
message("Generating violin plots...")
p_vln <- VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "PathResponse",
  ncol     = 3,
  pt.size  = 0
) & theme_bw()
ggsave(
  filename = file.path(OUT_FIG, "Bloc1_GSE207422_QC_violin.png"),
  plot     = p_vln,
  width = 12, height = 4, dpi = 300, bg = "white"
)
message("Saved: Bloc1_GSE207422_QC_violin.png")

# -----------------------------
# 8) Scatter plots
# -----------------------------
message("Generating scatter plots...")
p_sc1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme_bw()
p_sc2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme_bw()
p_sc  <- p_sc1 + p_sc2
ggsave(
  filename = file.path(OUT_FIG, "Bloc1_GSE207422_QC_scatter.png"),
  plot     = p_sc,
  width    = 10, height = 4, dpi = 300, bg = "white"
)
message("Saved: Bloc1_GSE207422_QC_scatter.png")

# -----------------------------
# 9) Save raw Seurat object
# -----------------------------
message("Saving raw Seurat object...")
saveRDS(seu, file.path(OUT_OBJ, "Bloc1_GSE207422_01_seu_raw.rds"))
message("Saved: Objects/Bloc1_GSE207422_01_seu_raw.rds")
message("DONE Bloc1 GSE207422 Script 01")