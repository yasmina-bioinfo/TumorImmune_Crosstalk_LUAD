#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Script 01: Load MTX matrix + filter LUAD cells
# Uses BPCells to avoid RAM overload
# Input:  Data/ (MTX + barcodes + genes + metadata)
# Output: Objects/01_LUAD_matrix_BPCells/
# ============================================================
suppressPackageStartupMessages({
  library(BPCells)
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR   <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD/Data"
OUT_DIR    <- file.path(DATA_DIR, "../Objects/01_LUAD_matrix_BPCells")

MTX_FILE      <- file.path(DATA_DIR, "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz")
BARCODES_FILE <- file.path(DATA_DIR, "GSE243013_barcodes.csv.gz")
GENES_FILE    <- file.path(DATA_DIR, "GSE243013_genes.csv.gz")
META_FILE     <- file.path(DATA_DIR, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")

# -----------------------------
# 1) Read metadata + extract LUAD barcodes
# -----------------------------
message("Reading metadata...")
meta <- fread(META_FILE)

# Extract LUAD cell barcodes only
luad_barcodes <- meta[cancer_type == "LUAD", cellID]
message("Number of LUAD cells: ", length(luad_barcodes))

# -----------------------------
# 2) Open MTX matrix with BPCells (reads from disk, no RAM load)
# -----------------------------
# NOTE: genes file has header 'geneSymbol'; must extract column explicitly
# NOTE: barcodes file has header 'barcode';same approach
# Without this, BPCells receives X+1 names for X rows/cols, leads to "error"

# NOTE: this MTX is cell x gene (not gene x cell)
# row_names = barcodes, col_names = genes
# Matrix will be transposed later when creating Seurat object

message("Opening MTX matrix with BPCells...")

barcode_names <- fread(BARCODES_FILE, header = TRUE)$barcode

gene_names <- fread(GENES_FILE, header = TRUE)$geneSymbol

mat <- import_matrix_market(
  mtx_path  = MTX_FILE,
  row_names = barcode_names, #cells = rows in MTX
  col_names = gene_names     # genes = columns in this MTX
)

message("Full matrix dim (genes x cells): ", paste(dim(mat), collapse = " x "))

# -----------------------------
# 3) Filter to LUAD cells only
# -----------------------------
# Transpose: BPCells loaded as cells x genes, we need genes x cells
mat <- t(mat)
# Convert to integer matrix for efficient compression
mat <- convert_matrix_type(mat, type = "uint32_t")

message("Transposed matrix dim (genes x cells): ", paste(dim(mat), collapse = " x "))

# Verify colnames are now barcodes
head(colnames(mat))
# NOTE: MTX file is stored as cells x genes; must transpose to genes x cells
# before filtering by barcode (colnames)

message("Filtering to LUAD cells...")
mat_luad <- mat[, luad_barcodes]
message("LUAD matrix dim (genes x cells): ", paste(dim(mat_luad), collapse = " x "))

# -----------------------------
# 4) Save LUAD matrix to disk (BPCells format)
# -----------------------------
message("Saving LUAD matrix to disk...")
write_matrix_dir(mat_luad, OUT_DIR)
message("Saved: Objects/01_LUAD_matrix_BPCells/")
message("DONE Script 01")