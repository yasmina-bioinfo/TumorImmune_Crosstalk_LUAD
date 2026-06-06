#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Script 02: Clinical metadata harmonization (LUAD only)
# Input:  Data/GSE243013_NSCLC_immune_scRNA_metadata.csv.gz
# Output: Results/Tables/02_clinical_metadata_LUAD.csv
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR    <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD/Data"
RESULTS_DIR <- file.path(DATA_DIR, "../Results/Tables")
META_FILE   <- file.path(DATA_DIR, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")

# -----------------------------
# 1) Read metadata
# -----------------------------
message("Reading metadata...")
meta <- fread(META_FILE)

# -----------------------------
# 2) Filter to LUAD only
# -----------------------------
meta_luad <- meta[cancer_type == "LUAD"]
message("LUAD cells: ", nrow(meta_luad))

# -----------------------------
# 3) Extract one row per patient (deduplicate by sampleID)
# -----------------------------
# Keep only patient-level columns (same value for all cells of a patient)
cols_patient <- c("sampleID", "pathological_response", "cancer_type",
                  "gender", "age", "anti-PD1_therapy", "radiological_response")

# Check all expected columns exist
missing_cols <- setdiff(cols_patient, colnames(meta_luad))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# Deduplicate: one row per patient
clinical_luad <- unique(meta_luad[, ..cols_patient, with = FALSE], by = "sampleID")
message("Number of unique LUAD patients: ", nrow(clinical_luad))

# -----------------------------
# 4) Inspect values for consistency
# -----------------------------
message("--- Consistency check ---")
message("pathological_response values: ", 
        paste(unique(clinical_luad$pathological_response), collapse = ", "))
message("cancer_type values: ", 
        paste(unique(clinical_luad$cancer_type), collapse = ", "))
message("gender values: ", 
        paste(unique(clinical_luad$gender), collapse = ", "))
message("age range: ", min(clinical_luad$age, na.rm = TRUE), 
        " - ", max(clinical_luad$age, na.rm = TRUE))
# NOTE: column name contains a hyphen — must use backticks `` ` `` to reference it in R
message("anti-PD1_therapy values: ", 
        paste(unique(clinical_luad$`anti-PD1_therapy`), collapse = ", "))
message("radiological_response values: ", 
        paste(unique(clinical_luad$radiological_response), collapse = ", "))

# Check for NA values
message("--- NA check ---")
print(colSums(is.na(clinical_luad)))

# NOTE: 'unknowm' is a typo in the original dataset (anti-PD1_therapy and 
# radiological_response). Kept as-is for traceability.
# NOTE: anti-PD1_therapy is heterogeneous (13 values including combinations)
# ;treat as complex covariable, not suitable for simple stratification.

# -----------------------------
# 5) Save clean clinical table
# -----------------------------
fwrite(clinical_luad, file.path(RESULTS_DIR, "02_clinical_metadata_LUAD.csv"))
message("Saved: Results/Tables/02_clinical_metadata_LUAD.csv")
message("DONE Script 02")