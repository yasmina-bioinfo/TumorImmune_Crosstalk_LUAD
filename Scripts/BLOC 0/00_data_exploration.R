#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Script 00: Metadata exploration (LUAD only)
# Output: Data/cancer_type_response_table.csv
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD/Data"

META_FILE <- file.path(DATA_DIR, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")

# -----------------------------
# 1) Read metadata
# -----------------------------
message("Reading metadata...")
meta <- fread(META_FILE)

# -----------------------------
# 2) Explore structure
# -----------------------------
message("Metadata columns: ", paste(colnames(meta), collapse = ", "))
message("Dimensions: ", paste(dim(meta), collapse = " x "))

# -----------------------------
# 3) Table cancer_type x pathological_response
# -----------------------------
message("Distribution cancer_type x pathological_response:")
print(table(meta$cancer_type, meta$pathological_response, useNA = "ifany"))

# -----------------------------
# 4) Count patients per response group (LUAD only)
# -----------------------------
meta_luad <- meta[cancer_type == "LUAD"]

patient_counts <- meta_luad[, .(n_patients = uniqueN(sampleID)), 
                            by = pathological_response]

message("Number of patients per response group (LUAD):")
print(patient_counts)

# Table 1 — cancer_type x pathological_response
write.csv(
  as.data.frame(table(meta$cancer_type, meta$pathological_response)),
  file.path(DATA_DIR, "../Results/Tables/cancer_type_response_table.csv"),
  row.names = FALSE
)

# Table 2 — patients per response group (LUAD only)
write.csv(
  patient_counts,
  file.path(DATA_DIR, "../Results/Tables/patient_counts_per_response_LUAD.csv"),
  row.names = FALSE
)

message("Table saved: cancer_type_response_table.csv")
message("DONE Script 00")