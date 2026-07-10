#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script (ad hoc): Individual HLA-D gene
# expression by TAM subtype, patient-level pseudobulk Wilcoxon
# Purpose: resolve whether individual HLA-D genes (HLA-DQB1,
# HLA-DRB3) engaging LAG3 in communication analysis follow the
# same direction as the aggregate CollecTRI RFX/CIITA regulon
# score, or diverge from it.
# Input:  Objects/Bloc4B_04_seu_TAMs_combined.rds
# Output: Results/Tables/BLOC5/Bloc5_HLA_D_individual_genes_wilcox.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5")

# 1) Load TAM object
seu_TAM <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_04_seu_TAMs_combined.rds"))

# 2) Apply short TAM labels (same mapping used throughout Bloc 5)
tam_short_labels <- c(
  "TAM_like_MRC1"                                            = "MRC1+ M2-like",
  "TAM_like_SPP1"                                            = "SPP1+ immunosuppressive",
  "TAM_like_resident_M2 (iron metabolism/anti-inflammatory)" = "Resident M2",
  "TAM_like_IFN (PD-L1+/IDO1+/CXCL9+)"                     = "IFN-stimulated",
  "TAM_like_monocyte (classical inflammatory)"               = "Monocyte-derived",
  "TAM_like_lipid (CCL18+/AKR+)"                            = "Lipid-associated",
  "TAM_like_stress (HSP-high/M1-like)"                      = "Stress-response",
  "TAM_like_regulatory (glucocorticoid-responsive)"          = "Regulatory",
  "TAM_like_M2 (SIGLEC8+/CCL18+)"                           = "M2-SIGLEC8+"
)
seu_TAM$cell_type <- unname(tam_short_labels[as.character(seu_TAM$combined_annotation)])

# 3) Restrict to the 3 subtypes involved in the LAG3 loops
subtypes_of_interest <- c("M2-SIGLEC8+", "MRC1+ M2-like", "SPP1+ immunosuppressive")
seu_sub <- subset(seu_TAM, subset = cell_type %in% subtypes_of_interest)

# 4) Genes of interest
genes_of_interest <- c("HLA-DQB1", "HLA-DRB3")
genes_of_interest <- intersect(genes_of_interest, rownames(seu_sub))
message("Genes found in object: ", paste(genes_of_interest, collapse = ", "))

# 5) Extract expression, build patient-level pseudobulk
DefaultAssay(seu_sub) <- "RNA"
expr_mat <- GetAssayData(seu_sub, layer = "data")[genes_of_interest, , drop = FALSE]

meta <- seu_sub@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, cell_type, Sample, PathResponse)

expr_df <- as.data.frame(t(as.matrix(expr_mat))) %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(meta, by = "cell_id")

# 6) Patient-level pseudobulk mean, per subtype x patient
pseudobulk <- expr_df %>%
  group_by(cell_type, Sample, PathResponse) %>%
  summarise(across(all_of(genes_of_interest), mean), .groups = "drop")

# 7) Wilcoxon MPR vs NMPR, within each subtype, for each gene
results <- expand.grid(subtype = subtypes_of_interest, gene = genes_of_interest, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    mpr_vals  = list(pseudobulk %>% filter(cell_type == subtype, PathResponse == "MPR") %>% pull(gene)),
    nmpr_vals = list(pseudobulk %>% filter(cell_type == subtype, PathResponse == "NMPR") %>% pull(gene)),
    n_MPR = length(mpr_vals), n_NMPR = length(nmpr_vals),
    median_MPR = ifelse(n_MPR > 0, median(mpr_vals), NA),
    median_NMPR = ifelse(n_NMPR > 0, median(nmpr_vals), NA),
    p_value = ifelse(n_MPR >= 3 & n_NMPR >= 3, wilcox.test(mpr_vals, nmpr_vals)$p.value, NA)
  ) %>%
  ungroup() %>%
  select(subtype, gene, n_MPR, n_NMPR, median_MPR, median_NMPR, p_value) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

fwrite(results, file.path(OUT_TAB, "Bloc5_HLA_D_individual_genes_wilcox.csv"))
message("Saved: Bloc5_HLA_D_individual_genes_wilcox.csv")
print(results)