#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 09: Global network intensity test
# Purpose: instead of testing hundreds of individual ligand-receptor
# pairs separately (diluted by multiple-testing correction), compute
# ONE composite communication intensity score per patient (sum of
# interaction scores across all prioritized CD8<->TAM pairs), then
# test this single summary score between conditions.
# Approach inspired by published CellChat "network intensity"
# comparisons (e.g., macrophage communication intensity in ESCC,
# Cancer Immunol Immunother 2024) rather than per-pair testing.
# Input:  Objects/Bloc5_01_seu_TME_GSE207422.rds (from corrected Script 01)
#         Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_differential_wilcox.csv
# Output: Results/Tables/BLOC5/Bloc5_09_LIANA_GSE207422_network_intensity.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(Matrix)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5")

# 1) Load the already-tested prioritized interaction list (from corrected Script 01)
diff_results <- fread(file.path(OUT_TAB, "Bloc5_01_LIANA_GSE207422_differential_wilcox.csv"))
message("Prioritized CD8<->TAM interactions available: ", nrow(diff_results))

# 2) Load TME object (already has cell_type, Sample, PathResponse)
seu_TME <- readRDS(file.path(DATA_DIR, "Objects/Bloc5_01_seu_TME_GSE207422.rds"))

expr_mat <- GetAssayData(seu_TME, layer = "data")

meta_TME <- seu_TME@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, cell_type, Sample, PathResponse) %>%
  mutate(group_key = paste(cell_type, Sample, sep = "__"))

groups <- unique(meta_TME$group_key)
group_idx <- match(meta_TME$group_key, groups)
indicator <- Matrix::sparseMatrix(i = seq_along(group_idx), j = group_idx,
                                  x = 1, dims = c(nrow(meta_TME), length(groups)))
colnames(indicator) <- groups

group_sums   <- expr_mat %*% indicator
group_counts <- Matrix::colSums(indicator)
group_means  <- sweep(as.matrix(group_sums), 2, group_counts, "/")

group_response <- meta_TME %>% distinct(group_key, cell_type, Sample, PathResponse)
unique_samples <- group_response %>% distinct(Sample, PathResponse)

get_scores <- function(gene_complex, ctype) {
  genes <- strsplit(gene_complex, "_")[[1]]
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) return(rep(NA_real_, nrow(unique_samples)))
  
  keys <- paste(ctype, unique_samples$Sample, sep = "__")
  present <- keys %in% colnames(group_means)
  
  if (length(genes) == 1) {
    vals <- rep(NA_real_, length(keys))
    vals[present] <- group_means[genes, keys[present]]
  } else {
    sub_vals <- sapply(genes, function(g) {
      v <- rep(NA_real_, length(keys))
      v[present] <- group_means[g, keys[present]]
      v
    })
    vals <- exp(rowMeans(log1p(sub_vals))) - 1
  }
  vals
}

# 3) For each patient, compute the SUM of interaction scores across
#    ALL prioritized CD8<->TAM pairs (composite network intensity)
message("Computing per-patient composite network intensity...")

interaction_matrix <- matrix(NA_real_, nrow = nrow(unique_samples), ncol = nrow(diff_results))

for (i in seq_len(nrow(diff_results))) {
  row <- diff_results[i, ]
  lig_vals <- get_scores(row$ligand, row$source)
  rec_vals <- get_scores(row$receptor, row$target)
  interaction_matrix[, i] <- lig_vals * rec_vals
}

# Composite score per patient = sum across all interactions (NA treated as 0 contribution)
composite_score <- rowSums(interaction_matrix, na.rm = TRUE)

network_df <- data.frame(Sample = unique_samples$Sample,
                         PathResponse = unique_samples$PathResponse,
                         composite_intensity = composite_score,
                         n_interactions_contributing = rowSums(!is.na(interaction_matrix)))

# 4) Single Wilcoxon test on the composite score
mpr_vals  <- network_df$composite_intensity[network_df$PathResponse == "MPR"]
nmpr_vals <- network_df$composite_intensity[network_df$PathResponse == "NMPR"]

test <- wilcox.test(mpr_vals, nmpr_vals)
message("Global network intensity test: p = ", signif(test$p.value, 4),
        " (median MPR = ", round(median(mpr_vals), 2),
        ", median NMPR = ", round(median(nmpr_vals), 2), ")")

fwrite(network_df, file.path(OUT_TAB, "Bloc5_09_LIANA_GSE207422_network_intensity.csv"))
message("Saved: Bloc5_01b_LIANA_GSE207422_network_intensity.csv")