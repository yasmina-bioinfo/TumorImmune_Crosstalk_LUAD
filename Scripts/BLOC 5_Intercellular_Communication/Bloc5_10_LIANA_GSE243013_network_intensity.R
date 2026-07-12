#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc5 Script 10: Global network intensity test
# Purpose: same composite communication intensity approach as
# Script 10 for GSE207422 — one summary score per patient instead
# of hundreds of individual ligand-receptor tests, to check whether
# a global signal exists even when individual pairs do not survive
# multiple-testing correction.
# Input:  Objects/Bloc5_03_seu_TME_GSE243013.rds
#         Results/Tables/BLOC5/GSE243013/Bloc5_05_LIANA_GSE243013_differential_wilcox.csv
# Output: Results/Tables/BLOC5/GSE243013/Bloc5_10_LIANA_GSE243013_network_intensity.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(Matrix)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5/GSE243013")

# 1) Load the already-tested prioritized interaction list (from Script 05)
diff_results <- fread(file.path(OUT_TAB, "Bloc5_05_LIANA_GSE243013_differential_wilcox.csv"))
message("Prioritized CD8<->TAM interactions available: ", nrow(diff_results))

# 2) Load TME object
seu_TME <- readRDS(file.path(DATA_DIR, "Objects/Bloc5_03_seu_TME_GSE243013.rds"))

expr_mat <- GetAssayData(seu_TME, layer = "data")

meta_TME <- seu_TME@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, cell_type, sampleID, pathological_response) %>%
  mutate(group_key = paste(cell_type, sampleID, sep = "__"))

groups <- unique(meta_TME$group_key)
group_idx <- match(meta_TME$group_key, groups)
indicator <- Matrix::sparseMatrix(i = seq_along(group_idx), j = group_idx,
                                  x = 1, dims = c(nrow(meta_TME), length(groups)))
colnames(indicator) <- groups

group_sums   <- expr_mat %*% indicator
group_counts <- Matrix::colSums(indicator)
group_means  <- sweep(as.matrix(group_sums), 2, group_counts, "/")

group_response <- meta_TME %>% distinct(group_key, cell_type, sampleID, pathological_response)
unique_samples <- group_response %>% distinct(sampleID, pathological_response)

get_scores <- function(gene_complex, ctype) {
  genes <- strsplit(gene_complex, "_")[[1]]
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) return(rep(NA_real_, nrow(unique_samples)))
  
  keys <- paste(ctype, unique_samples$sampleID, sep = "__")
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

# 3) Composite network intensity per patient
message("Computing per-patient composite network intensity...")

interaction_matrix <- matrix(NA_real_, nrow = nrow(unique_samples), ncol = nrow(diff_results))

for (i in seq_len(nrow(diff_results))) {
  row <- diff_results[i, ]
  lig_vals <- get_scores(row$ligand, row$source)
  rec_vals <- get_scores(row$receptor, row$target)
  interaction_matrix[, i] <- lig_vals * rec_vals
}

composite_score <- rowSums(interaction_matrix, na.rm = TRUE)

network_df <- data.frame(sampleID = unique_samples$sampleID,
                         pathological_response = unique_samples$pathological_response,
                         composite_intensity = composite_score,
                         n_interactions_contributing = rowSums(!is.na(interaction_matrix)))

# 4) Wilcoxon test on composite score, MPR vs non-MPR
mpr_vals    <- network_df$composite_intensity[network_df$pathological_response == "MPR"]
nonmpr_vals <- network_df$composite_intensity[network_df$pathological_response == "non-MPR"]

test <- wilcox.test(mpr_vals, nonmpr_vals)
message("Global network intensity test: p = ", signif(test$p.value, 4),
        " (median MPR = ", round(median(mpr_vals), 2),
        ", median non-MPR = ", round(median(nonmpr_vals), 2), ")")

fwrite(network_df, file.path(OUT_TAB, "Bloc5_10_LIANA_GSE243013_network_intensity.csv"))
message("Saved: Bloc5_10_LIANA_GSE243013_network_intensity.csv")