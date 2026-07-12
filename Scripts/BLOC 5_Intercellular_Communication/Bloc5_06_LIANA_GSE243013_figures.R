#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc5 Script 06: LIANA+ figures — main axis, all TAM subtypes
# Produces comparative dotplots for CD8 <-> TAMs axis, all 7 TAM subtypes.
# Each interaction shown is also tested statistically (Wilcoxon,
# patient-level pseudobulk, BH-corrected). aggregate_rank <= 0.05 is
# used only to PRIORITIZE which interactions to display (not a p-value).
# NOTE: broader in scope than the cells-of-interest script that follows
# (all 7 TAM subtypes here vs. the CollecTRI-discriminant subset there);
# both carry the same statistical rigor.
# Input:  Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_MPR_aggregated.csv
#         Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_non_MPR_aggregated.csv
#         Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_pCR_aggregated.csv
#         Objects/Bloc5_03_seu_TME_GSE243013.rds
# Output: Results/Figures/BLOC5_Communication/GSE243013/
#         Results/Tables/BLOC5/GSE243013/Bloc5_06_LIANA_GSE243013_wilcox_[axis].csv
# Reference: Dimitrov et al., Nature Communications 2022 (LIANA+)
# ============================================================

suppressPackageStartupMessages({
  library(liana)
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(Matrix)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE243013")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5/GSE243013")

# ============================================================
# 1) Load per condition results + TME object (for testing)
# ============================================================
message("Loading LIANA per condition results...")
liana_MPR    <- fread(file.path(OUT_TAB, "Bloc5_03_LIANA_GSE243013_MPR_aggregated.csv"))
liana_nonMPR <- fread(file.path(OUT_TAB, "Bloc5_03_LIANA_GSE243013_non_MPR_aggregated.csv"))
liana_pCR    <- fread(file.path(OUT_TAB, "Bloc5_03_LIANA_GSE243013_pCR_aggregated.csv"))

message("Loading TME object for statistical testing...")
seu_TME <- readRDS(file.path(DATA_DIR, "Objects/Bloc5_03_seu_TME_GSE243013.rds"))

cd8_types <- c("CD8.TEX", "CD8.TPEX")

short_tam_labels <- c(
  "Proliferating TAMs (cycling/MKI67+)"                                = "Proliferating",
  "Inflammatory monocyte-derived TAMs (FCN1+/S100A8+)"                 = "Monocyte FCN1+",
  "Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like)"  = "Resident M2",
  "TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs)"        = "Lipid-associated",
  "IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+)"      = "IFN-stimulated",
  "Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high)"    = "Stress-response",
  "Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+)"           = "Classical Mono-derived"
)

# ============================================================
# 2) Precompute pseudobulk means ONCE for the whole script
#    (same sparse matrix approach as Script 05, deduplicated by
#    unique sampleID to avoid the pseudo-replication bug found there)
# ============================================================
message("Precomputing pseudobulk expression means (once for all figures)...")

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

get_scores <- function(gene_complex, ctype) {
  genes <- strsplit(gene_complex, "_")[[1]]
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) return(NULL)
  
  # Deduplicated unique patients (fix for the Script 05 pseudo-replication bug)
  unique_samples <- group_response %>% distinct(sampleID, pathological_response)
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
  data.frame(sampleID = unique_samples$sampleID, pathological_response = unique_samples$pathological_response,
             mean_expr = vals)
}

test_interaction <- function(src, tgt, lig, rec, cond1, cond2) {
  lig_score <- get_scores(lig, src)
  rec_score <- get_scores(rec, tgt)
  if (is.null(lig_score) || is.null(rec_score)) return(c(p_value = NA, median_1 = NA, median_2 = NA))
  
  merged <- inner_join(lig_score, rec_score, by = c("sampleID", "pathological_response"),
                       suffix = c("_lig", "_rec")) %>%
    filter(!is.na(mean_expr_lig), !is.na(mean_expr_rec)) %>%
    mutate(score = mean_expr_lig * mean_expr_rec)
  
  g1 <- merged$score[merged$pathological_response == cond1]
  g2 <- merged$score[merged$pathological_response == cond2]
  if (length(g1) < 3 | length(g2) < 3) return(c(p_value = NA, median_1 = NA, median_2 = NA))
  
  p <- wilcox.test(g1, g2)$p.value
  c(p_value = p, median_1 = median(g1), median_2 = median(g2))
}

# ============================================================
# 3) Comparison function — with Wilcoxon testing + empty-data guard
# ============================================================
plot_liana_comparison <- function(liana_cond1, liana_cond2,
                                  source_groups, target_groups,
                                  label_cond1 = "MPR", label_cond2 = "non-MPR",
                                  ntop = 5, title = "", axis_name = "axis") {
  
  process_cond <- function(df, label) {
    df %>%
      filter(source %in% source_groups, target %in% target_groups,
             aggregate_rank <= 0.05) %>%
      group_by(target) %>%
      arrange(aggregate_rank) %>%
      slice_head(n = ntop) %>%
      ungroup() %>%
      mutate(condition = label,
             interaction = paste0(ligand.complex, " - ", receptor.complex))
  }
  
  df_all <- bind_rows(process_cond(liana_cond1, label_cond1),
                      process_cond(liana_cond2, label_cond2))
  
  if (nrow(df_all) == 0) {
    warning("No interactions found for axis: ", axis_name, " — skipping.")
    return(NULL)
  }
  
  unique_ints <- df_all %>% distinct(source, target, ligand.complex, receptor.complex)
  message("  Testing ", nrow(unique_ints), " interactions for axis: ", axis_name)
  
  test_results <- t(mapply(test_interaction,
                           unique_ints$source, unique_ints$target,
                           unique_ints$ligand.complex, unique_ints$receptor.complex,
                           MoreArgs = list(cond1 = label_cond1, cond2 = label_cond2)))
  unique_ints$p_value  <- test_results[, "p_value"]
  unique_ints[[paste0("median_", label_cond1)]] <- test_results[, "median_1"]
  unique_ints[[paste0("median_", label_cond2)]] <- test_results[, "median_2"]
  unique_ints$p_adj <- p.adjust(unique_ints$p_value, method = "BH")
  
  fwrite(unique_ints, file.path(OUT_TAB, paste0("Bloc5_06_LIANA_GSE243013_wilcox_", axis_name, ".csv")))
  message("  Saved: Bloc5_06_LIANA_GSE243013_wilcox_", axis_name, ".csv (",
          sum(unique_ints$p_adj < 0.05, na.rm = TRUE), " significant of ", nrow(unique_ints), ")")
  
  df_all <- df_all %>%
    left_join(unique_ints %>% select(source, target, ligand.complex, receptor.complex, p_adj),
              by = c("source", "target", "ligand.complex", "receptor.complex"))
  
  df_all$source <- ifelse(df_all$source %in% names(short_tam_labels),
                          short_tam_labels[df_all$source], df_all$source)
  df_all$target <- ifelse(df_all$target %in% names(short_tam_labels),
                          short_tam_labels[df_all$target], df_all$target)
  
  df_all$condition <- factor(df_all$condition, levels = c(label_cond1, label_cond2))
  color_map <- setNames(c("#4393C3", "#D73027"), c(label_cond1, label_cond2))
  
  ggplot(df_all, aes(x = interaction, y = target,
                     size = sca.LRscore, color = condition,
                     alpha = ifelse(p_adj < 0.05, 1, 0.35))) +
    geom_point() +
    scale_color_manual(values = color_map) +
    scale_size_continuous(range = c(3, 10)) +
    scale_alpha_identity() +
    facet_wrap(~condition, ncol = 1) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 13, color = "black"),
          axis.text.y  = element_text(size = 14, color = "black"),
          strip.text   = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 13),
          plot.title   = element_text(size = 15, face = "bold", hjust = 0.5)) +
    labs(title = title, x = "Ligand - Receptor", y = "Target",
         size = "LR Score", color = "Condition",
         caption = "Faded points: Wilcoxon p_adj >= 0.05 (patient-level pseudobulk)")
}

# ============================================================
# 4) Generate figures — MPR vs non-MPR (preprint-relevant)
# ============================================================
message("Generating MPR vs non-MPR figures...")

p1 <- plot_liana_comparison(liana_MPR, liana_nonMPR, cd8_types, names(short_tam_labels),
                            label_cond1 = "MPR", label_cond2 = "non-MPR",
                            title = "CD8 → TAMs interactions — MPR vs non-MPR (GSE243013)",
                            axis_name = "CD8_to_TAMs_MPR_nonMPR")
if (!is.null(p1)) ggsave(file.path(OUT_FIG, "Bloc5_04_LIANA_CD8_TAMs_MPR_nonMPR.png"), p1, width = 14, height = 10, dpi = 300, bg = "white")

p2 <- plot_liana_comparison(liana_nonMPR, liana_MPR, names(short_tam_labels), cd8_types,
                            label_cond1 = "non-MPR", label_cond2 = "MPR",
                            title = "TAMs → CD8 interactions — MPR vs non-MPR (GSE243013)",
                            axis_name = "TAMs_to_CD8_MPR_nonMPR")
if (!is.null(p2)) ggsave(file.path(OUT_FIG, "Bloc5_04_LIANA_TAMs_CD8_MPR_nonMPR.png"), p2, width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================
# 5) Generate figures — MPR vs pCR (perspectives-relevant)
# ============================================================
message("Generating MPR vs pCR figures...")

p3 <- plot_liana_comparison(liana_MPR, liana_pCR, cd8_types, names(short_tam_labels),
                            label_cond1 = "MPR", label_cond2 = "pCR",
                            title = "CD8 → TAMs interactions — MPR vs pCR (GSE243013)",
                            axis_name = "CD8_to_TAMs_MPR_pCR")
if (!is.null(p3)) ggsave(file.path(OUT_FIG, "Bloc5_04_LIANA_CD8_TAMs_MPR_pCR.png"), p3, width = 14, height = 10, dpi = 300, bg = "white")

p4 <- plot_liana_comparison(liana_MPR, liana_pCR, names(short_tam_labels), cd8_types,
                            label_cond1 = "MPR", label_cond2 = "pCR",
                            title = "TAMs → CD8 interactions — MPR vs pCR (GSE243013)",
                            axis_name = "TAMs_to_CD8_MPR_pCR")
if (!is.null(p4)) ggsave(file.path(OUT_FIG, "Bloc5_04_LIANA_TAMs_CD8_MPR_pCR.png"), p4, width = 14, height = 10, dpi = 300, bg = "white")

message("DONE LIANA+ GSE243013 figures — with Wilcoxon testing")