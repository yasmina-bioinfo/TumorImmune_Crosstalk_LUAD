#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 02b: LIANA+ figures — preprint
# Filtered figures on cells of interest guided by Blocs 3 and 4
# CD8 priority: CD8.TEX, CD8.TPEX
# TAMs priority (based on CollecTRI significant TF count, natural
# drop-off after Resident M2): MRC1+ M2-like, SPP1+ immunosuppressive,
# IFN-stimulated, M2-SIGLEC8+, Resident M2
# Epithelial: Tumor + Normal + Ciliated (all kept)
# Each interaction shown is tested statistically (Wilcoxon,
# patient-level pseudobulk, BH-corrected). aggregate_rank <= 0.05
# is used only to PRIORITIZE which interactions to display.
# Input:  Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv
#         Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_NMPR_aggregated.csv
#         Objects/Bloc5_01_seu_TME_GSE207422.rds
# Output: Results/Figures/BLOC5_Communication/GSE207422/Preprint/
#         Results/Tables/BLOC5/Bloc5_02b_LIANA_GSE207422_wilcox_[axis].csv
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
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE207422/Preprint")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) Load per condition results + TME object (for testing)
# ============================================================
message("Loading LIANA per condition results...")
liana_MPR  <- fread(file.path(OUT_TAB, "Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv"))
liana_NMPR <- fread(file.path(OUT_TAB, "Bloc5_01_LIANA_GSE207422_NMPR_aggregated.csv"))

message("Loading TME object for statistical testing...")
seu_TME <- readRDS(file.path(DATA_DIR, "Objects/Bloc5_01_seu_TME_GSE207422.rds"))

# ============================================================
# 2) Cells of interest — guided by Blocs 3 and 4 CollecTRI findings
# ============================================================

cd8_poi <- c("CD8.TEX", "CD8.TPEX")

# TAMs priority: 5 discriminant subtypes (CollecTRI significant TF count,
# natural drop-off after Resident M2 — Stress-response, Lipid-associated,
# Regulatory, Monocyte-derived excluded)
tam_poi_raw <- c(
  "MRC1+ M2-like",
  "SPP1+ immunosuppressive",
  "IFN-stimulated",
  "M2-SIGLEC8+",
  "Resident M2"
)

epi_types <- c("tumor_MPR", "normal_MPR", "tumor_NMPR", "normal_NMPR", "Ciliated")

# ============================================================
# 3) Precompute pseudobulk means ONCE for the whole script
# ============================================================
message("Precomputing pseudobulk expression means (once for all figures)...")

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

get_scores <- function(gene_complex, ctype) {
  genes <- strsplit(gene_complex, "_")[[1]]
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) return(NULL)
  
  keys <- paste(ctype, group_response$Sample, sep = "__")
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
  data.frame(Sample = group_response$Sample, PathResponse = group_response$PathResponse,
             mean_expr = vals)
}

epi_pairs <- list(
  tumor  = c(mpr = "tumor_MPR",  nmpr = "tumor_NMPR"),
  normal = c(mpr = "normal_MPR", nmpr = "normal_NMPR")
)
epi_condition_encoded <- unlist(epi_pairs)

test_interaction <- function(src, tgt, lig, rec) {
  src_encoded <- src %in% epi_condition_encoded
  tgt_encoded <- tgt %in% epi_condition_encoded
  
  if (!src_encoded && !tgt_encoded) {
    lig_score <- get_scores(lig, src)
    rec_score <- get_scores(rec, tgt)
    if (is.null(lig_score) || is.null(rec_score)) return(c(NA, NA, NA))
    merged <- inner_join(lig_score, rec_score, by = c("Sample", "PathResponse"),
                         suffix = c("_lig", "_rec")) %>%
      filter(!is.na(mean_expr_lig), !is.na(mean_expr_rec)) %>%
      mutate(score = mean_expr_lig * mean_expr_rec)
    mpr  <- merged$score[merged$PathResponse == "MPR"]
    nmpr <- merged$score[merged$PathResponse == "NMPR"]
  } else {
    epi_side <- if (src_encoded) src else tgt
    epi_base <- names(epi_pairs)[sapply(epi_pairs, function(x) epi_side %in% x)]
    if (length(epi_base) == 0) return(c(NA, NA, NA))
    partner <- if (src_encoded) tgt else src
    
    score_for <- function(epi_group) {
      lig_ct <- if (src_encoded) partner else epi_group
      rec_ct <- if (src_encoded) epi_group else partner
      ls <- get_scores(lig, lig_ct); rs <- get_scores(rec, rec_ct)
      if (is.null(ls) || is.null(rs)) return(NULL)
      inner_join(ls, rs, by = c("Sample", "PathResponse"), suffix = c("_lig", "_rec")) %>%
        filter(!is.na(mean_expr_lig), !is.na(mean_expr_rec)) %>%
        mutate(score = mean_expr_lig * mean_expr_rec)
    }
    mpr_data  <- score_for(epi_pairs[[epi_base]]["mpr"])
    nmpr_data <- score_for(epi_pairs[[epi_base]]["nmpr"])
    if (is.null(mpr_data) || is.null(nmpr_data)) return(c(NA, NA, NA))
    mpr  <- mpr_data$score[mpr_data$PathResponse == "MPR"]
    nmpr <- nmpr_data$score[nmpr_data$PathResponse == "NMPR"]
  }
  
  if (length(mpr) < 3 | length(nmpr) < 3) return(c(NA, NA, NA))
  p <- wilcox.test(mpr, nmpr)$p.value
  c(p_value = p, median_MPR = median(mpr), median_NMPR = median(nmpr))
}

# ============================================================
# 4) Comparison function — with Wilcoxon testing + empty-data guard
# ============================================================
plot_liana_comparison <- function(liana_MPR, liana_NMPR,
                                  source_groups, target_groups,
                                  label_cond1 = "MPR", label_cond2 = "NMPR",
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
  
  df_all <- bind_rows(process_cond(liana_MPR, label_cond1),
                      process_cond(liana_NMPR, label_cond2))
  
  if (nrow(df_all) == 0) {
    warning("No interactions found for axis: ", axis_name, " — skipping.")
    return(NULL)
  }
  
  # Test using the ORIGINAL literal labels (tumor_MPR, tumor_NMPR, etc.)
  # — test_interaction() needs these exact names to know which arm to compare
  unique_ints <- df_all %>% distinct(source, target, ligand.complex, receptor.complex)
  message("  Testing ", nrow(unique_ints), " interactions for axis: ", axis_name)
  test_results <- t(mapply(test_interaction,
                           unique_ints$source, unique_ints$target,
                           unique_ints$ligand.complex, unique_ints$receptor.complex))
  unique_ints$p_value    <- test_results[, "p_value"]
  unique_ints$median_MPR <- test_results[, "median_MPR"]
  unique_ints$median_NMPR <- test_results[, "median_NMPR"]
  unique_ints$p_adj <- p.adjust(unique_ints$p_value, method = "BH")
  
  # NOW collapse tumor_MPR/tumor_NMPR -> "Tumor" and normal_MPR/normal_NMPR
  # -> "Normal", AFTER testing. Duplicate rows (same axis, tested twice
  # under two labels) now share identical results, so distinct() here
  # is safe and loses no information.
  unique_ints$source <- case_when(
    unique_ints$source %in% c("tumor_MPR", "tumor_NMPR") ~ "Tumor",
    unique_ints$source %in% c("normal_MPR", "normal_NMPR") ~ "Normal",
    TRUE ~ unique_ints$source)
  unique_ints$target <- case_when(
    unique_ints$target %in% c("tumor_MPR", "tumor_NMPR") ~ "Tumor",
    unique_ints$target %in% c("normal_MPR", "normal_NMPR") ~ "Normal",
    TRUE ~ unique_ints$target)
  unique_ints <- unique_ints %>% distinct(source, target, ligand.complex, receptor.complex, .keep_all = TRUE)
  
  fwrite(unique_ints, file.path(OUT_TAB, paste0("Bloc5_02b_LIANA_GSE207422_wilcox_", axis_name, ".csv")))
  message("  Saved: Bloc5_02b_LIANA_GSE207422_wilcox_", axis_name, ".csv (",
          sum(unique_ints$p_adj < 0.05, na.rm = TRUE), " significant of ", nrow(unique_ints), ")")
  
  # Apply the same collapse to df_all, so it can be joined to unique_ints
  df_all$source <- case_when(
    df_all$source %in% c("tumor_MPR", "tumor_NMPR") ~ "Tumor",
    df_all$source %in% c("normal_MPR", "normal_NMPR") ~ "Normal",
    TRUE ~ df_all$source)
  df_all$target <- case_when(
    df_all$target %in% c("tumor_MPR", "tumor_NMPR") ~ "Tumor",
    df_all$target %in% c("normal_MPR", "normal_NMPR") ~ "Normal",
    TRUE ~ df_all$target)
  
  df_all <- df_all %>%
    left_join(unique_ints %>% select(source, target, ligand.complex, receptor.complex, p_adj),
              by = c("source", "target", "ligand.complex", "receptor.complex"))
  
  # Cosmetic rename for figure display (TAM short labels + epithelial labels)
  df_all$source <- ifelse(df_all$source %in% names(short_tam_labels),
                          short_tam_labels[df_all$source], df_all$source)
  df_all$target <- ifelse(df_all$target %in% names(short_tam_labels),
                          short_tam_labels[df_all$target], df_all$target)
  
  df_all$source <- ifelse(df_all$source == "Tumor", "Tumor epithelial",
                          ifelse(df_all$source == "Normal", "Normal epithelial", df_all$source))
  df_all$target <- ifelse(df_all$target == "Tumor", "Tumor epithelial",
                          ifelse(df_all$target == "Normal", "Normal epithelial", df_all$target))
  
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
          axis.text.y  = element_text(size = 13, color = "black"),
          strip.text   = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 13),
          plot.title   = element_text(size = 15, face = "bold", hjust = 0.5)) +
    labs(title = title, x = "Ligand - Receptor", y = "Target",
         size = "LR Score", color = "Condition",
         caption = "Faded points: Wilcoxon p_adj >= 0.05 (patient-level pseudobulk)")
}

# ============================================================
# 5) Generate preprint figures — cells of interest only
# ============================================================
message("Generating preprint figures — cells of interest...")

p1 <- plot_liana_comparison(liana_MPR, liana_NMPR, cd8_poi, tam_poi_raw,
                            title = "CD8 → TAMs interactions — MPR vs NMPR (GSE207422)",
                            axis_name = "CD8_to_TAMs")
if (!is.null(p1)) ggsave(file.path(OUT_FIG, "Bloc5_02b_LIANA_CD8_TAMs_preprint.png"), p1, width = 14, height = 10, dpi = 300, bg = "white")

p2 <- plot_liana_comparison(liana_MPR, liana_NMPR, tam_poi_raw, cd8_poi,
                            title = "TAMs → CD8 interactions — MPR vs NMPR (GSE207422)",
                            axis_name = "TAMs_to_CD8")
if (!is.null(p2)) ggsave(file.path(OUT_FIG, "Bloc5_02b_LIANA_TAMs_CD8_preprint.png"), p2, width = 14, height = 10, dpi = 300, bg = "white")

p3 <- plot_liana_comparison(liana_MPR, liana_NMPR, cd8_poi, epi_types,
                            title = "CD8 → Epithelial interactions — MPR vs NMPR (GSE207422)",
                            axis_name = "CD8_to_Epithelial")
if (!is.null(p3)) ggsave(file.path(OUT_FIG, "Bloc5_02b_LIANA_CD8_Epithelial_preprint.png"), p3, width = 14, height = 10, dpi = 300, bg = "white")

p4 <- plot_liana_comparison(liana_MPR, liana_NMPR, epi_types, cd8_poi,
                            title = "Epithelial → CD8 interactions — MPR vs NMPR (GSE207422)",
                            axis_name = "Epithelial_to_CD8")
if (!is.null(p4)) ggsave(file.path(OUT_FIG, "Bloc5_02b_LIANA_Epithelial_CD8_preprint.png"), p4, width = 14, height = 10, dpi = 300, bg = "white")

p5 <- plot_liana_comparison(liana_MPR, liana_NMPR, tam_poi_raw, epi_types,
                            title = "TAMs → Epithelial interactions — MPR vs NMPR (GSE207422)",
                            axis_name = "TAMs_to_Epithelial")
if (!is.null(p5)) ggsave(file.path(OUT_FIG, "Bloc5_02b_LIANA_TAMs_Epithelial_preprint.png"), p5, width = 14, height = 10, dpi = 300, bg = "white")

p6 <- plot_liana_comparison(liana_MPR, liana_NMPR, epi_types, tam_poi_raw,
                            title = "Epithelial → TAMs interactions — MPR vs NMPR (GSE207422)",
                            axis_name = "Epithelial_to_TAMs")
if (!is.null(p6)) ggsave(file.path(OUT_FIG, "Bloc5_02b_LIANA_Epithelial_TAMs_preprint.png"), p6, width = 14, height = 10, dpi = 300, bg = "white")

message("DONE — LIANA+ GSE207422 preprint figures — cells of interest")