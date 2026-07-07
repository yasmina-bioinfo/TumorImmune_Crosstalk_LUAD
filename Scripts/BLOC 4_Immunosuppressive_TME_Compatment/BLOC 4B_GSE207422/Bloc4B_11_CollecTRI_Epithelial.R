#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 11: CollecTRI TF activity on epithelial cells
# Infers transcription factor activity using decoupleR + CollecTRI
# Compares TF programs across 5 final groups:
# tumor_MPR, normal_MPR, tumor_NMPR, normal_NMPR, Ciliated
# TF selection: between-group variance, consistent with CD8/TAM methodology
# Input:  Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/CollecTRI_Epithelial/
#         Results/Tables/Bloc4B_11_CollecTRI_Epithelial_TF_activity.csv
#         Results/Tables/Bloc4B_11_CollecTRI_Epithelial_wilcox_results.csv
# Reference: Müller-Dott et al., Nucleic Acids Research 2023
#            Badia-i-Mompel et al., Bioinformatics Advances 2022
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(decoupleR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(pheatmap)
  library(patchwork)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/CollecTRI_Epithelial")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# 1) Load epithelial object with final groups
message("Loading epithelial object...")
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds"))
seu_Epi <- subset(seu_Epi, cells = colnames(seu_Epi)[!is.na(seu_Epi$final_group)])
message("Cells with final_group: ", ncol(seu_Epi))
print(table(seu_Epi$final_group))
DefaultAssay(seu_Epi) <- "RNA"

# 2) Load CollecTRI network
message("Loading CollecTRI network...")
net <- read.csv(file.path(DATA_DIR, "Data/collectri_network.csv"))
message("CollecTRI: ", nrow(net), " interactions, ",
        length(unique(net$source)), " TFs")

# 3) Extract normalized count matrix
message("Extracting count matrix...")
mat <- as(GetAssayData(seu_Epi, layer = "data"), "dgCMatrix")

# 4) Run ULM
message("Running decoupleR ULM...")
tf_acts <- run_ulm(mat     = mat,
                   net     = net,
                   .source = "source",
                   .target = "target",
                   .mor    = "mor",
                   minsize = 5)
message("TF activity computed for ", length(unique(tf_acts$source)), " TFs")

saveRDS(tf_acts, file.path(DATA_DIR, "Objects/Bloc4B_11_tf_acts_raw.rds"))
message("Saved: Objects/Bloc4B_11_tf_acts_raw.rds")

tf_scores <- tf_acts %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score)

meta <- seu_Epi@meta.data %>%
  dplyr::select(final_group, PathResponse) %>%
  tibble::rownames_to_column("condition")

# 5) Compute mean TF activity per final_group, for ALL TFs
message("Computing mean TF activity per final_group (all TFs)...")

tf_summary_all <- tf_scores %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(final_group)) %>%
  group_by(source, final_group) %>%
  summarise(mean_activity = mean(score), .groups = "drop")

# 6) Select TFs by BETWEEN-GROUP variance (variance of final_group means)
message("Selecting TFs by between-group variance...")

tf_var_between <- tf_summary_all %>%
  group_by(source) %>%
  summarise(variance = var(mean_activity), .groups = "drop") %>%
  arrange(desc(variance))

top20_tfs <- tf_var_between$source[1:20]
key_tfs   <- tf_var_between$source[1:6]
message("Top 20 TFs (between-group variance): ", paste(top20_tfs, collapse = ", "))
message("Top 6 TFs (between-group variance): ", paste(key_tfs, collapse = ", "))

tf_summary <- tf_summary_all %>%
  filter(source %in% top20_tfs)

fwrite(as.data.frame(tf_summary),
       file.path(OUT_TAB, "Bloc4B_11_CollecTRI_Epithelial_TF_activity.csv"))
message("Saved: CSV")

# 7) Wilcoxon tests — 3 comparisons, BH-corrected within each,
#    consistent with the UCell Hallmark comparisons for this compartment
message("Running Wilcoxon tests for the 3 comparisons...")

run_comparison <- function(scores, meta, filter_col, filter_vals, group_col, level_a, level_b, label) {
  full <- scores %>% left_join(meta, by = "condition") %>%
    filter(!!rlang::sym(filter_col) %in% filter_vals, source %in% top20_tfs)
  results <- lapply(top20_tfs, function(tf) {
    a <- full$score[full[[group_col]] == level_a & full$source == tf]
    b <- full$score[full[[group_col]] == level_b & full$source == tf]
    if (length(a) < 3 | length(b) < 3) return(NULL)
    test <- wilcox.test(a, b)
    data.frame(comparison = label, source = tf, p_value = test$p.value,
               median_A = median(a), median_B = median(b))
  }) %>% bind_rows()
  names(results)[names(results) == "median_A"] <- paste0("median_", level_a)
  names(results)[names(results) == "median_B"] <- paste0("median_", level_b)
  results
}

res_tumor    <- run_comparison(tf_scores, meta, "final_group", c("tumor_MPR","tumor_NMPR"),
                               "final_group", "tumor_MPR", "tumor_NMPR", "tumor_MPR_vs_NMPR")
res_normal   <- run_comparison(tf_scores, meta, "final_group", c("normal_MPR","normal_NMPR"),
                               "final_group", "normal_MPR", "normal_NMPR", "normal_MPR_vs_NMPR")
res_ciliated <- run_comparison(tf_scores, meta, "final_group", "Ciliated",
                               "PathResponse", "MPR", "NMPR", "Ciliated_MPR_vs_NMPR")

wilcox_results <- bind_rows(res_tumor, res_normal, res_ciliated) %>%
  group_by(comparison) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(comparison, p_adj)

fwrite(wilcox_results,
       file.path(OUT_TAB, "Bloc4B_11_CollecTRI_Epithelial_wilcox_results.csv"))
message("Saved: Bloc4B_11_CollecTRI_Epithelial_wilcox_results.csv")

for (comp in unique(wilcox_results$comparison)) {
  n_sig <- sum(wilcox_results$comparison == comp & wilcox_results$p_adj < 0.05, na.rm = TRUE)
  message(comp, ": ", n_sig, "/", length(top20_tfs), " TFs significant (p_adj < 0.05)")
}

# 8) Heatmap, 20 TFs x 5 groups (descriptive, row-scaled)
message("Generating heatmap...")

green_palette <- colorRampPalette(c("white", "#006400"))(100)

tf_heatmap <- tf_summary %>%
  pivot_wider(names_from = final_group, values_from = mean_activity) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

col_order <- c("tumor_MPR", "normal_MPR", "tumor_NMPR", "normal_NMPR", "Ciliated")
col_order <- col_order[col_order %in% colnames(tf_heatmap)]
tf_heatmap <- tf_heatmap[, col_order]

png(file.path(OUT_FIG, "Bloc4B_11_CollecTRI_Epithelial_heatmap.png"),
    width = 12, height = 10, units = "in", res = 300)
pheatmap(tf_heatmap,
         scale        = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color        = green_palette,
         main         = "TF activity — Epithelial cells GSE207422",
         fontsize_row = 11,
         fontsize_col = 13,
         angle_col    = 45,
         border_color = NA)
dev.off()
message("Saved: Heatmap")

# 9) Violin plots: top 6 TFs by between-group variance
message("Generating violin plots...")

group_colors <- c(
  "tumor_MPR"   = "#E63946",
  "normal_MPR"  = "#457B9D",
  "tumor_NMPR"  = "#F4A261",
  "normal_NMPR" = "#2A9D8F",
  "Ciliated"    = "#ADB5BD"
)

tf_key <- tf_scores %>%
  filter(source %in% key_tfs) %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(final_group))

tf_key_tumor <- tf_key %>% filter(final_group %in% c("tumor_MPR", "tumor_NMPR")) %>%
  mutate(final_group = factor(final_group, levels = c("tumor_MPR", "tumor_NMPR")))

plot_list_tumor <- lapply(key_tfs, function(tf) {
  tf_key_tumor %>% filter(source == tf) %>%
    ggplot(aes(x = final_group, y = score, fill = final_group)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = group_colors) +
    theme_bw() +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
    labs(title = tf, y = "TF activity (ULM score)")
})
p_tumor <- wrap_plots(plot_list_tumor, ncol = 3)
ggsave(file.path(OUT_FIG, "Bloc4B_11_CollecTRI_Epithelial_top6_violin_tumor.png"),
       p_tumor, width = 14, height = 8, dpi = 300, bg = "white")

tf_key_normal <- tf_key %>% filter(final_group %in% c("normal_MPR", "normal_NMPR", "Ciliated")) %>%
  mutate(final_group = factor(final_group, levels = c("normal_MPR", "normal_NMPR", "Ciliated")))

plot_list_normal <- lapply(key_tfs, function(tf) {
  tf_key_normal %>% filter(source == tf) %>%
    ggplot(aes(x = final_group, y = score, fill = final_group)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = group_colors) +
    theme_bw() +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
    labs(title = tf, y = "TF activity (ULM score)")
})
p_normal <- wrap_plots(plot_list_normal, ncol = 3)
ggsave(file.path(OUT_FIG, "Bloc4B_11_CollecTRI_Epithelial_top6_violin_normal.png"),
       p_normal, width = 14, height = 8, dpi = 300, bg = "white")

# 10) Save object
saveRDS(seu_Epi, file.path(DATA_DIR, "Objects/Bloc4B_11_seu_Epithelial_TF.rds"))
message("Saved: Objects/Bloc4B_11_seu_Epithelial_TF.rds")

message("DONE CollecTRI Epithelial")