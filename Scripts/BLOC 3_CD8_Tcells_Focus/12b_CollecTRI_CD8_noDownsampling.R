#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 12: CollecTRI TF activity on CD8 T cells
# Infers transcription factor activity using decoupleR + CollecTRI
# Tests H1: intrinsic TF differences between MPR, non-MPR, pCR
# CollecTRI network loaded from local CSV (OmnipathR server issues in WSL)
# Final: CD8.TEX, CD8.TPEX only — no downsampling
# CD8.EM excluded (secondary to narrative, RAM constraint)
# pCR retained in metadata and in the exploratory 3-condition heatmap only —
# excluded from the preprint figures and from the statistical test (Wilcoxon),
# reserved for the Perspectives section (see pseudobulk analyses, Script 06/13)
# Input:  Objects/Bloc3_08_seu_CD8_ProjecTILs.rds
# Output: Results/Figures/CD8/Exploratory/Bloc3_CollecTRI_heatmap_MPR_nonMPR_pCR.png
#         Results/Figures/CD8/Exploratory/Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR_pCR.png
#         Results/Figures/CD8/Preprint/Bloc3_CollecTRI_heatmap_MPR_nonMPR.png
#         Results/Figures/CD8/Preprint/Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR.png
#         Results/Tables/Bloc3_CollecTRI_TF_activity.csv
#         Results/Tables/Bloc3_CollecTRI_TF_wilcox_bystate.csv
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

# Paths
DATA_DIR            <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
#DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ              <- file.path(DATA_DIR, "Objects/Bloc3_08_seu_CD8_ProjecTILs.rds")
OUT_FIG_PREPRINT    <- file.path(DATA_DIR, "Results/Figures/CD8/Preprint")
OUT_FIG_EXPLORATORY <- file.path(DATA_DIR, "Results/Figures/CD8/Exploratory")
dir.create(OUT_FIG_PREPRINT, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_EXPLORATORY, recursive = TRUE, showWarnings = FALSE)
OUT_TAB             <- file.path(DATA_DIR, "Results/Tables")

# 1) Load CD8 ProjecTILs object
message("Loading CD8 ProjecTILs object...")
seu_CD8 <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_CD8))
DefaultAssay(seu_CD8) <- "RNA"

# 2) Load CollecTRI network
message("Loading CollecTRI network...")
net <- read.csv(file.path(DATA_DIR, "Data/collectri_network.csv"))
message("CollecTRI: ", nrow(net), " interactions, ",
        length(unique(net$source)), " TFs")

seu_sub <- subset(seu_CD8, subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX"))

# 3) Extract normalized count matrix
message("Extracting count matrix...")
mat <- as(GetAssayData(seu_sub, layer = "data"), "CsparseMatrix")

# 4) Run ULM (Univariate Linear Model) for TF activity inference
message("Running decoupleR ULM...")
tf_acts <- run_ulm(mat     = mat,
                   net     = net,
                   .source = "source",
                   .target = "target",
                   .mor    = "mor",
                   minsize = 5)
message("TF activity computed for ", length(unique(tf_acts$source)), " TFs")

# Save raw per-cell TF scores in case a future adjustment is needed
# without re-running decoupleR (the most expensive step)
saveRDS(tf_acts, file.path(DATA_DIR, "Objects/Bloc3_12_tf_acts_raw.rds"))
message("Saved: Objects/Bloc3_12_tf_acts_raw.rds")

tf_scores <- tf_acts %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score)

meta <- seu_sub@meta.data %>%
  select(response  = pathological_response,
         cd8_state = functional.cluster) %>%
  tibble::rownames_to_column("condition")

# 5) Compute mean TF activity per CD8 state x response group — for ALL TFs
#    (done BEFORE selection, so selection is based on between-group variance
#    rather than per-cell variance — aligned with the project's goal of
#    comparing response groups)
message("Computing mean TF activity per CD8 state x response group (all TFs)...")

tf_summary_all <- tf_scores %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(cd8_state)) %>%
  group_by(source, cd8_state, response) %>%
  summarise(mean_activity = mean(score), .groups = "drop")

# 6) Select TFs by BETWEEN-GROUP variance (variance of state x response means)
message("Selecting TFs by between-group variance...")

tf_var_between <- tf_summary_all %>%
  group_by(source) %>%
  summarise(variance = var(mean_activity), .groups = "drop") %>%
  arrange(desc(variance))

top20_tfs <- tf_var_between$source[1:20]
key_tfs_present <- tf_var_between$source[1:6]
message("Top 20 TFs (between-group variance): ", paste(top20_tfs, collapse = ", "))
message("Top 6 TFs (between-group variance): ", paste(key_tfs_present, collapse = ", "))

tf_summary <- tf_summary_all %>%
  filter(source %in% top20_tfs)

fwrite(as.data.frame(tf_summary),
       file.path(OUT_TAB, "Bloc3_CollecTRI_TF_activity.csv"))
message("Saved: Bloc3_CollecTRI_TF_activity.csv")

# 6b) Wilcoxon test per CD8 state (TEX, TPEX), MPR vs non-MPR, on raw per-cell
#     scores — BH-corrected within state. Restricted to the Top 20 TFs and to
#     MPR/non-MPR (pCR reserved for perspectives, not part of the main text).
message("Running Wilcoxon tests MPR vs non-MPR, per CD8 state, on raw TF scores...")

tf_scores_meta <- tf_scores %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(cd8_state),
         response %in% c("MPR", "non-MPR"),
         source %in% top20_tfs)

wilcox_bystate <- lapply(c("CD8.TEX", "CD8.TPEX"), function(state) {
  lapply(top20_tfs, function(tf) {
    mpr    <- tf_scores_meta$score[tf_scores_meta$cd8_state == state &
                                     tf_scores_meta$response  == "MPR" &
                                     tf_scores_meta$source    == tf]
    nonmpr <- tf_scores_meta$score[tf_scores_meta$cd8_state == state &
                                     tf_scores_meta$response  == "non-MPR" &
                                     tf_scores_meta$source    == tf]
    if (length(mpr) < 3 | length(nonmpr) < 3) return(NULL)
    test <- wilcox.test(mpr, nonmpr)
    data.frame(cd8_state     = state,
               source        = tf,
               p_value       = test$p.value,
               median_MPR    = median(mpr),
               median_nonMPR = median(nonmpr))
  }) %>% bind_rows()
}) %>% bind_rows()

wilcox_bystate <- wilcox_bystate %>%
  group_by(cd8_state) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(cd8_state, p_adj)

fwrite(wilcox_bystate,
       file.path(OUT_TAB, "Bloc3_CollecTRI_TF_wilcox_bystate.csv"))
message("Saved: Bloc3_CollecTRI_TF_wilcox_bystate.csv")
print(wilcox_bystate)

# 7) Build heatmap matrix — EXPLORATORY, 3 conditions (MPR/non-MPR/pCR)
#    NOT part of the preprint — kept on disk only, per project convention
tf_heatmap <- tf_summary %>%
  mutate(col_label = paste0(cd8_state, "\n", response)) %>%
  select(source, col_label, mean_activity) %>%
  pivot_wider(names_from = col_label, values_from = mean_activity) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

cd8_states   <- c("CD8.TEX", "CD8.TPEX")
conditions   <- c("non-MPR", "MPR", "pCR")
col_order    <- as.vector(outer(cd8_states, conditions,
                                function(s, c) paste0(s, "\n", c)))
col_order    <- col_order[col_order %in% colnames(tf_heatmap)]
tf_heatmap   <- tf_heatmap[, col_order]

col_anno <- data.frame(
  Response = sub(".*\n", "", colnames(tf_heatmap)),
  row.names = colnames(tf_heatmap)
)
anno_colors <- list(
  Response = c("MPR"     = "#4393C3",
               "non-MPR" = "#D73027",
               "pCR"     = "#1A7A1A")
)

message("Generating exploratory heatmap (3 conditions)...")
green_palette <- colorRampPalette(c("white", "#006400"))(100)

png(file.path(OUT_FIG_EXPLORATORY, "Bloc3_CollecTRI_heatmap_MPR_nonMPR_pCR.png"),
    width = 16, height = 10, units = "in", res = 300)
pheatmap(tf_heatmap,
         scale            = "row",
         cluster_cols     = FALSE,
         cluster_rows     = TRUE,
         color            = green_palette,
         annotation_col   = col_anno,
         annotation_colors = anno_colors,
         main             = "TF activity (CollecTRI/decoupleR) — CD8 T cells GSE243013 [EXPLORATORY]",
         fontsize_row     = 11,
         fontsize_col     = 10,
         angle_col        = 45,
         border_color     = NA,
         gaps_col         = seq(3, ncol(tf_heatmap) - 3, by = 3))
dev.off()
message("Saved: Bloc3_CollecTRI_heatmap_MPR_nonMPR_pCR.png (exploratory only)")

# 8) Heatmap MPR vs non-MPR only — PREPRINT figure
tf_heatmap_2cond <- tf_summary %>%
  filter(response %in% c("MPR", "non-MPR")) %>%
  mutate(col_label = paste0(cd8_state, "\n", response)) %>%
  select(source, col_label, mean_activity) %>%
  pivot_wider(names_from = col_label, values_from = mean_activity) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

col_order_2cond <- as.vector(outer(c("CD8.TEX", "CD8.TPEX"),
                                   c("non-MPR", "MPR"),
                                   function(s, c) paste0(s, "\n", c)))
col_order_2cond <- col_order_2cond[col_order_2cond %in% colnames(tf_heatmap_2cond)]
tf_heatmap_2cond <- tf_heatmap_2cond[, col_order_2cond]

col_anno_2cond <- data.frame(
  Response = sub(".*\n", "", colnames(tf_heatmap_2cond)),
  row.names = colnames(tf_heatmap_2cond)
)

anno_colors_2cond <- list(
  Response = c("MPR"     = "#4393C3",
               "non-MPR" = "#D73027")
)

message("Generating preprint heatmap (MPR vs non-MPR)...")
png(file.path(OUT_FIG_PREPRINT, "Bloc3_CollecTRI_heatmap_MPR_nonMPR.png"),
    width = 10, height = 10, units = "in", res = 300)
pheatmap(tf_heatmap_2cond,
         scale             = "row",
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         color             = green_palette,
         annotation_col    = col_anno_2cond,
         annotation_colors = anno_colors_2cond,
         main              = "TF activity (CollecTRI/decoupleR) — CD8 T cells GSE243013",
         fontsize_row      = 11,
         fontsize_col      = 10,
         angle_col         = 45,
         border_color      = NA)
dev.off()
message("Saved: Bloc3_CollecTRI_heatmap_MPR_nonMPR.png")

# 9) Violin plot: Top 6 TFs by between-group variance — EXPLORATORY, 3 conditions
message("Generating exploratory violin plots (3 conditions)...")

tf_key <- tf_scores %>%
  filter(source %in% key_tfs_present) %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(cd8_state),
         cd8_state %in% c("CD8.TEX", "CD8.TPEX"))

plot_list <- lapply(key_tfs_present, function(tf) {
  tf_key %>%
    filter(source == tf) %>%
    ggplot(aes(x = response, y = score, fill = response)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = c("MPR"     = "#4393C3",
                                 "non-MPR" = "#D73027",
                                 "pCR"     = "#1A7A1A")) +
    facet_wrap(~cd8_state) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(title = tf, y = "TF activity (ULM score)")
})

p_key <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(OUT_FIG_EXPLORATORY, "Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR_pCR.png"),
       p_key, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved: Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR_pCR.png (exploratory only)")

# 9b) Violin plot MPR vs non-MPR only — PREPRINT figure
message("Generating preprint violin plots (MPR vs non-MPR)...")

tf_key_2cond <- tf_scores %>%
  filter(source %in% key_tfs_present) %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(cd8_state),
         cd8_state %in% c("CD8.TEX", "CD8.TPEX"),
         response %in% c("MPR", "non-MPR"))

plot_list_2cond <- lapply(key_tfs_present, function(tf) {
  tf_key_2cond %>%
    filter(source == tf) %>%
    ggplot(aes(x = response, y = score, fill = response)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = c("MPR"     = "#4393C3",
                                 "non-MPR" = "#D73027")) +
    facet_wrap(~cd8_state) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(title = tf, y = "TF activity (ULM score)")
})

p_key_2cond <- wrap_plots(plot_list_2cond, ncol = 2)
ggsave(file.path(OUT_FIG_PREPRINT, "Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR.png"),
       p_key_2cond, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved: Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR.png")

# 10) Save updated object
message("Saving updated CD8 object...")
saveRDS(seu_CD8, file.path(DATA_DIR, "Objects/Bloc3_12b_seu_CD8_TF.rds"))
message("Saved: Objects/Bloc3_12b_seu_CD8_TF.rds")
message("DONE Bloc3 Script 12")