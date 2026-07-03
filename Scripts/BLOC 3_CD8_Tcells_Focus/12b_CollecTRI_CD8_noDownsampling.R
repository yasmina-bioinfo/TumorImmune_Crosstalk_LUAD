#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 12: CollecTRI TF activity on CD8 T cells
# Infers transcription factor activity using decoupleR + CollecTRI
# Tests H1: intrinsic TF differences between MPR, non-MPR, pCR
# CollecTRI network loaded from local CSV (OmnipathR server issues in WSL)
# Final: CD8.TEX, CD8.TPEX only — no downsampling
# CD8.EM excluded (secondary to narrative, RAM constraint)
# Input:  Objects/Bloc3_08_seu_CD8_ProjecTILs.rds
# Output: Results/Figures/CD8/Bloc3_CollecTRI_heatmap.png
#         Results/Figures/CD8/Bloc3_CollecTRI_key_TFs_violin.png
#         Results/Tables/Bloc3_CollecTRI_TF_activity.csv
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
DATA_DIR         <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ           <- file.path(DATA_DIR, "Objects/Bloc3_08_seu_CD8_ProjecTILs.rds")
OUT_FIG_PREPRINT <- file.path(DATA_DIR, "Results/Figures/CD8/Preprint")
dir.create(OUT_FIG_PREPRINT, recursive = TRUE, showWarnings = FALSE)
OUT_TAB          <- file.path(DATA_DIR, "Results/Tables")

# 1) Load CD8 ProjecTILs object
message("Loading CD8 ProjecTILs object...")
seu_CD8 <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_CD8))
DefaultAssay(seu_CD8) <- "RNA"

# 2) Load CollecTRI network
# NOTE: network downloaded from OmnipathR via decoupleR in RStudio Windows
# and saved locally to avoid OmnipathR server issues in WSL
message("Loading CollecTRI network...")
net <- read.csv(file.path(DATA_DIR, "Data/collectri_network.csv"))
message("CollecTRI: ", nrow(net), " interactions, ",
        length(unique(net$source)), " TFs")

seu_sub <- subset(seu_CD8, subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX"))

# 3) Extract normalized count matrix
# NOTE: decoupleR requires a genes x cells matrix
# Use normalized RNA counts (not scaled) for run_ulm
message("Extracting count matrix...")
mat <- as(GetAssayData(seu_sub, layer = "data"), "CsparseMatrix")

# 4) Run ULM (Univariate Linear Model) for TF activity inference
# NOTE: run_ulm is the recommended method for CollecTRI
# minsize = 5: minimum number of targets per TF
message("Running decoupleR ULM...")
tf_acts <- run_ulm(mat     = mat,
                   net     = net,
                   .source = "source",
                   .target = "target",
                   .mor    = "mor",
                   minsize = 5)

message("TF activity computed for ", length(unique(tf_acts$source)), " TFs")

# 5) Select top 20 TFs by variance across cells
# NOTE: top 20 for clarity, covers most informative TFs
# Avoids confirmation bias by selecting by variance, not by hypothesis
message("Selecting top 20 TFs by variance...")
tf_scores <- tf_acts %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score)

tf_var <- tf_scores %>%
  group_by(source) %>%
  summarise(variance = var(score), .groups = "drop") %>%
  arrange(desc(variance))

top20_tfs <- tf_var$source[1:20]
message("Top 20 TFs: ", paste(top20_tfs, collapse = ", "))

# 6) Compute mean TF activity per CD8 state x response group
message("Computing mean TF activity per CD8 state x response group...")

meta <- seu_sub@meta.data %>%
  select(response  = pathological_response,
         cd8_state = functional.cluster) %>%
  tibble::rownames_to_column("condition")

tf_summary <- tf_scores %>%
  filter(source %in% top20_tfs) %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(cd8_state)) %>%
  group_by(source, cd8_state, response) %>%
  summarise(mean_activity = mean(score), .groups = "drop")

# Save summary table
fwrite(as.data.frame(tf_summary),
       file.path(OUT_TAB, "Bloc3_CollecTRI_TF_activity.csv"))
message("Saved: Bloc3_CollecTRI_TF_activity.csv")

# 7) Build heatmap matrix
# Rows = TFs, Columns = CD8 state x condition
# Format: CD8state_condition (e.g. CD8.TEX_MPR)
tf_heatmap <- tf_summary %>%
  mutate(col_label = paste0(cd8_state, "\n", response)) %>%
  select(source, col_label, mean_activity) %>%
  pivot_wider(names_from = col_label, values_from = mean_activity) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# Order columns by CD8 state then condition
cd8_states   <- c("CD8.CM", "CD8.EM", "CD8.TEX", "CD8.TPEX",
                  "CD8.NaiveLike", "CD8.TEMRA", "CD8.MAIT")
conditions   <- c("non-MPR", "MPR", "pCR")
col_order    <- as.vector(outer(cd8_states, conditions,
                                function(s, c) paste0(s, "\n", c)))
col_order    <- col_order[col_order %in% colnames(tf_heatmap)]
tf_heatmap   <- tf_heatmap[, col_order]

# Annotation for columns: response group on top
col_anno <- data.frame(
  Response = sub(".*\n", "", colnames(tf_heatmap)),
  row.names = colnames(tf_heatmap)
)
anno_colors <- list(
  Response = c("MPR"     = "#4393C3",
               "non-MPR" = "#D73027",
               "pCR"     = "#1A7A1A")
)

# 8) Generate heatmap
# White to dark green palette, consistent with portfolio
message("Generating heatmap...")
green_palette <- colorRampPalette(c("white", "#006400"))(100)

png(file.path(OUT_FIG_PREPRINT, "Bloc3_CollecTRI_heatmap.png"),
    width = 16, height = 10, units = "in", res = 300)
pheatmap(tf_heatmap,
         scale            = "row",
         cluster_cols     = FALSE,
         cluster_rows     = TRUE,
         color            = green_palette,
         annotation_col   = col_anno,
         annotation_colors = anno_colors,
         main             = "TF activity (CollecTRI/decoupleR) — CD8 T cells GSE243013",
         fontsize_row     = 11,
         fontsize_col     = 10,
         angle_col        = 45,
         border_color     = NA,
         gaps_col         = seq(3, ncol(tf_heatmap) - 3, by = 3))
dev.off()
message("Saved: Bloc3_CollecTRI_heatmap.png")

# 8b) Heatmap MPR vs non-MPR only (preprint figure)
tf_heatmap_2cond <- tf_summary %>%
  filter(response %in% c("MPR", "non-MPR")) %>%
  mutate(col_label = paste0(cd8_state, "\n", response)) %>%
  select(source, col_label, mean_activity) %>%
  pivot_wider(names_from = col_label, values_from = mean_activity) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# Order columns
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

# 9) Violin plot: Top 6 TFs by variance (objective selection, no confirmation bias)
message("Generating violin plots for key TFs...")
key_tfs_present <- tf_var$source[1:6]
message("Top 10 TFs: ", paste(key_tfs_present, collapse = ", "))

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
ggsave(file.path(OUT_FIG_PREPRINT, "Bloc3_CollecTRI_key_TFs_violin.png"),
       p_key, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved: Bloc3_CollecTRI_key_TFs_violin.png")

# 9b) Violin plot MPR vs non-MPR only (preprint)
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