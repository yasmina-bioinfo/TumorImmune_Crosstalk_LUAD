#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 12: CollecTRI TF activity on CD8 T cells
# Infers transcription factor activity using decoupleR + CollecTRI
# Tests H1: intrinsic TF differences between MPR, non-MPR, pCR
# Compares with portfolio findings (STAT2 in non-MPR, ELK4 in MPR)
# NOTE: RAM constraint, downsampling strategy documented for reproducibility:
# Attempt 1: Full CD8 object (57,587 cells) → 14.6 GiB required → Killed
# Attempt 2: TEX + TPEX + EM + CM subset (56,899 cells) → 14.5 GiB → Killed (WSL)
# Attempt 3: TEX + TPEX + EM, downsample 16,000/state → ~33,496 cells → Killed (WSL)
# Attempt 4: TEX + TPEX + EM, downsample 10,000/state → ~21,821 cells
#            → RStudio Windows (more stable memory management) → SUCCESS
# Final: CD8.TEX, CD8.TPEX, CD8.EM, max 10,000 cells/state, set.seed(42)
# Full analysis on all 7 CD8 states requires a compute server with >32GB RAM
# CollecTRI network loaded from local CSV (OmnipathR server issues in WSL)
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
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_08_seu_CD8_ProjecTILs.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

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

# 2b) Subset to key CD8 populations for RAM efficiency
# NOTE: Full CD8 object (57,587 cells) requires 14.5 GiB for run_ulm
# Restricted to CD8.TEX, CD8.TPEX, CD8.EM, most relevant to narrative
# Downsampled to max 10,000 cells per state, set.seed(42) for reproducibility
seu_sub <- subset(seu_CD8, subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX", "CD8.EM"))
set.seed(42)
cells_keep <- seu_sub@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::group_by(functional.cluster) %>%
  dplyr::slice_sample(n = 10000) %>%
  dplyr::pull(barcode)
seu_sub <- subset(seu_sub, cells = cells_keep)
message("Downsampled: ", ncol(seu_sub), " cells")

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

png(file.path(OUT_FIG, "Bloc3_CollecTRI_heatmap.png"),
    width = 16, height = 10, units = "in", res = 300)
pheatmap(tf_heatmap,
         scale            = "row",
         cluster_cols     = FALSE,
         cluster_rows     = TRUE,
         color            = green_palette,
         annotation_col   = col_anno,
         annotation_colors = anno_colors,
         main             = "TF activity (CollecTRI/decoupleR) — CD8 T cells GSE243013 LUAD",
         fontsize_row     = 11,
         fontsize_col     = 10,
         angle_col        = 45,
         border_color     = NA,
         gaps_col         = seq(3, ncol(tf_heatmap) - 3, by = 3))
dev.off()
message("Saved: Bloc3_CollecTRI_heatmap.png")

# 9) Violin plot, key TFs from portfolio (STAT2, ELK4, ELK1, TBX21)
# Focus on CD8.TEX and CD8.TPEX, most relevant for narrative
message("Generating violin plots for key TFs...")

key_tfs <- c("STAT2", "ELK4", "ELK1", "TBX21", "IRF7", "STAT1")
key_tfs_present <- key_tfs[key_tfs %in% unique(tf_acts$source)]
message("Key TFs present: ", paste(key_tfs_present, collapse = ", "))

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
    scale_fill_manual(values = c("MPR"     = "#D73027",
                                 "non-MPR" = "#4393C3",
                                 "pCR"     = "#1A7A1A")) +
    facet_wrap(~cd8_state) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(title = tf, y = "TF activity (ULM score)")
})

p_key <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc3_CollecTRI_key_TFs_violin.png"),
       p_key, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved: Bloc3_CollecTRI_key_TFs_violin.png")

# 10) Save updated object
message("Saving updated CD8 object...")
saveRDS(seu_CD8, file.path(DATA_DIR, "Objects/Bloc3_12_seu_CD8_TF.rds"))
message("Saved: Objects/Bloc3_12_seu_CD8_TF.rds")
message("DONE Bloc3 Script 12")