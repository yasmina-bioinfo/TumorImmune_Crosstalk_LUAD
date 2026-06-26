#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 12b: CollecTRI figures MPR vs non-MPR
# Produces heatmap and violin plots from existing TF activity CSV
# pCR excluded from main preprint narrative
# Input:  Results/Tables/Bloc3_CollecTRI_TF_activity.csv
# Output: Results/Figures/CD8/Bloc3_CollecTRI_heatmap_MPR_nonMPR.png
#         Results/Figures/CD8/Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR.png
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(pheatmap)
  library(patchwork)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load existing TF activity results
message("Loading CollecTRI TF activity...")
tf_summary <- fread(file.path(OUT_TAB, "Bloc3_CollecTRI_TF_activity.csv"))

# 2) Filter MPR and non-MPR only
tf_summary <- tf_summary %>% filter(response %in% c("MPR", "non-MPR"))
message("Conditions retained: ", paste(unique(tf_summary$response), collapse=", "))

# 3) Top 20 TFs by variance
tf_var <- tf_summary %>%
  group_by(source) %>%
  summarise(variance = var(mean_activity), .groups = "drop") %>%
  arrange(desc(variance))

top20_tfs <- tf_var$source[1:20]
message("Top 20 TFs: ", paste(top20_tfs, collapse = ", "))

tf_summary <- tf_summary %>% filter(source %in% top20_tfs)

# 4) Heatmap
tf_heatmap <- tf_summary %>%
  mutate(col_label = paste0(cd8_state, "\n", response)) %>%
  select(source, col_label, mean_activity) %>%
  pivot_wider(names_from = col_label, values_from = mean_activity) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

cd8_states <- c("CD8.TEX", "CD8.TPEX", "CD8.EM")
conditions <- c("non-MPR", "MPR")
col_order  <- as.vector(outer(cd8_states, conditions,
                              function(s, c) paste0(s, "\n", c)))
col_order  <- col_order[col_order %in% colnames(tf_heatmap)]
tf_heatmap <- tf_heatmap[, col_order]

col_anno <- data.frame(
  Response  = sub(".*\n", "", colnames(tf_heatmap)),
  row.names = colnames(tf_heatmap)
)
anno_colors <- list(
  Response = c("MPR" = "#4393C3", "non-MPR" = "#D73027")
)

green_palette <- colorRampPalette(c("white", "#006400"))(100)

png(file.path(OUT_FIG, "Bloc3_CollecTRI_heatmap_MPR_nonMPR.png"),
    width = 12, height = 10, units = "in", res = 300)
pheatmap(tf_heatmap,
         scale             = "row",
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         color             = green_palette,
         annotation_col    = col_anno,
         annotation_colors = anno_colors,
         main              = "TF activity — CD8 T cells GSE243013 (MPR vs non-MPR)",
         fontsize_row      = 11,
         fontsize_col      = 11,
         angle_col         = 45,
         border_color      = NA,
         gaps_col          = seq(2, ncol(tf_heatmap) - 2, by = 2))
dev.off()
message("Saved: heatmap")

# 5) Violin plots top 6 TFs
key_tfs <- tf_var$source[1:6]

# Reload cell-level scores from CSV if available, otherwise use mean
# Here we use mean_activity per state x response as proxy
plot_list <- lapply(key_tfs, function(tf) {
  tf_summary %>%
    filter(source == tf,
           cd8_state %in% c("CD8.TEX", "CD8.TPEX")) %>%
    ggplot(aes(x = response, y = mean_activity, fill = response)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c("MPR" = "#4393C3", "non-MPR" = "#D73027")) +
    facet_wrap(~cd8_state) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(title = tf, y = "Mean TF activity")
})

p_key <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR.png"),
       p_key, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved: barplot TFs")

message("DONE Bloc3 Script 12b")