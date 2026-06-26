#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc4A Script 07b: CollecTRI heatmap TAMs MPR vs non-MPR
# pCR excluded from main preprint narrative
# Input:  Results/Tables/Bloc4A_CollecTRI_TF_activity.csv
# Output: Results/Figures/BLOC4A_TAMs/Bloc4A_CollecTRI_heatmap_MPR_nonMPR.png
# ============================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(pheatmap)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs")

# 1) Load and filter MPR vs non-MPR
tf_summary <- fread(file.path(DATA_DIR, "Results/Tables/Bloc4A_CollecTRI_TF_activity.csv"))
tf_summary <- tf_summary %>% filter(response %in% c("MPR", "non-MPR"))
message("Conditions: ", paste(unique(tf_summary$response), collapse=", "))

# 2) Top 20 TFs by variance
tf_var <- tf_summary %>%
  group_by(source) %>%
  summarise(variance = var(mean_activity), .groups = "drop") %>%
  arrange(desc(variance))

top20_tfs <- tf_var$source[1:20]
message("Top 20 TFs: ", paste(top20_tfs, collapse=", "))
tf_summary <- tf_summary %>% filter(source %in% top20_tfs)

# 3) Build heatmap matrix
tf_heatmap <- tf_summary %>%
  mutate(col_label = paste0(tam_state, "\n", response)) %>%
  select(source, col_label, mean_activity) %>%
  pivot_wider(names_from = col_label, values_from = mean_activity) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# Order columns
tam_states <- c("Resident M2", "LAMs", "Monocyte FCN1+",
                "Stress-response", "Proliferating",
                "IFN-stimulated", "Classical monocyte")
conditions <- c("non-MPR", "MPR")
col_order  <- as.vector(outer(tam_states, conditions,
                              function(s, c) paste0(s, "\n", c)))
col_order  <- col_order[col_order %in% colnames(tf_heatmap)]
tf_heatmap <- tf_heatmap[, col_order]

# 4) Column annotation
col_anno <- data.frame(
  Response  = sub(".*\n", "", colnames(tf_heatmap)),
  row.names = colnames(tf_heatmap)
)
anno_colors <- list(Response = c("MPR" = "#4393C3", "non-MPR" = "#D73027"))

# 5) Heatmap
green_palette <- colorRampPalette(c("white", "#006400"))(100)

png(file.path(OUT_FIG, "Bloc4A_CollecTRI_heatmap_MPR_nonMPR.png"),
    width = 14, height = 10, units = "in", res = 300)
pheatmap(tf_heatmap,
         scale             = "row",
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         color             = green_palette,
         annotation_col    = col_anno,
         annotation_colors = anno_colors,
         main              = "TF activity — TAMs GSE243013 (MPR vs non-MPR)",
         fontsize_row      = 11,
         fontsize_col      = 10,
         angle_col         = 45,
         border_color      = NA,
         gaps_col          = seq(2, ncol(tf_heatmap) - 2, by = 2))
dev.off()
message("Saved: Bloc4A_CollecTRI_heatmap_MPR_nonMPR.png")

message("DONE Bloc4A Script 07b")