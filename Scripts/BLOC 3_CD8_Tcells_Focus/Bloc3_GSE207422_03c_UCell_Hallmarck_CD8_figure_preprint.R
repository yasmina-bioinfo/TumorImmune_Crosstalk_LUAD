#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc3 Script 03c: Hallmark UCell preprint figures
# Biologically filtered discriminant signatures
# Main figure : heatmap 21 signatures with p-values
# Input:  Objects/Bloc3_GSE207422_02b_seu_CD8_Hallmark.rds
# Output: Results/Figures/CD8/Preprint/
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
})

DATA_DIR         <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG_PREPRINT <- file.path(DATA_DIR, "Results/Figures/CD8/Preprint")
dir.create(OUT_FIG_PREPRINT, recursive = TRUE, showWarnings = FALSE)

# 1) Load object with Hallmark UCell scores
message("Loading CD8 object with Hallmark scores...")
seu_CD8 <- readRDS(file.path(DATA_DIR,
                             "Objects/Bloc3_GSE207422_02b_seu_CD8_Hallmark.rds"))

# 2) Subset TEX and TPEX, MPR vs NMPR only
seu_sub <- subset(seu_CD8,
                  subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX") &
                    PathResponse %in% c("MPR", "NMPR"))

meta_sub <- seu_sub@meta.data %>%
  select(response  = PathResponse,
         cd8_state = functional.cluster,
         starts_with("HALLMARK"))

# 3) Biologically filtered signatures — 22 signatures
sig_filtered <- c(
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_ALLOGRAFT_REJECTION",       
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_HYPOXIA",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

# 4) Calculate p-values per signature
pval_row <- sapply(sig_filtered, function(sig) {
  mpr  <- seu_sub@meta.data[seu_sub@meta.data$PathResponse == "MPR",  sig]
  nmpr <- seu_sub@meta.data[seu_sub@meta.data$PathResponse == "NMPR", sig]
  mpr  <- mpr[!is.na(mpr)]
  nmpr <- nmpr[!is.na(nmpr)]
  wilcox.test(mpr, nmpr)$p.value
})

# Labels avec astérisques
labels_row_sig <- paste0(gsub("HALLMARK_", "", sig_filtered), " ",
                         ifelse(pval_row < 0.001, "***",
                                ifelse(pval_row < 0.01,  "**",
                                       ifelse(pval_row < 0.05,  "*", "ns"))))

# 5) Build heatmap matrix
mat_heatmap <- meta_sub %>%
  select(response, cd8_state, all_of(sig_filtered)) %>%
  mutate(group = paste0(cd8_state, "\n", response)) %>%
  select(-response, -cd8_state) %>%
  pivot_longer(cols = all_of(sig_filtered),
               names_to = "signature",
               values_to = "score") %>%
  group_by(group, signature) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(signature = gsub("HALLMARK_", "", signature),
         group = factor(group, levels = c(
           "CD8.TEX\nMPR", "CD8.TEX\nNMPR",
           "CD8.TPEX\nMPR", "CD8.TPEX\nNMPR"))) %>%
  pivot_wider(names_from = group, values_from = mean_score) %>%
  tibble::column_to_rownames("signature") %>%
  as.matrix()

# Column annotation
col_anno <- data.frame(
  Response = c("MPR", "NMPR", "MPR", "NMPR"),
  row.names = colnames(mat_heatmap))

anno_colors <- list(
  Response = c("MPR" = "#4393C3", "NMPR" = "#D73027"))

green_palette <- colorRampPalette(c("white", "#006400"))(100)

# 6) Heatmap — main preprint figure
message("Generating main heatmap with p-values...")

png(file.path(OUT_FIG_PREPRINT, "Bloc3_GSE207422_Hallmark_CD8_heatmap_main.png"),
    width = 10, height = 10, units = "in", res = 300)
pheatmap(mat_heatmap,
         scale             = "row",
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         color             = green_palette,
         annotation_col    = col_anno,
         annotation_colors = anno_colors,
         labels_row        = labels_row_sig,
         main              = "Hallmark UCell scores — CD8 T cells GSE207422\n*** p<0.001  ** p<0.01  * p<0.05  ns not significant",
         fontsize_row      = 10,
         fontsize_col      = 10,
         angle_col         = 45,
         border_color      = NA)
dev.off()
message("Saved: main heatmap with 22 signatures and p-values")

message("DONE Bloc3 Script 03c")