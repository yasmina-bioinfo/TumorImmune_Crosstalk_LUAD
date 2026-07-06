#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 10c: Hallmark UCell preprint figures
# Biologically filtered discriminant signatures 
# Figure principale : heatmap 21 signatures
# Input:  Objects/Bloc3_GSE243013_12c_seu_CD8_Hallmark.rds
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
                             "Objects/Bloc3_GSE243013_12c_seu_CD8_Hallmark.rds"))

# 2) Subset TEX and TPEX, MPR vs non-MPR only
seu_sub <- subset(seu_CD8,
                  subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX") &
                    pathological_response %in% c("MPR", "non-MPR"))

meta_sub <- seu_sub@meta.data %>%
  select(response  = pathological_response,
         cd8_state = functional.cluster,
         starts_with("HALLMARK"))

# 3) Biologically filtered signatures — 21 signatures
sig_filtered <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_HYPOXIA",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_COMPLEMENT"
)

# 4) Build heatmap matrix (no significance labels — significance is reported
# separately in Results text using Script 10b's state-specific _bysubtype table)

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
           "CD8.TEX\nMPR", "CD8.TEX\nnon-MPR",
           "CD8.TPEX\nMPR", "CD8.TPEX\nnon-MPR"))) %>%
  pivot_wider(names_from = group, values_from = mean_score) %>%
  tibble::column_to_rownames("signature") %>%
  as.matrix()

labels_row_clean <- gsub("HALLMARK_", "", sig_filtered)

# Column annotation
col_anno <- data.frame(
  Response = c("MPR", "non-MPR", "MPR", "non-MPR"),
  row.names = colnames(mat_heatmap))

anno_colors <- list(
  Response = c("MPR" = "#4393C3", "non-MPR" = "#D73027"))

green_palette <- colorRampPalette(c("white", "#006400"))(100)

png(file.path(OUT_FIG_PREPRINT, "Bloc3_GSE243013_Hallmark_CD8_heatmap_main.png"),
    width = 10, height = 10, units = "in", res = 300)
pheatmap(mat_heatmap,
         scale             = "row",
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         color             = green_palette,
         annotation_col    = col_anno,
         annotation_colors = anno_colors,
         labels_row        = labels_row_clean,
         main = "Hallmark UCell scores — CD8 T cells GSE243013",
         fontsize_row      = 10,
         fontsize_col      = 10,
         angle_col         = 45,
         border_color      = NA)
dev.off()

message("Saved: main heatmap 21 signatures with p-values")

message("DONE Bloc3 Script 10c")