#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc4A Script 05d: UCell Hallmark TAMs preprint figures
# Biologically filtered discriminant signatures — heatmap with p-values
# Input:  Objects/GSE243013_seu_TAMs_Hallmark.rds
# Output: Results/Figures/BLOC4A_TAMs/Preprint/
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

DATA_DIR         <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG_PREPRINT <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs/Preprint")
dir.create(OUT_FIG_PREPRINT, recursive = TRUE, showWarnings = FALSE)

# 1) Load object
message("Loading TAMs Hallmark object...")
seu_TAM <- readRDS(file.path(DATA_DIR,
                             "Objects/GSE243013_seu_TAMs_Hallmark.rds"))

# 2) Define filtered signatures
sig_filtered <- c(
  "HALLMARK_COMPLEMENT",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_G2M_CHECKPOINT"
)

# 3) Calculate p-values MPR vs non-MPR
pval_row <- sapply(sig_filtered, function(sig) {
  mpr  <- seu_TAM@meta.data[seu_TAM@meta.data$pathological_response == "MPR",  sig]
  nmpr <- seu_TAM@meta.data[seu_TAM@meta.data$pathological_response == "non-MPR", sig]
  mpr  <- mpr[!is.na(mpr)]
  nmpr <- nmpr[!is.na(nmpr)]
  wilcox.test(mpr, nmpr)$p.value
})

# Labels with asterisks
labels_row_sig <- paste0(gsub("HALLMARK_", "", sig_filtered), " ",
                         ifelse(pval_row < 0.001, "***",
                                ifelse(pval_row < 0.01,  "**",
                                       ifelse(pval_row < 0.05,  "*", "ns"))))

# 4) Build heatmap matrix
meta_sub <- seu_TAM@meta.data %>%
  select(response = pathological_response,
         all_of(sig_filtered))

mat_heatmap <- meta_sub %>%
  pivot_longer(cols = all_of(sig_filtered),
               names_to = "signature",
               values_to = "score") %>%
  group_by(response, signature) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(signature = gsub("HALLMARK_", "", signature),
         response = factor(response, levels = c("MPR", "non-MPR"))) %>%
  pivot_wider(names_from = response, values_from = mean_score) %>%
  tibble::column_to_rownames("signature") %>%
  as.matrix()

# Column annotation
col_anno <- data.frame(
  Response = c("MPR", "non-MPR"),
  row.names = colnames(mat_heatmap))

anno_colors <- list(
  Response = c("MPR" = "#4393C3", "non-MPR" = "#D73027"))

green_palette <- colorRampPalette(c("white", "#006400"))(100)

# 5) Heatmap
message("Generating heatmap...")

# 5a) Heatmaps global TAMs
png(file.path(OUT_FIG_PREPRINT, "GSE243013_Hallmark_TAMs_heatmap_main.png"),
    width = 8, height = 10, units = "in", res = 300)
pheatmap(mat_heatmap,
         scale             = "row",
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         color             = green_palette,
         annotation_col    = col_anno,
         annotation_colors = anno_colors,
         labels_row        = labels_row_sig,
         main              = "Hallmark UCell scores — TAMs GSE243013\n*** p<0.001  ** p<0.01  * p<0.05  ns not significant",
         fontsize_row      = 10,
         fontsize_col      = 11,
         angle_col         = 45,
         border_color      = NA)
dev.off()

message("Saved: GSE243013_Hallmark_TAMs_heatmap_main.png")

# 5b) Heatmap UCell Hallmark by subtype

# Mapping full names to short names
short_names <- c(
  "Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like)" = "Resident M2",
  "TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs)" = "LAMs",
  "Inflammatory monocyte-derived TAMs (FCN1+/S100A8+)" = "Monocyte FCN1+",
  "Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high)" = "Stress-response",
  "Proliferating TAMs (cycling/MKI67+)" = "Proliferating",
  "IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+)" = "IFN-stimulated",
  "Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+)" = "Classical-Mono"
)

meta_subtype <- seu_TAM@meta.data %>%
  select(subtype  = final_annotation,
         response = pathological_response,
         all_of(sig_filtered)) %>%
  mutate(subtype = short_names[subtype])

mat_subtype <- meta_subtype %>%
  pivot_longer(cols = all_of(sig_filtered),
               names_to = "signature",
               values_to = "score") %>%
  group_by(subtype, response, signature) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(signature = gsub("HALLMARK_", "", signature),
         group = paste0(subtype, "\n", response)) %>%
  select(signature, group, mean_score) %>%
  pivot_wider(names_from = group, values_from = mean_score) %>%
  tibble::column_to_rownames("signature") %>%
  as.matrix()

# Order columns
tam_states <- c("Resident M2", "LAMs", "Monocyte FCN1+",
                "Stress-response", "Proliferating",
                "IFN-stimulated", "Classical-Mono")
conditions <- c("non-MPR", "MPR")

col_order <- as.vector(outer(tam_states, conditions,
                             function(s, c) paste0(s, "\n", c)))
col_order <- col_order[col_order %in% colnames(mat_subtype)]
mat_subtype <- mat_subtype[, col_order]

# Column annotation
col_anno <- data.frame(
  Response  = sub(".*\n", "", colnames(mat_subtype)),
  row.names = colnames(mat_subtype))

anno_colors <- list(
  Response = c("MPR" = "#4393C3", "non-MPR" = "#D73027"))

png(file.path(OUT_FIG_PREPRINT,
              "GSE243013_Hallmark_TAMs_heatmap_bysubtype.png"),
    width = 16, height = 10, units = "in", res = 300)
pheatmap(mat_subtype,
         scale             = "row",
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         color             = green_palette,
         annotation_col    = col_anno,
         annotation_colors = anno_colors,
         main              = "Hallmark UCell scores — TAMs by subtype GSE243013",
         fontsize_row      = 10,
         fontsize_col      = 11,
         angle_col         = 45,
         border_color      = NA,
         gaps_col          = 7)
dev.off()
message("Saved: heatmap by subtype")

message("DONE Script 05d")