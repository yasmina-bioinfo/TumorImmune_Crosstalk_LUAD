#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 05c: UCell Hallmark TAMs preprint figures
# Biologically filtered discriminant signatures — heatmap with p-values
# Input:  Objects/GSE207422_seu_TAMs_Hallmark.rds
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Preprint/
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

DATA_DIR         <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG_PREPRINT <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/Preprint")
dir.create(OUT_FIG_PREPRINT, recursive = TRUE, showWarnings = FALSE)

# 1) Load object
message("Loading TAMs Hallmark object...")
seu_TAM <- readRDS(file.path(DATA_DIR,
                             "Objects/GSE207422_seu_TAMs_Hallmark.rds"))

# 2) Define filtered signatures
sig_filtered <- c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_NOTCH_SIGNALING"
)

# 3) Calculate p-values MPR vs NMPR
pval_row <- sapply(sig_filtered, function(sig) {
  mpr  <- seu_TAM@meta.data[seu_TAM@meta.data$PathResponse == "MPR",  sig]
  nmpr <- seu_TAM@meta.data[seu_TAM@meta.data$PathResponse == "NMPR", sig]
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
  select(response = PathResponse,
         all_of(sig_filtered))

mat_heatmap <- meta_sub %>%
  pivot_longer(cols = all_of(sig_filtered),
               names_to = "signature",
               values_to = "score") %>%
  group_by(response, signature) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(signature = gsub("HALLMARK_", "", signature),
         response = factor(response, levels = c("MPR", "NMPR"))) %>%
  pivot_wider(names_from = response, values_from = mean_score) %>%
  tibble::column_to_rownames("signature") %>%
  as.matrix()

# Column annotation
col_anno <- data.frame(
  Response = c("MPR", "NMPR"),
  row.names = colnames(mat_heatmap))

anno_colors <- list(
  Response = c("MPR" = "#4393C3", "NMPR" = "#D73027"))

green_palette <- colorRampPalette(c("white", "#006400"))(100)

# 5) Heatmap
message("Generating heatmap...")
png(file.path(OUT_FIG_PREPRINT, "GSE207422_Hallmark_TAMs_heatmap_main.png"),
    width = 8, height = 10, units = "in", res = 300)
pheatmap(mat_heatmap,
         scale             = "row",
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         color             = green_palette,
         annotation_col    = col_anno,
         annotation_colors = anno_colors,
         labels_row        = labels_row_sig,
         main              = "Hallmark UCell scores — TAMs GSE207422\n*** p<0.001  ** p<0.01  * p<0.05  ns not significant",
         fontsize_row      = 10,
         fontsize_col      = 11,
         angle_col         = 45,
         border_color      = NA)
dev.off()
message("Saved: GSE207422_Hallmark_TAMs_heatmap_main.png")

message("DONE Script 05c")