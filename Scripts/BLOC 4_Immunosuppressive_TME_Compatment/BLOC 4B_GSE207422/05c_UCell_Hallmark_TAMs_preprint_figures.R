#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 05c: UCell Hallmark TAMs preprint figures
# Biologically filtered discriminant signatures — heatmap with p-values
# Input:  Objects/GSE207422_seu_TAMs_Hallmark.rds
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Preprint/
#         Results/Tables/BLOC4B_Epithelial_TAMs/
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(data.table)   # FIX: required for fwrite()
})

DATA_DIR         <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG_PREPRINT <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/Preprint")
OUT_TAB          <- file.path(DATA_DIR, "Results/Tables/BLOC4B_Epithelial_TAMs")  # FIX: was undefined
dir.create(OUT_FIG_PREPRINT, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB,          recursive = TRUE, showWarnings = FALSE)

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

# 3) Calculate p-values MPR vs NMPR (global, all TAM subtypes pooled)
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

# 4) Build heatmap matrix (global, all TAM subtypes pooled)
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

# 5a) Global Heatmap
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

# 5b) Heatmap by subtype
# Mapping full names to short names
short_names_207 <- c(
  "TAM_like_IFN (PD-L1+/IDO1+/CXCL9+)"                       = "IFN-stimulated",
  "TAM_like_lipid (CCL18+/AKR+)"                              = "Lipid-associated",
  "TAM_like_M2 (SIGLEC8+/CCL18+)"                             = "M2-SIGLEC8+",
  "TAM_like_monocyte (classical inflammatory)"                 = "Monocyte-derived",
  "TAM_like_MRC1"                                              = "MRC1+ M2-like",
  "TAM_like_regulatory (glucocorticoid-responsive)"            = "Regulatory",
  "TAM_like_resident_M2 (iron metabolism/anti-inflammatory)"   = "Resident M2",
  "TAM_like_SPP1"                                              = "SPP1+",
  "TAM_like_stress (HSP-high/M1-like)"                         = "Stress-response"
)

# FIX: seu_TAM_207 -> seu_TAM (object was never defined under that name)
meta_subtype_207 <- seu_TAM@meta.data %>%
  select(subtype  = combined_annotation,
         response = PathResponse,
         all_of(sig_filtered)) %>%                 # FIX: sig_filtered_207 -> sig_filtered
  mutate(subtype = short_names_207[subtype])

# FIX: guard against silent NA mapping if a subtype label doesn't match
# short_names_207 exactly (typo, extra space, unmapped category, etc.)
if (any(is.na(meta_subtype_207$subtype))) {
  bad_labels <- unique(seu_TAM@meta.data$combined_annotation[
    is.na(short_names_207[seu_TAM@meta.data$combined_annotation])
  ])
  stop("Unmapped subtype label(s) in combined_annotation: ",
       paste(bad_labels, collapse = ", "))
}

mat_subtype_207 <- meta_subtype_207 %>%
  pivot_longer(cols = all_of(sig_filtered),        # FIX: sig_filtered_207 -> sig_filtered
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

# Wilcoxon by subtype
# NOTE: p.adjust (BH) is applied WITHIN each subtype (23 tests per subtype),
# not across the full 9 x 23 = 207 tests. This answers "does this signature
# discriminate MPR/NMPR within this specific subtype", not a global question.
# Document this choice explicitly in Methods so it's clear this differs from
# the global (all-subtypes-pooled) Wilcoxon table used elsewhere.
wilcox_subtype_207 <- lapply(unique(meta_subtype_207$subtype), function(st) {
  lapply(sig_filtered, function(sig) {              # FIX: sig_filtered_207 -> sig_filtered
    mpr  <- meta_subtype_207[meta_subtype_207$subtype == st &
                               meta_subtype_207$response == "MPR", sig]
    nmpr <- meta_subtype_207[meta_subtype_207$subtype == st &
                               meta_subtype_207$response == "NMPR", sig]
    mpr  <- mpr[!is.na(mpr)]
    nmpr <- nmpr[!is.na(nmpr)]
    if (length(mpr) < 3 | length(nmpr) < 3) return(NULL)
    test <- wilcox.test(mpr, nmpr)
    data.frame(subtype       = st,
               signature     = sig,
               p_value       = test$p.value,
               median_MPR    = median(mpr),
               median_NMPR   = median(nmpr))
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)

wilcox_subtype_207 <- wilcox_subtype_207 %>%
  group_by(subtype) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(subtype, p_adj)

fwrite(wilcox_subtype_207,
       file.path(OUT_TAB,
                 "GSE207422_UCell_Hallmark_TAMs_wilcox_bysubtype.csv"))
message("Saved: GSE207422 wilcox by subtype CSV")

# Figure

# Order columns
tam_states_207 <- c("Resident M2", "MRC1+ M2-like", "SPP1+",
                    "Monocyte-derived", "Lipid-associated",
                    "Stress-response", "IFN-stimulated",
                    "Regulatory", "M2-SIGLEC8+")
conditions_207 <- c("NMPR", "MPR")

col_order_207 <- as.vector(outer(tam_states_207, conditions_207,
                                 function(s, c) paste0(s, "\n", c)))
col_order_207 <- col_order_207[col_order_207 %in% colnames(mat_subtype_207)]
mat_subtype_207 <- mat_subtype_207[, col_order_207]

# Column annotation
col_anno_207 <- data.frame(
  Response  = sub(".*\n", "", colnames(mat_subtype_207)),
  row.names = colnames(mat_subtype_207))

anno_colors_207 <- list(
  Response = c("MPR" = "#4393C3", "NMPR" = "#D73027"))

png(file.path(OUT_FIG_PREPRINT,
              "GSE207422_Hallmark_TAMs_heatmap_bysubtype.png"),
    width = 16, height = 10, units = "in", res = 300)
pheatmap(mat_subtype_207,
         scale             = "row",
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         color             = green_palette,
         annotation_col    = col_anno_207,
         annotation_colors = anno_colors_207,
         main              = "Hallmark UCell scores — TAMs by subtype GSE207422",
         fontsize_row      = 10,
         fontsize_col      = 11,
         angle_col         = 45,
         border_color      = NA,
         gaps_col          = 9)
dev.off()
message("Saved: GSE207422_Hallmark_TAMs_heatmap_bysubtype.png")

message("DONE Script 05c")