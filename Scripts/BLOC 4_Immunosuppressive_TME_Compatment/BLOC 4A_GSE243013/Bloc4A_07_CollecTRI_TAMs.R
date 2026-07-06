#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc4A Script 07: CollecTRI TF activity on TAMs
# Infers transcription factor activity using decoupleR + CollecTRI
# Compares TF programs across TAM subtypes and response groups (MPR/non-MPR/pCR)
# Input:  Objects/Bloc4A_05_seu_TAMs_UCell.rds
# Output: Results/Figures/BLOC4A_TAMs/Bloc4A_CollecTRI_heatmap.png
#         Results/Figures/BLOC4A_TAMs/Bloc4A_CollecTRI_key_TFs_violin.png
#         Results/Tables/Bloc4A_CollecTRI_TF_activity.csv
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
DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
#DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ  <- file.path(DATA_DIR, "Objects/Bloc4A_05_seu_TAMs_UCell.rds")
OUT_FIG <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load TAMs UCell object
message("Loading TAMs UCell object...")
seu_TAMs <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_TAMs))
DefaultAssay(seu_TAMs) <- "RNA"

# Short labels for TAM subtypes , consistent with Bloc4A UCell script
# Added BEFORE downsampling so seu_sub inherits tam_short from seu_TAMs
short_labels <- c(
  "Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like)" = "Resident M2",
  "TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs)"        = "LAMs",
  "Inflammatory monocyte-derived TAMs (FCN1+/S100A8+)"                 = "Monocyte FCN1+",
  "Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high)"    = "Stress-response",
  "Proliferating TAMs (cycling/MKI67+)"                                = "Proliferating",
  "IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+)"      = "IFN-stimulated",
  "Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+)"           = "Classical-Mono"
)
seu_TAMs$tam_short <- unname(short_labels[as.character(seu_TAMs$final_annotation)])

# 2) Load CollecTRI network
message("Loading CollecTRI network...")
net <- read.csv(file.path(DATA_DIR, "Data/collectri_network.csv"))
message("CollecTRI: ", nrow(net), " interactions, ",
        length(unique(net$source)), " TFs")

seu_sub <- seu_TAMs
message("Using full TAMs population: ", ncol(seu_sub), " cells")

# 3) Extract normalized count matrix
# NOTE: TAMs object uses BPCells storage, requires explicit conversion
message("Extracting count matrix...")
mat <- GetAssayData(seu_sub, layer = "data")
mat <- as(mat, "dgCMatrix")

# 4) Run ULM (Univariate Linear Model) for TF activity inference
message("Running decoupleR ULM...")
tf_acts <- run_ulm(mat     = mat,
                   net     = net,
                   .source = "source",
                   .target = "target",
                   .mor    = "mor",
                   minsize = 5)
message("TF activity computed for ", length(unique(tf_acts$source)), " TFs")

saveRDS(tf_acts, file.path(DATA_DIR, "Objects/Bloc4A_07_tf_acts_raw.rds"))
message("Saved: Objects/Bloc4A_07_tf_acts_raw.rds")

# 5) Compute mean TF activity per TAM subtype x response group — for ALL TFs
message("Computing mean TF activity per TAM subtype x response group (all TFs)...")

tf_scores <- tf_acts %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score)

meta <- seu_sub@meta.data %>%
  select(response  = pathological_response,
         tam_state = tam_short) %>%
  tibble::rownames_to_column("condition")

tf_summary_all <- tf_scores %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(tam_state)) %>%
  group_by(source, tam_state, response) %>%
  summarise(mean_activity = mean(score), .groups = "drop")

# 6) Select TFs by BETWEEN-GROUP variance (variance of subtype x response means),
#    not per-cell variance, aligned with the project's goal of comparing
#    response groups across TAM subtypes.
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
       file.path(OUT_TAB, "Bloc4A_CollecTRI_TF_activity.csv"))
message("Saved: Bloc4A_CollecTRI_TF_activity.csv")

# 6b) Wilcoxon test per TAM subtype, MPR vs non-MPR, on raw per-cell scores
#     BH-corrected within subtype. Restricted to the Top 20 TFs and to
#     MPR/non-MPR (pCR reserved for perspectives, not part of the main text).
message("Running Wilcoxon tests MPR vs non-MPR, per TAM subtype, on raw TF scores...")

tf_scores_meta <- tf_scores %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(tam_state),
         response %in% c("MPR", "non-MPR"),
         source %in% top20_tfs)

tam_states_all <- unique(tf_scores_meta$tam_state)

wilcox_bysubtype <- lapply(tam_states_all, function(state) {
  lapply(top20_tfs, function(tf) {
    mpr    <- tf_scores_meta$score[tf_scores_meta$tam_state == state &
                                     tf_scores_meta$response  == "MPR" &
                                     tf_scores_meta$source    == tf]
    nonmpr <- tf_scores_meta$score[tf_scores_meta$tam_state == state &
                                     tf_scores_meta$response  == "non-MPR" &
                                     tf_scores_meta$source    == tf]
    if (length(mpr) < 3 | length(nonmpr) < 3) return(NULL)
    test <- wilcox.test(mpr, nonmpr)
    data.frame(tam_state     = state,
               source        = tf,
               p_value       = test$p.value,
               median_MPR    = median(mpr),
               median_nonMPR = median(nonmpr))
  }) %>% bind_rows()
}) %>% bind_rows()

wilcox_bysubtype <- wilcox_bysubtype %>%
  group_by(tam_state) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(tam_state, p_adj)

fwrite(wilcox_bysubtype,
       file.path(OUT_TAB, "Bloc4A_CollecTRI_TF_wilcox_bysubtype.csv"))
message("Saved: Bloc4A_CollecTRI_TF_wilcox_bysubtype.csv")
print(wilcox_bysubtype)

# 7) Build heatmap matrix, one per response group
green_palette <- colorRampPalette(c("white", "#006400"))(100)

for (resp in c("MPR", "non-MPR", "pCR")) {
  
  tf_heatmap_resp <- tf_summary %>%
    filter(response == resp) %>%
    select(source, tam_state, mean_activity) %>%
    pivot_wider(names_from = tam_state, values_from = mean_activity) %>%
    tibble::column_to_rownames("source") %>%
    as.matrix()
  
  col_order <- c("IFN-stimulated", "Resident M2", "Monocyte FCN1+",
                 "Stress-response", "Proliferating", "Classical-Mono", "LAMs")
  col_order <- col_order[col_order %in% colnames(tf_heatmap_resp)]
  tf_heatmap_resp <- tf_heatmap_resp[, col_order]
  
  png(file.path(OUT_FIG, paste0("Bloc4A_CollecTRI_heatmap_", resp, ".png")),
      width = 12, height = 10, units = "in", res = 300, type = "cairo")
  pheatmap(tf_heatmap_resp,
           scale        = "row",
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           color        = green_palette,
           main         = paste0("TF activity — TAMs GSE243013 LUAD — ", resp),
           fontsize_row = 11,
           fontsize_col = 13,
           angle_col    = 45,
           border_color = NA)
  dev.off()
  message("Saved: Bloc4A_CollecTRI_heatmap_", resp, ".png")
}

# 8) Violin plots
message("Generating violin plots for key TFs...")

tf_key <- tf_scores %>%
  filter(source %in% key_tfs_present) %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(tam_state))

plot_list <- lapply(key_tfs_present, function(tf) {
  tf_key %>%
    filter(source == tf) %>%
    ggplot(aes(x = response, y = score, fill = response)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = c("MPR"     = "#4393C3",
                                 "non-MPR" = "#D73027",
                                 "pCR"     = "#1A7A1A")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 11)) +
    labs(title = tf, y = "TF activity (ULM score)")
})

p_key <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc4A_CollecTRI_key_TFs_violin.png"),
       p_key, width = 14, height = 16, dpi = 300, bg = "white")
message("Saved: Bloc4A_CollecTRI_key_TFs_violin.png")

# 10) Save updated TAMs object
message("Saving updated TAMs object...")
saveRDS(seu_TAMs, file.path(DATA_DIR, "Objects/Bloc4A_07_seu_TAMs_TF.rds"))
message("Saved: Objects/Bloc4A_07_seu_TAMs_TF.rds")
message("DONE Bloc4A Script 07")