#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 06: CollecTRI TF activity on TAMs
# Infers transcription factor activity using decoupleR + CollecTRI
# Compares TF programs across TAM subtypes and response groups (MPR/NMPR)
# Input:  Objects/Bloc4B_04_seu_TAMs_combined.rds
# Output: Results/Figures/BLOC4B_TAMs/Bloc4B_GSE207422_CollecTRI_heatmap_MPR.png
#         Results/Figures/BLOC4B_TAMs/Bloc4B_GSE207422_CollecTRI_heatmap_NMPR.png
#         Results/Figures/BLOC4B_TAMs/Bloc4B_GSE207422_CollecTRI_key_TFs_violin.png
#         Results/Tables/Bloc4B_GSE207422_CollecTRI_TF_activity.csv
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
#DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ  <- file.path(DATA_DIR, "Objects/Bloc4B_04_seu_TAMs_combined.rds")
OUT_FIG <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB <- file.path(DATA_DIR, "Results/Tables")

# 1) Load TAMs combined object
message("Loading TAMs combined object...")
seu_TAMs <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_TAMs))
DefaultAssay(seu_TAMs) <- "RNA"

# Short labels for TAM subtypes — consistent with Bloc4B UCell script
# Added BEFORE downsampling so seu_sub inherits tam_short from seu_TAMs
short_labels <- c(
  "TAM_like_MRC1"                                            = "MRC1+ M2-like",
  "TAM_like_SPP1"                                            = "SPP1+ immunosuppressive",
  "TAM_like_resident_M2 (iron metabolism/anti-inflammatory)" = "Resident M2",
  "TAM_like_IFN (PD-L1+/IDO1+/CXCL9+)"                     = "IFN-stimulated",
  "TAM_like_monocyte (classical inflammatory)"               = "Monocyte-derived",
  "TAM_like_lipid (CCL18+/AKR+)"                            = "Lipid-associated",
  "TAM_like_stress (HSP-high/M1-like)"                      = "Stress-response",
  "TAM_like_regulatory (glucocorticoid-responsive)"          = "Regulatory",
  "TAM_like_M2 (SIGLEC8+/CCL18+)"                           = "M2-SIGLEC8+"
)
seu_TAMs$tam_short <- unname(short_labels[as.character(seu_TAMs$combined_annotation)])

# 2) Load CollecTRI network
message("Loading CollecTRI network...")
net <- read.csv(file.path(DATA_DIR, "Data/collectri_network.csv"))
message("CollecTRI: ", nrow(net), " interactions, ",
        length(unique(net$source)), " TFs")

seu_sub <- seu_TAMs
message("Using full TAMs population: ", ncol(seu_sub), " cells")

# 3) Extract normalized count matrix
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

saveRDS(tf_acts, file.path(DATA_DIR, "Objects/Bloc4B_06_tf_acts_raw.rds"))
message("Saved: Objects/Bloc4B_06_tf_acts_raw.rds")

# 5) Compute mean TF activity per TAM subtype x response group, for ALL TFs
#    (done BEFORE any TF selection, so that selection itself can be based on
#    between-group variance rather than per-cell variance)
message("Computing mean TF activity per TAM subtype x response group (all TFs)...")

meta <- seu_sub@meta.data %>%
  select(response  = PathResponse,
         tam_state = tam_short) %>%
  tibble::rownames_to_column("condition")

tf_scores <- tf_acts %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score)

tf_summary_all <- tf_scores %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(tam_state)) %>%
  group_by(source, tam_state, response) %>%
  summarise(mean_activity = mean(score), .groups = "drop")

# 6) Select TFs by BETWEEN-GROUP variance (variance of subtype x response means),
#    not per-cell variance, this is the criterion aligned with the project's
#    goal of comparing MPR vs NMPR response groups across TAM subtypes.
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
       file.path(OUT_TAB, "Bloc4B_GSE207422_CollecTRI_TF_activity.csv"))
message("Saved: Bloc4B_GSE207422_CollecTRI_TF_activity.csv")

# 6b) Wilcoxon test per TAM subtype, MPR vs NMPR, on raw per-cell scores
#     BH-corrected within subtype.
message("Running Wilcoxon tests MPR vs NMPR, per TAM subtype, on raw TF scores...")

tf_scores_meta <- tf_scores %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(response), !is.na(tam_state),
         source %in% top20_tfs)

tam_states_all <- unique(tf_scores_meta$tam_state)

wilcox_bysubtype <- lapply(tam_states_all, function(state) {
  lapply(top20_tfs, function(tf) {
    mpr  <- tf_scores_meta$score[tf_scores_meta$tam_state == state &
                                   tf_scores_meta$response  == "MPR" &
                                   tf_scores_meta$source    == tf]
    nmpr <- tf_scores_meta$score[tf_scores_meta$tam_state == state &
                                   tf_scores_meta$response  == "NMPR" &
                                   tf_scores_meta$source    == tf]
    if (length(mpr) < 3 | length(nmpr) < 3) return(NULL)
    test <- wilcox.test(mpr, nmpr)
    data.frame(tam_state  = state,
               source     = tf,
               p_value    = test$p.value,
               median_MPR = median(mpr),
               median_NMPR = median(nmpr))
  }) %>% bind_rows()
}) %>% bind_rows()

wilcox_bysubtype <- wilcox_bysubtype %>%
  group_by(tam_state) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(tam_state, p_adj)

fwrite(wilcox_bysubtype,
       file.path(OUT_TAB, "Bloc4B_GSE207422_CollecTRI_TF_wilcox_bysubtype.csv"))
message("Saved: Bloc4B_GSE207422_CollecTRI_TF_wilcox_bysubtype.csv")
print(wilcox_bysubtype)


# 7) Build heatmap matrix, one per response group
green_palette <- colorRampPalette(c("white", "#006400"))(100)

for (resp in c("MPR", "NMPR")) {
  
  tf_heatmap_resp <- tf_summary %>%
    filter(response == resp) %>%
    select(source, tam_state, mean_activity) %>%
    pivot_wider(names_from = tam_state, values_from = mean_activity) %>%
    tibble::column_to_rownames("source") %>%
    as.matrix()
  
  col_order <- c("IFN-stimulated", "Resident M2", "MRC1+ M2-like",
                 "SPP1+ immunosuppressive", "Monocyte-derived",
                 "Lipid-associated", "Stress-response",
                 "Regulatory", "M2-SIGLEC8+")
  col_order <- col_order[col_order %in% colnames(tf_heatmap_resp)]
  tf_heatmap_resp <- tf_heatmap_resp[, col_order]
  
  png(file.path(OUT_FIG, paste0("Bloc4B_GSE207422_CollecTRI_heatmap_", resp, ".png")),
      width = 14, height = 10, units = "in", res = 300)
  pheatmap(tf_heatmap_resp,
           scale        = "row",
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           color        = green_palette,
           main         = paste0("TF activity — TAMs GSE207422 NSCLC — ", resp),
           fontsize_row = 11,
           fontsize_col = 13,
           angle_col    = 45,
           border_color = NA)
  dev.off()
  message("Saved: Bloc4B_GSE207422_CollecTRI_heatmap_", resp, ".png")
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
    scale_fill_manual(values = c("MPR" = "#4393C3", "NMPR" = "#D73027")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(title = tf, y = "TF activity (ULM score)")
})

p_key <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc4B_GSE207422_CollecTRI_key_TFs_violin.png"),
       p_key, width = 14, height = 16, dpi = 300, bg = "white")
message("Saved: Bloc4B_GSE207422_CollecTRI_key_TFs_violin.png")