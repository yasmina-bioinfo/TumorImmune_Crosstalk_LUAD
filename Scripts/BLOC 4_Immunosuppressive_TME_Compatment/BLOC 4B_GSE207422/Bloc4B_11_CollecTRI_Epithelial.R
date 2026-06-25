#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 11: CollecTRI TF activity on epithelial cells
# Infers transcription factor activity using decoupleR + CollecTRI
# Compares TF programs across 5 final groups:
# tumor_MPR, normal_MPR, tumor_NMPR, normal_NMPR, Ciliated
# Input:  Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/CollecTRI_Epithelial/
#         Results/Tables/Bloc4B_11_CollecTRI_Epithelial_TF_activity.csv
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

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/CollecTRI_Epithelial")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# 1) Load epithelial object with final groups
message("Loading epithelial object...")
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds"))
message("Cells: ", ncol(seu_Epi))
DefaultAssay(seu_Epi) <- "RNA"

# Keep only cells with final_group assigned
seu_Epi <- subset(seu_Epi, cells = colnames(seu_Epi)[!is.na(seu_Epi$final_group)])
message("Cells with final_group: ", ncol(seu_Epi))
print(table(seu_Epi$final_group))

# 2) Load CollecTRI network
message("Loading CollecTRI network...")
net <- read.csv(file.path(DATA_DIR, "Data/collectri_network.csv"))
message("CollecTRI: ", nrow(net), " interactions, ",
        length(unique(net$source)), " TFs")

# 3) Extract normalized count matrix
message("Extracting count matrix...")
mat <- GetAssayData(seu_Epi, layer = "data")
mat <- as(mat, "dgCMatrix")

# 4) Run ULM
message("Running decoupleR ULM...")
tf_acts <- run_ulm(mat     = mat,
                   net     = net,
                   .source = "source",
                   .target = "target",
                   .mor    = "mor",
                   minsize = 5)
message("TF activity computed for ", length(unique(tf_acts$source)), " TFs")

# 5) Select top 20 TFs by variance
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

# 6) Compute mean TF activity per final_group
message("Computing mean TF activity per final_group...")

meta <- seu_Epi@meta.data %>%
  dplyr::select(final_group) %>%
  tibble::rownames_to_column("condition")

tf_summary <- tf_scores %>%
  filter(source %in% top20_tfs) %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(final_group)) %>%
  group_by(source, final_group) %>%
  summarise(mean_activity = mean(score), .groups = "drop")

fwrite(as.data.frame(tf_summary),
       file.path(OUT_TAB, "Bloc4B_11_CollecTRI_Epithelial_TF_activity.csv"))
message("Saved: CSV")

# 7) Heatmap, 20 TFs x 5 groups
message("Generating heatmap...")

green_palette <- colorRampPalette(c("white", "#006400"))(100)

tf_heatmap <- tf_summary %>%
  pivot_wider(names_from = final_group, values_from = mean_activity) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# Order columns
col_order <- c("tumor_MPR", "normal_MPR", "tumor_NMPR", "normal_NMPR", "Ciliated")
col_order <- col_order[col_order %in% colnames(tf_heatmap)]
tf_heatmap <- tf_heatmap[, col_order]

png(file.path(OUT_FIG, "Bloc4B_11_CollecTRI_Epithelial_heatmap.png"),
    width = 12, height = 10, units = "in", res = 300)
pheatmap(tf_heatmap,
         scale        = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color        = green_palette,
         main         = "TF activity — Epithelial cells GSE207422",
         fontsize_row = 11,
         fontsize_col = 13,
         angle_col    = 45,
         border_color = NA)
dev.off()
message("Saved: Heatmap")

# 8) Violin plots : top 6 TFs by variance
message("Generating violin plots...")

key_tfs <- tf_var$source[1:6]
message("Top 6 TFs: ", paste(key_tfs, collapse = ", "))

group_colors <- c(
  "tumor_MPR"   = "#E63946",
  "normal_MPR"  = "#457B9D",
  "tumor_NMPR"  = "#F4A261",
  "normal_NMPR" = "#2A9D8F",
  "Ciliated"    = "#ADB5BD"
)

tf_key <- tf_scores %>%
  filter(source %in% key_tfs) %>%
  left_join(meta, by = "condition") %>%
  filter(!is.na(final_group))

tf_key$final_group <- factor(tf_key$final_group,
                             levels = c("tumor_MPR", "normal_MPR",
                                        "tumor_NMPR", "normal_NMPR",
                                        "Ciliated"))

# Figure 1  Tumor groups
tf_key_tumor <- tf_key %>% filter(final_group %in% c("tumor_MPR", "tumor_NMPR"))
tf_key_tumor$final_group <- factor(tf_key_tumor$final_group,
                                   levels = c("tumor_MPR", "tumor_NMPR"))

plot_list_tumor <- lapply(key_tfs, function(tf) {
  tf_key_tumor %>%
    filter(source == tf) %>%
    ggplot(aes(x = final_group, y = score, fill = final_group)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = group_colors) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
    labs(title = tf, y = "TF activity (ULM score)")
})

p_tumor <- wrap_plots(plot_list_tumor, ncol = 3)
ggsave(file.path(OUT_FIG, "Bloc4B_11_CollecTRI_Epithelial_top6_violin_tumor.png"),
       p_tumor, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: Violin tumor")

# Figure 2  Normal groups
tf_key_normal <- tf_key %>% filter(final_group %in% c("normal_MPR", "normal_NMPR", "Ciliated"))
tf_key_normal$final_group <- factor(tf_key_normal$final_group,
                                    levels = c("normal_MPR", "normal_NMPR", "Ciliated"))

plot_list_normal <- lapply(key_tfs, function(tf) {
  tf_key_normal %>%
    filter(source == tf) %>%
    ggplot(aes(x = final_group, y = score, fill = final_group)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = group_colors) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
    labs(title = tf, y = "TF activity (ULM score)")
})

p_normal <- wrap_plots(plot_list_normal, ncol = 3)
ggsave(file.path(OUT_FIG, "Bloc4B_11_CollecTRI_Epithelial_top6_violin_normal.png"),
       p_normal, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: Violin normal")


# 10) Save object
saveRDS(seu_Epi, file.path(DATA_DIR, "Objects/Bloc4B_11_seu_Epithelial_TF.rds"))
message("Saved: Objects/Bloc4B_11_seu_Epithelial_TF.rds")

message("DONE CollecTRI Epithelial")