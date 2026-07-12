#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc2 Script 05: TME visualization with final annotation
# UMAP colored by final annotation + barplot proportions by response
# Input:  Objects/Bloc2_GSE207422_03_seu_sctype.rds
#         Results/Tables/Bloc2_GSE207422_consensus_annotation.csv
# Output: Results/Figures/Annotations_TME_GSE207422/Bloc2_GSE207422_UMAP_final_annotation.png
#         Results/Figures/Annotations_TME_GSE207422/Bloc2_GSE207422_UMAP_final_split_response.png
#         Results/Figures/Annotations_TME_GSE207422/Bloc2_GSE207422_Barplot_proportions.png
#         Results/Tables/Bloc2_GSE207422_chisq_test.csv
#         Results/Tables/Bloc2_GSE207422_fisher_posthoc.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(scales)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc2_GSE207422_03_seu_sctype.rds")
ANNOT    <- file.path(DATA_DIR, "Results/Tables/Bloc2_GSE207422_consensus_annotation.csv")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/Annotations_TME_GSE207422")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

message("Loading Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu))

message("Mapping final annotation...")
annot <- fread(ANNOT)
annot$cluster <- as.integer(annot$cluster)

annot_map <- setNames(annot$final_annotation, annot$cluster)
test <- annot_map[as.character(seu$seurat_clusters)]
names(test) <- colnames(seu)
seu$final_annotation <- test

message("Final annotation distribution:")
print(table(seu$final_annotation, useNA = "ifany"))

# Color palette, matching GSE207422's 19 final cell types, organized by lineage
tme_colors <- c(
  "CD8 T cells (effector/cytotoxic)"                      = "#f4a582",
  "CD8 T cells (exhausted/cytotoxic)"                      = "#D73027",
  "CD4 T cells (naive/memory)"                             = "#2171b5",
  "Tregs (tumor-infiltrating)"                             = "#08519c",
  "T cells (mixed CD4/CD8, low purity)"                    = "#6baed6",
  "NK cells (cytotoxic)"                                   = "#e6550d",
  "Proliferating cells (cycling)"                          = "#9970AB",
  "B cells"                                                = "#33a02c",
  "Plasma cells"                                           = "#006d2c",
  "Neutrophils"                                             = "#fdae6b",
  "Monocytes (inflammatory)"                               = "#a63603",
  "Alveolar macrophages (tissue-resident)"                 = "#7f2704",
  "Dendritic cells (cDC2)"                                 = "#BCBD22",
  "Plasmacytoid dendritic cells (pDCs)"                    = "#17BECF",
  "Mast cells"                                             = "#fdbf6f",
  "Tumor epithelial cells (squamous)"                      = "#993404",
  "Epithelial cells (AT2, normal)"                         = "#41ab5d",
  "Ciliated epithelial cells"                              = "#78c679",
  "Stromal cells (mesenchymal/endothelial, unresolved)"    = "#969696"
)

seu$final_annotation <- factor(seu$final_annotation, levels = names(tme_colors))

message("Generating UMAP final annotation...")

png(file.path(OUT_FIG, "Bloc2_GSE207422_UMAP_final_annotation.png"),
    width = 12, height = 8, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "final_annotation",
              cols      = tme_colors,
              label     = FALSE,
              pt.size   = 0.1,
              raster    = FALSE) +
        theme_classic() +
        labs(title = "TME global annotation — GSE207422 NSCLC") +
        theme(plot.title       = element_text(size = 14, face = "bold"),
              legend.text      = element_text(size = 9),
              legend.key.size  = unit(0.5, "cm"),
              legend.title     = element_blank()))
dev.off()
message("Saved: Bloc2_GSE207422_UMAP_final_annotation.png")

png(file.path(OUT_FIG, "Bloc2_GSE207422_UMAP_final_split_response.png"),
    width = 14, height = 6, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "final_annotation",
              split.by  = "PathResponse",
              cols      = tme_colors,
              label     = FALSE,
              pt.size   = 0.1,
              raster    = FALSE) +
        theme_classic() +
        labs(title = "TME annotation split by pathological response") +
        theme(plot.title      = element_text(size = 13, face = "bold"),
              strip.text      = element_text(size = 12, face = "bold"),
              legend.text     = element_text(size = 8),
              legend.key.size = unit(0.4, "cm"),
              legend.title    = element_blank()))
dev.off()
message("Saved: Bloc2_GSE207422_UMAP_final_split_response.png")

message("Computing proportions...")

df <- seu@meta.data %>%
  filter(!is.na(final_annotation), !is.na(PathResponse)) %>%
  transmute(response  = PathResponse,
            cell_type = final_annotation)

df_prop <- df %>%
  count(response, cell_type, name = "n") %>%
  group_by(response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_bar <- ggplot(df_prop,
                aes(x = response, y = prop, fill = cell_type)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = tme_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.1)) +
  ylab("Cell-type proportion") +
  theme_classic() +
  theme(axis.title.x  = element_blank(),
        axis.text.x   = element_text(size = 13, face = "bold"),
        axis.text.y   = element_text(size = 11),
        axis.title.y  = element_text(size = 12),
        legend.position = "right",
        legend.text   = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        legend.title  = element_blank()) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Bloc2_GSE207422_Barplot_proportions.png"),
       p_bar, width = 8, height = 7, dpi = 300, bg = "white")
message("Saved: Bloc2_GSE207422_Barplot_proportions.png")

message("Running statistical tests...")

cont_table <- table(seu$final_annotation, seu$PathResponse)

chisq_res <- chisq.test(cont_table)
message("Chi-2 test p-value: ", chisq_res$p.value)

chisq_summary <- data.frame(
  statistic = chisq_res$statistic,
  df        = chisq_res$parameter,
  p.value   = chisq_res$p.value
)
fwrite(chisq_summary, file.path(OUT_TAB, "Bloc2_GSE207422_chisq_test.csv"))
message("Saved: Bloc2_GSE207422_chisq_test.csv")

# NOTE: only MPR vs NMPR comparison here — GSE207422's corrected cohort
# has no pCR group (excluded, see Progress log)
fisher_test <- fisher.test(cont_table, simulate.p.value = TRUE, B = 10000)
fisher_df <- data.frame(
  comparison = "MPR vs NMPR",
  p.value    = fisher_test$p.value
)
print(fisher_df)
fwrite(fisher_df, file.path(OUT_TAB, "Bloc2_GSE207422_fisher_posthoc.csv"))
message("Saved: Bloc2_GSE207422_fisher_posthoc.csv")

message("Saving updated Seurat object...")
saveRDS(seu, file.path(DATA_DIR, "Objects/Bloc2_GSE207422_05_seu_final_annotated.rds"))
message("Saved: Objects/Bloc2_GSE207422_05_seu_final_annotated.rds")
message("DONE Bloc2 GSE207422 Script 05")