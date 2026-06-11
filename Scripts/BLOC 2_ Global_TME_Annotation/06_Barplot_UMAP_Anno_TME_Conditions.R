#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc2 Script 06: TME visualization with final annotation
# UMAP colored by final annotation + barplot proportions by response
# Input:  Objects/Bloc2_03_seu_sctype.rds
#         Results/Tables/Bloc2_consensus_annotation.csv
# Output: Results/Figures/Annotations/Bloc2_UMAP_final_annotation.png
#         Results/Figures/Annotations/Bloc2_UMAP_final_split_response.png
#         Results/Figures/Annotations/Bloc2_Barplot_proportions.png
#         Results/Tables/Bloc2_chisq_test.csv
#         Results/Tables/Bloc2_fisher_posthoc.csv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(scales)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc2_03_seu_sctype.rds")
ANNOT    <- file.path(DATA_DIR, "Results/Tables/Bloc2_consensus_annotation.csv")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/Annotations")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load Seurat object
message("Loading Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu))

# 2) Load consensus annotation and map final_annotation to Seurat object
message("Mapping final annotation...")
annot <- fread(ANNOT)

# NOTE: cluster column must match seurat_clusters type
annot$cluster <- as.integer(annot$cluster)

# Map final_annotation per cell
# NOTE: same barcode naming issue as sctype — must reassign names to cell barcodes
# Direct assignment fails — must use intermediate vector with barcode names
annot_map <- setNames(annot$final_annotation, annot$cluster)
test <- annot_map[as.character(seu$seurat_clusters)]
names(test) <- colnames(seu)
seu$final_annotation <- test

message("Final annotation distribution:")
print(table(seu$final_annotation, useNA = "ifany"))

# 3) Define color palette
# Adapted from project portfolio palette
# Colors organized by cell lineage for visual coherence
tme_colors <- c(
  # T cells — blues
  "CD4 T cells (naive/memory activated)"  = "#2171b5",
  "CD4 T cells (naive/resting)"           = "#a6cee3",
  "CD4 Tfh/exhausted CD4"                 = "#6baed6",
  "Tregs (tumor-infiltrating)"            = "#08519c",
  "CD8 T cells (exhausted/cytotoxic TRM)" = "#D73027",
  "CD8 T cells (effector memory)"         = "#f4a582",
  "NKT / gamma-delta T cells"             = "#fc8d59",
  "Proliferating cells (cycling)"         = "#9970AB",
  "IFN-stimulated cells (ISG-high)"       = "#c994c7",
  
  # NK — teal
  "NK cells (cytotoxic)"                  = "#e6550d",
  
  # B cells — greens
  "B cells"                               = "#33a02c",
  "B cells (naive/transitional)"          = "#b2df8a",
  "B cells (memory)"                      = "#74c476",
  "Plasma cells"                          = "#006d2c",
  
  # Myeloid — oranges
  "Monocytes (inflammatory)"              = "#e6550d",
  "Monocytes (classical CD14+)"           = "#fdae6b",
  "TAMs (M2-like/immunosuppressive)"      = "#a63603",
  "Dendritic cells (mregDC)"              = "#BCBD22",
  "Plasmacytoid Dendritic cells (pDCs)"   = "#17BECF",
  "Mast cells"                            = "#fdbf6f"
)

# Ensure factor order matches palette
seu$final_annotation <- factor(seu$final_annotation, levels = names(tme_colors))

# 4) UMAP : final annotation
message("Generating UMAP final annotation...")

png(file.path(OUT_FIG, "Bloc2_UMAP_final_annotation.png"),
    width = 12, height = 8, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "final_annotation",
              cols      = tme_colors,
              label     = FALSE,
              pt.size   = 0.1,
              raster    = FALSE) +
        theme_classic() +
        labs(title = "TME global annotation — GSE243013 LUAD") +
        theme(plot.title       = element_text(size = 14, face = "bold"),
              legend.text      = element_text(size = 9),
              legend.key.size  = unit(0.5, "cm"),
              legend.title     = element_blank()))
dev.off()
message("Saved: Bloc2_UMAP_final_annotation.png")

# UMAP split by pathological response
png(file.path(OUT_FIG, "Bloc2_UMAP_final_split_response.png"),
    width = 18, height = 6, units = "in", res = 300)
print(DimPlot(seu,
              reduction = "umap",
              group.by  = "final_annotation",
              split.by  = "pathological_response",
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
message("Saved: Bloc2_UMAP_final_split_response.png")

# 5) Barplot : proportions by response group
message("Computing proportions...")

df <- seu@meta.data %>%
  filter(!is.na(final_annotation), !is.na(pathological_response)) %>%
  transmute(response  = pathological_response,
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
  annotate("text",
           x        = 2,
           y        = 1.08,
           label    = "Chi-2 p < 2.2e-16",
           size     = 3.5,
           fontface = "italic") +
  ylab("Cell-type proportion") +
  labs(caption = "Note: ~5% cells excluded (clusters 5, 9, 13, 18 pending Azimuth confirmation)") +
  theme_classic() +
  theme(axis.title.x  = element_blank(),
        axis.text.x   = element_text(size = 13, face = "bold"),
        axis.text.y   = element_text(size = 11),
        axis.title.y  = element_text(size = 12),
        legend.position = "right",
        legend.text   = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        legend.title  = element_blank(),
        plot.caption  = element_text(size = 8, face = "italic", hjust = 0)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Bloc2_Barplot_proportions.png"),
       p_bar, width = 8, height = 7, dpi = 300, bg = "white")
message("Saved.")
message("Saved: Bloc2_Barplot_proportions.png")

# 6) Statistical tests
message("Running statistical tests...")

# Contingency table: cell types x response groups
cont_table <- table(seu$final_annotation, seu$pathological_response)

# Chi-2 global test
chisq_res <- chisq.test(cont_table)
message("Chi-2 test p-value: ", chisq_res$p.value)

chisq_summary <- data.frame(
  statistic = chisq_res$statistic,
  df        = chisq_res$parameter,
  p.value   = chisq_res$p.value
)
fwrite(chisq_summary, file.path(OUT_TAB, "Bloc2_chisq_test.csv"))
message("Saved: Bloc2_chisq_test.csv")

# Fisher post-hoc: pairwise comparisons between response groups
# NOTE: Bonferroni correction applied for multiple testing
comparisons <- list(
  c("MPR", "non-MPR"),
  c("pCR", "non-MPR"),
  c("MPR", "pCR")
)

fisher_results <- lapply(comparisons, function(pair) {
  sub_table <- cont_table[, pair]
  test <- fisher.test(sub_table, simulate.p.value = TRUE, B = 10000)
  data.frame(
    comparison = paste(pair, collapse = " vs "),
    p.value    = test$p.value,
    p.bonf     = min(test$p.value * length(comparisons), 1)
  )
})

fisher_df <- do.call(rbind, fisher_results)
print(fisher_df)
fwrite(fisher_df, file.path(OUT_TAB, "Bloc2_fisher_posthoc.csv"))
message("Saved: Bloc2_fisher_posthoc.csv")

# 7) Save updated Seurat object with final_annotation
message("Saving updated Seurat object...")
saveRDS(seu, file.path(DATA_DIR, "Objects/Bloc2_06_seu_final_annotated.rds"))
message("Saved: Objects/Bloc2_06_seu_final_annotated.rds")
message("DONE Bloc2 Script 06")