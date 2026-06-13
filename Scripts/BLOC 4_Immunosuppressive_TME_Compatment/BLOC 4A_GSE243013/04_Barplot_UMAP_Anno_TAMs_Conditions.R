#!/usr/bin/env Rscript
# ============================================================
# GSE243013  Bloc4A Script 04: TAMs final annotation + UMAP + Barplot
# 8 TAM subclusters identified from cluster 7 (M2-like) of global TME
# Cluster 6 excluded (lymphocyte contamination)
# Input:  Objects/Bloc4A_02_seu_TAMs_clustered.rds
# Output: Results/Figures/BLOC4A_TAMs/Bloc4A_UMAP_TAMs_annotation.png
#         Results/Figures/BLOC4A_TAMs/Bloc4A_UMAP_TAMs_split_response.png
#         Results/Figures/BLOC4A_TAMs/Bloc4A_Barplot_TAMs_proportions.png
#         Results/Tables/BLOC4A/Bloc4A_TAMs_chisq_test.csv
#         Results/Tables/BLOC4A/Bloc4A_TAMs_fisher_posthoc.csv
#         Objects/Bloc4A_04_seu_TAMs_annotated.rds
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
IN_OBJ  <- file.path(DATA_DIR, "Objects/Bloc4A_02_seu_TAM_clustered.rds")
OUT_FIG <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs")
OUT_TAB <- file.path(DATA_DIR, "Results/Tables/BLOC4A")

# 1) Load Seurat object
message("Loading Seurat object...")
seu <- readRDS(IN_OBJ)
seu_TAM <- seu
rm(seu); gc()
message("TAMs: ", ncol(seu_TAM))

# 2) Load consensus annotation and map final_annotation to Seurat object
message("Mapping final annotation...")
# Manual annotation, hardcoded from top20 markers analysis
# NOTE: cluster 6 excluded (lymphocyte contamination)
tam_annotation <- data.frame(
  cluster = 0:7,
  final_annotation = c(
    "Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like)",
    "TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs)",
    "Inflammatory monocyte-derived TAMs (FCN1+/S100A8+)",
    "Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high)",
    "Proliferating TAMs (cycling/MKI67+)",
    "IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+)",
    "EXCLUDE — lymphocyte contamination",
    "Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+)"
  )
)

# Map annotation to Seurat object
annot_map <- setNames(tam_annotation$final_annotation, tam_annotation$cluster)
test <- annot_map[as.character(seu_TAM$seurat_clusters)]
names(test) <- colnames(seu_TAM)
seu_TAM$final_annotation <- test

message("Annotation distribution:")
print(table(seu_TAM$final_annotation, useNA = "ifany"))


# 3) Define color palette
# Colors organized by cell lineage for visual coherence
tam_colors <- c(
  "Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like)" = "#D73027",  # rouge saturé
  "TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs)"        = "#A50026",  # rouge foncé saturé
  "Inflammatory monocyte-derived TAMs (FCN1+/S100A8+)"                 = "#FDAE61",  # orange pâle
  "Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high)"    = "#FEE090",  # jaune pâle
  "Proliferating TAMs (cycling/MKI67+)"                                = "#ABD9E9",  # bleu pâle
  "IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+)"      = "#4575B4",  # bleu saturé
  "Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+)"           = "#FDAE61"   # orange pâle
)

# Ensure factor order matches palette
seu_TAM$final_annotation <- factor(seu_TAM$final_annotation, levels = names(tam_colors))

# 4) UMAP : final annotation
message("Generating UMAP final annotation...")
seu_TAM <- subset(seu_TAM, 
                   subset = final_annotation != "EXCLUDE — lymphocyte contamination")
png(file.path(OUT_FIG, "Bloc4A_UMAP_TAMs_annotation.png"),
    width = 18, height = 10, units = "in", res = 300)
print(DimPlot(seu_TAM,
              reduction = "umap",
              group.by  = "final_annotation",
              cols      = tam_colors,
              label     = FALSE,
              pt.size   = 0.1,
              raster    = FALSE) +
        theme_classic() +
        labs(title = "TAM subtypes — GSE243013 LUAD") +
        theme(plot.title       = element_text(size = 14, face = "bold"),
              legend.text      = element_text(size = 14),
              legend.key.size  = unit(0.6, "cm"),
              legend.title     = element_blank()))
dev.off()
message("Saved: Bloc4A_UMAP_TAMs_annotation.png")

# UMAP split by pathological response
png(file.path(OUT_FIG, "Bloc4A_UMAP_TAMs_split_response.png"),
    width = 18, height = 10, units = "in", res = 300)
print(DimPlot(seu_TAM,
              reduction = "umap",
              group.by  = "final_annotation",
              split.by  = "pathological_response",
              cols      = tam_colors,
              label     = FALSE,
              pt.size   = 0.1,
              raster    = FALSE) +
        theme_classic() +
        labs(title = "TAM subtypes split by pathological response") +
        theme(plot.title      = element_text(size = 13, face = "bold"),
              strip.text      = element_text(size = 12, face = "bold"),
              legend.text     = element_text(size = 14),
              legend.key.size = unit(0.7, "cm"),
              legend.title    = element_blank()))
dev.off()
message("Saved: Bloc4A_UMAP_TAMs_split_response.png")

# 5) Barplot : proportions by response group
message("Computing proportions...")

df <- seu_TAM@meta.data %>%
  filter(!is.na(final_annotation), 
         !is.na(pathological_response),
         final_annotation != "EXCLUDE — lymphocyte contamination") %>%
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
  scale_fill_manual(values = tam_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.1)) +
  annotate("text",
           x        = 2,
           y        = 1.08,
           label    = "Chi-2 p < 2.2e-16",
           size     = 3.5,
           fontface = "italic") +
  ylab("TAM subtype proportion") +
  theme_classic() +
  theme(axis.title.x  = element_blank(),
        axis.text.x   = element_text(size = 13, face = "bold"),
        axis.text.y   = element_text(size = 11),
        axis.title.y  = element_text(size = 12),
        legend.position = "right",
        legend.text   = element_text(size = 12),
        legend.key.size = unit(0.6, "cm"),
        legend.title  = element_blank(),
        plot.caption  = element_text(size = 8, face = "italic", hjust = 0)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Bloc4A_Barplot_TAMs_proportions.png"),
       p_bar, width = 8, height = 7, dpi = 300, bg = "white")
message("Saved.")
message("Saved: Bloc4A_Barplot_TAMs_proportions.png")

# 6) Statistical tests
message("Running statistical tests...")

# Contingency table: cell types x response groups
cont_table <- table(seu_TAM$final_annotation, seu_TAM$pathological_response)

# Chi-2 global test
chisq_res <- chisq.test(cont_table)
message("Chi-2 test p-value: ", chisq_res$p.value)

chisq_summary <- data.frame(
  statistic = chisq_res$statistic,
  df        = chisq_res$parameter,
  p.value   = chisq_res$p.value
)
fwrite(chisq_summary, file.path(OUT_TAB, "Bloc4A_TAMs_chisq_test.csv"))
message("Saved: Bloc4A_TAMs_chisq_test.csv")

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
fwrite(fisher_df, file.path(OUT_TAB, "Bloc4A_TAMs_fisher_posthoc.csv"))
message("Saved: Bloc4A_TAMs_fisher_posthoc.csv")

# 7) Save updated Seurat object with final_annotation
message("Saving updated Seurat object...")
saveRDS(seu_TAM, file.path(DATA_DIR, "Objects/Bloc4A_04_seu_TAMs_annotated.rds"))
message("Saved: Objects/Bloc4A_04_seu_TAMs_annotated.rds")
message("DONE Bloc2 Script 06")